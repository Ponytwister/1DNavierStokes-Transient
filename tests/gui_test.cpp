#include <main_window.h>
#include <controls_dialog.h>
#include <model_controls.h>
#include <QCheckBox>
#include <QComboBox>
#include <gtest/gtest.h>
#include <QApplication>
#include <QElapsedTimer>
#include <QFile>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QPixmap>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <thread>

namespace {
template<class T> T* widget(QWidget& window, const char* name) {
    auto* result = window.findChild<T*>(name);
    if (!result) throw std::runtime_error(name);
    return result;
}
bool until(const std::function<bool()>& ready) {
    QElapsedTimer clock; clock.start();
    while (!ready() && clock.elapsed() < 15000) {
        QApplication::processEvents(); std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return ready();
}
struct Inputs {
    QTemporaryDir directory;
    QString database;
    Inputs() {
        if (!directory.isValid()) throw std::runtime_error("Temporary directory failed");
        database = directory.filePath("inputs.db");
        QFile fixture(QString::fromUtf8(TSENSOR_FIXTURE_DIR) + "/workflow.sql");
        if (!fixture.open(QIODevice::ReadOnly)) throw std::runtime_error("Fixture missing");
        execute(fixture.readAll());
    }
    double execute(const QByteArray& query) {
        sqlite3* raw = nullptr;
        if (sqlite3_open(database.toUtf8().constData(), &raw) != SQLITE_OK) throw std::runtime_error("Open failed");
        std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
        double value = 0;
        auto callback = [](void* ptr, int, char** values, char**) {
            if (values[0]) *static_cast<double*>(ptr) = std::stod(values[0]);
            return 0;
        };
        if (sqlite3_exec(db.get(), query.constData(), callback, &value, nullptr) != SQLITE_OK)
            throw std::runtime_error(sqlite3_errmsg(db.get()));
        return value;
    }
    void choose(MainWindow& window) {
        widget<QLineEdit>(window, "databasePath")->setText(database);
        widget<QLineEdit>(window, "outputPath")->setText(directory.filePath("reports"));
    }
};

int controlRow(QTableWidget* table, const QString& name) {
    for (int i = 0; i < table->rowCount(); ++i)
        if (table->item(i, 0)->text() == name) return i;
    throw std::runtime_error("Control row missing");
}
model_controls::Row& control(std::vector<model_controls::Row>& rows, const QString& name) {
    for (auto& row : rows) if (row.name == name) return row;
    throw std::runtime_error("Control missing");
}

TEST(ModelControls, SupportsRealColumnNamesNullsAndPreservesOtherData)
{
    Inputs input;
    input.execute("CREATE TABLE controls_real(Parameter TEXT NOT NULL, Setting ANY);"
                  "INSERT INTO controls_real SELECT * FROM model_controls;"
                  "DROP TABLE model_controls; ALTER TABLE controls_real RENAME TO model_controls;"
                  "INSERT INTO model_controls VALUES('run_solver',NULL),('save_normalized_profiles','false')");
    auto original = model_controls::load(input.database);
    EXPECT_EQ(original.nameColumn, "Parameter");
    EXPECT_FALSE(control(original.rows, "run_solver").value.has_value());
    auto edited = original.rows;
    control(edited, "run_solver").value = "false";
    control(edited, "max_iterations").value = "12";
    model_controls::save(input.database, original, edited);
    auto loaded = model_controls::load(input.database);
    EXPECT_EQ(loaded.rows, edited);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    control(edited, "run_solver").value.reset();
    model_controls::save(input.database, loaded, edited);
    loaded = model_controls::load(input.database);
    EXPECT_FALSE(control(loaded.rows, "run_solver").value.has_value());
}

TEST(ModelControls, ValidationConflictAndTransactionRollback)
{
    Inputs input;
    auto original = model_controls::load(input.database);
    auto edited = original.rows;
    control(edited, "width resolution (X)").value = "0";
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_THROW(model_controls::validate({"convergence_epsx", "nan"}), std::runtime_error);
    EXPECT_THROW(model_controls::validate({"max_iterations", "3abc"}), std::runtime_error);
    EXPECT_THROW(model_controls::validate({"run_solver", "yes"}), std::runtime_error);
    EXPECT_NO_THROW(model_controls::validate({"convergence_epsx", "1e-12"}));
    edited = original.rows;
    control(edited, "debug_level").value = "4";
    control(edited, "max_iterations").value = "12";
    input.execute("CREATE TRIGGER reject_control BEFORE UPDATE ON model_controls "
                  "WHEN NEW.criterion='max_iterations' BEGIN SELECT RAISE(ABORT,'test failure'); END");
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original); // Earlier debug_level update rolled back.
    input.execute("DROP TRIGGER reject_control; UPDATE model_controls SET value='5' WHERE criterion='debug_level'");
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='debug_level'"), 5);
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 3);
    EXPECT_THROW(model_controls::load(input.directory.filePath("missing.db")), std::runtime_error);
    EXPECT_FALSE(QFile::exists(input.directory.filePath("missing.db")));
}

TEST(ModelControls, UnknownRowsAndAmbiguousNamesCannotBeEdited)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('future_control','keep me')");
    auto original = model_controls::load(input.database);
    auto edited = original.rows;
    control(edited, "future_control").value = "changed";
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    input.execute("CREATE TABLE duplicate_controls(Parameter TEXT, Setting ANY);"
                  "INSERT INTO duplicate_controls VALUES('debug_level',1),('debug_level',2);"
                  "DROP TABLE model_controls; ALTER TABLE duplicate_controls RENAME TO model_controls");
    EXPECT_THROW(model_controls::load(input.database), std::runtime_error);
}

TEST(Gui, ControlsCancelValidationAndSaveBeforeRun)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('run_solver','true'),('save_normalized_profiles','false')");
    {
        ControlsDialog dialog(input.database); dialog.show();
        auto* save = widget<QPushButton>(dialog, "saveControlsButton");
        ASSERT_TRUE(until([&] { return save->isEnabled(); }));
        auto* table = widget<QTableWidget>(dialog, "controlsTable");
        qobject_cast<QLineEdit*>(table->cellWidget(controlRow(table, "max_iterations"), 1))->setText("15");
        EXPECT_FALSE(table->cellWidget(controlRow(table, "save_normalized_profiles"), 1)->isEnabled());
        if (const auto capture = qEnvironmentVariable("NAVIER_CONTROLS_CAPTURE"); !capture.isEmpty()) {
            QApplication::processEvents(); EXPECT_TRUE(dialog.grab().save(capture));
        }
        dialog.reject();
        EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 3);
    }
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(widget<QPushButton>(window, "exportButton")->isEnabled());
    QTimer::singleShot(0, &window, [&] {
        auto* dialog = dynamic_cast<ControlsDialog*>(QApplication::activeModalWidget());
        if (!dialog) { ADD_FAILURE() << "Controls dialog missing"; return; }
        auto* save = widget<QPushButton>(*dialog, "saveControlsButton");
        if (!until([&] { return save->isEnabled(); })) { ADD_FAILURE() << "Loading timed out"; dialog->reject(); return; }
        EXPECT_FALSE(run->isEnabled());
        auto* table = widget<QTableWidget>(*dialog, "controlsTable");
        auto* iterations = qobject_cast<QLineEdit*>(table->cellWidget(controlRow(table, "max_iterations"), 1));
        iterations->setText("bad"); save->click();
        EXPECT_TRUE(widget<QLabel>(*dialog, "controlsStatus")->text().contains("max_iterations"));
        EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 3);
        iterations->setText("15");
        qobject_cast<QComboBox*>(table->cellWidget(controlRow(table, "run_solver"), 1))->setCurrentText("false");
        save->click();
        EXPECT_FALSE(save->isEnabled());
    });
    widget<QPushButton>(window, "modelControlsButton")->click();
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 15);
    EXPECT_FALSE(widget<QPushButton>(window, "exportButton")->isEnabled());
    EXPECT_EQ(widget<QTableWidget>(window, "parameterTable")->rowCount(), 0);
    run->click();
    EXPECT_FALSE(widget<QPushButton>(window, "modelControlsButton")->isEnabled());
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Optimizer iterations: Not run"));
}

TEST(Gui, RunExportAndExplicitSaves)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    auto* exportButton = widget<QPushButton>(window, "exportButton");
    EXPECT_FALSE(exportButton->isEnabled());
    run->click();
    EXPECT_FALSE(widget<QLineEdit>(window, "databasePath")->isEnabled());
    EXPECT_FALSE(run->isEnabled());
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(exportButton->isEnabled());
    EXPECT_EQ(widget<QTableWidget>(window, "parameterTable")->rowCount(), 1);
    EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Termination code:"));
    if (const auto capture = qEnvironmentVariable("NAVIER_GUI_CAPTURE"); !capture.isEmpty()) {
        QApplication::processEvents();
        EXPECT_TRUE(window.grab().save(capture));
    }
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    exportButton->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("reports/uniform,.txt")));
    widget<QPushButton>(window, "profilesButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 36);
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    QTimer::singleShot(0, [] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget()))
            box->button(QMessageBox::Yes)->click();
    });
    widget<QPushButton>(window, "inputsButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), 1);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    // Changing databases must not leave actions attached to an old session.
    widget<QLineEdit>(window, "databasePath")->setText("another.db");
    EXPECT_FALSE(exportButton->isEnabled());
    EXPECT_EQ(widget<QTableWidget>(window, "parameterTable")->rowCount(), 0);
}

TEST(Gui, CancelAndCloseDuringCalculation)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='100000' WHERE criterion='length/time resolution (Z)'");
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    bool heartbeat = false;
    run->click();
    QTimer::singleShot(0, &window, [&] {
        heartbeat = true;
        EXPECT_FALSE(run->isEnabled());
        widget<QPushButton>(window, "cancelButton")->click();
        // Programmatic changes can refresh controls even while inputs are locked.
        widget<QLineEdit>(window, "outputPath")->setText(input.directory.filePath("cancelled"));
        EXPECT_FALSE(widget<QPushButton>(window, "cancelButton")->isEnabled());
    });
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(heartbeat);
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("cancelled"));
    EXPECT_FALSE(widget<QPushButton>(window, "exportButton")->isEnabled());
    run->click();
    EXPECT_FALSE(window.close()); // Close deferred while the worker stops.
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
}

TEST(Gui, FailureCanBeCorrectedAndRetried)
{
    Inputs input;
    input.execute("UPDATE alglib_input SET SCALE=0");
    MainWindow window; input.choose(window);
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_FALSE(widget<QPushButton>(window, "inputsButton")->isEnabled());
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("zero passed to SCALE"));
    input.execute("UPDATE alglib_input SET SCALE=1");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(widget<QPushButton>(window, "exportButton")->isEnabled());
    widget<QLineEdit>(window, "outputPath")->setText(input.database); // File, not directory.
    widget<QPushButton>(window, "exportButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("export", Qt::CaseInsensitive));
    widget<QLineEdit>(window, "outputPath")->setText(input.directory.filePath("retry"));
    widget<QPushButton>(window, "exportButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("retry/uniform,.txt")));
}
TEST(Gui, DisabledProfileSaveAndCloseDuringExport)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='false' WHERE criterion='save_model_profiles'");
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_FALSE(widget<QPushButton>(window, "profilesButton")->isEnabled());
    ASSERT_TRUE(widget<QPushButton>(window, "exportButton")->isEnabled());
    widget<QPushButton>(window, "exportButton")->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("reports/uniform,.txt")));
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
}

TEST(Gui, FailedExportAbortsCloseAndAllowsRetry)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    auto* exportButton = widget<QPushButton>(window, "exportButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    widget<QLineEdit>(window, "outputPath")->setText(input.database);
    exportButton->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return run->isEnabled() || !window.isVisible(); }));
    ASSERT_TRUE(window.isVisible());
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("export", Qt::CaseInsensitive));
    ASSERT_TRUE(exportButton->isEnabled());
    widget<QLineEdit>(window, "outputPath")->setText(input.directory.filePath("retry"));
    exportButton->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("retry/uniform,.txt")));
}
} // namespace

int main(int argc, char** argv) {
    QApplication application(argc, argv);
    QApplication::setQuitOnLastWindowClosed(false);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
