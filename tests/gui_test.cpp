#include <main_window.h>
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
template<class T> T* widget(MainWindow& window, const char* name) {
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
