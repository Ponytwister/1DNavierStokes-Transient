#include <main_window.h>
#include <controls_dialog.h>
#include <model_controls.h>
#include <setup_file.h>
#include <QCheckBox>
#include <QComboBox>
#include <QListWidget>
#include <QToolButton>
#include <QMenu>
#include <QAction>
#include <QMenuBar>
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
#include <QTabWidget>
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
        execute("UPDATE model_controls SET value='3' WHERE criterion='debug_level'");
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

model_controls::Row& control(std::vector<model_controls::Row>& rows, const QString& name) {
    for (auto& row : rows) if (row.name == name) return row;
    throw std::runtime_error("Control missing");
}

TEST(Gui, FileMenuContainsSetupAndPathControls)
{
    MainWindow window;
    auto* file = window.menuBar()->actions().front()->menu();
    ASSERT_NE(file, nullptr);
    EXPECT_TRUE(file->actions().contains(widget<QAction>(window, "openSetupAction")));
    EXPECT_TRUE(file->actions().contains(widget<QAction>(window, "saveSetupAction")));
    EXPECT_TRUE(widget<QAction>(window, "openSetupAction")->isEnabled());
    EXPECT_FALSE(widget<QAction>(window, "saveSetupAction")->isEnabled());
    EXPECT_EQ(window.centralWidget()->findChild<QLineEdit*>("databasePath"), nullptr);
    EXPECT_EQ(window.centralWidget()->findChild<QLineEdit*>("outputPath"), nullptr);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    EXPECT_EQ(tabs->tabText(1), "Model controls");
    EXPECT_FALSE(tabs->isTabEnabled(1));
    Inputs input; input.choose(window);
    EXPECT_TRUE(tabs->isTabEnabled(1));
    EXPECT_TRUE(widget<QAction>(window, "saveSetupAction")->isEnabled());
}

TEST(SetupFile, RoundTripAndTransactionalRestore)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('future_control',NULL)");
    const auto original = model_controls::load(input.database);
    const auto filename = input.directory.filePath(QString::fromUtf8("saved setup ü.navier.json"));
    setup_file::Setup setup{input.database, input.directory.filePath("output with spaces"), original.rows};
    setup_file::save(filename, setup);
    const auto loaded = setup_file::load(filename);
    EXPECT_EQ(loaded.database, setup.database);
    EXPECT_EQ(loaded.outputDirectory, setup.outputDirectory);
    EXPECT_EQ(loaded.controls, original.rows);
    input.execute("UPDATE model_controls SET value='7' WHERE criterion='max_iterations'");
    const auto changed = model_controls::load(input.database);
    setup_file::restore(loaded, changed);
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    EXPECT_THROW(setup_file::restore(loaded, changed), std::runtime_error);
    auto incompatible = loaded;
    incompatible.controls.pop_back();
    EXPECT_THROW(setup_file::restore(incompatible, original), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_THROW(setup_file::save(input.database, setup), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    setup.controls.push_back(setup.controls.front());
    EXPECT_THROW(setup_file::save(filename, setup), std::runtime_error);
    EXPECT_EQ(setup_file::load(filename).controls, original.rows);
}

TEST(SetupFile, ValidatesFileAndResolvesRelativePaths)
{
    QTemporaryDir directory;
    const auto filename = directory.filePath("setup.json");
    auto write = [&](const QByteArray& bytes) {
        QFile file(filename); ASSERT_TRUE(file.open(QIODevice::WriteOnly)); file.write(bytes);
    };
    write(R"({"format":"navier-setup","version":1,"database":"inputs.db","outputDirectory":"reports","controls":[{"name":"max_iterations","value":"3"}]})");
    const auto loaded = setup_file::load(filename);
    EXPECT_EQ(loaded.database, directory.filePath("inputs.db"));
    EXPECT_EQ(loaded.outputDirectory, directory.filePath("reports"));
    for (const auto& bytes : {
        QByteArray("{"),
        QByteArray(R"({"format":"navier-setup","version":2})"),
        QByteArray(R"({"format":"navier-setup","version":1,"database":"a","outputDirectory":"b","controls":[{"name":"max_iterations","value":"-1"}]})"),
        QByteArray(R"({"format":"navier-setup","version":1,"database":"a","outputDirectory":"b","controls":[{"name":"max_iterations"}]})"),
        QByteArray(R"({"format":"navier-setup","version":1,"database":"a","outputDirectory":"b","controls":[{"name":"x","value":null},{"name":"x","value":""}]})")}) {
        write(bytes);
        EXPECT_THROW(setup_file::load(filename), std::runtime_error);
    }
    EXPECT_THROW(setup_file::load(directory.filePath("missing.json")), std::runtime_error);
    EXPECT_THROW(setup_file::save(directory.filePath("missing/setup.json"), loaded), std::runtime_error);
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
    control(edited, "universal_solve_for").value.reset();
    model_controls::save(input.database, loaded, edited);
    loaded = model_controls::load(input.database);
    EXPECT_FALSE(control(loaded.rows, "universal_solve_for").value.has_value());
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

TEST(ModelControls, RequiredValuesAndStrictNumericBounds)
{
    for (const auto* value : {"", "-1", "7", "1.5", "nan"})
        EXPECT_THROW(model_controls::validate({"debug_level", value}), std::runtime_error);
    for (const auto* value : {"0", "6"})
        EXPECT_NO_THROW(model_controls::validate({"debug_level", value}));
    for (const auto* value : {"", "0", "-1e-9", "1e-3", "1", "nan", "inf"})
        EXPECT_THROW(model_controls::validate({"convergence_epsx", value}), std::runtime_error);
    EXPECT_NO_THROW(model_controls::validate({"convergence_epsx", "9.99e-4"}));
    for (const auto* name : {"debug_level", "experiment_name", "run_solver", "convergence_epsx"})
        EXPECT_THROW(model_controls::validate({name, std::nullopt}), std::runtime_error);
    EXPECT_NO_THROW(model_controls::validate({"universal_solve_for", std::nullopt}));
    EXPECT_NO_THROW(model_controls::validate({"universal_solve_for", ""}));
}

TEST(Gui, InvalidUnchangedFieldsBlockBothActions)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value=NULL WHERE criterion='experiment_name'");
    ControlsDialog dialog(input.database); dialog.show();
    auto* use = widget<QPushButton>(dialog, "useControlsButton");
    auto* save = widget<QPushButton>(dialog, "saveControlsButton");
    ASSERT_TRUE(until([&] { return use->isEnabled(); }));
    const auto original = model_controls::load(input.database);
    for (auto* button : {use, save}) {
        button->click();
        EXPECT_TRUE(dialog.isVisible());
        EXPECT_TRUE(widget<QLabel>(dialog, "controlsStatus")->text().contains("experiment_name"));
        EXPECT_EQ(model_controls::load(input.database), original);
    }
    widget<QListWidget>(dialog, "experimentChoices")->item(0)->setCheckState(Qt::Checked);
    auto* parameters = widget<QListWidget>(dialog, "parameterChoices");
    for (int i = 0; i < parameters->count(); ++i) parameters->item(i)->setCheckState(Qt::Unchecked);
    use->click();
    EXPECT_EQ(dialog.result(), QDialog::Accepted);
    auto values = dialog.values().rows;
    EXPECT_FALSE(control(values, "universal_solve_for").value);
}

TEST(Gui, ExperimentChecklistSupportsMultipleSelectionsAndExplicitPersistence)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME) VALUES('second'),('third')");
    const auto original = model_controls::load(input.database);
    auto choose = [](ControlsDialog& dialog, const QString& name, Qt::CheckState state) {
        auto* list = widget<QListWidget>(dialog, "experimentChoices");
        const auto items = list->findItems(name, Qt::MatchExactly);
        ASSERT_EQ(items.size(), 1);
        items.front()->setCheckState(state);
    };
    {
        ControlsDialog dialog(input.database); dialog.show();
        ASSERT_TRUE(until([&] { return widget<QPushButton>(dialog, "useControlsButton")->isEnabled(); }));
        EXPECT_EQ(widget<QListWidget>(dialog, "experimentChoices")->count(), 3);
        EXPECT_EQ(dialog.findChild<QLineEdit*>("experiment_name"), nullptr);
        choose(dialog, "second", Qt::Checked);
        dialog.reject();
        EXPECT_EQ(model_controls::load(input.database), original);
    }
    ControlsDialog dialog(input.database); dialog.show();
    ASSERT_TRUE(until([&] { return widget<QPushButton>(dialog, "useControlsButton")->isEnabled(); }));
    choose(dialog, "second", Qt::Checked);
    EXPECT_TRUE(widget<QToolButton>(dialog, "experiment_name")->text().contains("second"));
    if (const auto capture = qEnvironmentVariable("NAVIER_EXPERIMENTS_CAPTURE"); !capture.isEmpty()) {
        auto* picker = widget<QToolButton>(dialog, "experiment_name");
        picker->menu()->popup(picker->mapToGlobal(QPoint(0, picker->height())));
        QApplication::processEvents();
        EXPECT_TRUE(picker->menu()->grab().save(capture));
        picker->menu()->hide();
    }
    widget<QPushButton>(dialog, "useControlsButton")->click();
    ASSERT_EQ(dialog.result(), QDialog::Accepted);
    auto memory = dialog.values();
    EXPECT_EQ(control(memory.rows, "experiment_name").value, "uniform second");
    EXPECT_EQ(model_controls::load(input.database), original);
    ControlsDialog reopened(input.database, nullptr, memory); reopened.show();
    auto* save = widget<QPushButton>(reopened, "saveControlsButton");
    ASSERT_TRUE(until([&] { return save->isEnabled(); }));
    choose(reopened, "uniform", Qt::Unchecked);
    choose(reopened, "second", Qt::Unchecked);
    for (auto* button : {save, widget<QPushButton>(reopened, "useControlsButton")}) {
        button->click(); EXPECT_TRUE(reopened.isVisible());
        EXPECT_TRUE(widget<QLabel>(reopened, "controlsStatus")->text().contains("experiment_name"));
    }
    choose(reopened, "second", Qt::Checked);
    choose(reopened, "third", Qt::Checked);
    save->click();
    ASSERT_TRUE(until([&] { return reopened.result() == QDialog::Accepted; }));
    auto stored = model_controls::load(input.database);
    EXPECT_EQ(control(stored.rows, "experiment_name").value, "second third");
}

TEST(Gui, MissingSavedExperimentMustBeDeselected)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='missing uniform' WHERE criterion='experiment_name'");
    ControlsDialog dialog(input.database); dialog.show();
    auto* use = widget<QPushButton>(dialog, "useControlsButton");
    ASSERT_TRUE(until([&] { return use->isEnabled(); }));
    use->click(); EXPECT_TRUE(dialog.isVisible());
    EXPECT_TRUE(widget<QLabel>(dialog, "controlsStatus")->text().contains("missing"));
    widget<QListWidget>(dialog, "experimentChoices")->findItems("missing", Qt::MatchExactly).front()->setCheckState(Qt::Unchecked);
    use->click(); EXPECT_EQ(dialog.result(), QDialog::Accepted);
}

TEST(ModelControls, ParameterChoicesDependOnAnyExperimentsReactionCount)
{
    Inputs input;
    const QStringList base{"p1", "kon1", "keq1", "left_edge", "width", "QE1"};
    EXPECT_EQ(model_controls::solvableParameters(input.database), base);
    input.execute("INSERT INTO experiments(NAME,REACTIONS) VALUES('other',NULL)");
    EXPECT_EQ(model_controls::solvableParameters(input.database), base);
    input.execute("UPDATE experiments SET REACTIONS='  first  ' WHERE NAME='other'");
    EXPECT_EQ(model_controls::solvableParameters(input.database), base);
    input.execute("UPDATE experiments SET REACTIONS='first second' WHERE NAME='other'");
    auto expanded = base; expanded.append({"p2", "kon2", "keq2", "QE2"});
    EXPECT_EQ(model_controls::solvableParameters(input.database), expanded);
}

TEST(Gui, ParameterChecklistSupportsMemoryDefaultsAndEmptySelection)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME,REACTIONS) VALUES('other','first second')");
    const auto original = model_controls::load(input.database);
    ControlsDialog dialog(input.database); dialog.show();
    auto* use = widget<QPushButton>(dialog, "useControlsButton");
    ASSERT_TRUE(until([&] { return use->isEnabled(); }));
    EXPECT_EQ(dialog.findChild<QLineEdit*>("universal_solve_for"), nullptr);
    auto* choices = widget<QListWidget>(dialog, "parameterChoices");
    ASSERT_EQ(choices->count(), 10);
    choices->findItems("kon2", Qt::MatchExactly).front()->setCheckState(Qt::Checked);
    use->click(); ASSERT_EQ(dialog.result(), QDialog::Accepted);
    auto memory = dialog.values();
    EXPECT_EQ(control(memory.rows, "universal_solve_for").value, "keq1 kon2");
    EXPECT_EQ(model_controls::load(input.database), original);
    ControlsDialog reopened(input.database, nullptr, memory); reopened.show();
    auto* save = widget<QPushButton>(reopened, "saveControlsButton");
    ASSERT_TRUE(until([&] { return save->isEnabled(); }));
    auto* storedChoices = widget<QListWidget>(reopened, "parameterChoices");
    EXPECT_EQ(storedChoices->findItems("kon2", Qt::MatchExactly).front()->checkState(), Qt::Checked);
    save->click(); ASSERT_TRUE(until([&] { return reopened.result() == QDialog::Accepted; }));
    auto stored = model_controls::load(input.database);
    EXPECT_EQ(control(stored.rows, "universal_solve_for").value, "keq1 kon2");
    ControlsDialog empty(input.database); empty.show();
    auto* emptySave = widget<QPushButton>(empty, "saveControlsButton");
    ASSERT_TRUE(until([&] { return emptySave->isEnabled(); }));
    auto* emptyChoices = widget<QListWidget>(empty, "parameterChoices");
    for (int i = 0; i < emptyChoices->count(); ++i) emptyChoices->item(i)->setCheckState(Qt::Unchecked);
    emptySave->click(); ASSERT_TRUE(until([&] { return empty.result() == QDialog::Accepted; }));
    stored = model_controls::load(input.database);
    EXPECT_FALSE(control(stored.rows, "universal_solve_for").value);
}

TEST(Gui, ControlsCancelValidationAndSaveBeforeRun)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('run_solver','true'),('save_normalized_profiles','false'),"
                  "('exp_left_padding','0'),('exp_right_padding','0'),('disable_reactions','false'),"
                  "('disable_reverse_reactions','false'),('use_alglib_init_values','true'),('convergence_epsx','1e-12')");
    {
        ControlsDialog dialog(input.database); dialog.show();
        auto* save = widget<QPushButton>(dialog, "saveControlsButton");
        ASSERT_TRUE(until([&] { return save->isEnabled(); }));

        widget<QLineEdit>(dialog, "max_iterations")->setText("15");
        EXPECT_EQ(dialog.findChild<QCheckBox*>("save_normalized_profiles"), nullptr);
        EXPECT_EQ(dialog.findChild<QCheckBox*>("save_model_profiles"), nullptr);
        EXPECT_EQ(dialog.findChild<QTableWidget*>(), nullptr);
        EXPECT_EQ(dialog.findChildren<QCheckBox*>().size(), 4);
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
        auto* dialog = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
        if (!dialog) { ADD_FAILURE() << "Controls dialog missing"; return; }
        auto* save = widget<QPushButton>(*dialog, "saveControlsButton");
        if (!until([&] { return save->isEnabled(); })) { ADD_FAILURE() << "Loading timed out"; dialog->reject(); return; }
        EXPECT_FALSE(run->isEnabled());
        EXPECT_FALSE(widget<QAction>(window, "openSetupAction")->isEnabled());
        EXPECT_FALSE(widget<QAction>(window, "saveSetupAction")->isEnabled());

        auto* iterations = widget<QLineEdit>(*dialog, "max_iterations");
        iterations->setText("bad"); save->click();
        EXPECT_TRUE(widget<QLabel>(*dialog, "controlsStatus")->text().contains("max_iterations"));
        EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 3);
        iterations->setText("15");
        widget<QCheckBox>(*dialog, "run_solver")->setChecked(false);
        save->click();
        EXPECT_FALSE(save->isEnabled());
    });
    widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(1);
    ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "runButton")->isEnabled(); }));
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 15);
    EXPECT_FALSE(widget<QPushButton>(window, "exportButton")->isEnabled());
    EXPECT_EQ(widget<QTableWidget>(window, "parameterTable")->rowCount(), 0);
    run->click();
    EXPECT_FALSE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(1));
    EXPECT_FALSE(widget<QAction>(window, "openSetupAction")->isEnabled());
    EXPECT_FALSE(widget<QAction>(window, "saveSetupAction")->isEnabled());
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(widget<QAction>(window, "openSetupAction")->isEnabled());
    EXPECT_TRUE(widget<QAction>(window, "saveSetupAction")->isEnabled());
    EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Optimizer iterations: Not run"));
}

TEST(ModelControls, WorkerUsesSnapshotBeforeLoadingExperimentInputs)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='false' WHERE criterion='save_model_profiles'");
    const auto defaults = model_controls::load(input.database);
    auto rows = defaults.rows;
    control(rows, "max_iterations").value = "19";
    control(rows, "width resolution (X)").value = "7";
    control(rows, "save_model_profiles").value.reset();
    tsensor_workflow::control_values values;
    for (const auto& row : rows)
        values.emplace_back(row.name.toStdString(), row.value ? std::optional<std::string>(row.value->toStdString()) : std::nullopt);
    values.emplace_back("run_solver", "false");
    tsensor_workflow::background_runner runner;
    runner.start(std::filesystem::path(input.database.toStdWString()), {}, values);
    runner.wait();
    auto result = runner.take_result();
    if (result.failure) std::rethrow_exception(result.failure);
    ASSERT_NE(result.session, nullptr);
    EXPECT_EQ(result.session->parameters().max_iterations, 19);
    EXPECT_EQ(result.session->parameters().X, 7);
    EXPECT_TRUE(result.session->parameters().save_model_profiles);
    EXPECT_FALSE(result.result->optimizer_ran);
    EXPECT_EQ(model_controls::load(input.database), defaults);
}

TEST(Gui, MemoryControlsSurviveReopenAndRunWithoutChangingDefaults)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('run_solver','true')");
    const auto defaults = model_controls::load(input.database);
    MainWindow window; input.choose(window); window.show();
    auto editControls = [&](bool cancel) {
        QTimer::singleShot(0, &window, [&] {
            auto* dialog = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
            ASSERT_NE(dialog, nullptr);
            EXPECT_FALSE(dialog->isWindow());
            EXPECT_EQ(QApplication::activeModalWidget(), nullptr);
            auto* use = widget<QPushButton>(*dialog, "useControlsButton");
            ASSERT_TRUE(until([&] { return use->isEnabled(); }));

            auto* solver = widget<QCheckBox>(*dialog, "run_solver");
            if (cancel) {
                EXPECT_FALSE(solver->isChecked());
                solver->setChecked(true); dialog->reject();
            } else {
                auto* iterations = widget<QLineEdit>(*dialog, "max_iterations");
                iterations->setText("bad"); use->click();
                EXPECT_TRUE(dialog->isVisible());
                iterations->setText("9"); solver->setChecked(false);
                auto* tabs = widget<QTabWidget>(window, "mainTabs");
                tabs->setCurrentIndex(0);
                EXPECT_FALSE(widget<QPushButton>(window, "runButton")->isEnabled());
                tabs->setCurrentIndex(1);
                EXPECT_EQ(iterations->text(), "9");
                if (const auto capture = qEnvironmentVariable("NAVIER_CONTROLS_TAB_CAPTURE"); !capture.isEmpty()) {
                    QApplication::processEvents(); EXPECT_TRUE(window.grab().save(capture));
                }
                use->click();
            }
        });
        widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(1);
        ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "runButton")->isEnabled(); }));
    };
    editControls(false);
    editControls(true);
    auto* run = widget<QPushButton>(window, "runButton");
    for (int i = 0; i < 2; ++i) {
        run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
        EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Optimizer iterations: Not run"));
        EXPECT_EQ(model_controls::load(input.database), defaults);
    }
    // Updating defaults after a memory-only edit must compare with database values.
    QTimer::singleShot(0, &window, [&] {
        auto* dialog = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
        ASSERT_NE(dialog, nullptr);
        auto* save = widget<QPushButton>(*dialog, "saveControlsButton");
        ASSERT_TRUE(until([&] { return save->isEnabled(); }));
        save->click();
    });
    widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(1);
    ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "runButton")->isEnabled(); }));
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 9);
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
TEST(Gui, ExplicitProfileSaveIgnoresLegacyFlagAndCloseDuringExport)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='false' WHERE criterion='save_model_profiles'");
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(widget<QPushButton>(window, "profilesButton")->isEnabled());
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    widget<QPushButton>(window, "profilesButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_GT(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    EXPECT_EQ(input.execute("SELECT value='false' FROM model_controls WHERE criterion='save_model_profiles'"), 1);
    ASSERT_TRUE(widget<QPushButton>(window, "exportButton")->isEnabled());
    widget<QPushButton>(window, "exportButton")->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("reports/uniform,.txt")));
    EXPECT_GT(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
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
