#include "main_window.h"
#include "experiments_tab.h"
#include "database_table_tab.h"
#include "controls_dialog.h"
#include "setup_file.h"
#include "report_tab.h"
#include "profiles_tab.h"
#include <QActionGroup>
#include <QApplication>
#include <QCloseEvent>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QDir>
#include <QDialog>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QMenuBar>
#include <QMetaObject>
#include <QKeySequence>
#include <QSettings>
#include <QWidgetAction>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QTabWidget>
#include <QTimer>
#include <QVBoxLayout>
#include <sqlite3.h>
#include <stdexcept>

using namespace tsensor_workflow;
namespace {
std::filesystem::path path(const QString& value) { return std::filesystem::path(value.toStdWString()); }
QString text(const std::string& value) { return QString::fromUtf8(value.data(), static_cast<int>(value.size())); }
QString terminationLabel(const std::optional<alglib::ae_int_t>& code) {
    if (!code) return "Not run";
    switch (*code) {
    case -8: return "-8 — Numerical failure (NaN/Inf in the function or Jacobian)";
    case -5: return "-5 — Inappropriate optimizer configuration";
    case -3: return "-3 — Inconsistent constraints";
    case 2: return "2 — Relative step reached the EpsX tolerance";
    case 5: return "5 — Maximum iterations reached";
    case 7: return "7 — Further improvement is impossible under the stopping conditions";
    case 8: return "8 — Stopped by a user termination request";
    default: return QString::number(*code) + " — Unrecognized termination code";
    }
}
QString legacyReportFilename(const parameters_t& parameters) {
    QString filename;
    for (const auto& run : parameters.experiment_runs)
        filename += text(run.name) + ",";
    return filename + ".txt";
}
}

MainWindow::MainWindow()
{
    setWindowTitle("Navier — Transient Model");
    resize(900, 720);
    auto* central = new QWidget(this);
    setCentralWidget(central);
    auto* outer = new QVBoxLayout(central);
    tabs_ = new QTabWidget; tabs_->setObjectName("mainTabs");
    outer->addWidget(tabs_);
    auto* results = new QWidget;
    tabs_->addTab(results, "Run and results");
    controlsPage_ = new QWidget;
    new QVBoxLayout(controlsPage_);
    tabs_->addTab(controlsPage_, "Model controls");
    experimentsPage_ = new ExperimentsTab;
    tabs_->addTab(experimentsPage_, "Experiments");
    reactionsPage_ = new DatabaseTableTab(DatabaseTableTab::Table::reactions);
    tabs_->addTab(reactionsPage_, "Reactions");
    speciesPage_ = new DatabaseTableTab(DatabaseTableTab::Table::species);
    tabs_->addTab(speciesPage_, "Species");
    alglibPage_ = new DatabaseTableTab(DatabaseTableTab::Table::alglib);
    tabs_->addTab(alglibPage_, "Variables");
    rawProfilesPage_ = new DatabaseTableTab(DatabaseTableTab::Table::raw_profile);
    tabs_->addTab(rawProfilesPage_, "Raw profiles");
    reportPage_ = new ReportTab;
    reportPage_->generateRequested = [this] { prepareReport(); };
    tabs_->addTab(reportPage_, "Report");
    profilesPage_ = new ProfilesTab;
    profilesPage_->singleRunRequested = [this](const SingleProfileRunRequest& request) {
        tsensor_workflow::single_profile_selection selection;
        selection.experiment_name = request.experiment.toStdString();
        selection.dataset_name = request.dataset.toStdString();
        for (const auto& variable : request.variables) selection.variables.push_back(variable.toStdString());
        selection.initial_values = request.initialValues;
        selection.run_solver = request.runSolver;
        startRun(std::move(selection));
    };
    tabs_->addTab(profilesPage_, "Profiles");
    auto selectedExperiments = [this] {
        const auto snapshot = modelControls_ ? *modelControls_ : model_controls::load(database_->text().trimmed());
        for (const auto& row : snapshot.rows)
            if (row.name == "experiment_name") return row.value.value_or(QString{}).split(' ', Qt::SkipEmptyParts);
        return QStringList{};
    };
    experimentsPage_->experimentFilter()->selectedNames = selectedExperiments;
    rawProfilesPage_->experimentFilter()->selectedNames = selectedExperiments;
    alglibPage_->experimentFilter()->selectedNames = [this] {
        const auto database = database_->text().trimmed();
        return model_controls::selectedVariables(database,
            modelControls_ ? *modelControls_ : model_controls::load(database));
    };
    reactionsPage_->experimentFilter()->selectedNames = [this, selectedExperiments] {
        return model_controls::experimentReferences(database_->text().trimmed(), selectedExperiments(),
            model_controls::ExperimentReferences::reactions);
    };
    speciesPage_->experimentFilter()->selectedNames = [this, selectedExperiments] {
        return model_controls::experimentReferences(database_->text().trimmed(), selectedExperiments(),
            model_controls::ExperimentReferences::species);
    };
    auto* layout = new QVBoxLayout(results);
    auto* title = new QLabel("Navier transient model");
    auto font = title->font(); font.setPointSize(18); title->setFont(font);
    layout->addWidget(title);
    auto* file = menuBar()->addMenu("&File");
    newDatabase_ = file->addAction("&New database..."); newDatabase_->setObjectName("newDatabaseAction");
    openSetup_ = file->addAction("&Open setup..."); openSetup_->setObjectName("openSetupAction");
    saveSetup_ = file->addAction("&Save setup..."); saveSetup_->setObjectName("saveSetupAction");
    file->addSeparator();
    auto makePath = [&](const QString& label, const char* name, QLineEdit*& edit, QPushButton*& browse) {
        auto* menu = file->addMenu(label);
        auto* panel = new QWidget;
        auto* row = new QHBoxLayout(panel);
        edit = new QLineEdit; edit->setObjectName(name);
        edit->setMinimumWidth(360);
        browse = new QPushButton("Browse...");
        row->addWidget(edit); row->addWidget(browse);
        auto* action = new QWidgetAction(menu);
        action->setDefaultWidget(panel); menu->addAction(action);
    };
    makePath("&Database", "databasePath", database_, browseDatabase_);
    makePath("&Output directory", "outputPath", output_, browseOutput_);
    connect(openSetup_, &QAction::triggered, this, [this] { openSetup(); });
    connect(saveSetup_, &QAction::triggered, this, [this] { saveSetup(); });
    connect(newDatabase_, &QAction::triggered, this, [this] { createDatabase(); });

    auto* edit = menuBar()->addMenu("&Edit");
    auto addEditAction = [this, edit](const QString& label, QKeySequence::StandardKey key, const char* method) {
        auto* action = edit->addAction(label);
        action->setShortcut(QKeySequence(key));
        connect(action, &QAction::triggered, this, [method] {
            if (auto* focus = QApplication::focusWidget())
                QMetaObject::invokeMethod(focus, method, Qt::DirectConnection);
        });
        return action;
    };
    addEditAction("Undo", QKeySequence::Undo, "undo");
    addEditAction("Redo", QKeySequence::Redo, "redo");
    edit->addSeparator();
    addEditAction("Cut", QKeySequence::Cut, "cut");
    addEditAction("Copy", QKeySequence::Copy, "copy");
    addEditAction("Paste", QKeySequence::Paste, "paste");
    edit->addSeparator();
    addEditAction("Select All", QKeySequence::SelectAll, "selectAll");

    auto* view = menuBar()->addMenu("&View");
    auto* navigate = view->addMenu("Go to tab");
    auto* tabActions = new QActionGroup(this);
    tabActions->setExclusive(true);
    for (int index = 0; index < tabs_->count(); ++index) {
        auto* action = navigate->addAction(tabs_->tabText(index));
        action->setCheckable(true);
        action->setChecked(index == tabs_->currentIndex());
        tabActions->addAction(action);
        connect(action, &QAction::triggered, this, [this, index] { tabs_->setCurrentIndex(index); });
    }
    connect(tabs_, &QTabWidget::currentChanged, this, [this, tabActions](int index) {
        if (index >= 0 && index < tabActions->actions().size())
            tabActions->actions().at(index)->setChecked(true);
    });
    view->addSeparator();
    auto* showLog = view->addAction("Show progress log");
    showLog->setCheckable(true);
    showLog->setChecked(true);
    connect(showLog, &QAction::toggled, this, [this](bool visible) { log_->setVisible(visible); });

    auto* help = menuBar()->addMenu("&Help");
    auto* guide = help->addAction("Using Navier...");
    connect(guide, &QAction::triggered, this, [this] {
        QMessageBox::information(this, "Using Navier",
            "Choose a database and output directory in File. Configure model controls and experiment data, then select Run and results and click Run. "
            "Use the Report tab to prepare reports. Results and fitted inputs are written only when you choose an explicit save action.\n\n"
            "Aliases link model values to variables selected for solving. Enter an alias such as #kon1 in a reaction, species, or raw-profile value field, "
            "then select kon1 in the applicable solve-for section and provide its bounds in Variables. You can pair a numeric fallback with the alias, "
            "for example 23.7#kon1. When kon1 is selected, the solver uses its variable value; when it is not selected, the model uses 23.7 as a fixed value. "
            "Saving fitted parameters updates the numeric value and keeps the alias.\n\n"
            "For reaction coefficients and exponents, a leading minus sign negates the alias. For example, -#p1 means negative p1, while -23.7#p1 "
            "means negative p1 when selected and -23.7 when not selected. Alias-only entries such as -#p1 need the alias selected; the numeric part in "
            "-23.7#p1 supplies the fallback when it is not selected.");
    });
    auto* about = help->addAction("About Navier");
    connect(about, &QAction::triggered, this, [this] {
        QMessageBox::about(this, "About Navier",
            "<b>Navier — Transient Model</b><p>A transient 1D Navier–Stokes model and parameter fitting application.</p>");
    });

    auto* note = new QLabel("Running may create solution records in the selected database. Results and fitted inputs are saved only when you choose to save.");
    note->setWordWrap(true); layout->addWidget(note);
    auto* controls = new QHBoxLayout;
    run_ = new QPushButton("Run"); run_->setObjectName("runButton");
    pause_ = new QPushButton("Pause"); pause_->setObjectName("pauseButton");
    cancel_ = new QPushButton("Cancel"); cancel_->setObjectName("cancelButton");
    controls->addWidget(run_); controls->addWidget(pause_); controls->addWidget(cancel_); controls->addStretch();
    layout->addLayout(controls);
    status_ = new QLabel("Choose a database to begin."); status_->setObjectName("runStatus");
    status_->setTextFormat(Qt::PlainText); status_->setWordWrap(true); layout->addWidget(status_);
    activity_ = new QProgressBar; activity_->setTextVisible(false); layout->addWidget(activity_);
    summary_ = new QLabel; summary_->setObjectName("resultSummary"); summary_->setWordWrap(true);
    resultDatabase_ = new QLabel; resultDatabase_->setTextFormat(Qt::PlainText); resultDatabase_->setWordWrap(true);
    layout->addWidget(summary_); layout->addWidget(resultDatabase_);
    values_ = new QTableWidget(0, 4); values_->setObjectName("parameterTable");
    values_->setHorizontalHeaderLabels({"Source", "Parameter", "Initial value", "Value"});
    values_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    values_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    layout->addWidget(values_, 1);
    auto* saves = new QHBoxLayout;
    auto* prepare = new QPushButton("Prepare report..."); prepare->setObjectName("prepareRunReportButton");
    connect(prepare, &QPushButton::clicked, this, [this] { tabs_->setCurrentWidget(reportPage_); prepareReport(); });
    export_ = new QPushButton("Export legacy report"); export_->setObjectName("exportButton");
    profiles_ = new QPushButton("Save profiles"); profiles_->setObjectName("profilesButton");
    inputs_ = new QPushButton("Save fitted inputs…"); inputs_->setObjectName("inputsButton");
    saves->addWidget(prepare); saves->addWidget(export_); saves->addWidget(profiles_); saves->addWidget(inputs_);
    layout->addLayout(saves);
    log_ = new QPlainTextEdit; log_->setReadOnly(true); log_->setMaximumBlockCount(500);
    log_->setObjectName("progressLog"); layout->addWidget(log_, 1);
    connect(browseDatabase_, &QPushButton::clicked, this, [this] {
        auto selected = QFileDialog::getOpenFileName(this, "Choose experiment database", database_->text(), "SQLite databases (*.db *.sqlite *.sqlite3);;All files (*)");
        if (!selected.isEmpty()) database_->setText(selected);
    });
    connect(browseOutput_, &QPushButton::clicked, this, [this] {
        auto selected = QFileDialog::getExistingDirectory(this, "Choose output directory", output_->text());
        if (!selected.isEmpty()) output_->setText(selected);
    });
    connect(database_, &QLineEdit::textChanged, this, [this] {
        if (work_ != Work::idle) return;
        discardControlsEditor(); presetFilename_.clear();
        if (tabs_->currentIndex() == 1) tabs_->setCurrentIndex(0);
        modelControls_.reset(); experimentsPage_->clear();
        reactionsPage_->clear(); speciesPage_->clear(); alglibPage_->clear(); rawProfilesPage_->clear(); clearResult();
        // Each table owns its in-memory rows for the selected database.
        reloadTables();
        profilesPage_->setDatabase(database_->text().trimmed());
    });
    connect(output_, &QLineEdit::textChanged, this, [this] { updateControls(); });
    connect(run_, &QPushButton::clicked, this, [this] { startRun(); });
    connect(pause_, &QPushButton::clicked, this, [this] {
        if (runner_.request_pause()) {
            pause_->setEnabled(false);
            setStatus("Pause requested. Finishing the current optimizer iteration…");
        }
    });
    connect(tabs_, &QTabWidget::currentChanged, this, [this](int index) {
        if (tabs_->widget(index) == profilesPage_ && !controlsEditor_ && work_ == Work::idle && !closing_) {
            const QFileInfo input(database_->text().trimmed());
            if (input.isFile()) {
                try {
                    loadControls();
                    profilesPage_->setSingleRunContext(input.absoluteFilePath(), *modelControls_);
                } catch (...) { reportFailure(std::current_exception()); }
            }
        }
        if (index != 1 || controlsEditor_ || work_ != Work::idle || closing_) return;
        const QFileInfo input(database_->text().trimmed());
        if (!input.isFile()) {
            tabs_->setCurrentIndex(0);
            setStatus("Choose an existing database file.");
            return;
        }

        auto* editor = new ControlsDialog(input.absoluteFilePath(), controlsPage_, modelControls_);
        controlsEditor_ = editor;
        editor->configurePresetButton(presetFilename_.isEmpty() ? "Save preset..." : "Update preset", [this] { savePreset(false); });
        editor->setWindowFlags(Qt::Widget);
        editor->setObjectName("modelControlsEditor");
        controlsPage_->layout()->addWidget(editor);
        connect(editor, &QDialog::finished, this, [this, editor](int result) {
            if (result == QDialog::Accepted) {
                modelControls_ = editor->values();
                profilesPage_->setSingleRunContext(database_->text().trimmed(), *modelControls_);
                applyTableFilters();
                clearResult();
            }
            discardControlsEditor();
            tabs_->setCurrentIndex(0);
            updateControls();
        });
        editor->show();
    });
    connect(cancel_, &QPushButton::clicked, this, [this] {
        runner_.request_cancel(); cancelling_ = true; updateControls(); setStatus("Cancellation requested. Waiting for the calculation to stop…");
    });
    connect(export_, &QPushButton::clicked, this, [this] { save(operation::export_results); });
    connect(profiles_, &QPushButton::clicked, this, [this] { save(operation::save_model_profiles); });
    connect(inputs_, &QPushButton::clicked, this, [this] {
        if (QMessageBox::question(this, "Save fitted inputs", "Replace fitted initial inputs in\n" + activeDatabase_ + "?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
            save(operation::save_fitted_parameters);
    });
    experimentsPage_->saved = [this] { applyTableFilters(); clearResult(); setStatus("Experiment saved. Run again to calculate results."); };
    auto referenceSaved = [this] { clearResult(); setStatus("Reference data saved. Run again to calculate results."); };
    reactionsPage_->saved = referenceSaved;
    speciesPage_->saved = referenceSaved;
    alglibPage_->saved = referenceSaved;
    rawProfilesPage_->saved = [this] { applyTableFilters(); clearResult(); setStatus("Raw profile saved. Run again to calculate results."); };
    timer_ = new QTimer(this);
    connect(timer_, &QTimer::timeout, this, [this] { poll(); });
    timer_->start(50);
    profilesPage_->setDatabase(database_->text().trimmed());
    updateControls();
}

bool MainWindow::openReportFile(const QString& filename) {
    const auto suffix = QFileInfo(filename).suffix().toLower();
    if (suffix != "txt" && suffix != "csv" && suffix != "tsv") return false;
    tabs_->setCurrentWidget(reportPage_);
    return reportPage_->loadFile(filename);
}

void MainWindow::reloadTables() {
    const auto database = database_->text().trimmed();
    experimentsPage_->load(database);
    reactionsPage_->load(database);
    speciesPage_->load(database);
    alglibPage_->load(database);
    rawProfilesPage_->load(database);
}

void MainWindow::applyTableFilters() {
    alglibPage_->experimentFilter()->apply();
    experimentsPage_->experimentFilter()->apply();
    rawProfilesPage_->experimentFilter()->apply();
    reactionsPage_->experimentFilter()->apply();
    speciesPage_->experimentFilter()->apply();
}

void MainWindow::setStatus(const QString& value) {
    status_->setText(value.size() > 240 ? value.left(240) + "... (see log)" : value);
    log_->appendPlainText(value);
}

void MainWindow::discardControlsEditor() {
    if (!controlsEditor_) return;
    controlsPage_->layout()->removeWidget(controlsEditor_);
    controlsEditor_->hide();
    controlsEditor_->deleteLater();
    controlsEditor_ = nullptr;
}

void MainWindow::saveSetup() { savePreset(true); }

void MainWindow::savePreset(bool saveAs) {
    if (work_ != Work::idle || closing_) return;
    try {
        loadControls();
        const auto snapshot = controlsEditor_ ? controlsEditor_->draft() : *modelControls_;
        const auto filename = saveAs || presetFilename_.isEmpty()
            ? QFileDialog::getSaveFileName(this, "Save preset", presetFilename_, "Navier setup (*.navier.json)")
            : presetFilename_;
        if (filename.isEmpty()) return;
        setup_file::save(filename, {database_->text().trimmed(), output_->text().trimmed(), snapshot.rows, profilesPage_->settings()});
        presetFilename_ = filename;
        if (controlsEditor_) controlsEditor_->configurePresetButton("Update preset", [this] { savePreset(false); });
        setStatus("Preset saved: " + filename);
        if (controlsEditor_) controlsEditor_->showStatus("Preset saved: " + filename);
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::openSetup() {
    if (work_ != Work::idle || closing_) return;
    QSettings settings("Navier", "Navier");
    auto startingDirectory = settings.value("lastPresetDirectory").toString();
    if (!QDir(startingDirectory).exists()) startingDirectory.clear();
    const auto filename = QFileDialog::getOpenFileName(this, "Open setup", startingDirectory,
        "Navier setup (*.navier.json);;JSON files (*.json)");
    if (filename.isEmpty()) return;
    try {
        const auto setup = setup_file::load(filename);
        const auto original = model_controls::load(setup.database);
        if (original.rows.size() != setup.controls.size()) throw std::runtime_error("Setup controls do not match this database.");
        for (std::size_t i = 0; i < original.rows.size(); ++i) {
            if (original.rows[i].name != setup.controls[i].name ||
                (model_controls::kind(original.rows[i].name) == model_controls::Kind::unknown && original.rows[i] != setup.controls[i]))
                throw std::runtime_error("Setup controls do not match this database.");
        }
        discardControlsEditor();
        clearResult();
        database_->setText(setup.database); output_->setText(setup.outputDirectory);
        profilesPage_->setSettings(setup.profiles);
        modelControls_ = original; modelControls_->rows = setup.controls;
        profilesPage_->setSingleRunContext(setup.database, *modelControls_);
        presetFilename_ = filename;
        tabs_->setCurrentIndex(0);
        experimentsPage_->experimentFilter()->apply();
        rawProfilesPage_->experimentFilter()->apply();
        reactionsPage_->experimentFilter()->apply();
        speciesPage_->experimentFilter()->apply();
        alglibPage_->experimentFilter()->apply();
        settings.setValue("lastPresetDirectory", QFileInfo(filename).absolutePath());
        settings.sync();
        setStatus("Setup opened: " + filename);
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::createDatabase() {
    if (work_ != Work::idle || closing_) return;

    QFileDialog dialog(this, "Create experiment database");
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setNameFilter("SQLite databases (*.db *.sqlite *.sqlite3);;All files (*)");
    dialog.setDefaultSuffix("db");
    dialog.setOption(QFileDialog::DontConfirmOverwrite);
    dialog.selectFile("new-experiment.db");
    if (dialog.exec() != QDialog::Accepted || dialog.selectedFiles().isEmpty()) return;

    const auto filename = dialog.selectedFiles().front();
    QFile schemaFile(":/schema/blank_database.sql");
    if (!schemaFile.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(this, "Cannot create database", "The built-in database schema is unavailable.");
        return;
    }
    const auto schema = schemaFile.readAll();
    if (schemaFile.error() != QFileDevice::NoError) {
        QMessageBox::critical(this, "Cannot create database", "The built-in database schema could not be read.");
        return;
    }

    QFile databaseFile(filename);
    if (!databaseFile.open(QIODevice::WriteOnly | QIODevice::NewOnly)) {
        QMessageBox::warning(this, "Cannot create database",
            "The destination could not be created. Choose a new path; existing files are never replaced.");
        return;
    }
    databaseFile.close();

    sqlite3* database = nullptr;
    const auto encodedFilename = filename.toUtf8();
    int rc = sqlite3_open_v2(encodedFilename.constData(), &database, SQLITE_OPEN_READWRITE, nullptr);
    char* sqlError = nullptr;
    if (rc == SQLITE_OK)
        rc = sqlite3_exec(database, schema.constData(), nullptr, nullptr, &sqlError);

    QString error;
    if (rc != SQLITE_OK) {
        error = sqlError ? QString::fromUtf8(sqlError)
                         : QString::fromUtf8(database ? sqlite3_errmsg(database) : sqlite3_errstr(rc));
    }
    if (sqlError) sqlite3_free(sqlError);
    if (database) {
        const int closeResult = sqlite3_close_v2(database);
        if (rc == SQLITE_OK && closeResult != SQLITE_OK) {
            rc = closeResult;
            error = "The database was created but could not be closed cleanly.";
        }
    }
    if (rc != SQLITE_OK) {
        QFile::remove(filename);
        QMessageBox::critical(this, "Cannot create database", error);
        return;
    }

    discardControlsEditor();
    clearResult();
    database_->setText(QDir::toNativeSeparators(QFileInfo(filename).absoluteFilePath()));
    tabs_->setCurrentIndex(0);
    setStatus("Empty experiment database created: " + filename);
    QMessageBox::information(this, "Database created",
        "The empty database is ready. Add experiment controls and input data before running the model.");
}

void MainWindow::clearResult()
{
    reportPage_->clearGenerated();
    session_.reset(); activeDatabase_.clear(); values_->setRowCount(0); summary_->clear(); resultDatabase_->clear();
    originalResidualSum_.reset();
    latestResidualSum_.reset();
    status_->setText("Choose a database, then run the calculation.");
    updateControls();
}

void MainWindow::updateControls()
{
    const bool idle = work_ == Work::idle && !closing_;
    newDatabase_->setEnabled(idle);
    openSetup_->setEnabled(idle);
    saveSetup_->setEnabled(idle && !database_->text().trimmed().isEmpty() && !output_->text().trimmed().isEmpty());
    database_->setEnabled(idle); browseDatabase_->setEnabled(idle);
    output_->setEnabled(idle); browseOutput_->setEnabled(idle);
    run_->setEnabled(idle && !database_->text().trimmed().isEmpty());
    profilesPage_->setSingleRunEnabled(idle && !database_->text().trimmed().isEmpty());
    tabs_->setTabEnabled(1, idle && !database_->text().trimmed().isEmpty());
    controlsPage_->setEnabled(work_ == Work::idle && !closing_);
    tabs_->setTabEnabled(2, work_ == Work::idle && !closing_);
    experimentsPage_->setEnabled(work_ == Work::idle && !closing_);
    tabs_->setTabEnabled(3, work_ == Work::idle && !closing_);
    tabs_->setTabEnabled(4, work_ == Work::idle && !closing_);
    tabs_->setTabEnabled(5, work_ == Work::idle && !closing_);
    tabs_->setTabEnabled(6, work_ == Work::idle && !closing_);
    rawProfilesPage_->setEnabled(work_ == Work::idle && !closing_);
    reactionsPage_->setEnabled(work_ == Work::idle && !closing_);
    speciesPage_->setEnabled(work_ == Work::idle && !closing_);
    alglibPage_->setEnabled(work_ == Work::idle && !closing_);
    reactionsPage_->setEditingEnabled(idle);
    speciesPage_->setEditingEnabled(idle);
    alglibPage_->setEditingEnabled(idle);
    rawProfilesPage_->setEditingEnabled(idle);
    cancel_->setEnabled(work_ == Work::solve && !closing_ && !cancelling_);
    pause_->setEnabled(work_ == Work::solve && runner_.status() == background_state::running && !closing_);
    export_->setEnabled(idle && session_ && !output_->text().trimmed().isEmpty());
    reportPage_->setRunAvailable(idle && session_);
    findChild<QPushButton*>("prepareRunReportButton")->setEnabled(idle && session_);
    profiles_->setEnabled(idle && session_);
    profiles_->setToolTip("Save result profiles to the run's database.");
    inputs_->setEnabled(idle && session_);
    activity_->setRange(0, idle ? 1 : 0); activity_->setValue(0);
    activity_->setVisible(work_ != Work::idle);
}

void MainWindow::loadControls() {
    if (!modelControls_) modelControls_ = model_controls::load(database_->text().trimmed());
}

void MainWindow::startRun(std::optional<tsensor_workflow::single_profile_selection> selection)
{
    if (work_ != Work::idle || closing_) return;
    const QFileInfo input(database_->text().trimmed());
    if (!input.isFile()) { setStatus("Choose an existing database file."); return; }
    try {
        loadControls();
        if (controlsEditor_) {
            const auto draft = controlsEditor_->draft(false);
            if (draft.rows != modelControls_->rows) {
                if (QMessageBox::question(this, "Accept model controls",
                    "Accept the model control changes and run?", QMessageBox::Yes | QMessageBox::Cancel,
                    QMessageBox::Cancel) != QMessageBox::Yes) return;
                modelControls_ = controlsEditor_->draft();
                experimentsPage_->experimentFilter()->apply();
                rawProfilesPage_->experimentFilter()->apply();
                reactionsPage_->experimentFilter()->apply();
                speciesPage_->experimentFilter()->apply();
                alglibPage_->experimentFilter()->apply();
            }
        }
    } catch (...) { reportFailure(std::current_exception()); return; }
    clearResult(); log_->clear(); activeDatabase_ = input.absoluteFilePath();
    try {
        loadControls();
        control_values controls;
        bool hasRunSolverControl = false;
        for (const auto& row : modelControls_->rows) {
            // Desktop profile saving is an explicit button action, never a model control.
            if (row.name != "save_model_profiles" && row.name != "save_normalized_profiles") {
                auto value = row.value ? std::optional<std::string>(row.value->toStdString()) : std::nullopt;
                if (selection && row.name == "experiment_name") value = selection->experiment_name;
                if (selection && row.name == "run_solver") {
                    value = selection->run_solver ? "true" : "false";
                    hasRunSolverControl = true;
                }
                controls.emplace_back(row.name.toStdString(), std::move(value));
            }
        }
        if (selection && !hasRunSolverControl)
            controls.emplace_back("run_solver", selection->run_solver ? "true" : "false");
        runner_.start(path(activeDatabase_), {}, std::move(controls), true, std::move(selection));
        profilesPage_->beginLiveRun();
        cancelling_ = false; work_ = Work::solve; setStatus("Running…"); updateControls();
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::reportFailure(std::exception_ptr failure)
{
    try { std::rethrow_exception(failure); }
    catch (const workflow_error& error) { setStatus(text(operation_name(error.action)) + ": " + text(error.what())); }
    catch (const std::exception& error) { setStatus("Operation failed: " + text(error.what())); }
    catch (...) { setStatus("Operation failed with an unknown error."); }
    if (controlsEditor_) controlsEditor_->showStatus(status_->text());
}

void MainWindow::save(operation requested)
{
    if (work_ != Work::idle || !session_ || closing_) return;
    std::filesystem::path outputFile;
    if (requested == operation::export_results) {
        const auto outputDirectory = output_->text().trimmed();
        if (outputDirectory.isEmpty()) return;
        const auto suggestedPath = QDir(outputDirectory).filePath(legacyReportFilename(session_->parameters()));
        auto selected = QFileDialog::getSaveFileName(this, "Export legacy report", suggestedPath,
            "Text files (*.txt);;All files (*)");
        if (selected.isEmpty()) return;
        if (QFileInfo(selected).suffix().isEmpty()) selected += ".txt";
        outputFile = path(selected);
    }
    try {
        action_ = std::async(std::launch::async, [session = session_.get(), requested, outputFile] {
            ActionResult result;
            try {
                if (requested == operation::export_results) {
                    auto written = session->export_results_as(outputFile);
                    result.reportFile = QString::fromStdWString(written.wstring());
                    result.message = "Exported report: " + QString::fromStdWString(written.wstring());
                } else if (requested == operation::save_model_profiles) {
                    session->save_model_profiles(); result.reloadTables = true; result.message = "Profiles saved. Open Profiles to plot the saved model profiles.";
                } else {
                    session->save_fitted_parameters(); result.reloadTables = true; result.message = "Fitted inputs saved.";
                }
            } catch (...) { result.error = std::current_exception(); }
            return result;
        });
        work_ = Work::save; setStatus("Saving…"); updateControls();
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::prepareReport()
{
    if (work_ != Work::idle || !session_ || closing_) return;
    try {
        action_ = std::async(std::launch::async, [session = session_.get()] {
            ActionResult result;
            try {
                result.reportText = text(session->generate_report());
                result.message = "Report prepared from run results. Choose data and formatting in Report, then export.";
            } catch (...) { result.error = std::current_exception(); }
            return result;
        });
        work_ = Work::save; setStatus("Preparing report from current results…"); updateControls();
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::showParameters(const std::vector<parameter_value>& parameters, bool completed)
{
    values_->setRowCount(static_cast<int>(parameters.size()));
    for (int i = 0; i < values_->rowCount(); ++i) {
        const auto& parameter = parameters[i];
        values_->setItem(i, 0, new QTableWidgetItem(text(parameter.source)));
        values_->setItem(i, 1, new QTableWidgetItem(text(parameter.name)));
        values_->setItem(i, 2, new QTableWidgetItem(QString::number(parameter.initial_value, 'g', 15)));
        values_->setItem(i, 3, new QTableWidgetItem(completed ? QString::number(parameter.value, 'g', 15) : QString{}));
    }
}

void MainWindow::poll()
{
    if (work_ == Work::solve) {
        // Observe completion before draining so terminal progress is not lost.
        const auto state = runner_.status();
        for (const auto& event : runner_.drain_events()) {
            if (event.kind == event_kind::parameters_initialized)
                showParameters(event.parameters, false);
            else if (event.kind == event_kind::evaluation || event.kind == event_kind::optimizer_progress) {
                if (event.kind == event_kind::optimizer_progress) {
                    showParameters(event.parameters, true);
                    if (!originalResidualSum_) originalResidualSum_ = event.sum_squared_residuals;
                    latestResidualSum_ = event.sum_squared_residuals;
                }
                QString summary = "Completed model evaluations: " + QString::number(event.evaluations.value_or(0));
                if (originalResidualSum_)
                    summary += " · Original sum of squared residuals: " + QString::number(*originalResidualSum_, 'g', 15);
                if (latestResidualSum_)
                    summary += " · Current sum of squared residuals: " + QString::number(*latestResidualSum_, 'g', 15);
                summary_->setText(summary);
            }
            else if (event.kind != event_kind::message || event.detail_level >= 3)
                log_->appendPlainText(text(event.message));
        }
        if (auto liveProfiles = runner_.take_live_profiles())
            profilesPage_->setLiveProfiles(std::move(*liveProfiles));
        if (state != background_state::running) {
            auto outcome = runner_.take_result(); work_ = Work::idle;
            if (state == background_state::cancelled) { summary_->clear(); setStatus("Run cancelled. No results were saved."); }
            else if (outcome.failure) { summary_->clear(); reportFailure(outcome.failure); }
            else {
                session_ = std::move(outcome.session);
                const auto& result = *outcome.result;
                summary_->setText(QString("Model evaluations: %1 · Optimizer iterations: %2 · Termination: %3")
                    .arg(result.residual_evaluations)
                    .arg(result.optimizer_iterations ? QString::number(*result.optimizer_iterations) : "Not run")
                    .arg(terminationLabel(result.termination_type)));
                resultDatabase_->setText("Results from: " + activeDatabase_);
                if (result.sum_squared_residuals)
                    summary_->setText(summary_->text() + " · Final sum of squared residuals: " +
                                      QString::number(*result.sum_squared_residuals, 'g', 15));
                originalResidualSum_ = result.original_sum_squared_residuals;
                if (originalResidualSum_)
                    summary_->setText(summary_->text() + " · Original sum of squared residuals: " +
                                      QString::number(*originalResidualSum_, 'g', 15));
                showParameters(result.parameters, true);
                const bool residualSumIncreased = result.original_sum_squared_residuals &&
                    result.sum_squared_residuals &&
                    *result.sum_squared_residuals > *result.original_sum_squared_residuals;
                if (result.paused)
                    setStatus("Paused at an optimizer iteration. Review or export the partial results, save profiles, or save fitted inputs for a later run.");
                else if (residualSumIncreased)
                    setStatus("Run finished with a higher residual sum than it started with.");
                else
                    setStatus("Run finished. Review the termination reason and results before saving.");
                if (residualSumIncreased && !closing_) {
                    const auto original = QString::number(*result.original_sum_squared_residuals, 'g', 15);
                    const auto final = QString::number(*result.sum_squared_residuals, 'g', 15);
                    QMessageBox::warning(this, "Residual sum increased",
                        "The final sum of squared residuals is greater than the original value.\n\n"
                        "Original: " + original + "\nFinal: " + final);
                }
            }
            updateControls();
        }
    } else if (work_ == Work::save && action_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        auto result = action_.get(); work_ = Work::idle;
        if (result.error) {
            // Keep errors and results visible for retry after a deferred close.
            closing_ = false;
            reportFailure(result.error);
        } else {
            if (result.reloadTables) {
                reloadTables();
                profilesPage_->setDatabase(activeDatabase_);
            }
            if (!result.reportFile.isEmpty()) {
                reportPage_->loadFile(result.reportFile);
                if (!closing_) tabs_->setCurrentWidget(reportPage_);
            }
            if (!result.reportText.isEmpty()) {
                reportPage_->loadGenerated(result.reportText, activeDatabase_, QDir(output_->text()).filePath("report"));
                if (!closing_) tabs_->setCurrentWidget(reportPage_);
            }
            setStatus(result.message);
        }
        updateControls();
    }
    if (closing_ && work_ == Work::idle) close();
}

void MainWindow::closeEvent(QCloseEvent* event)
{
    if (work_ == Work::idle) { timer_->stop(); event->accept(); return; }
    closing_ = true;
    if (work_ == Work::solve) runner_.request_cancel();
    setStatus(work_ == Work::solve ? "Stopping the calculation before closing…" : "Finishing the save before closing…");
    updateControls(); event->ignore();
}
