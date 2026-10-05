#include "main_window.h"
#include "experiments_tab.h"
#include "database_table_tab.h"
#include "controls_dialog.h"
#include "setup_file.h"
#include <QCloseEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QMenuBar>
#include <QWidgetAction>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QTableWidget>
#include <QTabWidget>
#include <QTimer>
#include <QVBoxLayout>

using namespace tsensor_workflow;
namespace {
std::filesystem::path path(const QString& value) { return std::filesystem::path(value.toStdWString()); }
QString text(const std::string& value) { return QString::fromUtf8(value.data(), static_cast<int>(value.size())); }
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
    auto selectedExperiments = [this] {
        const auto snapshot = modelControls_ ? *modelControls_ : model_controls::load(database_->text().trimmed());
        for (const auto& row : snapshot.rows)
            if (row.name == "experiment_name") return row.value.value_or(QString{}).split(' ', Qt::SkipEmptyParts);
        return QStringList{};
    };
    experimentsPage_->experimentFilter()->selectedNames = selectedExperiments;
    rawProfilesPage_->experimentFilter()->selectedNames = selectedExperiments;
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
    auto* note = new QLabel("Running may create solution records in the selected database. Results and fitted inputs are saved only when you choose to save.");
    note->setWordWrap(true); layout->addWidget(note);
    auto* controls = new QHBoxLayout;
    run_ = new QPushButton("Run"); run_->setObjectName("runButton");
    cancel_ = new QPushButton("Cancel"); cancel_->setObjectName("cancelButton");
    controls->addWidget(run_); controls->addWidget(cancel_); controls->addStretch();
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
    export_ = new QPushButton("Export report"); export_->setObjectName("exportButton");
    profiles_ = new QPushButton("Save profiles"); profiles_->setObjectName("profilesButton");
    inputs_ = new QPushButton("Save fitted inputs…"); inputs_->setObjectName("inputsButton");
    saves->addWidget(export_); saves->addWidget(profiles_); saves->addWidget(inputs_);
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
        const auto database = database_->text().trimmed();
        if (tabs_->currentIndex() == 2) experimentsPage_->load(database);
        if (tabs_->currentIndex() == 3) reactionsPage_->load(database);
        if (tabs_->currentIndex() == 4) speciesPage_->load(database);
        if (tabs_->currentIndex() == 5) alglibPage_->load(database);
        if (tabs_->currentIndex() == 6) rawProfilesPage_->load(database);
    });
    connect(output_, &QLineEdit::textChanged, this, [this] { updateControls(); });
    connect(run_, &QPushButton::clicked, this, [this] { startRun(); });
    connect(tabs_, &QTabWidget::currentChanged, this, [this](int index) {
        if (work_ == Work::idle && !closing_) {
            const auto database = database_->text().trimmed();
            if (index == 2) experimentsPage_->load(database);
            if (index == 3) reactionsPage_->load(database);
            if (index == 4) speciesPage_->load(database);
            if (index == 5) alglibPage_->load(database);
            if (index == 6) rawProfilesPage_->load(database);
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
    experimentsPage_->saved = [this] { clearResult(); setStatus("Experiment saved. Run again to calculate results."); };
    auto referenceSaved = [this] { clearResult(); setStatus("Reference data saved. Run again to calculate results."); };
    reactionsPage_->saved = referenceSaved;
    speciesPage_->saved = referenceSaved;
    alglibPage_->saved = referenceSaved;
    rawProfilesPage_->saved = [this] { clearResult(); setStatus("Raw profile saved. Run again to calculate results."); };
    timer_ = new QTimer(this);
    connect(timer_, &QTimer::timeout, this, [this] { poll(); });
    timer_->start(50);
    updateControls();
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
        setup_file::save(filename, {database_->text().trimmed(), output_->text().trimmed(), snapshot.rows});
        presetFilename_ = filename;
        if (controlsEditor_) controlsEditor_->configurePresetButton("Update preset", [this] { savePreset(false); });
        setStatus("Preset saved: " + filename);
        if (controlsEditor_) controlsEditor_->showStatus("Preset saved: " + filename);
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::openSetup() {
    if (work_ != Work::idle || closing_) return;
    const auto filename = QFileDialog::getOpenFileName(this, "Open setup", {}, "Navier setup (*.navier.json);;JSON files (*.json)");
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
        modelControls_ = original; modelControls_->rows = setup.controls;
        presetFilename_ = filename;
        tabs_->setCurrentIndex(0);
        experimentsPage_->experimentFilter()->apply();
        rawProfilesPage_->experimentFilter()->apply();
        reactionsPage_->experimentFilter()->apply();
        speciesPage_->experimentFilter()->apply();
        setStatus("Setup opened: " + filename);
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::clearResult()
{
    session_.reset(); activeDatabase_.clear(); values_->setRowCount(0); summary_->clear(); resultDatabase_->clear();
    status_->setText("Choose a database, then run the calculation.");
    updateControls();
}

void MainWindow::updateControls()
{
    const bool idle = work_ == Work::idle && !closing_;
    openSetup_->setEnabled(idle);
    saveSetup_->setEnabled(idle && !database_->text().trimmed().isEmpty() && !output_->text().trimmed().isEmpty());
    database_->setEnabled(idle); browseDatabase_->setEnabled(idle);
    output_->setEnabled(idle); browseOutput_->setEnabled(idle);
    run_->setEnabled(idle && !database_->text().trimmed().isEmpty());
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
    export_->setEnabled(idle && session_ && !output_->text().trimmed().isEmpty());
    profiles_->setEnabled(idle && session_);
    profiles_->setToolTip("Save result profiles to the run's database.");
    inputs_->setEnabled(idle && session_);
    activity_->setRange(0, idle ? 1 : 0); activity_->setValue(0);
    activity_->setVisible(work_ != Work::idle);
}

void MainWindow::loadControls() {
    if (!modelControls_) modelControls_ = model_controls::load(database_->text().trimmed());
}

void MainWindow::startRun()
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
            }
        }
    } catch (...) { reportFailure(std::current_exception()); return; }
    clearResult(); log_->clear(); activeDatabase_ = input.absoluteFilePath();
    try {
        loadControls();
        control_values controls;
        for (const auto& row : modelControls_->rows)
            // Desktop profile saving is an explicit button action, never a model control.
            if (row.name != "save_model_profiles" && row.name != "save_normalized_profiles")
                controls.emplace_back(row.name.toStdString(), row.value ? std::optional<std::string>(row.value->toStdString()) : std::nullopt);
        runner_.start(path(activeDatabase_), {}, std::move(controls));
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
    auto directory = path(output_->text().trimmed());
    if (requested == operation::export_results && directory.empty()) return;
    try {
        action_ = std::async(std::launch::async, [session = session_.get(), requested, directory] {
            ActionResult result;
            try {
                if (requested == operation::export_results) {
                    auto written = session->export_results(directory);
                    result.message = "Exported report: " + QString::fromStdWString(written.wstring());
                } else if (requested == operation::save_model_profiles) {
                    session->save_model_profiles(); result.message = "Profiles saved.";
                } else {
                    session->save_fitted_parameters(); result.message = "Fitted inputs saved.";
                }
            } catch (...) { result.error = std::current_exception(); }
            return result;
        });
        work_ = Work::save; setStatus("Saving…"); updateControls();
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
            else if (event.kind == event_kind::evaluation)
                summary_->setText("Completed model evaluations: " + QString::number(event.evaluations.value_or(0)));
            else if (event.kind != event_kind::message || event.detail_level >= 3)
                log_->appendPlainText(text(event.message));
        }
        if (state != background_state::running) {
            auto outcome = runner_.take_result(); work_ = Work::idle;
            if (state == background_state::cancelled) { summary_->clear(); setStatus("Run cancelled. No results were saved."); }
            else if (outcome.failure) { summary_->clear(); reportFailure(outcome.failure); }
            else {
                session_ = std::move(outcome.session);
                const auto& result = *outcome.result;
                summary_->setText(QString("Model evaluations: %1 · Optimizer iterations: %2 · Termination code: %3")
                    .arg(result.residual_evaluations)
                    .arg(result.optimizer_iterations ? QString::number(*result.optimizer_iterations) : "Not run")
                    .arg(result.termination_type ? QString::number(*result.termination_type) : "Not run"));
                resultDatabase_->setText("Results from: " + activeDatabase_);
                showParameters(result.parameters, true);
                setStatus("Run finished. Review the termination code and results before saving.");
            }
            updateControls();
        }
    } else if (work_ == Work::save && action_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        auto result = action_.get(); work_ = Work::idle;
        if (result.error) {
            // Keep errors and results visible for retry after a deferred close.
            closing_ = false;
            reportFailure(result.error);
        } else setStatus(result.message);
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
