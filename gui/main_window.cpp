#include "main_window.h"
#include <QCloseEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QTableWidget>
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
    auto* layout = new QVBoxLayout(central);
    auto* title = new QLabel("Navier transient model");
    auto font = title->font(); font.setPointSize(18); title->setFont(font);
    layout->addWidget(title);
    auto* form = new QFormLayout;
    auto makePath = [&](const QString& label, const char* name, QLineEdit*& edit, QPushButton*& browse) {
        auto* row = new QHBoxLayout;
        edit = new QLineEdit; edit->setObjectName(name);
        browse = new QPushButton("Browse…");
        row->addWidget(edit); row->addWidget(browse); form->addRow(label, row);
    };
    makePath("Database", "databasePath", database_, browseDatabase_);
    makePath("Output directory", "outputPath", output_, browseOutput_);
    layout->addLayout(form);
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
    values_ = new QTableWidget(0, 3); values_->setObjectName("parameterTable");
    values_->setHorizontalHeaderLabels({"Source", "Parameter", "Value"});
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
    connect(database_, &QLineEdit::textChanged, this, [this] { if (work_ == Work::idle) clearResult(); });
    connect(output_, &QLineEdit::textChanged, this, [this] { updateControls(); });
    connect(run_, &QPushButton::clicked, this, [this] { startRun(); });
    connect(cancel_, &QPushButton::clicked, this, [this] {
        runner_.request_cancel(); cancel_->setEnabled(false); setStatus("Cancellation requested. Waiting for the calculation to stop…");
    });
    connect(export_, &QPushButton::clicked, this, [this] { save(operation::export_results); });
    connect(profiles_, &QPushButton::clicked, this, [this] { save(operation::save_model_profiles); });
    connect(inputs_, &QPushButton::clicked, this, [this] {
        if (QMessageBox::question(this, "Save fitted inputs", "Replace fitted initial inputs in\n" + activeDatabase_ + "?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
            save(operation::save_fitted_parameters);
    });
    timer_ = new QTimer(this);
    connect(timer_, &QTimer::timeout, this, [this] { poll(); });
    timer_->start(50);
    updateControls();
}

void MainWindow::setStatus(const QString& value) {
    status_->setText(value.size() > 240 ? value.left(240) + "... (see log)" : value);
    log_->appendPlainText(value);
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
    database_->setEnabled(idle); browseDatabase_->setEnabled(idle);
    output_->setEnabled(idle); browseOutput_->setEnabled(idle);
    run_->setEnabled(idle && !database_->text().trimmed().isEmpty());
    cancel_->setEnabled(work_ == Work::solve && !closing_);
    export_->setEnabled(idle && session_ && !output_->text().trimmed().isEmpty());
    profiles_->setEnabled(idle && session_ && session_->parameters().save_model_profiles);
    profiles_->setToolTip(idle && session_ && !session_->parameters().save_model_profiles
        ? "Profile saving is disabled in this database's settings." : "Save result profiles to the run's database.");
    inputs_->setEnabled(idle && session_);
    activity_->setRange(0, idle ? 1 : 0); activity_->setValue(0);
    activity_->setVisible(work_ != Work::idle);
}

void MainWindow::startRun()
{
    if (work_ != Work::idle || closing_) return;
    const QFileInfo input(database_->text().trimmed());
    if (!input.isFile()) { setStatus("Choose an existing database file."); return; }
    clearResult(); log_->clear(); activeDatabase_ = input.absoluteFilePath();
    try {
        runner_.start(path(activeDatabase_));
        work_ = Work::solve; setStatus("Running…"); updateControls();
    } catch (...) { reportFailure(std::current_exception()); }
}

void MainWindow::reportFailure(std::exception_ptr failure)
{
    try { std::rethrow_exception(failure); }
    catch (const workflow_error& error) { setStatus(text(operation_name(error.action)) + ": " + text(error.what())); }
    catch (const std::exception& error) { setStatus("Operation failed: " + text(error.what())); }
    catch (...) { setStatus("Operation failed with an unknown error."); }
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

void MainWindow::poll()
{
    if (work_ == Work::solve) {
        for (const auto& event : runner_.drain_events()) {
            if (event.kind == event_kind::evaluation)
                summary_->setText("Completed model evaluations: " + QString::number(event.evaluations.value_or(0)));
            else if (event.kind != event_kind::message || event.detail_level >= 3)
                log_->appendPlainText(text(event.message));
        }
        const auto state = runner_.status();
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
                values_->setRowCount(static_cast<int>(result.parameters.size()));
                for (int i = 0; i < values_->rowCount(); ++i) {
                    const auto& parameter = result.parameters[i];
                    values_->setItem(i, 0, new QTableWidgetItem(text(parameter.source)));
                    values_->setItem(i, 1, new QTableWidgetItem(text(parameter.name)));
                    values_->setItem(i, 2, new QTableWidgetItem(QString::number(parameter.value, 'g', 15)));
                }
                setStatus("Run finished. Review the termination code and results before saving.");
            }
            updateControls();
        }
    } else if (work_ == Work::save && action_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        auto result = action_.get(); work_ = Work::idle;
        if (result.error) reportFailure(result.error); else setStatus(result.message);
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
