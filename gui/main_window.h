#pragma once
#include <QMainWindow>
#include <background_runner.h>
#include <future>
#include "model_controls.h"

class ExperimentsTab;
class DatabaseTableTab;
class QAction;
class QLineEdit;
class QPushButton;
class QLabel;
class QPlainTextEdit;
class QTableWidget;
class QProgressBar;
class QTimer;
class QTabWidget;

class MainWindow : public QMainWindow {
public:
    MainWindow();
protected:
    void closeEvent(QCloseEvent* event) override;
private:
    enum class Work { idle, solve, save };
    struct ActionResult { QString message; std::exception_ptr error; };
    void startRun();
    void saveSetup();
    void openSetup();
    void poll();
    void updateControls();
    void clearResult();
    void save(tsensor_workflow::operation action);
    void reportFailure(std::exception_ptr error);
    void setStatus(const QString& text);
    Work work_ = Work::idle;
    bool closing_ = false;
    bool cancelling_ = false;
    bool editingControls_ = false;
    QString activeDatabase_;
    std::optional<model_controls::Snapshot> modelControls_;
    void loadControls();
    tsensor_workflow::background_runner runner_;
    std::unique_ptr<tsensor_workflow::run_session> session_;
    // Destroy/join pending work before destroying the session it borrows.
    std::future<ActionResult> action_;
    QLineEdit *database_, *output_;
    QAction *openSetup_, *saveSetup_;
    QPushButton *browseDatabase_, *browseOutput_, *run_, *cancel_, *export_, *profiles_, *inputs_;
    QLabel *status_, *summary_, *resultDatabase_;
    QPlainTextEdit* log_;
    QTableWidget* values_;
    QProgressBar* activity_;
    QTimer* timer_;
    QTabWidget* tabs_;
    QWidget* controlsPage_;
    ExperimentsTab* experimentsPage_;
    DatabaseTableTab *reactionsPage_, *speciesPage_, *alglibPage_;
};
