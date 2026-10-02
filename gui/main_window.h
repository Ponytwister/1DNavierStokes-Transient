#pragma once
#include <QMainWindow>
#include <background_runner.h>
#include <future>

class QLineEdit;
class QPushButton;
class QLabel;
class QPlainTextEdit;
class QTableWidget;
class QProgressBar;
class QTimer;

class MainWindow : public QMainWindow {
public:
    MainWindow();
protected:
    void closeEvent(QCloseEvent* event) override;
private:
    enum class Work { idle, solve, save };
    struct ActionResult { QString message; std::exception_ptr error; };
    void startRun();
    void poll();
    void updateControls();
    void clearResult();
    void save(tsensor_workflow::operation action);
    void reportFailure(std::exception_ptr error);
    void setStatus(const QString& text);
    Work work_ = Work::idle;
    bool closing_ = false;
    bool cancelling_ = false;
    QString activeDatabase_;
    tsensor_workflow::background_runner runner_;
    std::unique_ptr<tsensor_workflow::run_session> session_;
    // Destroy/join pending work before destroying the session it borrows.
    std::future<ActionResult> action_;
    QLineEdit *database_, *output_;
    QPushButton *browseDatabase_, *browseOutput_, *run_, *cancel_, *export_, *profiles_, *inputs_;
    QLabel *status_, *summary_, *resultDatabase_;
    QPlainTextEdit* log_;
    QTableWidget* values_;
    QProgressBar* activity_;
    QTimer* timer_;
};
