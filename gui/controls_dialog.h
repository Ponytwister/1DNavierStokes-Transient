#pragma once
#include "model_controls.h"
#include <QDialog>
#include <future>

class QLabel;
class QPushButton;
class QTableWidget;
class ControlsDialog : public QDialog {
public:
    explicit ControlsDialog(const QString& database, QWidget* parent = nullptr);
    void reject() override;
private:
    void populate();
    void save();
    void poll();
    void busy(bool value);
    QString database_;
    model_controls::Snapshot original_;
    bool saving_ = false;
    std::future<model_controls::Snapshot> pending_;
    QLabel* status_;
    QTableWidget* table_;
    QPushButton *save_, *cancel_;
};
