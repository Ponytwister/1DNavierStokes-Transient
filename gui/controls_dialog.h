#pragma once
#include "model_controls.h"
#include <QDialog>
#include <future>

class QLabel;
class QPushButton;
class QFormLayout;
class ControlsDialog : public QDialog {
public:
    explicit ControlsDialog(const QString& database, QWidget* parent = nullptr,
                            std::optional<model_controls::Snapshot> current = std::nullopt);
    const model_controls::Snapshot& values() const { return current_; }
    bool updatedDefault() const { return saving_; }
    void reject() override;
private:
    void populate();
    void save(bool persist);
    void poll();
    void busy(bool value);
    QString database_;
    model_controls::Snapshot original_, current_;
    bool saving_ = false;
    std::future<model_controls::Snapshot> pending_;
    QLabel* status_;
    QWidget* fields_;
    QFormLayout* form_;
    std::vector<QWidget*> editors_;
    QPushButton *save_, *cancel_, *use_;
};
