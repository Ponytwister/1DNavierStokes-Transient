#pragma once
#include "model_controls.h"
#include <QDialog>
#include <future>
#include <functional>

class QLabel;
class QPushButton;
class QFormLayout;
class ControlsDialog : public QDialog {
public:
    explicit ControlsDialog(const QString& database, QWidget* parent = nullptr,
                            std::optional<model_controls::Snapshot> current = std::nullopt);
    const model_controls::Snapshot& values() const { return current_; }
    model_controls::Snapshot draft(bool validate = true) const;
    bool ready() const { return !pending_.valid() && !editors_.empty(); }
    void showStatus(const QString& message);
    void configurePresetButton(const QString& label, std::function<void()> savePreset);
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
    struct Loaded { model_controls::Snapshot controls; QStringList experiments; QStringList parameters; };
    std::future<Loaded> pending_;
    QStringList experiments_, parameters_;
    QLabel* status_;
    QWidget* fields_;
    QFormLayout* form_;
    std::vector<QWidget*> editors_;
    QPushButton *save_, *cancel_, *use_;
};
