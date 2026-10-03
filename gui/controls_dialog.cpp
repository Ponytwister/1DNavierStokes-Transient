#include "controls_dialog.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QTimer>
#include <QVBoxLayout>

ControlsDialog::ControlsDialog(const QString& database, QWidget* parent, std::optional<model_controls::Snapshot> current)
    : QDialog(parent), database_(database)
{
    setWindowTitle("Model controls"); setObjectName("modelControlsDialog"); resize(780, 600);
    auto* layout = new QVBoxLayout(this);
    auto* description = new QLabel("Edit controls for the next run. Update default saves to the database. Use values applies only in memory; Cancel discards edits.\nDatabase: " + database);
    description->setTextFormat(Qt::PlainText); description->setWordWrap(true); layout->addWidget(description);
    table_ = new QTableWidget(0, 3); table_->setObjectName("controlsTable");
    table_->setHorizontalHeaderLabels({"Control", "Value", "NULL (skip)"});
    table_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    table_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    layout->addWidget(table_);
    status_ = new QLabel("Loading controls..."); status_->setObjectName("controlsStatus");
    status_->setTextFormat(Qt::PlainText); status_->setWordWrap(true); layout->addWidget(status_);
    auto* buttons = new QDialogButtonBox;
    save_ = buttons->addButton("Update default", QDialogButtonBox::AcceptRole); save_->setObjectName("saveControlsButton");
    use_ = buttons->addButton("Use values", QDialogButtonBox::AcceptRole); use_->setObjectName("useControlsButton");
    use_->setDefault(true);
    cancel_ = buttons->addButton(QDialogButtonBox::Cancel); cancel_->setObjectName("cancelControlsButton");
    layout->addWidget(buttons);
    connect(save_, &QPushButton::clicked, this, [this] { save(true); });
    connect(use_, &QPushButton::clicked, this, [this] { save(false); });
    connect(cancel_, &QPushButton::clicked, this, &ControlsDialog::reject);
    auto* timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, [this] { poll(); }); timer->start(50);
    if (current) current_ = *current;
    busy(true);
    try { pending_ = std::async(std::launch::async, [database] { return model_controls::load(database); }); }
    catch (const std::exception& error) { busy(false); save_->setEnabled(false); use_->setEnabled(false); status_->setText(QString::fromUtf8(error.what())); }
}
void ControlsDialog::busy(bool value) {
    use_->setEnabled(!value); table_->setEnabled(!value); save_->setEnabled(!value); cancel_->setEnabled(!value);
}
void ControlsDialog::reject() {
    // A pending transaction must finish before its dialog and result are discarded.
    if (!pending_.valid()) QDialog::reject();
}
void ControlsDialog::populate() {
    table_->setRowCount(static_cast<int>(current_.rows.size()));
    for (int i = 0; i < table_->rowCount(); ++i) {
        const auto& row = current_.rows[i];
        const auto type = model_controls::kind(row.name);
        const bool supported = type != model_controls::Kind::unknown;
        auto* name = new QTableWidgetItem(row.name);
        name->setFlags(name->flags() & ~Qt::ItemIsEditable);
        name->setToolTip(model_controls::help(row.name)); table_->setItem(i, 0, name);
        QWidget* editor;
        if (type == model_controls::Kind::boolean || type == model_controls::Kind::scatter) {
            auto* combo = new QComboBox;
            combo->addItems(type == model_controls::Kind::boolean ? QStringList{"true", "false"} : QStringList{"none", "NS_ND"});
            if (row.value && combo->findText(*row.value) < 0) combo->addItem(*row.value);
            combo->setCurrentText(row.value.value_or(combo->itemText(0))); editor = combo;
        } else {
            auto* edit = new QLineEdit(row.value.value_or(QString{})); editor = edit;
        }
        editor->setToolTip(model_controls::help(row.name));
        editor->setEnabled(supported && row.value.has_value()); table_->setCellWidget(i, 1, editor);
        auto* null = new QCheckBox; null->setChecked(!row.value); null->setEnabled(supported);
        null->setToolTip("SQL NULL: the existing model loader skips this row. It does not supply a required value.");
        connect(null, &QCheckBox::toggled, editor, [editor, supported](bool checked) { editor->setEnabled(supported && !checked); });
        table_->setCellWidget(i, 2, null);
    }
}
void ControlsDialog::save(bool persist) {
    if (pending_.valid()) return;
    auto edited = current_.rows;
    for (int i = 0; i < table_->rowCount(); ++i) {
        if (model_controls::kind(edited[i].name) == model_controls::Kind::unknown) continue;
        if (qobject_cast<QCheckBox*>(table_->cellWidget(i, 2))->isChecked()) edited[i].value.reset();
        else if (auto* combo = qobject_cast<QComboBox*>(table_->cellWidget(i, 1))) edited[i].value = combo->currentText();
        else edited[i].value = qobject_cast<QLineEdit*>(table_->cellWidget(i, 1))->text();
    }
    try {
        for (std::size_t i = 0; i < edited.size(); ++i)
            if (edited[i] != current_.rows[i]) model_controls::validate(edited[i]);
        if (!persist) { current_.rows = std::move(edited); accept(); return; }
        pending_ = std::async(std::launch::async, [database = database_, original = original_, edited] {
            model_controls::save(database, original, edited);
            auto result = original; result.rows = edited; return result;
        });
        saving_ = true; busy(true); status_->setText("Saving model controls...");
    } catch (const std::exception& error) { status_->setText(QString::fromUtf8(error.what())); }
}
void ControlsDialog::poll() {
    if (!pending_.valid() || pending_.wait_for(std::chrono::seconds(0)) != std::future_status::ready) return;
    try {
        auto result = pending_.get();
        busy(false);
        if (saving_) { current_ = std::move(result); accept(); return; }
        original_ = result; if (current_.nameColumn.isEmpty()) current_ = std::move(result); populate();
        save_->setEnabled(!original_.rows.empty()); use_->setEnabled(!original_.rows.empty());
        status_->setText("Hover over a control for help. Unrecognized controls are read-only.");
    } catch (const std::exception& error) {
        busy(false); save_->setEnabled(saving_); use_->setEnabled(saving_); saving_ = false;
        status_->setText(QString::fromUtf8(error.what()));
    }
}
