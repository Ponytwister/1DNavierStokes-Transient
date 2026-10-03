#include "controls_dialog.h"
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QToolButton>
#include <QMenu>
#include <QListWidget>
#include <QWidgetAction>
#include <QRegularExpression>
#include <stdexcept>

#include <QTimer>
#include <QVBoxLayout>

namespace {
class ChecklistPicker : public QToolButton {
public:
    ChecklistPicker(const QStringList& available, const QString& selected, const QString& field,
                    const QString& choicesName, const QString& emptyText)
        : field_(field), emptyText_(emptyText) {
        setPopupMode(QToolButton::InstantPopup);
        setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
        setArrowType(Qt::DownArrow);
        setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        initial_ = selected.split(' ', Qt::SkipEmptyParts);
        auto* menu = new QMenu(this);
        auto* action = new QWidgetAction(menu);
        list_ = new QListWidget;
        list_->setObjectName(choicesName);
        list_->setSelectionMode(QAbstractItemView::NoSelection);
        list_->setMinimumSize(360, 220);
        auto names = available;
        for (const auto& name : initial_) if (!names.contains(name)) names.push_back(name);
        for (const auto& name : names) {
            const bool supported = available.contains(name) && !name.contains(QRegularExpression("\\s"));
            auto* item = new QListWidgetItem(name, list_);
            item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
            item->setData(Qt::UserRole, supported);
            item->setCheckState(initial_.contains(name) ? Qt::Checked : Qt::Unchecked);
            if (!supported) item->setToolTip(available.contains(name)
                ? "Names containing whitespace cannot be represented by the model's space-separated input."
                : "This saved selection is no longer available. Uncheck it before applying.");
        }
        action->setDefaultWidget(list_); menu->addAction(action); setMenu(menu);
        connect(list_, &QListWidget::itemChanged, this, [this] { updateText(); });
        updateText();
    }
    QString value() const {
        QStringList checked;
        for (int i = 0; i < list_->count(); ++i) {
            const auto* item = list_->item(i);
            if (item->checkState() != Qt::Checked) continue;
            if (!item->data(Qt::UserRole).toBool())
                throw std::runtime_error((field_ + ": unavailable or unsupported name: " + item->text()).toStdString());
            checked.push_back(item->text());
        }
        // Retain existing selection order; append newly selected names.
        QStringList ordered;
        for (const auto& name : initial_) if (checked.removeOne(name)) ordered.push_back(name);
        ordered.append(checked);
        return ordered.join(' ');
    }
private:
    void updateText() {
        QStringList names;
        for (int i = 0; i < list_->count(); ++i)
            if (list_->item(i)->checkState() == Qt::Checked) names.push_back(list_->item(i)->text());
        setText(names.isEmpty() ? emptyText_ : names.join(", "));
    }
    QListWidget* list_;
    QStringList initial_;
    QString field_, emptyText_;
};
}

ControlsDialog::ControlsDialog(const QString& database, QWidget* parent, std::optional<model_controls::Snapshot> current)
    : QDialog(parent), database_(database)
{
    setWindowTitle("Model controls"); setObjectName("modelControlsDialog"); resize(780, 780);
    auto* layout = new QVBoxLayout(this);
    auto* description = new QLabel("Edit controls for the next run. Update default saves to the database. Use values applies only in memory; Cancel discards edits.\nDatabase: " + database);
    description->setTextFormat(Qt::PlainText); description->setWordWrap(true); layout->addWidget(description);
    auto* scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    fields_ = new QWidget;
    form_ = new QFormLayout(fields_);
    form_->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    form_->setVerticalSpacing(12);
    scroll->setWidget(fields_); layout->addWidget(scroll);
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
    try { pending_ = std::async(std::launch::async, [database] { return Loaded{model_controls::load(database), model_controls::experimentNames(database), model_controls::solvableParameters(database)}; }); }
    catch (const std::exception& error) { busy(false); save_->setEnabled(false); use_->setEnabled(false); status_->setText(QString::fromUtf8(error.what())); }
}
void ControlsDialog::busy(bool value) {
    use_->setEnabled(!value); fields_->setEnabled(!value); save_->setEnabled(!value); cancel_->setEnabled(!value);
}
void ControlsDialog::reject() {
    // A pending transaction must finish before its dialog and result are discarded.
    if (!pending_.valid()) QDialog::reject();
}
void ControlsDialog::populate() {
    for (const auto& row : current_.rows) {
        const auto type = model_controls::kind(row.name);
        if (type == model_controls::Kind::unknown) { editors_.push_back(nullptr); continue; }
        QWidget* editor;
        if (row.name == "experiment_name") {
            editor = new ChecklistPicker(experiments_, row.value.value_or(QString{}), row.name,
                "experimentChoices", "Select experiments...");
        } else if (row.name == "universal_solve_for") {
            editor = new ChecklistPicker(parameters_, row.value.value_or(QString{}), row.name,
                "parameterChoices", "None (optional)");
        } else if (type == model_controls::Kind::boolean) {
            auto* check = new QCheckBox;
            check->setChecked(row.value == std::optional<QString>("true"));
            check->setProperty("invalidValue", !row.value || (*row.value != "true" && *row.value != "false"));
            connect(check, &QCheckBox::toggled, check, [check] { check->setProperty("invalidValue", false); });
            editor = check;
        } else if (type == model_controls::Kind::scatter) {
            auto* combo = new QComboBox;
            combo->addItems({"none", "NS_ND"});
            combo->setCurrentIndex(row.value ? combo->findText(*row.value) : -1);
            editor = combo;
        } else {
            auto* edit = new QLineEdit(row.value.value_or(QString{}));
            editor = edit;
        }
        editor->setObjectName(row.name);
        editor->setToolTip(model_controls::help(row.name));
        form_->addRow(row.name, editor);
        editors_.push_back(editor);
    }
}
void ControlsDialog::save(bool persist) {
    if (pending_.valid()) return;
    auto edited = current_.rows;
    try {
        for (std::size_t i = 0; i < editors_.size(); ++i) {
            auto* editor = editors_[i];
            if (!editor) continue;
            if (auto* picker = dynamic_cast<ChecklistPicker*>(editor)) {
                const auto value = picker->value();
                edited[i].value = edited[i].name == "universal_solve_for" && value.isEmpty()
                    ? std::nullopt : std::optional<QString>(value);
            }
            else if (auto* check = qobject_cast<QCheckBox*>(editor)) {
                if (check->property("invalidValue").toBool()) edited[i].value.reset();
                else edited[i].value = check->isChecked() ? "true" : "false";
            } else if (auto* combo = qobject_cast<QComboBox*>(editor)) edited[i].value = combo->currentText();
            else {
                const auto text = qobject_cast<QLineEdit*>(editor)->text();
                edited[i].value = edited[i].name == "universal_solve_for" && text.trimmed().isEmpty()
                    ? std::nullopt : std::optional<QString>(text);
            }
        }
        for (std::size_t i = 0; i < edited.size(); ++i)
            if (editors_[i]) model_controls::validate(edited[i]);
        if (!persist) { current_.rows = std::move(edited); accept(); return; }
        pending_ = std::async(std::launch::async, [database = database_, original = original_, edited] {
            model_controls::save(database, original, edited);
            auto result = original; result.rows = edited; return Loaded{result, {}, {}};
        });
        saving_ = true; busy(true); status_->setText("Saving model controls...");
    } catch (const std::exception& error) { status_->setText(QString::fromUtf8(error.what())); }
}
void ControlsDialog::poll() {
    if (!pending_.valid() || pending_.wait_for(std::chrono::seconds(0)) != std::future_status::ready) return;
    try {
        auto result = pending_.get();
        busy(false);
        if (saving_) { current_ = std::move(result.controls); accept(); return; }
        original_ = result.controls;
        if (current_.nameColumn.isEmpty()) current_ = std::move(result.controls);
        experiments_ = std::move(result.experiments);
        parameters_ = std::move(result.parameters); populate();
        save_->setEnabled(!original_.rows.empty()); use_->setEnabled(!original_.rows.empty());
        status_->setText("All values are required except universal_solve_for. Hover over a control for help.");
    } catch (const std::exception& error) {
        busy(false); save_->setEnabled(saving_); use_->setEnabled(saving_); saving_ = false;
        status_->setText(QString::fromUtf8(error.what()));
    }
}
