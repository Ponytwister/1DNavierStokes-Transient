#pragma once
#include <QToolButton>
#include <QMenu>
#include <QListWidget>
#include <QWidgetAction>
#include <QRegularExpression>
#include <stdexcept>

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
