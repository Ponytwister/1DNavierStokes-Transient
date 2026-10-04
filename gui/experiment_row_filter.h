#pragma once
#include <QCheckBox>
#include <QLabel>
#include <QSet>
#include <QStringList>
#include <QTableWidget>
#include <QVBoxLayout>
#include <functional>
#include <stdexcept>

// Hides rows without removing them, preserving unsaved experiment edits.
class ExperimentRowFilter : public QWidget {
public:
    ExperimentRowFilter(QTableWidget* table, const QString& name) : table_(table) {
        auto* layout = new QVBoxLayout(this);
        layout->setContentsMargins(0, 0, 0, 0);
        enabled_ = new QCheckBox("Only selected experiments");
        enabled_->setObjectName(name + "SelectedOnly");
        enabled_->setToolTip("Show rows whose NAME is selected in Model controls.");
        enabled_->setChecked(true);
        error_ = new QLabel;
        error_->setObjectName(name + "FilterStatus");
        error_->setTextFormat(Qt::PlainText); error_->setWordWrap(true);
        layout->addWidget(enabled_); layout->addWidget(error_);
        error_->hide();
        connect(enabled_, &QCheckBox::toggled, this, [this] { apply(); });
    }
    std::function<QStringList()> selectedNames;
    void apply() {
        error_->clear(); error_->hide();
        QSet<QString> names;
        if (enabled_->isChecked() && table_->rowCount() > 0) {
            try {
                if (!selectedNames) throw std::runtime_error("No model controls selection available.");
                const auto selected = selectedNames();
                names = QSet<QString>(selected.begin(), selected.end());
            } catch (const std::exception& error) {
                error_->setText("Cannot load experiment selection: " + QString::fromUtf8(error.what())
                    + " Uncheck Only selected experiments to view all rows.");
                error_->show();
            }
        }
        int nameColumn = -1;
        for (int col = 0; col < table_->columnCount(); ++col)
            if (table_->horizontalHeaderItem(col)->text() == "NAME") nameColumn = col;
        for (int row = 0; row < table_->rowCount(); ++row) {
            const auto* item = nameColumn < 0 ? nullptr : table_->item(row, nameColumn);
            const bool selected = item && !item->data(Qt::UserRole + 1).toBool() && names.contains(item->text());
            table_->setRowHidden(row, enabled_->isChecked() && !selected);
        }
    }
private:
    QTableWidget* table_;
    QCheckBox* enabled_;
    QLabel* error_;
};
