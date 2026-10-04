#pragma once
#include <QWidget>
#include "experiment_row_filter.h"
#include <functional>
#include <QStringList>

class QLabel;
class QPushButton;
class QTableWidget;

// Browse reference tables without loading the model or retaining a connection.
class DatabaseTableTab : public QWidget {
public:
    enum class Table { reactions, species, alglib, raw_profile };
    explicit DatabaseTableTab(Table table, QWidget* parent = nullptr);
    void load(QString database);
    void clear();
    ExperimentRowFilter* experimentFilter() const { return filter_; }
    void setEditingEnabled(bool enabled);
    std::function<void()> saved;
private:
    ExperimentRowFilter* filter_ = nullptr;
    void editRow(bool adding);
    void updateButtons();
    bool editingEnabled_ = true;
    bool loaded_ = false;
    QStringList columns_;
    Table tableKind_;
    QString database_;
    QString tableName_;
    QLabel* status_;
    QTableWidget* table_;
    QPushButton *refresh_, *add_, *modify_;
};
