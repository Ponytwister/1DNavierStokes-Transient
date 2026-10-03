#pragma once
#include <QWidget>
class QTableWidget;
class QLabel;

// Reads experiment metadata only; opening this tab never initializes a model.
class ExperimentsTab : public QWidget {
public:
    explicit ExperimentsTab(QWidget* parent = nullptr);
    void load(const QString& database);
    void clear();
private:
    QTableWidget* table_;
    QLabel* status_;
};
