#pragma once
#include <QWidget>
#include "experiment_row_filter.h"
#include <functional>
#include <array>
class QTableWidget;
class QLabel;
class QPushButton;

// Reads experiment metadata only; opening this tab never initializes a model.
class ExperimentsTab : public QWidget {
public:
    explicit ExperimentsTab(QWidget* parent = nullptr);
    void load(const QString& database);
    void clear();
    ExperimentRowFilter* experimentFilter() const { return filter_; }
    bool dirty() const { return dirty_; }
    std::function<void()> changed;
    std::function<void()> saved;
private:
    ExperimentRowFilter* filter_ = nullptr;
    void saveDimensions();
    void setDirty(bool value);
    QString database_;
    bool dirty_ = false;
    std::array<int, 3> dimensionColumns_{-1, -1, -1};
    int nameColumn_ = -1;
    QPushButton *save_, *reload_;
    QTableWidget* table_;
    QLabel* status_;
};
