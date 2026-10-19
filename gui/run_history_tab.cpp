#include "run_history_tab.h"
#include <QHeaderView>
#include <QAbstractItemView>
#include <QVariant>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <sqlite3.h>
#include <memory>

RunHistoryTab::RunHistoryTab(QWidget* parent) : QWidget(parent) {
    auto* layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("Run attempts and the reaction, species, and model-control values captured for each run."));
    runs_ = new QTableWidget(0, 5, this);
    runs_->setObjectName("runHistoryTable");
    runs_->setHorizontalHeaderLabels({"Run", "Completed (UTC)", "Experiments", "Status", "Evaluations"});
    runs_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    runs_->setSelectionBehavior(QAbstractItemView::SelectRows);
    runs_->setSelectionMode(QAbstractItemView::SingleSelection);
    runs_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    layout->addWidget(runs_, 1);
    details_ = new QPlainTextEdit(this);
    details_->setObjectName("runHistoryDetails");
    details_->setReadOnly(true);
    layout->addWidget(details_, 1);
    auto* actions = new QHBoxLayout;
    auto* refreshButton = new QPushButton("Refresh", this);
    load_ = new QPushButton("Load model controls", this);
    load_->setObjectName("loadRunHistoryControls");
    load_->setEnabled(false);
    actions->addWidget(refreshButton); actions->addWidget(load_); actions->addStretch();
    layout->addLayout(actions);
    connect(refreshButton, &QPushButton::clicked, this, [this] { refresh(); });
    connect(load_, &QPushButton::clicked, this, [this] { if (loadRequested) loadRequested(); });
    connect(runs_, &QTableWidget::itemSelectionChanged, this, [this] {
        selectedRun_ = 0; controls_.clear(); details_->clear(); load_->setEnabled(false);
        const auto selected = runs_->selectedItems();
        if (selected.isEmpty()) return;
        selectedRun_ = runs_->item(selected.front()->row(), 0)->data(Qt::UserRole).toLongLong();
        sqlite3* raw = nullptr;
        if (sqlite3_open_v2(database_.toUtf8().constData(), &raw, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) { if (raw) sqlite3_close(raw); return; }
        std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
        sqlite3_stmt* statement = nullptr;
        if (sqlite3_prepare_v2(raw, "SELECT CATEGORY,ITEM,FIELD,VALUE FROM run_history_values WHERE RUN_ID=? ORDER BY CATEGORY,ITEM,FIELD", -1, &statement, nullptr) != SQLITE_OK) return;
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> query(statement, sqlite3_finalize);
        sqlite3_bind_int64(statement, 1, selectedRun_);
        QString text;
        while (sqlite3_step(statement) == SQLITE_ROW) {
            auto col = [statement](int i) { return QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(statement, i)), sqlite3_column_bytes(statement, i)); };
            const auto category = col(0), item = col(1), field = col(2), value = col(3);
            const auto storedValue = sqlite3_column_type(statement, 3) == SQLITE_NULL
                ? std::optional<std::string>{} : std::optional<std::string>{value.toStdString()};
            text += category + " · " + item + " · " + field + " = " + (storedValue ? value : QString("NULL")) + "\n";
            if (category == "control") controls_.emplace_back(item.toStdString(), storedValue);
        }
        details_->setPlainText(text);
        load_->setEnabled(!controls_.empty());
    });
}

void RunHistoryTab::setDatabase(const QString& database) { database_ = database; refresh(); }

std::vector<std::pair<std::string, std::optional<std::string>>> RunHistoryTab::selectedControls() const { return controls_; }

void RunHistoryTab::refresh() {
    runs_->setRowCount(0); details_->clear(); controls_.clear(); selectedRun_ = 0; load_->setEnabled(false);
    sqlite3* raw = nullptr;
    if (database_.isEmpty() || sqlite3_open_v2(database_.toUtf8().constData(), &raw, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) { if (raw) sqlite3_close(raw); return; }
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    sqlite3_stmt* statement = nullptr;
    if (sqlite3_prepare_v2(raw, "SELECT ID,COMPLETED_AT_UTC,EXPERIMENTS,STATUS,RESIDUAL_EVALUATIONS FROM run_history ORDER BY ID DESC", -1, &statement, nullptr) != SQLITE_OK) return;
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> query(statement, sqlite3_finalize);
    auto cell = [statement](int index) { return QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(statement, index)), sqlite3_column_bytes(statement, index)); };
    while (sqlite3_step(statement) == SQLITE_ROW) {
        const int row = runs_->rowCount(); runs_->insertRow(row);
        auto* id = new QTableWidgetItem(QString::number(sqlite3_column_int64(statement, 0)));
        id->setData(Qt::UserRole, QVariant::fromValue<qlonglong>(sqlite3_column_int64(statement, 0)));
        runs_->setItem(row,0,id);
        for (int col = 1; col < 5; ++col) runs_->setItem(row,col,new QTableWidgetItem(cell(col)));
    }
}
