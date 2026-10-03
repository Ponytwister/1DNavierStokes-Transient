#include "experiments_tab.h"
#include <QHeaderView>
#include <QLabel>
#include <QTableWidget>
#include <QVBoxLayout>
#include <sqlite3.h>
#include <memory>
#include <stdexcept>

ExperimentsTab::ExperimentsTab(QWidget* parent) : QWidget(parent) {
    auto* layout = new QVBoxLayout(this);
    status_ = new QLabel("Choose a database from the File menu.");
    status_->setObjectName("experimentsStatus");
    status_->setTextFormat(Qt::PlainText); status_->setWordWrap(true);
    layout->addWidget(status_);
    table_ = new QTableWidget;
    table_->setObjectName("experimentsTable");
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    layout->addWidget(table_);
}

void ExperimentsTab::clear() {
    table_->clear(); table_->setRowCount(0); table_->setColumnCount(0);
    status_->setText("Choose a database from the File menu.");
}

void ExperimentsTab::load(const QString& database) {
    clear();
    if (database.trimmed().isEmpty()) return;
    try {
        sqlite3* raw = nullptr;
        const int opened = sqlite3_open_v2(database.toUtf8().constData(), &raw, SQLITE_OPEN_READONLY, nullptr);
        std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)> db(raw, sqlite3_close_v2);
        auto check = [&](int rc, int expected) {
            if (rc != expected) throw std::runtime_error(sqlite3_errmsg(db.get()));
        };
        check(opened, SQLITE_OK);
        sqlite3_busy_timeout(db.get(), 1000);
        sqlite3_stmt* query = nullptr;
        const int prepared = sqlite3_prepare_v2(db.get(), "SELECT * FROM experiments ORDER BY NAME", -1, &query, nullptr);
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(query, sqlite3_finalize);
        check(prepared, SQLITE_OK);
        table_->setColumnCount(sqlite3_column_count(query));
        QStringList headers;
        for (int col = 0; col < table_->columnCount(); ++col)
            headers << QString::fromUtf8(sqlite3_column_name(query, col));
        table_->setHorizontalHeaderLabels(headers);
        int rc;
        while ((rc = sqlite3_step(query)) == SQLITE_ROW) {
            const int row = table_->rowCount(); table_->insertRow(row);
            for (int col = 0; col < table_->columnCount(); ++col) {
                const bool null = sqlite3_column_type(query, col) == SQLITE_NULL;
                const auto value = null ? QString("NULL") : QString::fromUtf8(
                    reinterpret_cast<const char*>(sqlite3_column_text(query, col)), sqlite3_column_bytes(query, col));
                auto* item = new QTableWidgetItem(value);
                item->setToolTip(null ? "SQL NULL (no value)" : value);
                table_->setItem(row, col, item);
            }
        }
        check(rc, SQLITE_DONE);
        status_->setText(QString("%1 experiments in %2").arg(table_->rowCount()).arg(database));
    } catch (const std::exception& error) {
        clear(); status_->setText("Cannot load experiments: " + QString::fromUtf8(error.what()));
    }
}
