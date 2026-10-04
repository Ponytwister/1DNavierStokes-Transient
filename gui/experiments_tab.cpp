#include "experiments_tab.h"
#include <QHeaderView>
#include <QPushButton>
#include <QSignalBlocker>
#include <channel_dimensions.h>
#include <set>
#include <charconv>
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
    filter_ = new ExperimentRowFilter(table_, "experiments");
    layout->addWidget(filter_);
    layout->addWidget(table_);
    auto* buttons = new QHBoxLayout;
    save_ = new QPushButton("Save channel dimensions"); save_->setObjectName("saveDimensionsButton");
    reload_ = new QPushButton("Reload / discard edits"); reload_->setObjectName("reloadExperimentsButton");
    buttons->addWidget(save_); buttons->addWidget(reload_); buttons->addStretch(); layout->addLayout(buttons);
    save_->setEnabled(false);
    connect(table_, &QTableWidget::itemChanged, this, [this] { setDirty(true); });
    connect(save_, &QPushButton::clicked, this, [this] { saveDimensions(); });
    connect(reload_, &QPushButton::clicked, this, [this] { const auto database = database_; clear(); load(database); });
}

void ExperimentsTab::clear() {
    const QSignalBlocker blocker(table_);
    dimensionColumns_ = {-1, -1, -1}; nameColumn_ = -1;
    setDirty(false);
    table_->clear(); table_->setRowCount(0); table_->setColumnCount(0);
    filter_->apply();
    status_->setText("Choose a database from the File menu.");
}

void ExperimentsTab::load(const QString& database) {
    if (dirty_ && database == database_) { filter_->apply(); return; }
    clear();
    database_ = database;
    const QSignalBlocker blocker(table_);
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
        nameColumn_ = headers.indexOf("NAME");
        const QStringList dimensions = {"CHANNEL_WIDTH", "CHANNEL_HEIGHT", "CHANNEL_LENGTH"};
        const QStringList labels = {"Channel Width (m)", "Channel Height (m)", "Channel Length (m)"};
        bool editable = nameColumn_ >= 0;
        for (int i = 0; i < 3; ++i) {
            dimensionColumns_[i] = headers.indexOf(dimensions[i]);
            editable = editable && dimensionColumns_[i] >= 0;
            if (dimensionColumns_[i] >= 0) headers[dimensionColumns_[i]] = labels[i];
        }
        table_->setHorizontalHeaderLabels(headers);
        table_->setEditTriggers(editable ? QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed
                                        : QAbstractItemView::NoEditTriggers);
        int rc;
        while ((rc = sqlite3_step(query)) == SQLITE_ROW) {
            const int row = table_->rowCount(); table_->insertRow(row);
            for (int col = 0; col < table_->columnCount(); ++col) {
                const bool null = sqlite3_column_type(query, col) == SQLITE_NULL;
                const auto value = null ? QString("NULL") : QString::fromUtf8(
                    reinterpret_cast<const char*>(sqlite3_column_text(query, col)), sqlite3_column_bytes(query, col));
                auto* item = new QTableWidgetItem(value);
                const bool dimension = std::find(dimensionColumns_.begin(), dimensionColumns_.end(), col) != dimensionColumns_.end();
                item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | (editable && dimension ? Qt::ItemIsEditable : Qt::NoItemFlags));
                if (dimension) {
                    const auto number = sqlite3_column_double(query, col);
                    const int type = sqlite3_column_type(query, col);
                    if ((type != SQLITE_FLOAT && type != SQLITE_INTEGER) || !std::isfinite(number) || number <= 0)
                        throw std::runtime_error("Channel dimensions must be finite positive numeric values.");
                    char buffer[64];
                    const auto formatted = std::to_chars(buffer, buffer + sizeof(buffer), number);
                    if (formatted.ec != std::errc{}) throw std::runtime_error("Cannot format channel dimension.");
                    item->setText(QString::fromLatin1(buffer, static_cast<int>(formatted.ptr - buffer)));
                    item->setData(Qt::UserRole, number);
                }
                item->setData(Qt::UserRole + 1, null);
                item->setToolTip(null ? "SQL NULL (no value)" : value);
                table_->setItem(row, col, item);
            }
        }
        check(rc, SQLITE_DONE);
        filter_->apply();
        // Keep physical dimensions visible beside NAME even in wide legacy tables.
        if (editable) {
            table_->horizontalHeader()->moveSection(table_->horizontalHeader()->visualIndex(nameColumn_), 0);
            for (int i = 0; i < 3; ++i)
                table_->horizontalHeader()->moveSection(table_->horizontalHeader()->visualIndex(dimensionColumns_[i]), i + 1);
        }
        status_->setText(QString("%1 experiments in %2. %3").arg(table_->rowCount()).arg(database).arg(editable
            ? "Double-click channel dimensions to edit (meters), then Save. Other metadata is read-only."
            : "Legacy dimensions: 5e-4 x 4e-5 x 0.025 m. Apply migration 001 to edit channel dimensions."));
    } catch (const std::exception& error) {
        clear(); status_->setText("Cannot load experiments: " + QString::fromUtf8(error.what()));
    }
}

void ExperimentsTab::setDirty(bool value) {
    dirty_ = value; save_->setEnabled(value);
    if (changed) changed();
}

void ExperimentsTab::saveDimensions() {
    if (!dirty_) return;
    try {
        // Validate the entire edit before beginning a transaction.
        std::vector<std::array<double, 3>> dimensions;
        std::set<QString> names;
        for (int row = 0; row < table_->rowCount(); ++row) {
            const auto name = table_->item(row, nameColumn_)->text();
            if (!names.insert(name).second) throw std::runtime_error("Duplicate experiment names cannot be edited.");
            std::array<double, 3> values;
            for (int i = 0; i < 3; ++i) {
                bool ok = false;
                values[i] = table_->item(row, dimensionColumns_[i])->text().toDouble(&ok);
                if (!ok) throw std::invalid_argument("Channel dimensions must be numbers in meters.");
            }
            const channel_dimensions valid(values[0], values[1], values[2]);
            dimensions.push_back(values);
        }
        sqlite3* raw = nullptr;
        const int opened = sqlite3_open_v2(database_.toUtf8().constData(), &raw, SQLITE_OPEN_READWRITE, nullptr);
        std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)> db(raw, sqlite3_close_v2);
        auto check = [&](int rc, int expected = SQLITE_OK) {
            if (rc != expected) throw std::runtime_error(sqlite3_errmsg(db.get()));
        };
        check(opened); sqlite3_busy_timeout(db.get(), 1000);
        check(sqlite3_exec(db.get(), "BEGIN IMMEDIATE", nullptr, nullptr, nullptr));
        try {
            sqlite3_stmt* query = nullptr;
            const int prepared = sqlite3_prepare_v2(db.get(),
                "UPDATE experiments SET CHANNEL_WIDTH=?, CHANNEL_HEIGHT=?, CHANNEL_LENGTH=? "
                "WHERE NAME=? AND CHANNEL_WIDTH IS ? AND CHANNEL_HEIGHT IS ? AND CHANNEL_LENGTH IS ?", -1, &query, nullptr);
            std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(query, sqlite3_finalize);
            check(prepared);
            for (int row = 0; row < table_->rowCount(); ++row) {
                bool modified = false;
                for (int i = 0; i < 3; ++i) {
                    const double original = table_->item(row, dimensionColumns_[i])->data(Qt::UserRole).toDouble();
                    modified = modified || dimensions[row][i] != original;
                    check(sqlite3_bind_double(query, i + 1, dimensions[row][i]));
                    check(sqlite3_bind_double(query, i + 5, original));
                }
                if (modified) {
                    const auto name = table_->item(row, nameColumn_)->text().toUtf8();
                    check(sqlite3_bind_text(query, 4, name.constData(), name.size(), SQLITE_TRANSIENT));
                    check(sqlite3_step(query), SQLITE_DONE);
                    if (sqlite3_changes(db.get()) != 1)
                        throw std::runtime_error("Experiment dimensions changed in the database. Reload before saving.");
                }
                check(sqlite3_reset(query));
            }
            check(sqlite3_exec(db.get(), "COMMIT", nullptr, nullptr, nullptr));
        } catch (...) {
            sqlite3_exec(db.get(), "ROLLBACK", nullptr, nullptr, nullptr); throw;
        }
        setDirty(false);
        load(database_);
        status_->setText("Channel dimensions saved in meters. Run again to calculate results.");
        if (saved) saved();
    } catch (const std::exception& error) {
        status_->setText("Cannot save channel dimensions: " + QString::fromUtf8(error.what()));
    }
}
