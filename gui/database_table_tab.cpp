#include "database_table_tab.h"
#include <QHeaderView>
#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QLineEdit>
#include <QScrollArea>
#include <QRegularExpression>
#include <QVariant>
#include <cmath>
#include <vector>
#include <QLabel>
#include <QPushButton>
#include <QStringList>
#include <QTableWidget>
#include <QVBoxLayout>
#include <sqlite3.h>
#include <memory>
#include <stdexcept>

DatabaseTableTab::DatabaseTableTab(Table table, QWidget* parent)
    : QWidget(parent), tableKind_(table), tableName_(table == Table::reactions ? "reactions" : table == Table::species ? "species" : table == Table::alglib ? "alglib_input" : "raw_profile") {
    auto* layout = new QVBoxLayout(this);
    status_ = new QLabel;
    status_->setObjectName(tableName_ + "Status");
    status_->setTextFormat(Qt::PlainText); status_->setWordWrap(true);
    layout->addWidget(status_);
    table_ = new QTableWidget;
    table_->setObjectName(tableName_ + "Table");
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    if (tableKind_ == Table::raw_profile) {
        filter_ = new ExperimentRowFilter(table_, "raw_profile");
        layout->addWidget(filter_);
    }
    layout->addWidget(table_, 1);
    auto* buttons = new QHBoxLayout;
    refresh_ = new QPushButton("Refresh");
    refresh_->setObjectName(tableName_ + "RefreshButton");
    add_ = new QPushButton("Add..."); add_->setObjectName(tableName_ + "AddButton");
    modify_ = new QPushButton("Modify..."); modify_->setObjectName(tableName_ + "ModifyButton");
    buttons->addWidget(add_); buttons->addWidget(modify_);
    buttons->addWidget(refresh_); buttons->addStretch(); layout->addLayout(buttons);
    connect(refresh_, &QPushButton::clicked, this, [this] { load(database_); });
    connect(add_, &QPushButton::clicked, this, [this] { editRow(true); });
    connect(modify_, &QPushButton::clicked, this, [this] { editRow(false); });
    connect(table_, &QTableWidget::itemSelectionChanged, this, [this] { updateButtons(); });
    if (tableKind_ == Table::raw_profile) { add_->hide(); modify_->hide(); }
    clear();
}

void DatabaseTableTab::clear() {
    loaded_ = false; columns_.clear();
    database_.clear();
    table_->clear(); table_->setRowCount(0); table_->setColumnCount(0);
    refresh_->setEnabled(false); updateButtons();
    if (filter_) filter_->apply();
    status_->setText("Choose a database from the File menu.");
}

void DatabaseTableTab::load(QString database) {
    clear();
    database_ = database;
    if (database.isEmpty()) return;
    refresh_->setEnabled(true);
    try {
        sqlite3* raw = nullptr;
        const int opened = sqlite3_open_v2(database.toUtf8().constData(), &raw, SQLITE_OPEN_READONLY, nullptr);
        std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)> db(raw, sqlite3_close_v2);
        auto check = [&](int rc, int expected) {
            if (rc != expected) throw std::runtime_error(sqlite3_errmsg(db.get()));
        };
        check(opened, SQLITE_OK);
        sqlite3_busy_timeout(db.get(), 1000);
        // Fixed queries: no user-controlled table or column names enter SQL.
        const char* sql = tableKind_ == Table::reactions
            ? "SELECT * FROM reactions ORDER BY REACTION_NAME"
            : tableKind_ == Table::species ? "SELECT * FROM species ORDER BY SPECIES_NAME"
            : tableKind_ == Table::alglib ? "SELECT * FROM alglib_input ORDER BY VARIABLE"
            : "SELECT * FROM raw_profile ORDER BY NAME";
        sqlite3_stmt* query = nullptr;
        const int prepared = sqlite3_prepare_v2(db.get(), sql, -1, &query, nullptr);
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(query, sqlite3_finalize);
        check(prepared, SQLITE_OK);
        table_->setColumnCount(sqlite3_column_count(query));
        QStringList headers;
        for (int col = 0; col < table_->columnCount(); ++col)
            headers << QString::fromUtf8(sqlite3_column_name(query, col));
        columns_ = headers;
        table_->setHorizontalHeaderLabels(headers);
        int rc;
        while ((rc = sqlite3_step(query)) == SQLITE_ROW) {
            const int row = table_->rowCount(); table_->insertRow(row);
            for (int col = 0; col < table_->columnCount(); ++col) {
                const bool null = sqlite3_column_type(query, col) == SQLITE_NULL;
                const auto value = null ? QString("NULL") : QString::fromUtf8(
                    reinterpret_cast<const char*>(sqlite3_column_text(query, col)), sqlite3_column_bytes(query, col));
                auto* item = new QTableWidgetItem(value);
                // Retain SQLite storage types for exact optimistic concurrency checks.
                QVariant original;
                switch (sqlite3_column_type(query, col)) {
                case SQLITE_INTEGER: original = QVariant::fromValue<qlonglong>(sqlite3_column_int64(query, col)); break;
                case SQLITE_FLOAT: original = sqlite3_column_double(query, col); break;
                case SQLITE_TEXT: original = value; break;
                case SQLITE_BLOB: original = QByteArray(static_cast<const char*>(sqlite3_column_blob(query, col)), sqlite3_column_bytes(query, col)); break;
                }
                item->setData(Qt::UserRole, original);
                item->setData(Qt::UserRole + 1, null);
                item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
                item->setToolTip(null ? "SQL NULL (no value)" : value);
                table_->setItem(row, col, item);
            }
        }
        check(rc, SQLITE_DONE);
        loaded_ = true; updateButtons();
        if (filter_) filter_->apply();
        status_->setText(QString("%1 %2 rows in %3. Use Add or select a row and choose Modify.")
            .arg(table_->rowCount()).arg(tableName_).arg(database));
        if (tableKind_ == Table::raw_profile)
            status_->setText(QString("%1 raw profile rows in %2. Read-only.").arg(table_->rowCount()).arg(database));
    } catch (const std::exception& error) {
        loaded_ = false; columns_.clear(); updateButtons();
        table_->clear(); table_->setRowCount(0); table_->setColumnCount(0);
        status_->setText("Cannot load " + tableName_ + ": " + QString::fromUtf8(error.what()));
    }
}

namespace {
QString quoteIdentifier(QString name) { return '"' + name.replace('"', "\"\"") + '"'; }
QStringList editableColumns(DatabaseTableTab::Table table) {
    if (table == DatabaseTableTab::Table::alglib)
        return {"VARIABLE", "INITIAL VALUE", "LOWER BOUND", "UPPER BOUND", "SCALE"};
    return table == DatabaseTableTab::Table::reactions ? QStringList{"REACTION_NAME", "SPECIES", "COEFFICIENTS", "Ks", "EXPONENTS"}
                     : QStringList{"SPECIES_NAME", "SPECIES_TYPE", "DIFFUSION_RATE", "QE", "PARTICLE_DIAMETER", "PARTICLE_DENSITY", "MOLECULAR_WEIGHT"};
}
void bindValue(sqlite3* db, sqlite3_stmt* query, int index, const QVariant& value) {
    int rc;
    if (!value.isValid()) rc = sqlite3_bind_null(query, index);
    else if (value.metaType().id() == QMetaType::LongLong) rc = sqlite3_bind_int64(query, index, value.toLongLong());
    else if (value.metaType().id() == QMetaType::Double) rc = sqlite3_bind_double(query, index, value.toDouble());
    else if (value.metaType().id() == QMetaType::QByteArray) {
        const auto bytes = value.toByteArray();
        rc = sqlite3_bind_blob(query, index, bytes.constData(), bytes.size(), SQLITE_TRANSIENT);
    } else {
        const auto bytes = value.toString().toUtf8();
        rc = sqlite3_bind_text(query, index, bytes.constData(), bytes.size(), SQLITE_TRANSIENT);
    }
    if (rc != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
}
}

void DatabaseTableTab::setEditingEnabled(bool enabled) {
    editingEnabled_ = enabled; updateButtons();
}

void DatabaseTableTab::updateButtons() {
    const auto expected = editableColumns(tableKind_);
    bool supported = loaded_ && editingEnabled_ && tableKind_ != Table::raw_profile;
    for (const auto& column : expected) supported = supported && columns_.contains(column);
    add_->setEnabled(supported);
    modify_->setEnabled(supported && !table_->selectedItems().isEmpty());
}

void DatabaseTableTab::editRow(bool adding) {
    if (tableKind_ == Table::raw_profile || !loaded_ || !editingEnabled_ || (!adding && table_->selectedItems().isEmpty())) return;
    const bool reactions = tableKind_ == Table::reactions;
    const bool alglib = tableKind_ == Table::alglib;
    const auto fields = editableColumns(tableKind_);
    for (const auto& field : fields) if (!columns_.contains(field)) return;
    const auto database = database_;
    const auto columns = columns_;
    const int row = table_->currentRow();
    std::vector<QVariant> original;
    if (!adding) for (int col = 0; col < table_->columnCount(); ++col)
        original.push_back(table_->item(row, col)->data(Qt::UserRole));

    QDialog dialog(this); dialog.setObjectName(tableName_ + "Editor");
    dialog.setWindowTitle((adding ? "Add " : "Modify ") + QString(reactions ? "reaction" : alglib ? "ALGLIB input" : "species"));
    dialog.resize(660, 440);
    auto* layout = new QVBoxLayout(&dialog);
    auto* note = new QLabel("Values use the database's existing units. Names cannot be changed when modifying a row, to preserve references.");
    note->setWordWrap(true); layout->addWidget(note);
    auto* scroll = new QScrollArea; scroll->setWidgetResizable(true);
    auto* panel = new QWidget; auto* form = new QFormLayout(panel);
    scroll->setWidget(panel); layout->addWidget(scroll);
    std::vector<QLineEdit*> edits;
    std::vector<QCheckBox*> nulls;
    for (const auto& field : fields) {
        const int col = columns.indexOf(field);
        const QVariant value = adding ? QVariant{} : original[col];
        auto* line = new QLineEdit(value.toString()); line->setObjectName(field);
        auto* null = new QCheckBox("NULL"); null->setObjectName(field + "Null");
        null->setChecked(!value.isValid());
        // Reaction lists and names are required; optional species numeric fields retain SQL NULL.
        const bool nullable = tableKind_ == Table::species && field != "SPECIES_NAME" && field != "SPECIES_TYPE";
        if (!nullable) null->setChecked(false);
        null->setVisible(nullable);
        line->setEnabled(!null->isChecked());
        connect(null, &QCheckBox::toggled, line, [line](bool checked) { line->setEnabled(!checked); });
        if (!adding && field == fields.front()) line->setReadOnly(true);
        auto* fieldRow = new QHBoxLayout; fieldRow->addWidget(line); fieldRow->addWidget(null);
        form->addRow(field, fieldRow); edits.push_back(line); nulls.push_back(null);
    }
    auto* error = new QLabel; error->setObjectName("referenceEditorStatus");
    error->setTextFormat(Qt::PlainText); error->setWordWrap(true); layout->addWidget(error);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Save | QDialogButtonBox::Cancel);
    buttons->button(QDialogButtonBox::Save)->setObjectName("saveReferenceButton");
    buttons->button(QDialogButtonBox::Cancel)->setObjectName("cancelReferenceButton");
    layout->addWidget(buttons);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    bool wrote = false;
    connect(buttons, &QDialogButtonBox::accepted, &dialog, [&] {
        try {
            std::vector<QVariant> values;
            for (int i = 0; i < fields.size(); ++i) {
                const QString text = edits[i]->text();
                QVariant value = nulls[i]->isChecked() ? QVariant{} : QVariant(text);
                if (value.isValid() && !reactions && i >= (alglib ? 1 : 2)) {
                    bool ok = false; const double number = text.toDouble(&ok);
                    if (!ok || !std::isfinite(number)) throw std::invalid_argument((fields[i] + (alglib ? ": enter a finite number." : ": enter a finite number or select NULL.")).toStdString());
                    value = number;
                }
                if (!adding) {
                    const auto& old = original[columns.indexOf(fields[i])];
                    if (text == old.toString() && nulls[i]->isChecked() == !old.isValid()) value = old;
                }
                values.push_back(value);
            }
            const auto name = values[0].toString();
            if (name.isEmpty() || name.contains(QRegularExpression("[\\s']")))
                throw std::invalid_argument("Names must be nonempty and contain no whitespace or apostrophes.");
            QStringList species;
            if (reactions) {
                species = values[1].toString().split(' ', Qt::SkipEmptyParts);
                if (species.isEmpty()) throw std::invalid_argument("SPECIES must contain at least one species name.");
                for (int i = 2; i < 5; ++i) {
                    const auto tokens = values[i].toString().split(' ', Qt::SkipEmptyParts);
                    const int expected = i == 3 ? 2 : species.size();
                    if (tokens.size() != expected) throw std::invalid_argument((fields[i] + ": incorrect number of values.").toStdString());
                    for (const auto& token : tokens) {
                        bool ok = false; const double number = token.toDouble(&ok);
                        if (!ok || !std::isfinite(number)) throw std::invalid_argument((fields[i] + ": values must be finite numbers.").toStdString());
                    }
                    values[i] = tokens.join(' ');
                }
                values[1] = species.join(' ');
            } else if (alglib) {
                const double initial = values[1].toDouble();
                const double lower = values[2].toDouble(), upper = values[3].toDouble();
                if (lower > upper) throw std::invalid_argument("LOWER BOUND must not exceed UPPER BOUND.");
                if (initial < lower || initial > upper) throw std::invalid_argument("INITIAL VALUE must lie within the bounds.");
                if (values[4].toDouble() == 0) throw std::invalid_argument("SCALE must be nonzero.");
            } else if (values[1].toString() != "molecule" && values[1].toString() != "particle")
                throw std::invalid_argument("SPECIES_TYPE must be molecule or particle.");

            std::vector<int> changedFields;
            for (int i = adding ? 0 : 1; i < fields.size(); ++i)
                if (adding || values[i] != original[columns.indexOf(fields[i])]) changedFields.push_back(i);
            if (!adding && changedFields.empty()) { dialog.accept(); return; }
            sqlite3* raw = nullptr;
            const int opened = sqlite3_open_v2(database.toUtf8().constData(), &raw, SQLITE_OPEN_READWRITE, nullptr);
            std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)> db(raw, sqlite3_close_v2);
            auto check = [&](int rc, int expected = SQLITE_OK) {
                if (rc != expected) throw std::runtime_error(sqlite3_errmsg(db.get()));
            };
            check(opened); sqlite3_busy_timeout(db.get(), 1000);
            check(sqlite3_exec(db.get(), "PRAGMA foreign_keys=ON", nullptr, nullptr, nullptr));
            check(sqlite3_exec(db.get(), "BEGIN IMMEDIATE", nullptr, nullptr, nullptr));
            try {
                using Statement = std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)>;
                auto prepare = [&](const QString& sql) {
                    sqlite3_stmt* query = nullptr;
                    const int rc = sqlite3_prepare_v2(db.get(), sql.toUtf8().constData(), -1, &query, nullptr);
                    Statement result(query, sqlite3_finalize); check(rc); return result;
                };
                auto identity = prepare("SELECT count(*) FROM " + quoteIdentifier(tableName_) + " WHERE " + quoteIdentifier(fields.front()) + "=?");
                bindValue(db.get(), identity.get(), 1, values[0]); check(sqlite3_step(identity.get()), SQLITE_ROW);
                if (sqlite3_column_int(identity.get(), 0) != (adding ? 0 : 1))
                    throw std::runtime_error(adding ? "That name already exists." : "The selected row is missing or ambiguous. Refresh before modifying.");
                identity.reset();
                for (const auto& specie : species) {
                    auto lookup = prepare("SELECT count(*) FROM species WHERE SPECIES_NAME=?");
                    bindValue(db.get(), lookup.get(), 1, specie); check(sqlite3_step(lookup.get()), SQLITE_ROW);
                    if (sqlite3_column_int(lookup.get(), 0) != 1) throw std::invalid_argument(("Unknown or ambiguous species: " + specie).toStdString());
                }
                QStringList assignments, placeholders, predicates;
                for (const int fieldIndex : changedFields) {
                    const auto& field = fields[fieldIndex];
                    assignments << (adding ? quoteIdentifier(field) : quoteIdentifier(field) + "=?");
                    placeholders << "?";
                }
                QString sql;
                if (adding) sql = "INSERT INTO " + quoteIdentifier(tableName_) + " (" + assignments.join(',') + ") VALUES (" + placeholders.join(',') + ")";
                else {
                    for (const auto& column : columns) predicates << quoteIdentifier(column) + " IS ?";
                    sql = "UPDATE " + quoteIdentifier(tableName_) + " SET " + assignments.join(',') + " WHERE " + predicates.join(" AND ");
                }
                auto statement = prepare(sql);
                int index = 1;
                for (const int fieldIndex : changedFields) bindValue(db.get(), statement.get(), index++, values[fieldIndex]);
                if (!adding) for (const auto& value : original) bindValue(db.get(), statement.get(), index++, value);
                check(sqlite3_step(statement.get()), SQLITE_DONE);
                if (sqlite3_changes(db.get()) != 1) throw std::runtime_error("The row changed in the database. Cancel and Refresh before modifying.");
                statement.reset();
                check(sqlite3_exec(db.get(), "COMMIT", nullptr, nullptr, nullptr));
            } catch (...) { sqlite3_exec(db.get(), "ROLLBACK", nullptr, nullptr, nullptr); throw; }
            wrote = true; dialog.accept();
        } catch (const std::exception& failure) { error->setText(QString::fromUtf8(failure.what())); }
    });
    if (dialog.exec() == QDialog::Accepted && wrote) {
        load(database);
        if (saved) saved();
    }
}
