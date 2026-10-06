#include "model_controls.h"
#include <sqlite3.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace model_controls {
namespace {
[[noreturn]] void fail(const QString& message) { throw std::runtime_error(message.toUtf8().constData()); }
using Database = std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)>;
using Statement = std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)>;
Database open(const QString& path, bool writable) {
    sqlite3* raw = nullptr;
    const int rc = sqlite3_open_v2(path.toUtf8().constData(), &raw,
        writable ? SQLITE_OPEN_READWRITE : SQLITE_OPEN_READONLY, nullptr);
    Database db(raw, sqlite3_close_v2);
    if (rc != SQLITE_OK) fail("Cannot open controls database: " + QString::fromUtf8(sqlite3_errmsg(raw)));
    sqlite3_busy_timeout(db.get(), 1000);
    return db;
}
Statement prepare(sqlite3* db, const QString& sql) {
    sqlite3_stmt* raw = nullptr;
    const int rc = sqlite3_prepare_v2(db, sql.toUtf8().constData(), -1, &raw, nullptr);
    Statement statement(raw, sqlite3_finalize);
    if (rc != SQLITE_OK) fail(QString::fromUtf8(sqlite3_errmsg(db)));
    return statement;
}
void execute(sqlite3* db, const char* sql) {
    if (sqlite3_exec(db, sql, nullptr, nullptr, nullptr) != SQLITE_OK)
        fail(QString::fromUtf8(sqlite3_errmsg(db)));
}
QString column(sqlite3_stmt* statement, int index) {
    return QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(statement, index)),
                             sqlite3_column_bytes(statement, index));
}
QString identifier(QString value) { return '"' + value.replace('"', "\"\"") + '"'; }
Snapshot read(sqlite3* db) {
    auto statement = prepare(db, "SELECT * FROM model_controls");
    if (sqlite3_column_count(statement.get()) < 2) fail("model_controls must have at least two columns.");
    Snapshot result{QString::fromUtf8(sqlite3_column_name(statement.get(), 0)),
                    QString::fromUtf8(sqlite3_column_name(statement.get(), 1)), {}};
    int rc;
    while ((rc = sqlite3_step(statement.get())) == SQLITE_ROW) {
        if (sqlite3_column_type(statement.get(), 0) == SQLITE_NULL) fail("A model control has no name.");
        Row row{column(statement.get(), 0), std::nullopt};
        if (sqlite3_column_type(statement.get(), 1) != SQLITE_NULL) row.value = column(statement.get(), 1);
        result.rows.push_back(std::move(row));
    }
    if (rc != SQLITE_DONE) fail(QString::fromUtf8(sqlite3_errmsg(db)));
    std::sort(result.rows.begin(), result.rows.end(), [](const Row& a, const Row& b) { return a.name < b.name; });
    for (std::size_t i = 1; i < result.rows.size(); ++i)
        if (result.rows[i-1].name == result.rows[i].name) fail("Duplicate model control: " + result.rows[i].name);
    return result;
}
void bind(sqlite3* db, sqlite3_stmt* statement, int index, const std::optional<QString>& value) {
    const auto bytes = value ? value->toUtf8() : QByteArray{};
    int rc = value ? sqlite3_bind_text(statement, index, bytes.constData(), bytes.size(), SQLITE_TRANSIENT)
                   : sqlite3_bind_null(statement, index);
    if (rc != SQLITE_OK) fail(QString::fromUtf8(sqlite3_errmsg(db)));
}
}

Kind kind(const QString& name) {
    if (name == "width resolution (X)" || name == "length/time resolution (Z)") return Kind::positive_integer;
    if (name == "exp_left_padding" || name == "exp_right_padding" || name == "max_iterations" || name == "debug_level") return Kind::nonnegative_integer;
    if (name == "run_solver" || name == "disable_reactions" || name == "disable_reverse_reactions" ||
        name == "use_alglib_init_values") return Kind::boolean;
    if (name == "convergence_epsx") return Kind::tolerance;
    if (name == "experiment_name" || name == "universal_solve_for") return Kind::text;
    if (name == "scatter_correction_type") return Kind::scatter;
    return Kind::unknown;
}
QString help(const QString& name) {
    if (name == "max_iterations") return "Nonnegative integer; 0 leaves the optimizer iteration limit unset.";
    if (name == "experiment_name") return "Select one or more experiments from the database.";
    if (name == "universal_solve_for") return "Select global parameters to solve for. No selections means NULL; parameter links are unchanged.";
    if (name == "convergence_epsx") return "Numeric tolerance strictly greater than 0 and less than 1e-3. Scientific notation is accepted.";
    if (name == "debug_level") return "Integer from 0 to 6; larger values suppress more messages.";
    if (name == "run_solver") return "true enables parameter fitting; false evaluates the model without fitting.";
    if (name == "scatter_correction_type") return "Supported correction methods: none and NS_ND.";
    if (kind(name) == Kind::positive_integer) return "Positive integer resolution, within the model's integer range.";
    if (kind(name) == Kind::nonnegative_integer) return "Nonnegative integer padding, in profile samples.";
    if (kind(name) == Kind::boolean) return "Choose true or false.";
    return "Not recognized by the current model loader; preserved without editing.";
}
void validate(const Row& row) {
    if (kind(row.name) == Kind::unknown) fail("Unsupported model control: " + row.name);
    if (!row.value || row.value->trimmed().isEmpty()) {
        if (row.name == "universal_solve_for") return;
        fail(row.name + ": a value is required.");
    }
    const auto& value = *row.value;
    bool ok = true;
    switch (kind(row.name)) {
    case Kind::positive_integer:
    case Kind::nonnegative_integer: {
        const int number = value.toInt(&ok);
        ok = ok && number >= (kind(row.name) == Kind::positive_integer ? 1 : 0);
        if (row.name == "debug_level") ok = ok && number <= 6;
        break;
    }
    case Kind::tolerance: {
        const double number = value.toDouble(&ok);
        ok = ok && std::isfinite(number) && number > 0 && number < 1e-3;
        break;
    }
    case Kind::boolean: ok = value == "true" || value == "false"; break;
    case Kind::scatter: ok = value == "none" || value == "NS_ND"; break;
    case Kind::text:
        ok = !value.contains(QChar::Null) && !value.contains('\n') && !value.contains('\r');
        if (row.name == "experiment_name") ok = ok && !value.trimmed().isEmpty();
        break;
    default: break;
    }
    if (!ok) fail(row.name + ": " + help(row.name));
}
QStringList experimentNames(const QString& database) {
    auto db = open(database, false);
    auto statement = prepare(db.get(), "SELECT DISTINCT NAME FROM experiments WHERE NAME IS NOT NULL ORDER BY NAME");
    QStringList names;
    int rc;
    while ((rc = sqlite3_step(statement.get())) == SQLITE_ROW) {
        const auto name = column(statement.get(), 0);
        if (!name.isEmpty()) names.push_back(name);
    }
    if (rc != SQLITE_DONE) fail(QString::fromUtf8(sqlite3_errmsg(db.get())));
    return names;
}
QStringList experimentReferences(const QString& database, const QStringList& selected,
                                 ExperimentReferences references) {
    if (selected.isEmpty()) return {};
    auto db = open(database, false);
    auto statement = prepare(db.get(), references == ExperimentReferences::reactions
        ? "SELECT NAME, REACTIONS FROM experiments" : "SELECT NAME, SPECIES FROM experiments");
    QStringList names;
    int rc;
    while ((rc = sqlite3_step(statement.get())) == SQLITE_ROW) {
        if (sqlite3_column_type(statement.get(), 0) == SQLITE_NULL ||
            !selected.contains(column(statement.get(), 0)) ||
            sqlite3_column_type(statement.get(), 1) == SQLITE_NULL) continue;
        names.append(column(statement.get(), 1).split(' ', Qt::SkipEmptyParts));
    }
    if (rc != SQLITE_DONE) fail(QString::fromUtf8(sqlite3_errmsg(db.get())));
    names.removeDuplicates();
    return names;
}
QStringList solvableParameters(const QString& database) {
    QStringList names{"p1", "kon1", "keq1", "left_edge", "width", "QE1"};
    auto db = open(database, false);
    auto statement = prepare(db.get(), "SELECT REACTIONS FROM experiments WHERE REACTIONS IS NOT NULL");
    int rc;
    while ((rc = sqlite3_step(statement.get())) == SQLITE_ROW) {
        // Experiment reactions use the model's space-separated name format.
        if (column(statement.get(), 0).split(' ', Qt::SkipEmptyParts).size() > 1) {
            names.append({"p2", "kon2", "keq2", "QE2"});
            return names;
        }
    }
    if (rc != SQLITE_DONE) fail(QString::fromUtf8(sqlite3_errmsg(db.get())));
    return names;
}
QStringList selectedVariables(const QString& database, const Snapshot& controls) {
    QStringList names, experiments;
    for (const auto& row : controls.rows) {
        const auto tokens = row.value.value_or(QString{}).split(' ', Qt::SkipEmptyParts);
        if (row.name == "universal_solve_for") names.append(tokens);
        if (row.name == "experiment_name") experiments = tokens;
    }
    auto db = open(database, false);
    for (const auto* table : {"experiments", "raw_profile"}) {
        const bool profiles = QString(table) == "raw_profile";
        auto statement = prepare(db.get(), "SELECT * FROM " + QString(table));
        int name = -1, variables = -1, omit = -1;
        for (int i = 0; i < sqlite3_column_count(statement.get()); ++i) {
            const QString field = QString::fromUtf8(sqlite3_column_name(statement.get(), i));
            if (field == "NAME") name = i;
            if (field == (profiles ? "INDEPENDENT_PARAMETERS_TO_SOLVE_FOR" : "PARAMETERS_TO_SOLVE_FOR")) variables = i;
            if (field == "OMIT") omit = i;
        }
        // Older databases can omit the optional parameter selection columns.
        if (name < 0 || variables < 0) continue;
        int rc;
        while ((rc = sqlite3_step(statement.get())) == SQLITE_ROW) {
            if (!experiments.contains(column(statement.get(), name))) continue;
            // Match the model's space-padded, case-insensitive true flag.
            auto omitted = omit < 0 ? QString{} : column(statement.get(), omit);
            while (omitted.startsWith(' ')) omitted.remove(0, 1);
            while (omitted.endsWith(' ')) omitted.chop(1);
            if (profiles && omitted.compare("true", Qt::CaseInsensitive) == 0) continue;
            names.append(column(statement.get(), variables).split(' ', Qt::SkipEmptyParts));
        }
        if (rc != SQLITE_DONE) fail(QString::fromUtf8(sqlite3_errmsg(db.get())));
    }
    names.removeDuplicates();
    return names;
}
Snapshot load(const QString& database) { auto db = open(database, false); return read(db.get()); }
void save(const QString& database, const Snapshot& original, const std::vector<Row>& edited) {
    if (edited.size() != original.rows.size()) fail("Control rows cannot be added or removed here.");
    for (std::size_t i = 0; i < edited.size(); ++i) {
        if (edited[i].name != original.rows[i].name) fail("Control names cannot be changed.");
        if (kind(edited[i].name) != Kind::unknown || edited[i] != original.rows[i]) validate(edited[i]);
    }
    auto db = open(database, true);
    execute(db.get(), "BEGIN IMMEDIATE");
    try {
        if (read(db.get()) != original) fail("Model controls changed in the database. Close and reopen the editor before saving.");
        auto update = prepare(db.get(), "UPDATE model_controls SET " + identifier(original.valueColumn) +
            "=? WHERE " + identifier(original.nameColumn) + "=?");
        for (std::size_t i = 0; i < edited.size(); ++i) {
            if (edited[i] == original.rows[i]) continue;
            sqlite3_reset(update.get());
            bind(db.get(), update.get(), 1, edited[i].value);
            bind(db.get(), update.get(), 2, edited[i].name);
            if (sqlite3_step(update.get()) != SQLITE_DONE) fail(QString::fromUtf8(sqlite3_errmsg(db.get())));
            if (sqlite3_changes(db.get()) != 1) fail("Expected to update exactly one model control.");
        }
        execute(db.get(), "COMMIT");
    } catch (...) { sqlite3_exec(db.get(), "ROLLBACK", nullptr, nullptr, nullptr); throw; }
}
}
