#pragma once
#include <QString>
#include <QStringList>
#include <optional>
#include <vector>

namespace model_controls {
enum class Kind { unknown, boolean, positive_integer, nonnegative_integer, tolerance, text, scatter };
struct Row {
    QString name;
    std::optional<QString> value;
    bool operator==(const Row&) const = default;
};
struct Snapshot {
    QString nameColumn, valueColumn;
    std::vector<Row> rows;
    bool operator==(const Snapshot&) const = default;
};
Kind kind(const QString& name);
QString help(const QString& name);
void validate(const Row& row);
// Opens existing databases only. Reading never loads a model or creates records.
Snapshot load(const QString& database);
enum class ExperimentReferences { reactions, species };
// Union of exact, space-separated references in the selected experiment rows.
QStringList experimentReferences(const QString& database, const QStringList& selected,
                                 ExperimentReferences references);
QStringList experimentNames(const QString& database);
QStringList solvableParameters(const QString& database);
// Compare against the original snapshot under a write transaction. Only changed
// values are written; any conflict or failure rolls the entire transaction back.
void save(const QString& database, const Snapshot& original, const std::vector<Row>& edited);
}
