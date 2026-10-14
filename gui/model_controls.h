#pragma once
#include <QString>
#include <QStringList>
#include <optional>
#include <map>
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
QStringList variableNames(const QString& database);
QStringList solvableParameters(const QString& database);
// Selected global/experiment variables and non-omitted profile-local variables.
QStringList selectedVariables(const QString& database, const Snapshot& controls);
QStringList profileDatasets(const QString& database, const QString& experiment);
QStringList applicableProfileVariables(const QString& database, const Snapshot& controls,
                                      const QString& experiment, const QString& dataset);
std::optional<double> profileVariableInitialValue(const QString& database, const Snapshot& controls,
                                                 const QString& experiment, const QString& dataset,
                                                 const QString& variable);
// Compare against the original snapshot under a write transaction. Only changed
// values are written; any conflict or failure rolls the entire transaction back.
void save(const QString& database, const Snapshot& original, const std::vector<Row>& edited);
}
