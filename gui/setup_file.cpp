#include "setup_file.h"
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QSet>
#include <algorithm>
#include <stdexcept>

namespace setup_file {
namespace {
[[noreturn]] void fail(const QString& message) { throw std::runtime_error(message.toUtf8().constData()); }
void validate(const Setup& setup) {
    for (const auto& path : {setup.database, setup.outputDirectory})
        if (path.trimmed().isEmpty() || path.contains(QChar::Null)) fail("Setup paths must not be empty or contain NUL.");
    QSet<QString> names;
    for (const auto& row : setup.controls) {
        if (row.name.isEmpty() || names.contains(row.name)) fail("Empty or duplicate model control name.");
        names.insert(row.name);
        if (model_controls::kind(row.name) != model_controls::Kind::unknown) model_controls::validate(row);
    }
}
}
void save(const QString& filename, const Setup& setup) {
    validate(setup);
    if (QFileInfo(filename).absoluteFilePath() == QFileInfo(setup.database).absoluteFilePath() ||
        (!QFileInfo(filename).canonicalFilePath().isEmpty() &&
         QFileInfo(filename).canonicalFilePath() == QFileInfo(setup.database).canonicalFilePath()))
        fail("The setup file cannot replace the experiment database.");
    QJsonArray controls;
    for (const auto& row : setup.controls)
        controls.append(QJsonObject{{"name", row.name}, {"value", row.value ? QJsonValue(*row.value) : QJsonValue(QJsonValue::Null)}});
    const auto bytes = QJsonDocument(QJsonObject{{"format", "navier-setup"}, {"version", 1},
        {"database", QFileInfo(setup.database).absoluteFilePath()},
        {"outputDirectory", QFileInfo(setup.outputDirectory).absoluteFilePath()}, {"controls", controls}}).toJson();
    QSaveFile file(filename);
    if (!file.open(QIODevice::WriteOnly)) fail(file.errorString());
    if (file.write(bytes) != bytes.size()) fail(file.errorString());
    if (!file.commit()) fail(file.errorString());
}
Setup load(const QString& filename) {
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) fail(file.errorString());
    QJsonParseError error;
    const auto document = QJsonDocument::fromJson(file.readAll(), &error);
    if (error.error != QJsonParseError::NoError || !document.isObject()) fail("Invalid setup JSON: " + error.errorString());
    const auto object = document.object();
    if (object["format"] != QJsonValue("navier-setup") || object["version"] != QJsonValue(1))
        fail("Unsupported setup format or version.");
    if (!object["database"].isString() || !object["outputDirectory"].isString() || !object["controls"].isArray())
        fail("Setup requires database, outputDirectory, and controls fields.");
    Setup setup{object["database"].toString(), object["outputDirectory"].toString(), {}};
    for (const auto value : object["controls"].toArray()) {
        if (!value.isObject()) fail("Invalid control row.");
        const auto row = value.toObject();
        if (!row["name"].isString() || (!row["value"].isString() && !row["value"].isNull())) fail("Invalid control name or value.");
        setup.controls.push_back({row["name"].toString(), row["value"].isNull() ? std::nullopt : std::optional<QString>(row["value"].toString())});
    }
    validate(setup);
    const auto directory = QFileInfo(filename).absoluteDir();
    setup.database = QDir::cleanPath(directory.absoluteFilePath(setup.database));
    setup.outputDirectory = QDir::cleanPath(directory.absoluteFilePath(setup.outputDirectory));
    std::sort(setup.controls.begin(), setup.controls.end(), [](const auto& a, const auto& b) { return a.name < b.name; });
    return setup;
}
void restore(const Setup& setup, const model_controls::Snapshot& original) {
    validate(setup);
    auto rows = setup.controls;
    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) { return a.name < b.name; });
    model_controls::save(setup.database, original, rows);
}
}
