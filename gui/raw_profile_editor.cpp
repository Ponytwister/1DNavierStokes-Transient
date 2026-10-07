#include "raw_profile_editor.h"
#include "checklist_picker.h"
#include "model_controls.h"
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QVBoxLayout>
#include <QCheckBox>
#include <QHBoxLayout>
#include <sqlite3.h>
#include <cmath>
#include <map>
#include <memory>

namespace {
using Database = std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)>;
using Statement = std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)>;
void check(sqlite3* db, int rc, int expected = SQLITE_OK) {
    if (rc != expected) throw std::runtime_error(sqlite3_errmsg(db));
}
QString quote(QString name) { return '"' + name.replace('"', "\"\"") + '"'; }
Statement prepare(sqlite3* db, const QString& sql) {
    sqlite3_stmt* raw = nullptr;
    const int rc = sqlite3_prepare_v2(db, sql.toUtf8().constData(), -1, &raw, nullptr);
    Statement result(raw, sqlite3_finalize); check(db, rc); return result;
}
void bind(sqlite3* db, sqlite3_stmt* query, int index, const QVariant& value) {
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
    check(db, rc);
}
QString modelDefaultUnits(const QString& database, const QString& experiment) {
    if (experiment.isEmpty()) return "mg/ml";
    try {
        sqlite3* raw = nullptr;
        const int opened = sqlite3_open_v2(database.toUtf8().constData(), &raw, SQLITE_OPEN_READONLY, nullptr);
        Database db(raw, sqlite3_close_v2); check(db.get(), opened);
        auto query = prepare(db.get(), "SELECT SPECIES, SPECIE_MODEL_CONC_UNITS FROM experiments WHERE NAME COLLATE BINARY=?");
        bind(db.get(), query.get(), 1, experiment);
        if (sqlite3_step(query.get()) != SQLITE_ROW) return "mg/ml";
        const auto species = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 0))).split(' ', Qt::SkipEmptyParts);
        const auto units = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 1))).split(' ', Qt::SkipEmptyParts);
        const auto unit = units.value(species.indexOf("FITC"));
        return unit.isEmpty() ? "mg/ml" : unit;
    } catch (...) { return "mg/ml"; }
}
}

bool editRawProfile(const QString& database, const QStringList& columns,
                    const std::vector<QVariant>& original, QWidget* parent) {
    const bool adding = original.empty();
    auto old = [&](const QString& name) -> QVariant {
        const int index = columns.indexOf(name);
        return adding || index < 0 ? QVariant{} : original.at(index);
    };
    const auto experiments = model_controls::experimentNames(database);
    const auto parameters = model_controls::solvableParameters(database);
    QDialog dialog(parent); dialog.setObjectName("raw_profileEditor");
    dialog.setWindowTitle(adding ? "Add raw profile" : "Modify raw profile");
    dialog.resize(760, 700);
    auto* layout = new QVBoxLayout(&dialog);
    auto* note = new QLabel("Enter intensity samples separated by spaces. NAME and WT_PERCENT identify a profile. Values use the database's existing units.");
    note->setWordWrap(true); layout->addWidget(note);
    auto* scroll = new QScrollArea; scroll->setWidgetResizable(true);
    auto* panel = new QWidget; auto* form = new QFormLayout(panel);
    form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    scroll->setWidget(panel); layout->addWidget(scroll);
    auto* name = new QComboBox; name->setObjectName("NAME"); name->addItems(experiments);
    name->setCurrentIndex(adding ? (experiments.isEmpty() ? -1 : 0) : name->findText(old("NAME").toString(), Qt::MatchExactly));
    form->addRow("NAME", name);
    std::map<QString, QLineEdit*> edits;
    QComboBox* entranceUnits = nullptr;
    QString defaultEntranceUnits;
    if (columns.contains("ENTRANCE_CONC_UNITS"))
        defaultEntranceUnits = modelDefaultUnits(database, adding ? name->currentText() : old("NAME").toString());
    auto addField = [&](const QString& field) {
        if (!columns.contains(field)) return;
        auto* line = new QLineEdit(old(field).toString()); line->setObjectName(field);
        if (adding && field == "CHANNEL_LEFT_EDGE") line->setText("0");
        edits[field] = line;
        if (field == "ENTRANCE_CONC" && columns.contains("ENTRANCE_CONC_UNITS")) {
            auto* row = new QHBoxLayout; row->addWidget(line);
            entranceUnits = new QComboBox; entranceUnits->setObjectName("ENTRANCE_CONC_UNITS");
            const QStringList supported{"umol", "mg/ml", "wt%", "g/ml"};
            for (const auto& unit : supported) entranceUnits->addItem(unit, unit);
            const auto savedUnits = old("ENTRANCE_CONC_UNITS").toString();
            const auto selectedUnits = savedUnits.isEmpty() ? defaultEntranceUnits : savedUnits;
            int index = entranceUnits->findData(selectedUnits);
            if (index < 0) { entranceUnits->addItem(selectedUnits + " (unsupported)", selectedUnits); index = entranceUnits->count() - 1; }
            entranceUnits->setCurrentIndex(index);
            row->addWidget(entranceUnits); form->addRow(field, row);
            if (!old("ENTRANCE_CONC_UNITS").isValid() || old("ENTRANCE_CONC_UNITS").toString().isEmpty()) {
                QObject::connect(name, &QComboBox::currentTextChanged, entranceUnits, [database, entranceUnits](const QString& experiment) {
                    const auto unit = modelDefaultUnits(database, experiment);
                    int index = entranceUnits->findData(unit);
                    if (index < 0) { entranceUnits->addItem(unit + " (unsupported)", unit); index = entranceUnits->count() - 1; }
                    entranceUnits->setCurrentIndex(index);
                });
            }
        } else form->addRow(field, line);
    };
    addField("WT_PERCENT");
    ChecklistPicker* picker = nullptr;
    const QString solveField = "INDEPENDENT_PARAMETERS_TO_SOLVE_FOR";
    if (columns.contains(solveField)) {
        picker = new ChecklistPicker(parameters, old(solveField).toString(), solveField, "rawProfileParameterChoices", "None (optional)");
        picker->setObjectName(solveField); form->addRow(solveField, picker);
    }
    auto* omit = new QComboBox; omit->setObjectName("OMIT");
    omit->addItem("NULL", QVariant{}); omit->addItem("false", "false"); omit->addItem("true", "true");
    if (old("OMIT").isValid()) {
        int index = omit->findText(old("OMIT").toString(), Qt::MatchExactly);
        if (index < 0) { omit->addItem(old("OMIT").toString(), old("OMIT")); index = omit->count() - 1; }
        omit->setCurrentIndex(index);
    }
    const int initialOmit = omit->currentIndex(); form->addRow("OMIT", omit);
    for (const auto& field : QStringList{"LEFT_EDGE", "WIDTH", "ENTRANCE_CONC", "CHANNEL_LEFT_EDGE", "CHANNEL_RIGHT_EDGE"})
        addField(field);
    const auto initialIntensity = old("INTENSITY_ARRAY").toString().replace('\t', ' ');
    auto* intensity = new QPlainTextEdit(initialIntensity); intensity->setObjectName("INTENSITY_ARRAY");
    intensity->setMinimumHeight(130); form->addRow("INTENSITY_ARRAY", intensity);
    auto* error = new QLabel; error->setObjectName("rawProfileEditorStatus");
    error->setTextFormat(Qt::PlainText); error->setWordWrap(true); layout->addWidget(error);
    if (!adding && name->currentIndex() < 0) error->setText("The saved experiment no longer exists. Choose a valid NAME before saving.");
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Save | QDialogButtonBox::Cancel);
    buttons->button(QDialogButtonBox::Save)->setObjectName("saveRawProfileButton");
    buttons->button(QDialogButtonBox::Cancel)->setObjectName("cancelRawProfileButton");
    layout->addWidget(buttons);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    bool wrote = false;
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, [&] {
        try {
            if (name->currentIndex() < 0) throw std::runtime_error("NAME: choose an existing experiment.");
            std::map<QString, QVariant> values;
            values["NAME"] = name->currentText();
            for (const auto& [field, line] : edits) {
                const auto text = line->text();
                if (text.trimmed().isEmpty()) { values[field] = QVariant{}; continue; }
                QVariant value = text;
                bool ok = false;
                if (field == "CHANNEL_LEFT_EDGE" || field == "CHANNEL_RIGHT_EDGE") {
                    const int number = text.toInt(&ok);
                    if (!ok || number < 0)
                        throw std::runtime_error((field + ": enter a non-negative integer.").toStdString());
                    value = QVariant::fromValue<qlonglong>(number);
                } else {
                    const double number = text.toDouble(&ok);
                    if (!ok || !std::isfinite(number)) throw std::runtime_error((field + ": enter a finite number.").toStdString());
                    if (field == "ENTRANCE_CONC" && number < 0)
                        throw std::runtime_error("ENTRANCE_CONC: enter a non-negative concentration.");
                    value = number;
                }
                values[field] = !adding && old(field).isValid() && text == old(field).toString() ? old(field) : value;
            }
            if (adding && columns.contains("ENTRANCE_CONC") && !values.at("ENTRANCE_CONC").isValid())
                throw std::runtime_error("ENTRANCE_CONC: enter a concentration for a new profile.");
            const auto samples = intensity->toPlainText().split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            if (samples.isEmpty()) throw std::runtime_error("INTENSITY_ARRAY: enter at least one numeric sample.");
            for (const auto& sample : samples) {
                bool ok = false; const double number = sample.toDouble(&ok);
                if (!ok || !std::isfinite(number)) throw std::runtime_error("INTENSITY_ARRAY: every sample must be a finite number.");
            }
            if (values.at("CHANNEL_RIGHT_EDGE").toLongLong() > samples.size())
                throw std::runtime_error("CHANNEL_RIGHT_EDGE must be no greater than the number of INTENSITY_ARRAY samples.");
            if (values.at("CHANNEL_RIGHT_EDGE").toLongLong() <= values.at("CHANNEL_LEFT_EDGE").toLongLong())
                throw std::runtime_error("CHANNEL_RIGHT_EDGE must be greater than CHANNEL_LEFT_EDGE.");
            // The existing solver consumes tab-separated samples. Preserve untouched storage exactly.
            values["INTENSITY_ARRAY"] = !adding && intensity->toPlainText() == initialIntensity
                ? old("INTENSITY_ARRAY") : QVariant(samples.join('\t'));
            values["OMIT"] = !adding && omit->currentIndex() == initialOmit ? old("OMIT") : omit->currentData();
            if (entranceUnits) {
                if (values.at("ENTRANCE_CONC").isValid()) {
                    const auto selected = entranceUnits->currentData().toString();
                    values["ENTRANCE_CONC_UNITS"] = !adding && selected == old("ENTRANCE_CONC_UNITS").toString()
                        ? old("ENTRANCE_CONC_UNITS") : QVariant(selected);
                } else values["ENTRANCE_CONC_UNITS"] = old("ENTRANCE_CONC_UNITS");
            }
            if (picker) {
                const auto selected = picker->value();
                values[solveField] = selected.isEmpty() ? QVariant{} : QVariant(selected);
                if (!adding && selected == old(solveField).toString()) values[solveField] = old(solveField);
            }
            QStringList changed;
            for (const auto& [field, value] : values)
                if (adding || value != old(field)) changed << field;
            if (changed.isEmpty()) { dialog.accept(); return; }
            sqlite3* raw = nullptr;
            const int opened = sqlite3_open_v2(database.toUtf8().constData(), &raw, SQLITE_OPEN_READWRITE, nullptr);
            Database db(raw, sqlite3_close_v2); check(db.get(), opened);
            sqlite3_busy_timeout(db.get(), 1000);
            check(db.get(), sqlite3_exec(db.get(), "PRAGMA foreign_keys=ON", nullptr, nullptr, nullptr));
            check(db.get(), sqlite3_exec(db.get(), "BEGIN IMMEDIATE", nullptr, nullptr, nullptr));
            try {
                auto experiment = prepare(db.get(), "SELECT count(*) FROM experiments WHERE NAME COLLATE BINARY=?");
                bind(db.get(), experiment.get(), 1, values.at("NAME")); check(db.get(), sqlite3_step(experiment.get()), SQLITE_ROW);
                if (sqlite3_column_int(experiment.get(), 0) != 1) throw std::runtime_error("NAME: experiment is missing or ambiguous. Cancel and reopen the editor.");
                experiment.reset();
                auto countKey = [&](const QVariant& keyName, const QVariant& keyWeight) {
                    auto query = prepare(db.get(), "SELECT count(*) FROM raw_profile WHERE NAME IS ? AND WT_PERCENT IS ?");
                    bind(db.get(), query.get(), 1, keyName); bind(db.get(), query.get(), 2, keyWeight);
                    check(db.get(), sqlite3_step(query.get()), SQLITE_ROW); return sqlite3_column_int(query.get(), 0);
                };
                if (!adding && countKey(old("NAME"), old("WT_PERCENT")) != 1)
                    throw std::runtime_error("The selected profile is missing or ambiguous. Cancel and Refresh.");
                if ((adding || values.at("NAME") != old("NAME") || values.at("WT_PERCENT") != old("WT_PERCENT")) &&
                    countKey(values.at("NAME"), values.at("WT_PERCENT")) != 0)
                    throw std::runtime_error("A profile with this NAME and WT_PERCENT already exists.");
                QStringList assignments, placeholders, predicates;
                for (const auto& field : changed) {
                    assignments << (adding ? quote(field) : quote(field) + "=?"); placeholders << "?";
                }
                QString sql;
                if (adding) sql = "INSERT INTO raw_profile (" + assignments.join(',') + ") VALUES (" + placeholders.join(',') + ")";
                else {
                    for (const auto& field : columns) predicates << quote(field) + " IS ?";
                    sql = "UPDATE raw_profile SET " + assignments.join(',') + " WHERE " + predicates.join(" AND ");
                }
                auto statement = prepare(db.get(), sql); int index = 1;
                for (const auto& field : changed) bind(db.get(), statement.get(), index++, values.at(field));
                if (!adding) for (const auto& value : original) bind(db.get(), statement.get(), index++, value);
                check(db.get(), sqlite3_step(statement.get()), SQLITE_DONE);
                if (sqlite3_changes(db.get()) != 1) throw std::runtime_error("The profile changed in the database. Cancel and Refresh before modifying.");
                statement.reset();
                check(db.get(), sqlite3_exec(db.get(), "COMMIT", nullptr, nullptr, nullptr));
            } catch (...) { sqlite3_exec(db.get(), "ROLLBACK", nullptr, nullptr, nullptr); throw; }
            wrote = true; dialog.accept();
        } catch (const std::exception& failure) { error->setText(QString::fromUtf8(failure.what())); }
    });
    dialog.exec();
    return wrote;
}
