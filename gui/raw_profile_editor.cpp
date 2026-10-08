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
#include <QGroupBox>
#include <QMap>
#include <QRegularExpression>
#include <sqlite3.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <optional>

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
struct SpeciesInfo { QString name; QString type; QString modelUnits; };
struct ExperimentCatalog { std::vector<SpeciesInfo> species; int entranceCount = 0; };
struct ConcentrationValue { int entranceNumber; QString species; double value; QString units; };
std::optional<double> scalarConcentration(const QVariant& value) {
    if (!value.isValid() || value.isNull()) return std::nullopt;
    bool ok = false;
    const double number = value.toDouble(&ok);
    if (!ok || !std::isfinite(number) || number < 0.0) return std::nullopt;
    return number;
}
ExperimentCatalog speciesCatalog(sqlite3* db, const QString& experiment) {
    ExperimentCatalog catalog;
    auto row = prepare(db, "SELECT SPECIES, SPECIE_MODEL_CONC_UNITS, ENTRANCE_FLOWRATE FROM experiments WHERE NAME COLLATE BINARY=?");
    bind(db, row.get(), 1, experiment);
    if (sqlite3_step(row.get()) != SQLITE_ROW) throw std::runtime_error("Cannot load species for the selected experiment.");
    const auto speciesText = sqlite3_column_type(row.get(), 0) == SQLITE_NULL ? QString{}
        : QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(row.get(), 0)));
    const auto unitsText = sqlite3_column_type(row.get(), 1) == SQLITE_NULL ? QString{}
        : QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(row.get(), 1)));
    const auto flowrateText = sqlite3_column_type(row.get(), 2) == SQLITE_NULL ? QString{}
        : QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(row.get(), 2)));
    const auto names = speciesText.split(' ', Qt::SkipEmptyParts);
    const auto units = unitsText.split(' ', Qt::SkipEmptyParts);
    catalog.entranceCount = flowrateText.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts).size();
    auto types = prepare(db, "SELECT SPECIES_NAME, SPECIES_TYPE FROM species");
    QMap<QString, QString> typeByName;
    int rc;
    while ((rc = sqlite3_step(types.get())) == SQLITE_ROW) {
        const auto name = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(types.get(), 0)));
        const auto type = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(types.get(), 1)));
        typeByName[name] = type;
    }
    if (rc != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
    for (int i = 0; i < names.size(); ++i)
        catalog.species.push_back({names[i], typeByName.value(names[i]), units.value(i)});
    return catalog;
}
const SpeciesInfo* speciesInfo(const ExperimentCatalog& catalog, const QString& species) {
    for (const auto& item : catalog.species) if (item.name == species) return &item;
    return nullptr;
}
QString defaultSpecies(const ExperimentCatalog& catalog) {
    if (speciesInfo(catalog, "FITC")) return "FITC";
    return catalog.species.empty() ? QString{} : catalog.species.front().name;
}
QString modelUnits(const ExperimentCatalog& catalog, const QString& species) {
    const auto* info = speciesInfo(catalog, species);
    return info ? info->modelUnits : QString{};
}
QStringList supportedUnits(const QString& type) {
    QStringList units{"umol", "mg/ml", "wt%", "g/ml"};
    if (type == "particle") units << "um2/ul" << "nm2/ul" << "mm2/nl";
    return units;
}
std::vector<ConcentrationValue> readConcentrations(sqlite3* db, const QVariant& name, const QVariant& weight) {
    sqlite3_stmt* tableRaw = nullptr;
    const int prepared = sqlite3_prepare_v2(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='raw_profile_entrance_concentrations'",
        -1, &tableRaw, nullptr);
    Statement table(tableRaw, sqlite3_finalize); check(db, prepared);
    const int found = sqlite3_step(table.get());
    if (found != SQLITE_ROW) {
        if (found != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
        throw std::runtime_error("Apply migrations 003 and 004 to edit per-species inlet concentrations.");
    }
    auto columns = prepare(db, "PRAGMA table_info(raw_profile_entrance_concentrations)");
    bool hasEntranceNumber = false;
    int columnRc;
    while ((columnRc = sqlite3_step(columns.get())) == SQLITE_ROW) {
        const auto* name = reinterpret_cast<const char*>(sqlite3_column_text(columns.get(), 1));
        if (name && QString::fromUtf8(name) == "ENTRANCE_NUMBER") hasEntranceNumber = true;
    }
    if (columnRc != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
    if (!hasEntranceNumber) throw std::runtime_error("Apply migration 004_raw_profile_entrance_numbers.sql to edit per-inlet concentrations.");
    auto query = prepare(db, "SELECT ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS FROM raw_profile_entrance_concentrations WHERE NAME IS ? AND WT_PERCENT IS ? ORDER BY ENTRANCE_NUMBER, SPECIES_NAME");
    bind(db, query.get(), 1, name); bind(db, query.get(), 2, weight);
    std::vector<ConcentrationValue> result;
    int rc;
    while ((rc = sqlite3_step(query.get())) == SQLITE_ROW) {
        const auto species = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 1)));
        const auto units = sqlite3_column_type(query.get(), 3) == SQLITE_NULL ? QString{}
            : QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 3)));
        result.push_back({sqlite3_column_int(query.get(), 0), species, sqlite3_column_double(query.get(), 2), units});
    }
    if (rc != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
    return result;
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
    sqlite3* catalogRaw = nullptr;
    const int catalogOpen = sqlite3_open_v2(database.toUtf8().constData(), &catalogRaw, SQLITE_OPEN_READONLY, nullptr);
    Database catalogDb(catalogRaw, sqlite3_close_v2); check(catalogDb.get(), catalogOpen);
    auto catalog = speciesCatalog(catalogDb.get(), adding ? name->currentText() : old("NAME").toString());
    std::vector<ConcentrationValue> originalConcentrations;
    if (!adding) {
        originalConcentrations = readConcentrations(catalogDb.get(), old("NAME"), old("WT_PERCENT"));
        for (const auto& concentration : originalConcentrations)
            if (concentration.entranceNumber > catalog.entranceCount)
                throw std::runtime_error("ENTRANCE_FLOWRATE has fewer inlets than this profile's saved concentration rows.");
        if (originalConcentrations.empty() && scalarConcentration(old("ENTRANCE_CONC"))) {
            const QString defaultUnits = modelUnits(catalog, "FITC");
            const QString savedUnits = old("ENTRANCE_CONC_UNITS").toString();
            originalConcentrations.push_back({1, "FITC", *scalarConcentration(old("ENTRANCE_CONC")), savedUnits.isEmpty() ? defaultUnits : savedUnits});
        }
    } else readConcentrations(catalogDb.get(), QVariant{}, QVariant{});
    auto addField = [&](const QString& field) {
        if (!columns.contains(field)) return;
        auto* line = new QLineEdit(old(field).toString()); line->setObjectName(field);
        if (adding && field == "CHANNEL_LEFT_EDGE") line->setText("0");
        edits[field] = line;
        form->addRow(field, line);
    };
    addField("WT_PERCENT");
    QLineEdit* entranceAlias = nullptr;
    if (columns.contains("ENTRANCE_CONC")) {
        const auto saved = old("ENTRANCE_CONC").toString();
        entranceAlias = new QLineEdit(saved.startsWith('#') ? saved : QString{});
        entranceAlias->setObjectName("entranceConcAlias");
        entranceAlias->setPlaceholderText("#variable_name (optional)");
        form->addRow("FITC concentration alias", entranceAlias);
    }
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
    for (const auto& field : QStringList{"LEFT_EDGE", "WIDTH", "CHANNEL_LEFT_EDGE", "CHANNEL_RIGHT_EDGE"})
        addField(field);
    struct ConcentrationWidgets { int entranceNumber; QString species; QLineEdit* value; QComboBox* units; };
    std::vector<ConcentrationWidgets> concentrationWidgets;
    auto* concentrationPanel = new QWidget;
    auto* concentrationLayout = new QVBoxLayout(concentrationPanel);
    concentrationLayout->setContentsMargins(0, 0, 0, 0);
    form->addRow("ENTRANCE_CONC", concentrationPanel);
    auto rebuildConcentrationRows = [&] {
        while (auto* item = concentrationLayout->takeAt(0)) {
            if (item->widget()) delete item->widget();
            delete item;
        }
        concentrationWidgets.clear();
        if (catalog.entranceCount < 1 || catalog.species.empty()) {
            concentrationLayout->addWidget(new QLabel("Select an experiment with species and one flowrate value per inlet.", concentrationPanel));
            return;
        }
        for (int entrance = 1; entrance <= catalog.entranceCount; ++entrance) {
            auto* group = new QGroupBox(QString("Entrance / inlet %1").arg(entrance), concentrationPanel);
            auto* groupLayout = new QFormLayout(group);
            auto* headings = new QWidget(group);
            auto* headingLayout = new QHBoxLayout(headings);
            headingLayout->setContentsMargins(0, 0, 0, 0);
            headingLayout->addWidget(new QLabel("Concentration", headings));
            headingLayout->addWidget(new QLabel("Units", headings));
            groupLayout->addRow(QString{}, headings);
            for (int speciesIndex = 0; speciesIndex < static_cast<int>(catalog.species.size()); ++speciesIndex) {
                const auto& specie = catalog.species.at(speciesIndex);
                const auto saved = std::find_if(originalConcentrations.cbegin(), originalConcentrations.cend(),
                    [&](const ConcentrationValue& value) { return value.entranceNumber == entrance && value.species == specie.name; });
                const QString initialValue = saved == originalConcentrations.cend()
                    ? QString{} : QString::number(saved->value, 'g', 17);
                const QString initialUnits = saved == originalConcentrations.cend() ? specie.modelUnits : saved->units;
                auto* row = new QWidget(group);
                auto* rowLayout = new QHBoxLayout(row); rowLayout->setContentsMargins(0, 0, 0, 0);
                auto* value = new QLineEdit(initialValue, row);
                value->setObjectName(QString("entranceConcValue_%1_%2").arg(entrance).arg(speciesIndex));
                auto* units = new QComboBox(row);
                units->setObjectName(QString("entranceConcUnits_%1_%2").arg(entrance).arg(speciesIndex));
                for (const auto& unit : supportedUnits(specie.type)) units->addItem(unit, unit);
                int unitIndex = units->findData(initialUnits);
                if (unitIndex < 0 && !initialUnits.isEmpty()) {
                    units->addItem(initialUnits + " (unsupported)", initialUnits);
                    unitIndex = units->count() - 1;
                }
                if (unitIndex < 0 && units->count() > 0) unitIndex = 0;
                units->setCurrentIndex(unitIndex);
                rowLayout->addWidget(value); rowLayout->addWidget(units);
                auto* speciesLabel = new QLabel(specie.name, group);
                speciesLabel->setObjectName(QString("entranceConcSpecies_%1_%2").arg(entrance).arg(speciesIndex));
                groupLayout->addRow(speciesLabel, row);
                concentrationWidgets.push_back({entrance, specie.name, value, units});
            }
            concentrationLayout->addWidget(group);
        }
    };
    rebuildConcentrationRows();
    const auto initialIntensity = old("INTENSITY_ARRAY").toString().replace('\t', ' ');
    auto* intensity = new QPlainTextEdit(initialIntensity); intensity->setObjectName("INTENSITY_ARRAY");
    intensity->setMinimumHeight(130); form->addRow("INTENSITY_ARRAY", intensity);
    auto* error = new QLabel; error->setObjectName("rawProfileEditorStatus");
    error->setTextFormat(Qt::PlainText); error->setWordWrap(true); layout->addWidget(error);
    if (!adding && name->currentIndex() < 0) error->setText("The saved experiment no longer exists. Choose a valid NAME before saving.");
    QObject::connect(name, &QComboBox::currentTextChanged, &dialog, [&](const QString& experiment) {
        try {
            catalog = speciesCatalog(catalogDb.get(), experiment);
            rebuildConcentrationRows();
        } catch (const std::exception& failure) {
            error->setText(QString::fromUtf8(failure.what()));
        }
    });

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
            std::vector<ConcentrationValue> editedConcentrations;
            QMap<QString, ConcentrationValue> concentrationByKey;
            bool allConcentrationsEntered = true;
            for (const auto& entry : concentrationWidgets) {
                const auto text = entry.value->text().trimmed();
                if (text.isEmpty()) { allConcentrationsEntered = false; continue; }
                bool ok = false;
                const double concentration = text.toDouble(&ok);
                if (!ok || !std::isfinite(concentration) || concentration < 0)
                    throw std::runtime_error("ENTRANCE_CONC: each value must be a finite, non-negative number.");
                const auto units = entry.units->currentData().toString();
                if (units.isEmpty()) throw std::runtime_error(("ENTRANCE_CONC: choose units for " + entry.species).toStdString());
                const ConcentrationValue value{entry.entranceNumber, entry.species, concentration, units};
                const QString key = QString::number(entry.entranceNumber) + QChar(0x1f) + entry.species;
                concentrationByKey.insert(key, value);
                editedConcentrations.push_back(value);
            }
            if ((adding || !old("INLET_COND_ID").isValid()) && !allConcentrationsEntered)
                throw std::runtime_error("ENTRANCE_CONC: enter a concentration for every species at every inlet.");
            if ((adding || !old("INLET_COND_ID").isValid()) && concentrationWidgets.empty())
                throw std::runtime_error("ENTRANCE_CONC: the selected experiment needs species and a flowrate value for every inlet.");
            auto findSpeciesConcentration = [](const std::vector<ConcentrationValue>& entries, int entrance, const QString& species) -> const ConcentrationValue* {
                for (const auto& entry : entries) if (entry.entranceNumber == entrance && entry.species == species) return &entry;
                return nullptr;
            };
            const auto fitc = findSpeciesConcentration(editedConcentrations, 1, "FITC");
            const auto aliasText = entranceAlias ? entranceAlias->text().trimmed() : QString{};
            if (!aliasText.isEmpty()) {
                if (!aliasText.startsWith('#') || aliasText.size() < 2 || aliasText.mid(1).contains(QRegularExpression("\\s")))
                    throw std::runtime_error("ENTRANCE_CONC alias: enter #variable_name with no spaces.");
                if (!model_controls::variableNames(database).contains(aliasText.mid(1)))
                    throw std::runtime_error("ENTRANCE_CONC alias must match a row in Variables.");
                if (!fitc || catalog.entranceCount < 1)
                    throw std::runtime_error("ENTRANCE_CONC alias requires FITC and at least one inlet.");
            }
            if (columns.contains("ENTRANCE_CONC")) values["ENTRANCE_CONC"] = !aliasText.isEmpty()
                ? QVariant(aliasText) : (fitc ? QVariant(fitc->value) : QVariant{});
            if (columns.contains("ENTRANCE_CONC_UNITS")) values["ENTRANCE_CONC_UNITS"] =
                fitc ? QVariant(fitc->units) : QVariant{};
            auto concentrationKey = [](int entrance, const QString& species) {
                return QString::number(entrance) + QChar(0x1f) + species;
            };
            QMap<QString, ConcentrationValue> originalByKey;
            for (const auto& entry : originalConcentrations) originalByKey[concentrationKey(entry.entranceNumber, entry.species)] = entry;
            bool concentrationsChanged = originalByKey.size() != concentrationByKey.size();
            if (!concentrationsChanged) {
                for (auto it = originalByKey.cbegin(); it != originalByKey.cend(); ++it) {
                    if (!concentrationByKey.contains(it.key()) ||
                        concentrationByKey.value(it.key()).value != it.value().value ||
                        concentrationByKey.value(it.key()).units != it.value().units) {
                        concentrationsChanged = true; break;
                    }
                }
            }
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
            if (picker) {
                const auto selected = picker->value();
                values[solveField] = selected.isEmpty() ? QVariant{} : QVariant(selected);
                if (!adding && selected == old(solveField).toString()) values[solveField] = old(solveField);
            }
            QStringList changed;
            for (const auto& [field, value] : values)
                if (adding || value != old(field)) changed << field;
            if (changed.isEmpty() && !concentrationsChanged) { dialog.accept(); return; }
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
                else if (!changed.isEmpty()) {
                    for (const auto& field : columns) predicates << quote(field) + " IS ?";
                    sql = "UPDATE raw_profile SET " + assignments.join(',') + " WHERE " + predicates.join(" AND ");
                }
                if (!sql.isEmpty()) {
                    auto statement = prepare(db.get(), sql); int index = 1;
                    for (const auto& field : changed) bind(db.get(), statement.get(), index++, values.at(field));
                    if (!adding) for (const auto& value : original) bind(db.get(), statement.get(), index++, value);
                    check(db.get(), sqlite3_step(statement.get()), SQLITE_DONE);
                    if (sqlite3_changes(db.get()) != 1) throw std::runtime_error("The profile changed in the database. Cancel and Refresh before modifying.");
                }
                auto removeOldConcentrations = prepare(db.get(), "DELETE FROM raw_profile_entrance_concentrations WHERE NAME IS ? AND WT_PERCENT IS ?");
                bind(db.get(), removeOldConcentrations.get(), 1, adding ? values.at("NAME") : old("NAME"));
                bind(db.get(), removeOldConcentrations.get(), 2, adding ? values.at("WT_PERCENT") : old("WT_PERCENT"));
                check(db.get(), sqlite3_step(removeOldConcentrations.get()), SQLITE_DONE);
                removeOldConcentrations.reset();
                auto insertConcentration = prepare(db.get(), "INSERT INTO raw_profile_entrance_concentrations(NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS) VALUES(?, ?, ?, ?, ?, ?)");
                for (const auto& concentration : editedConcentrations) {
                    sqlite3_reset(insertConcentration.get()); sqlite3_clear_bindings(insertConcentration.get());
                    bind(db.get(), insertConcentration.get(), 1, values.at("NAME"));
                    bind(db.get(), insertConcentration.get(), 2, values.at("WT_PERCENT"));
                    bind(db.get(), insertConcentration.get(), 3, static_cast<qlonglong>(concentration.entranceNumber));
                    bind(db.get(), insertConcentration.get(), 4, concentration.species);
                    bind(db.get(), insertConcentration.get(), 5, concentration.value);
                    bind(db.get(), insertConcentration.get(), 6, concentration.units);
                    check(db.get(), sqlite3_step(insertConcentration.get()), SQLITE_DONE);
                }
                insertConcentration.reset();
                check(db.get(), sqlite3_exec(db.get(), "COMMIT", nullptr, nullptr, nullptr));
            } catch (...) { sqlite3_exec(db.get(), "ROLLBACK", nullptr, nullptr, nullptr); throw; }
            wrote = true; dialog.accept();
        } catch (const std::exception& failure) { error->setText(QString::fromUtf8(failure.what())); }
    });
    dialog.exec();
    return wrote;
}
