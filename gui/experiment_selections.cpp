#include "experiment_selections.h"
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <memory>

namespace {
QMap<QString, QString> readPairs(sqlite3* db, const char* sql) {
    sqlite3_stmt* raw = nullptr;
    const int prepared = sqlite3_prepare_v2(db, sql, -1, &raw, nullptr);
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(raw, sqlite3_finalize);
    if (prepared != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
    QMap<QString, QString> result;
    int rc;
    while ((rc = sqlite3_step(raw)) == SQLITE_ROW) {
        auto text = [&](int i) { const auto* s = sqlite3_column_text(raw, i);
            return s ? QString::fromUtf8(reinterpret_cast<const char*>(s)) : QString{}; };
        if (result.contains(text(0))) throw std::runtime_error("Ambiguous species or reaction name: " + text(0).toStdString());
        result.insert(text(0), text(1));
    }
    if (rc != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
    return result;
}
QStringList supportedUnits(const QString& type) {
    QStringList result{"umol", "mg/ml", "wt%", "g/ml"};
    if (type == "particle") result << "um2/ul" << "nm2/ul" << "mm2/nl";
    return result;
}
void setSelected(ChecklistPicker* picker, const QString& name, bool selected) {
    auto* list = picker->findChild<QListWidget*>();
    for (int i = 0; i < list->count(); ++i)
        if (list->item(i)->text() == name) list->item(i)->setCheckState(selected ? Qt::Checked : Qt::Unchecked);
}
}

ExperimentSelections::ExperimentSelections(const QString& database, const QMap<QString, QVariant>& initial, QWidget* parent)
    : QObject(parent) {
    sqlite3* raw = nullptr;
    const int rc = sqlite3_open_v2(database.toUtf8().constData(), &raw, SQLITE_OPEN_READONLY, nullptr);
    std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)> db(raw, sqlite3_close_v2);
    if (rc != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(raw));
    types_ = readPairs(raw, "SELECT SPECIES_NAME, SPECIES_TYPE FROM species ORDER BY SPECIES_NAME");
    reactionSpecies_ = readPairs(raw, "SELECT REACTION_NAME, SPECIES FROM reactions ORDER BY REACTION_NAME");
    species_ = new ChecklistPicker(types_.keys(), initial.value("SPECIES").toString(), "Species", "experimentSpeciesChoices", "None");
    species_->setParent(parent);
    reactions_ = new ChecklistPicker(reactionSpecies_.keys(), initial.value("REACTIONS").toString(), "Reactions", "experimentReactionChoices", "None");
    reactions_->setParent(parent);
    fields_["SPECIES"] = species_; fields_["REACTIONS"] = reactions_;
    // Retain even unavailable saved species so opening the form never drops their units.
    const auto names = initial.value("SPECIES").toString().split(' ', Qt::SkipEmptyParts);
    for (const auto& field : {QString("SPECIE_INLET_CONC_UNITS"), QString("SPECIE_MODEL_CONC_UNITS")}) {
        units_[field] = {};
        auto* panel = new QWidget(parent); new QFormLayout(panel); fields_[field] = panel;
        const auto saved = initial.value(field).toString().split(' ', Qt::SkipEmptyParts);
        if (saved.size() > names.size()) throw std::runtime_error("Concentration unit list has more entries than Species. Correct the database list before editing.");
        for (int i = 0; i < names.size(); ++i) {
            auto* combo = new QComboBox(panel); combo->addItem("Select units", "");
            for (const auto& unit : supportedUnits(types_.value(names[i]))) combo->addItem(unit, unit);
            const auto unit = saved.value(i);
            int index = combo->findData(unit);
            if (index < 0) { combo->addItem(unit + " (unsupported)", unit); index = combo->count() - 1; }
            combo->setCurrentIndex(index); units_[field][names[i]] = combo;
        }
    }
    for (auto it = fields_.begin(); it != fields_.end(); ++it) {
        it.value()->setParent(parent); it.value()->setObjectName(it.key());
    }
    refreshUnits();
    connect(species_->findChild<QListWidget*>(), &QListWidget::itemChanged, this, [this] { refreshUnits(); });
}

QStringList ExperimentSelections::species() const {
    return species_->value().split(' ', Qt::SkipEmptyParts);
}

void ExperimentSelections::refreshUnits() {
    // Invalid saved names remain visible and are rejected on Save by ChecklistPicker.
    QStringList names;
    try { names = species(); } catch (...) {
        auto* list = species_->findChild<QListWidget*>();
        for (int i = 0; i < list->count(); ++i)
            if (list->item(i)->checkState() == Qt::Checked) names << list->item(i)->text();
    }
    for (auto it = units_.begin(); it != units_.end(); ++it) {
        auto* panel = fields_[it.key()]; auto* form = static_cast<QFormLayout*>(panel->layout());
        while (form->rowCount()) {
            auto row = form->takeRow(0);
            if (row.labelItem) { delete row.labelItem->widget(); delete row.labelItem; }
            if (row.fieldItem) { row.fieldItem->widget()->hide(); delete row.fieldItem; }
        }
        for (const auto& name : names) {
            auto*& combo = it.value()[name];
            if (!combo) {
                combo = new QComboBox(panel); combo->addItem("Select units", "");
                for (const auto& unit : supportedUnits(types_.value(name))) combo->addItem(unit, unit);
            }
            combo->setObjectName(it.key() + ":" + name);
            form->addRow(name, combo); combo->show();
        }
    }
}

QString ExperimentSelections::value(const QString& name) const {
    if (name == "SPECIES") return species_->value();
    if (name == "REACTIONS") return reactions_->value();
    QStringList result;
    for (const auto& specie : species()) {
        if (types_.value(specie) != "particle" && types_.value(specie) != "molecule")
            throw std::runtime_error(("Unsupported species type for " + specie).toStdString());
        const auto unit = units_.value(name).value(specie)->currentData().toString();
        if (!supportedUnits(types_.value(specie)).contains(unit))
            throw std::runtime_error(("Choose supported concentration units for " + specie + ".").toStdString());
        result << unit;
    }
    return result.join(' ');
}

bool ExperimentSelections::resolveReactions(QWidget* parent) {
    const auto selected = species();
    QStringList missing, offending;
    for (const auto& reaction : reactions_->value().split(' ', Qt::SkipEmptyParts)) {
        for (const auto& specie : reactionSpecies_[reaction].split(' ', Qt::SkipEmptyParts)) {
            if (!selected.contains(specie)) {
                if (!missing.contains(specie)) missing << specie;
                if (!offending.contains(reaction)) offending << reaction;
            }
        }
    }
    if (missing.isEmpty()) return true;
    QMessageBox prompt(QMessageBox::Warning, "Reaction species missing",
        "Selected reactions require these species: " + missing.join(", ") + ".\nAdd them and choose their units, remove the affected reactions, or cancel to continue editing.",
        QMessageBox::Cancel, parent);
    prompt.setObjectName("reactionSpeciesWarning");
    auto* add = prompt.addButton("Add species", QMessageBox::AcceptRole);
    auto* remove = prompt.addButton("Remove reactions", QMessageBox::DestructiveRole);
    for (const auto& name : missing) if (!types_.contains(name)) add->setEnabled(false);
    prompt.setDefaultButton(QMessageBox::Cancel); prompt.exec();
    if (prompt.clickedButton() == add) for (const auto& name : missing) setSelected(species_, name, true);
    if (prompt.clickedButton() == remove) for (const auto& name : offending) setSelected(reactions_, name, false);
    return false; // Review changed selections before saving; no partial database writes.
}

void ExperimentSelections::validateReferences(sqlite3* db) const {
    const auto types = readPairs(db, "SELECT SPECIES_NAME, SPECIES_TYPE FROM species");
    const auto reactions = readPairs(db, "SELECT REACTION_NAME, SPECIES FROM reactions");
    const auto selected = species();
    for (const auto& name : selected)
        if (!types.contains(name) || types[name] != types_.value(name))
            throw std::runtime_error("Species changed. Cancel and reopen the editor.");
    for (const auto& name : reactions_->value().split(' ', Qt::SkipEmptyParts)) {
        if (!reactions.contains(name) || reactions[name] != reactionSpecies_.value(name))
            throw std::runtime_error("Reaction changed. Cancel and reopen the editor.");
        for (const auto& specie : reactions[name].split(' ', Qt::SkipEmptyParts))
            if (!selected.contains(specie)) throw std::runtime_error("Reaction species are missing from Species.");
    }
}
