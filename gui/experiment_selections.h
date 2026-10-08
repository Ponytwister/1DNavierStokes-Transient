#pragma once
#include "checklist_picker.h"
#include <QComboBox>
#include <QFormLayout>
#include <QMap>
#include <QVariant>
#include <sqlite3.h>

// Owns the experiment's ordered selections and one concentration unit per species.
class ExperimentSelections : public QObject {
public:
    ExperimentSelections(const QString& database, const QMap<QString, QVariant>& initial, QWidget* parent);
    QWidget* field(const QString& name) const { return fields_.value(name); }
    QString value(const QString& name) const;
    bool resolveReactions(QWidget* parent);
    void validateReferences(sqlite3* db) const;
private:
    void refreshUnits();
    QStringList species() const;
    QMap<QString, QWidget*> fields_;
    QMap<QString, QString> types_, reactionSpecies_;
    QMap<QString, QMap<QString, QComboBox*>> units_;
    ChecklistPicker *species_, *reactions_;
    ChecklistPicker* parameters_ = nullptr;
};
