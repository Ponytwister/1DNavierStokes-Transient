#pragma once
#include <QStringList>
#include <QVariant>
#include <vector>
class QWidget;
// An empty original row means Add. Returns true only after a committed write.
bool editRawProfile(const QString& database, const QStringList& columns,
                    const std::vector<QVariant>& original, QWidget* parent);
