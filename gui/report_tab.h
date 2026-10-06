#pragma once
#include <QWidget>
#include <QStringList>

class QLabel;
class QPushButton;
class QTableView;

// Keep the export's textual numbers, repeated headers and blank rows intact.
namespace report_format {
using Rows = QList<QStringList>;
Rows parse(const QString& text);
QString tsv(const Rows& rows);
}

class ReportTab : public QWidget {
public:
    explicit ReportTab(QWidget* parent = nullptr);
    bool loadFile(const QString& filename);
private:
    QLabel* status_;
    QTableView* table_;
    QPushButton *copy_, *save_;
    report_format::Rows rows_;
    QString filename_;
};
