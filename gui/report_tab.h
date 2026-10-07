#pragma once
#include <QWidget>
#include <QStringList>
#include <functional>
#include <QMap>

class QLabel;
class QPushButton;
class QTableView;
class QListWidget;
class QCheckBox;
class QComboBox;
class QSpinBox;

// Keep the export's textual numbers, repeated headers and blank rows intact.
namespace report_format {
using Rows = QList<QStringList>;
Rows parse(const QString& text);
QString tsv(const Rows& rows);
QString profileType(const QString& label);
QStringList metricNames();
enum class MetadataLayout { stacked, separateBlock };
struct Options {
    QList<int> blocks;
    QStringList types;
    QList<int> columns;
    bool samples = true, headers = true, units = true, blankRows = true;
    char numberFormat = 0; // zero preserves the source text
    int decimals = 6;
    MetadataLayout metadataLayout = MetadataLayout::stacked;
    QList<int> metrics = {0, 1, 2, 3, 4, 5};
};
QStringList blockNames(const Rows& rows);
Rows select(const Rows& rows, const Options& options);
}

class ReportTab : public QWidget {
public:
    explicit ReportTab(QWidget* parent = nullptr);
    bool loadFile(const QString& filename);
    bool saveFile(const QString& filename);
    bool loadGenerated(const QString& text, const QString& source, const QString& suggestedFile);
    void setRunAvailable(bool available);
    void clearGenerated();
    std::function<void()> generateRequested;
private:
    bool loadSource(const QString& text, bool preserveSelections);
    void clearSource();
    void updatePreview();
    void rebuildColumns();
    QLabel* status_;
    QTableView* table_;
    QPushButton *copy_, *save_;
    QPushButton* generate_;
    report_format::Rows rows_;
    report_format::Rows selected_;
    QListWidget *blocks_, *types_, *columns_;
    QCheckBox *samples_, *headers_, *units_, *blankRows_;
    QComboBox* numberFormat_;
    QComboBox* metadataLayout_;
    QMap<int, Qt::CheckState> columnStates_;
    QSpinBox* decimals_;
    QString filename_;
    QString sourceLabel_;
    bool generated_ = false;
};
