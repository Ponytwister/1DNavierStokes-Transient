#pragma once
#include <QWidget>
#include <QStringList>
#include <functional>
#include <QMap>
#include <QJsonObject>

class QLabel;
class QPushButton;
class QTableView;
class QListWidget;
class QCheckBox;
class QComboBox;
class QSpinBox;
class QTableWidget;

// Keep the export's textual numbers, repeated headers and blank rows intact.
namespace report_format {
using Rows = QList<QStringList>;
Rows parse(const QString& text);
QString tsv(const Rows& rows);
QString profileType(const QString& label);
QStringList metricNames();
enum class MetadataLayout { stacked, separateBlock, mapped };
struct Options {
    QList<int> blocks;
    QStringList types;
    QList<int> columns;
    // Stable IDs 0-8 are source fields; 9-14 are named metadata values.
    QList<int> orderedColumns;
    // Metadata-column order, including unchecked fields. This controls the
    // stacked metadata units independently of separate-column visibility.
    QList<int> metadataFieldOrder;
    bool reorderRows = false; // Use block/type list order, retaining each block's axes.
    bool samples = true, headers = true, units = true, blankRows = true;
    bool legacyStyle = false;
    char numberFormat = 0; // zero preserves the source text
    int decimals = 6;
    MetadataLayout metadataLayout = MetadataLayout::stacked;
    QList<int> metrics = {0, 1, 2, 3, 4, 5};
    // Profile type -> stable metadata field ID (0-14; 7 is the stacked field).
    QMap<QString, int> profileMetadata;
};
QStringList blockNames(const Rows& rows);
Rows select(const Rows& rows, const Options& options);
QJsonObject legacyFormat();
}

class ReportTab : public QWidget {
public:
    explicit ReportTab(QWidget* parent = nullptr);
    bool loadFile(const QString& filename);
    bool saveFile(const QString& filename);
    bool loadGenerated(const QString& text, const QString& source, const QString& suggestedFile);
    bool saveFormatFile(const QString& filename);
    bool loadFormatFile(const QString& filename);
    QString formatFilename() const { return formatFilename_; }
    QJsonObject formatSettings() const;
    void setFormatSettings(const QJsonObject& settings);
    void applyLegacyFormat();
    void setRunAvailable(bool available);
    void clearGenerated();
    std::function<void()> generateRequested;
private:
    bool loadSource(const QString& text, bool preserveSelections);
    void applyFormatToControls(const QJsonObject& settings);
    void clearSource();
    void updatePreview();
    void rebuildColumns();
    void rebuildMetricTargets(const QMap<QString, int>& savedMappings = {});
    void updateMetricTargetStyles();
    QLabel* status_;
    QTableView* table_;
    QPushButton *copy_, *save_;
    QPushButton* generate_;
    report_format::Rows rows_;
    report_format::Rows selected_;
    QListWidget *blocks_, *types_, *columns_;
    QCheckBox *samples_, *headers_, *units_, *blankRows_;
    QComboBox* numberFormat_;
    QTableWidget* metricTargets_;
    QMap<int, Qt::CheckState> columnStates_;
    QList<int> columnOrder_;
    QSpinBox* decimals_;
    QString filename_;
    QString sourceLabel_;
    QString formatFilename_;
    QJsonObject pendingFormat_;
    bool generated_ = false;
    bool rowsReordered_ = false;
    bool legacyStyle_ = false;
};
