#pragma once

#include <QWidget>
#include <QVector>
#include <QStringList>
#include <QJsonObject>
#include <QColor>
#include <feedback.h>
#include <model_controls.h>
#include <functional>
#include <map>
#include <vector>

class QComboBox;
class QCheckBox;
class QDoubleSpinBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;

struct ProfileSeriesData {
    QString name;
    QVector<double> x;
    QVector<double> y;
    QColor lineColor;
    QColor fillColor;
    bool filled = false;
    int yAxis = 0;
};

struct SingleProfileRunRequest {
    QString experiment;
    QString dataset;
    QStringList variables;
    std::map<std::string, double> initialValues;
    bool runSolver = true;
};

class ProfilePlot : public QWidget {
public:
    explicit ProfilePlot(QWidget* parent = nullptr);
    void setData(const QString& title, const QString& xLabel, const QString& yLabel,
                 const QVector<ProfileSeriesData>& series);
    void setAppearance(bool grid, bool ticks, bool autoRangeX, bool autoRangeY, bool useSecondYAxis,
                       bool autoRangeY2, double xmin, double xmax, double ymin, double ymax, double y2min, double y2max,
                       int xTickCount, int yTickCount, int y2TickCount, const QString& xLabel, const QString& yLabel, const QString& y2Label,
                       const QString& title, double aspectRatio, int legendPosition,
                       int xNumberGap, int xTitleGap, int yNumberGap, int yTitleGap, int y2NumberGap, int y2TitleGap);
    QString displayTitle() const { return title_; }
protected:
    void paintEvent(QPaintEvent* event) override;
private:
    QString title_, xLabel_, yLabel_, y2Label_;
    QVector<ProfileSeriesData> series_;
    bool grid_ = true, ticks_ = true, autoRangeX_ = true, autoRangeY_ = true;
    bool useSecondYAxis_ = false, autoRangeY2_ = true;
    double xmin_ = 0, xmax_ = 1, ymin_ = 0, ymax_ = 1;
    double y2min_ = 0, y2max_ = 1;
    double aspectRatio_ = 1.6;
    int xTickCount_ = 5, yTickCount_ = 5, y2TickCount_ = 5;
    int legendPosition_ = 1, xNumberGap_ = 5, xTitleGap_ = 6, yNumberGap_ = 5, yTitleGap_ = 8, y2NumberGap_ = 5, y2TitleGap_ = 8;
};

class ProfilesTab : public QWidget {
public:
    explicit ProfilesTab(QWidget* parent = nullptr);
    void setDatabase(const QString& filename);
    void setSingleRunContext(const QString& filename, const model_controls::Snapshot& controls);
    void setSingleRunEnabled(bool enabled);
    std::function<void(const SingleProfileRunRequest&)> singleRunRequested;
    void beginLiveRun();
    void setLiveProfiles(std::vector<tsensor_workflow::live_profile_snapshot> profiles);
    QJsonObject settings() const;
    void setSettings(const QJsonObject& settings);
private:
    struct SeriesRow {
        QComboBox* dataset = nullptr;
        QComboBox* profile = nullptr;
        QComboBox* yAxis = nullptr;
        QPushButton* lineColorButton = nullptr;
        QPushButton* fillColorButton = nullptr;
        QCheckBox* fill = nullptr;
        QLineEdit* legend = nullptr;
        QColor lineColor, fillColor;
        QString lastProfile;
    };
    void loadDatabase();
    void loadReport(const QString& filename);
    void showLiveProfiles();
    void refreshSingleRunDatasets(const QString& preferredDataset = {});
    void refreshSingleRunVariables();
    void requestSingleRun();
    void updatePlot();
    void addSeriesRow(const QJsonObject& settings = {});
    void removeSeriesRow(SeriesRow* row);
    void refreshSeriesRows();
    void chooseLineColor(SeriesRow* row);
    void chooseFillColor(SeriesRow* row);
    QComboBox* source_;
    QTableWidget* seriesTable_;
    QPushButton* addSeriesButton_;
    QGroupBox* singleRunGroup_;
    QComboBox *singleRunExperiment_, *singleRunDataset_;
    QTableWidget* singleRunVariables_;
    QCheckBox* singleRunSolver_;
    QPushButton* singleRunButton_;
    QLabel* singleRunStatus_;
    model_controls::Snapshot singleRunControls_;
    QString singleRunDatabase_;
    QLabel* status_;
    QPushButton* openReport_;
    QPushButton *exportImage_;
    QLineEdit *title_, *xLabel_, *yLabel_, *y2Label_;
    QDoubleSpinBox *xMin_, *xMax_, *yMin_, *yMax_, *y2Min_, *y2Max_, *xNumberGap_, *xTitleGap_, *yNumberGap_, *yTitleGap_, *y2NumberGap_, *y2TitleGap_;
    QDoubleSpinBox* aspectRatio_;
    QSpinBox *xTicks_, *yTicks_, *y2Ticks_;
    QComboBox* legendPosition_;
    QCheckBox *autoRangeX_, *autoRangeY_, *autoRangeY2_, *useSecondYAxis_, *grid_, *ticks_;
    ProfilePlot* plot_;
    QString database_, reportFile_;
    QStringList datasetLabels_, datasetIds_, profiles_;
    QVector<SeriesRow*> seriesRows_;
    std::vector<tsensor_workflow::live_profile_snapshot> liveProfiles_;
};
