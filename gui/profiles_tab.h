#pragma once

#include <QWidget>
#include <QVector>
#include <QStringList>
#include <QJsonObject>
#include <QColor>

class QComboBox;
class QCheckBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;

class ProfilePlot : public QWidget {
public:
    explicit ProfilePlot(QWidget* parent = nullptr);
    void setData(const QString& title, const QString& xLabel, const QString& yLabel,
                 const QVector<double>& x, const QVector<double>& y);
    void setAppearance(const QColor& line, const QColor& fill, bool filled, bool grid,
                       bool ticks, bool autoRangeX, bool autoRangeY,
                       double xmin, double xmax, double ymin, double ymax,
                       int xTickCount, int yTickCount, const QString& xLabel, const QString& yLabel,
                       const QString& title, double aspectRatio);
    QString displayTitle() const { return title_; }
protected:
    void paintEvent(QPaintEvent* event) override;
private:
    QString title_, xLabel_, yLabel_;
    QVector<double> x_, y_;
    QColor lineColor_, fillColor_;
    bool filled_ = false, grid_ = true, ticks_ = true, autoRangeX_ = true, autoRangeY_ = true;
    double xmin_ = 0, xmax_ = 1, ymin_ = 0, ymax_ = 1;
    double aspectRatio_ = 1.6;
    int xTickCount_ = 5, yTickCount_ = 5;
};

class ProfilesTab : public QWidget {
public:
    explicit ProfilesTab(QWidget* parent = nullptr);
    void setDatabase(const QString& filename);
    QJsonObject settings() const;
    void setSettings(const QJsonObject& settings);
private:
    void loadDatabase();
    void loadReport(const QString& filename);
    void updatePlot();
    void chooseLineColor();
    void chooseFillColor();
    QComboBox *source_, *block_, *profile_;
    QLabel* status_;
    QPushButton* openReport_;
    QPushButton *exportImage_, *lineColorButton_, *fillColorButton_;
    QLineEdit *title_, *xLabel_, *yLabel_;
    QDoubleSpinBox *xMin_, *xMax_, *yMin_, *yMax_;
    QDoubleSpinBox* aspectRatio_;
    QSpinBox *xTicks_, *yTicks_;
    QCheckBox *autoRangeX_, *autoRangeY_, *grid_, *ticks_, *fill_;
    ProfilePlot* plot_;
    QString database_, reportFile_;
    QStringList blocks_, profiles_;
    QColor lineColor_{37, 99, 235}, fillColor_{37, 99, 235, 72};
};
