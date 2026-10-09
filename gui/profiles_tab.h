#pragma once

#include <QWidget>
#include <QVector>
#include <QStringList>

class QComboBox;
class QLabel;
class QPushButton;

class ProfilePlot : public QWidget {
public:
    explicit ProfilePlot(QWidget* parent = nullptr);
    void setData(const QString& title, const QString& xLabel, const QString& yLabel,
                 const QVector<double>& x, const QVector<double>& y);
protected:
    void paintEvent(QPaintEvent* event) override;
private:
    QString title_, xLabel_, yLabel_;
    QVector<double> x_, y_;
};

class ProfilesTab : public QWidget {
public:
    explicit ProfilesTab(QWidget* parent = nullptr);
    void setDatabase(const QString& filename);
private:
    void loadDatabase();
    void loadReport(const QString& filename);
    void updatePlot();
    QComboBox *source_, *block_, *profile_;
    QLabel* status_;
    QPushButton* openReport_;
    ProfilePlot* plot_;
    QString database_, reportFile_;
    QStringList blocks_, profiles_;
};
