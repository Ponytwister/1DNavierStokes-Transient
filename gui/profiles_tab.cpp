#include "profiles_tab.h"

#include "report_tab.h"
#include <QComboBox>
#include <QFileDialog>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QFile>
#include <QPainterPath>
#include <QColorDialog>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QLineEdit>
#include <QSpinBox>
#include <QTableWidget>
#include <QHeaderView>
#include <QGroupBox>
#include <QGridLayout>
#include <QScrollArea>
#include <QFontMetrics>
#include <QFileInfo>
#include <QMessageBox>
#include <QDateTime>
#include <QTimeZone>
#include <QRegularExpression>
#include <QImage>
#include <QJsonArray>
#include <QHash>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <sqlite3.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <iterator>

ProfilePlot::ProfilePlot(QWidget* parent) : QWidget(parent) { setMinimumSize(480, 320); }

void ProfilePlot::setData(const QString& title, const QString& xLabel, const QString& yLabel,
                          const QVector<ProfileSeriesData>& series)
{
    title_ = title; xLabel_ = xLabel; yLabel_ = yLabel; series_ = series; update();
}

void ProfilePlot::setAppearance(bool grid, bool ticks, bool autoRangeX, bool autoRangeY, bool useSecondYAxis,
                                bool autoRangeY2, double xmin, double xmax, double ymin, double ymax, double y2min, double y2max,
                                int xTickCount, int yTickCount, int y2TickCount, const QString& xLabel, const QString& yLabel, const QString& y2Label,
                                const QString& title, double aspectRatio, int legendPosition,
                                int xNumberGap, int xTitleGap, int yNumberGap, int yTitleGap, int y2NumberGap, int y2TitleGap)
{
    grid_ = grid; ticks_ = ticks;
    autoRangeX_ = autoRangeX; autoRangeY_ = autoRangeY;
    useSecondYAxis_ = useSecondYAxis; autoRangeY2_ = autoRangeY2;
    xmin_ = xmin; xmax_ = xmax; ymin_ = ymin; ymax_ = ymax;
    y2min_ = y2min; y2max_ = y2max;
    xTickCount_ = xTickCount; yTickCount_ = yTickCount; y2TickCount_ = y2TickCount; aspectRatio_ = aspectRatio;
    legendPosition_ = legendPosition;
    xNumberGap_ = xNumberGap; xTitleGap_ = xTitleGap; yNumberGap_ = yNumberGap; yTitleGap_ = yTitleGap; y2NumberGap_ = y2NumberGap; y2TitleGap_ = y2TitleGap;
    if (!xLabel.isEmpty()) xLabel_ = xLabel;
    if (!yLabel.isEmpty()) yLabel_ = yLabel;
    if (!y2Label.isEmpty()) y2Label_ = y2Label;
    if (!title.isEmpty()) title_ = title;
    update();
}

void ProfilePlot::paintEvent(QPaintEvent*)
{
    QPainter painter(this); painter.fillRect(rect(), palette().base()); painter.setRenderHint(QPainter::Antialiasing);
    const bool hasLegend = !series_.isEmpty();
    const qreal legendSpace = hasLegend ? 42 : 0;
    const qreal leftMargin = 104 + yNumberGap_ + yTitleGap_;
    const qreal rightMargin = useSecondYAxis_ ? 104 + y2NumberGap_ + y2TitleGap_ : 32;
    const qreal bottomMargin = 76 + xNumberGap_ + xTitleGap_ + legendSpace;
    QRectF area(leftMargin, 42, width() - leftMargin - rightMargin, height() - 42 - bottomMargin);
    if (hasLegend && legendPosition_ == 2) { area.adjust(100, 0, 0, 0); }
    if (hasLegend && legendPosition_ == 3) { area.adjust(0, 0, -100, 0); }
    if (hasLegend && legendPosition_ == 0) { area.adjust(0, 34, 0, 0); }
    if (aspectRatio_ > 0 && area.width() > 0 && area.height() > 0) {
        const double availableRatio = area.width() / area.height();
        if (availableRatio > aspectRatio_) {
            const double plotWidth = area.height() * aspectRatio_;
            area.setLeft(area.center().x() - plotWidth / 2);
            area.setWidth(plotWidth);
        } else {
            const double plotHeight = area.width() / aspectRatio_;
            area.setTop(area.center().y() - plotHeight / 2);
            area.setHeight(plotHeight);
        }
    }
    painter.setPen(palette().text().color()); painter.drawText(QRectF(8, 8, width() - 16, 24), Qt::AlignCenter, title_);
    if (series_.isEmpty() || area.width() <= 0 || area.height() <= 0) {
        painter.drawText(area, Qt::AlignCenter, "Select profile data to plot"); return;
    }
    double dataMinX = std::numeric_limits<double>::infinity();
    double dataMaxX = -std::numeric_limits<double>::infinity();
    double dataMinY = std::numeric_limits<double>::infinity(), dataMaxY = -std::numeric_limits<double>::infinity();
    double dataMinY2 = std::numeric_limits<double>::infinity(), dataMaxY2 = -std::numeric_limits<double>::infinity();
    for (const auto& series : series_) {
        const auto count = std::min(series.x.size(), series.y.size());
        for (qsizetype i = 0; i < count; ++i) {
            if (!std::isfinite(series.x[i]) || !std::isfinite(series.y[i])) continue;
            dataMinX = std::min(dataMinX, series.x[i]); dataMaxX = std::max(dataMaxX, series.x[i]);
            if (useSecondYAxis_ && series.yAxis == 1) { dataMinY2 = std::min(dataMinY2, series.y[i]); dataMaxY2 = std::max(dataMaxY2, series.y[i]); }
            else { dataMinY = std::min(dataMinY, series.y[i]); dataMaxY = std::max(dataMaxY, series.y[i]); }
        }
    }
    if (!std::isfinite(dataMinX)) {
        painter.drawText(area, Qt::AlignCenter, "Select profile data to plot"); return;
    }
    if (!std::isfinite(dataMinY)) { dataMinY = 0; dataMaxY = 1; }
    if (!std::isfinite(dataMinY2)) { dataMinY2 = 0; dataMaxY2 = 1; }
    double lowX = autoRangeX_ ? dataMinX : xmin_, highX = autoRangeX_ ? dataMaxX : xmax_;
    double lowY = autoRangeY_ ? dataMinY : ymin_, highY = autoRangeY_ ? dataMaxY : ymax_;
    double lowY2 = autoRangeY2_ ? dataMinY2 : y2min_, highY2 = autoRangeY2_ ? dataMaxY2 : y2max_;
    double dx = highX - lowX, dy = highY - lowY, dy2 = highY2 - lowY2;
    if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dy2)) return;
    if (dx == 0) dx = 1; if (dy == 0) dy = std::max(1.0, std::abs(lowY) * 0.1); if (dy2 == 0) dy2 = std::max(1.0, std::abs(lowY2) * 0.1);
    const double x0 = autoRangeX_ ? lowX - dx * .02 : lowX, x1 = autoRangeX_ ? highX + dx * .02 : highX;
    const double y0 = autoRangeY_ ? lowY - dy * .08 : lowY, y1 = autoRangeY_ ? highY + dy * .08 : highY;
    const double y20 = autoRangeY2_ ? lowY2 - dy2 * .08 : lowY2, y21 = autoRangeY2_ ? highY2 + dy2 * .08 : highY2;
    if (x1 <= x0 || y1 <= y0 || y21 <= y20) return;
    painter.setPen(QPen(palette().mid().color(), 1)); painter.drawRect(area);
    painter.setPen(palette().text().color());
    if (grid_) {
        painter.setPen(QPen(palette().mid().color().lighter(130), 1, Qt::DotLine));
        for (int i = 1; i < yTickCount_; ++i) { const qreal y = area.bottom() - area.height() * i / yTickCount_; painter.drawLine(QPointF(area.left(), y), QPointF(area.right(), y)); }
        for (int i = 1; i < xTickCount_; ++i) { const qreal x = area.left() + area.width() * i / xTickCount_; painter.drawLine(QPointF(x, area.top()), QPointF(x, area.bottom())); }
    }
    painter.setPen(palette().text().color());
    for (int i = 0; i <= yTickCount_; ++i) {
        const qreal y = area.bottom() - area.height() * i / yTickCount_;
        const double value = y0 + (y1 - y0) * i / yTickCount_;
        if (ticks_) painter.drawLine(QPointF(area.left() - 4, y), QPointF(area.left(), y));
        painter.drawText(QRectF(area.left() - 68 - yNumberGap_, y - 10, 62, 20), Qt::AlignRight | Qt::AlignVCenter, QString::number(value, 'g', 4));
    }
    if (useSecondYAxis_) {
        for (int i = 0; i <= y2TickCount_; ++i) {
            const qreal y = area.bottom() - area.height() * i / y2TickCount_;
            const double value = y20 + (y21 - y20) * i / y2TickCount_;
            if (ticks_) painter.drawLine(QPointF(area.right(), y), QPointF(area.right() + 4, y));
            painter.drawText(QRectF(area.right() + 6 + y2NumberGap_, y - 10, 62, 20), Qt::AlignLeft | Qt::AlignVCenter, QString::number(value, 'g', 4));
        }
    }
    for (int i = 0; i <= xTickCount_; ++i) {
        const qreal x = area.left() + area.width() * i / xTickCount_;
        const double value = x0 + (x1 - x0) * i / xTickCount_;
        if (ticks_) painter.drawLine(QPointF(x, area.bottom()), QPointF(x, area.bottom() + 4));
        painter.drawText(QRectF(x - 35, area.bottom() + 4 + xNumberGap_, 70, 16), Qt::AlignHCenter | Qt::AlignTop, QString::number(value, 'g', 4));
    }
    painter.drawText(QRectF(area.left(), area.bottom() + 20 + xNumberGap_ + xTitleGap_, area.width(), 20), Qt::AlignCenter, xLabel_);
    painter.save(); painter.translate(area.left() - 68 - yNumberGap_ - yTitleGap_ - 10, area.center().y()); painter.rotate(-90);
    painter.drawText(QRectF(-area.height()/2, -10, area.height(), 20), Qt::AlignCenter, yLabel_); painter.restore();
    if (useSecondYAxis_) {
        painter.save(); painter.translate(area.right() + 68 + y2NumberGap_ + y2TitleGap_ + 10, area.center().y()); painter.rotate(90);
        painter.drawText(QRectF(-area.height()/2, -10, area.height(), 20), Qt::AlignCenter, y2Label_); painter.restore();
    }
    painter.save(); painter.setClipRect(area);
    for (qsizetype seriesIndex = 0; seriesIndex < series_.size(); ++seriesIndex) {
        const auto& series = series_[seriesIndex];
        const bool onSecondAxis = useSecondYAxis_ && series.yAxis == 1;
        const double seriesY0 = onSecondAxis ? y20 : y0, seriesY1 = onSecondAxis ? y21 : y1;
        const double baseline = std::clamp(0.0, seriesY0, seriesY1);
        const qreal baseY = area.bottom() - (baseline - seriesY0) / (seriesY1 - seriesY0) * area.height();
        const auto count = std::min(series.x.size(), series.y.size());
        QPainterPath path;
        for (qsizetype i = 0; i < count; ++i) {
            if (!std::isfinite(series.x[i]) || !std::isfinite(series.y[i])) continue;
            const QPointF point(area.left() + (series.x[i] - x0) / (x1 - x0) * area.width(),
                                area.bottom() - (series.y[i] - seriesY0) / (seriesY1 - seriesY0) * area.height());
            if (path.elementCount() == 0) path.moveTo(point); else path.lineTo(point);
        }
        if (path.elementCount() < 2) continue;
        if (series.filled) {
            QPainterPath fillPath = path;
            const auto first = fillPath.elementAt(0), last = fillPath.elementAt(fillPath.elementCount() - 1);
            fillPath.lineTo(QPointF(last.x, baseY)); fillPath.lineTo(QPointF(first.x, baseY)); fillPath.closeSubpath();
            painter.fillPath(fillPath, series.fillColor);
        }
        painter.setBrush(Qt::NoBrush);
        painter.setPen(QPen(series.lineColor, 2));
        painter.drawPath(path);
    }
    painter.restore();

    if (hasLegend) {
        const QFontMetrics metrics(painter.font());
        qreal legendX = legendPosition_ == 2 ? 8 : area.left();
        qreal legendY = legendPosition_ == 0 ? 20 : area.bottom() + 46;
        const qreal rowHeight = metrics.height() + 2;
        for (qsizetype i = 0; i < series_.size(); ++i) {
            const qreal itemWidth = 20 + metrics.horizontalAdvance(series_[i].name) + 16;
            if (legendPosition_ == 2 || legendPosition_ == 3) {
                legendX = legendPosition_ == 2 ? 8 : width() - 92;
                legendY = 52 + i * rowHeight;
            } else if (legendX + itemWidth > area.right() && legendX > area.left()) {
                legendX = area.left(); legendY += rowHeight;
            }
            painter.setPen(QPen(series_[i].lineColor, 2));
            painter.drawLine(QPointF(legendX, legendY + metrics.ascent() / 2),
                             QPointF(legendX + 14, legendY + metrics.ascent() / 2));
            painter.setPen(palette().text().color());
            painter.drawText(QPointF(legendX + 19, legendY + metrics.ascent()), series_[i].name);
            legendX += itemWidth;
        }
    }
}

ProfilesTab::ProfilesTab(QWidget* parent) : QWidget(parent)
{
    auto* outer = new QVBoxLayout(this);
    auto* scrollArea = new QScrollArea(this); scrollArea->setWidgetResizable(true);
    scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    auto* content = new QWidget; auto* layout = new QVBoxLayout(content);
    outer->addWidget(scrollArea); scrollArea->setWidget(content);
    auto* controls = new QHBoxLayout;
    source_ = new QComboBox; source_->setObjectName("profilesSource"); source_->addItems({"Database", "Report file", "Live run"});
    openReport_ = new QPushButton("Open report..."); openReport_->setObjectName("openProfileReportButton");
    exportImage_ = new QPushButton("Export image..."); exportImage_->setObjectName("exportProfileImageButton");
    controls->addWidget(new QLabel("Source:")); controls->addWidget(source_); controls->addStretch(); controls->addWidget(openReport_); controls->addWidget(exportImage_); layout->addLayout(controls);
    singleRunGroup_ = new QGroupBox("Run one profile"); singleRunGroup_->setObjectName("singleProfileRunRegion");
    auto* singleRunLayout = new QGridLayout(singleRunGroup_);
    singleRunExperiment_ = new QComboBox; singleRunExperiment_->setObjectName("singleRunExperiment");
    singleRunDataset_ = new QComboBox; singleRunDataset_->setObjectName("singleRunDataset");
    singleRunSolver_ = new QCheckBox("Run solver (fit selected variables)"); singleRunSolver_->setObjectName("singleRunSolverToggle"); singleRunSolver_->setChecked(true);
    singleRunVariables_ = new QTableWidget(0, 3); singleRunVariables_->setObjectName("singleRunVariables");
    singleRunVariables_->setHorizontalHeaderLabels({"Variable", "Fit", "Starting value"});
    singleRunVariables_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    singleRunVariables_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    singleRunVariables_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
    singleRunVariables_->verticalHeader()->hide(); singleRunVariables_->setMinimumHeight(110); singleRunVariables_->setMaximumHeight(180);
    singleRunButton_ = new QPushButton("Run"); singleRunButton_->setObjectName("singleRunStartButton");
    singleRunPause_ = new QPushButton("Pause"); singleRunPause_->setObjectName("singleRunPauseButton");
    singleRunStop_ = new QPushButton("Stop"); singleRunStop_->setObjectName("singleRunStopButton");
    singleRunSaveInputs_ = new QPushButton("Save fitted inputs…"); singleRunSaveInputs_->setObjectName("singleRunSaveInputsButton");
    singleRunStatus_ = new QLabel("Choose an experiment and dataset to configure a profile run."); singleRunStatus_->setObjectName("singleRunStatus");
    singleRunLayout->addWidget(new QLabel("Experiment"), 0, 0); singleRunLayout->addWidget(singleRunExperiment_, 0, 1);
    singleRunLayout->addWidget(new QLabel("Dataset"), 0, 2); singleRunLayout->addWidget(singleRunDataset_, 0, 3);
    singleRunLayout->addWidget(singleRunSolver_, 1, 0, 1, 4);
    singleRunLayout->addWidget(singleRunVariables_, 2, 0, 1, 4);
    singleRunLayout->addWidget(singleRunStatus_, 3, 0, 1, 4);
    singleRunLayout->addWidget(singleRunButton_, 4, 0); singleRunLayout->addWidget(singleRunPause_, 4, 1);
    singleRunLayout->addWidget(singleRunStop_, 4, 2); singleRunLayout->addWidget(singleRunSaveInputs_, 4, 3);
    singleRunExperiment_->setEnabled(false); singleRunDataset_->setEnabled(false);
    singleRunSolver_->setEnabled(false); singleRunVariables_->setEnabled(false); singleRunButton_->setEnabled(false);
    singleRunPause_->setEnabled(false); singleRunStop_->setEnabled(false); singleRunSaveInputs_->setEnabled(false);
    auto spin = [](const char* name, double value) { auto* box = new QDoubleSpinBox; box->setObjectName(name); box->setDecimals(6); box->setRange(-1e12, 1e12); box->setSingleStep(.1); box->setValue(value); box->setKeyboardTracking(false); return box; };
    title_ = new QLineEdit; title_->setObjectName("profileGraphTitle"); title_->setPlaceholderText("Graph title (automatic)");
    xLabel_ = new QLineEdit; xLabel_->setObjectName("profileXLabel"); xLabel_->setPlaceholderText("X axis title (automatic)");
    yLabel_ = new QLineEdit; yLabel_->setObjectName("profileYLabel"); yLabel_->setPlaceholderText("Y axis title (automatic)");
    autoRangeX_ = new QCheckBox("Auto range"); autoRangeX_->setObjectName("profileAutoRangeX"); autoRangeX_->setChecked(true);
    autoRangeY_ = new QCheckBox("Auto range"); autoRangeY_->setObjectName("profileAutoRangeY"); autoRangeY_->setChecked(true);
    grid_ = new QCheckBox("Grid"); grid_->setObjectName("profileGrid"); grid_->setChecked(true);
    ticks_ = new QCheckBox("Tick marks"); ticks_->setObjectName("profileTicks"); ticks_->setChecked(true);
    xMin_ = spin("profileXMin", 0); xMax_ = spin("profileXMax", 1); yMin_ = spin("profileYMin", 0); yMax_ = spin("profileYMax", 1);
    xNumberGap_ = spin("profileXNumberGap", 5); xTitleGap_ = spin("profileXTitleGap", 6); yNumberGap_ = spin("profileYNumberGap", 5); yTitleGap_ = spin("profileYTitleGap", 8);
    for (auto* gap : {xNumberGap_, xTitleGap_, yNumberGap_, yTitleGap_}) { gap->setRange(0, 50); gap->setDecimals(0); gap->setSingleStep(1); }
    xTicks_ = new QSpinBox; xTicks_->setObjectName("profileXTickCount"); xTicks_->setRange(2, 12); xTicks_->setValue(5);
    yTicks_ = new QSpinBox; yTicks_->setObjectName("profileYTickCount"); yTicks_->setRange(2, 12); yTicks_->setValue(5);
    aspectRatio_ = spin("profileAspectRatio", 1.6); aspectRatio_->setRange(.25, 4); aspectRatio_->setDecimals(2); aspectRatio_->setSingleStep(.1);
    legendPosition_ = new QComboBox; legendPosition_->setObjectName("profileLegendPosition"); legendPosition_->addItems({"Top", "Bottom", "Left", "Right"}); legendPosition_->setCurrentIndex(1);

    auto* seriesGroup = new QGroupBox("Profiles to plot"); auto* seriesLayout = new QVBoxLayout(seriesGroup);
    seriesTable_ = new QTableWidget(0, 8); seriesTable_->setObjectName("profileSeriesTable");
    seriesTable_->setHorizontalHeaderLabels({"Dataset / saved version", "Profile", "Y axis", "Line", "Fill", "Fill under curve", "Legend name", ""});
    seriesTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch); seriesTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch); seriesTable_->horizontalHeader()->setSectionResizeMode(6, QHeaderView::Stretch); seriesTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents); seriesTable_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents); seriesTable_->horizontalHeader()->setSectionResizeMode(4, QHeaderView::ResizeToContents); seriesTable_->horizontalHeader()->setSectionResizeMode(5, QHeaderView::ResizeToContents); seriesTable_->horizontalHeader()->setSectionResizeMode(7, QHeaderView::ResizeToContents); seriesTable_->verticalHeader()->hide(); seriesTable_->setMinimumHeight(116);
    addSeriesButton_ = new QPushButton("Add profile"); addSeriesButton_->setObjectName("addProfileSeriesButton"); seriesLayout->addWidget(seriesTable_); seriesLayout->addWidget(addSeriesButton_, 0, Qt::AlignLeft);
    auto* graphGroup = new QGroupBox("Graph"); auto* graphLayout = new QGridLayout(graphGroup);
    graphLayout->addWidget(new QLabel("Graph title"), 0, 0); graphLayout->addWidget(title_, 0, 1, 1, 3);
    graphLayout->addWidget(new QLabel("Aspect W:H"), 1, 0); graphLayout->addWidget(aspectRatio_, 1, 1); graphLayout->addWidget(grid_, 1, 2); graphLayout->addWidget(ticks_, 1, 3);
    graphLayout->addWidget(new QLabel("Legend position"), 2, 0); graphLayout->addWidget(legendPosition_, 2, 1);
    auto* graphAndProfiles = new QHBoxLayout; graphAndProfiles->addWidget(seriesGroup, 3); graphAndProfiles->addWidget(graphGroup, 2); layout->addLayout(graphAndProfiles);
    auto makeAxis = [&](const QString& name, QLineEdit* label, QCheckBox* autoRange, QDoubleSpinBox* minimum, QDoubleSpinBox* maximum, QSpinBox* tickCount, QDoubleSpinBox* numberGap, QDoubleSpinBox* titleGap) {
        auto* group = new QGroupBox(name); auto* box = new QGridLayout(group); box->addWidget(label, 0, 0, 1, 4); box->addWidget(autoRange, 1, 0, 1, 4);
        box->addWidget(new QLabel("Minimum"), 2, 0); box->addWidget(minimum, 2, 1); box->addWidget(new QLabel("Maximum"), 2, 2); box->addWidget(maximum, 2, 3);
        box->addWidget(new QLabel("Axis ticks"), 3, 0); box->addWidget(tickCount, 3, 1); box->addWidget(new QLabel("Number gap"), 3, 2); box->addWidget(numberGap, 3, 3);
        box->addWidget(new QLabel("Title gap"), 4, 2); box->addWidget(titleGap, 4, 3); return group;
    };
    useSecondYAxis_ = new QCheckBox("Enable right Y axis"); useSecondYAxis_->setObjectName("profileEnableSecondYAxis");
    auto* axes = new QHBoxLayout; axes->addWidget(makeAxis("X axis", xLabel_, autoRangeX_, xMin_, xMax_, xTicks_, xNumberGap_, xTitleGap_)); axes->addWidget(makeAxis("Left Y axis", yLabel_, autoRangeY_, yMin_, yMax_, yTicks_, yNumberGap_, yTitleGap_));
    y2Label_ = new QLineEdit; y2Label_->setObjectName("profileY2Label"); y2Label_->setPlaceholderText("Right Y axis title");
    autoRangeY2_ = new QCheckBox("Auto range"); autoRangeY2_->setObjectName("profileAutoRangeY2"); autoRangeY2_->setChecked(true);
    y2Min_ = spin("profileY2Min", 0); y2Max_ = spin("profileY2Max", 1); y2Ticks_ = new QSpinBox; y2Ticks_->setObjectName("profileY2TickCount"); y2Ticks_->setRange(2, 12); y2Ticks_->setValue(5);
    y2NumberGap_ = spin("profileY2NumberGap", 5); y2TitleGap_ = spin("profileY2TitleGap", 8);
    for (auto* gap : {y2NumberGap_, y2TitleGap_}) { gap->setRange(0, 50); gap->setDecimals(0); gap->setSingleStep(1); }
    auto* y2Group = new QGroupBox("Right Y axis"); auto* y2Layout = new QVBoxLayout(y2Group); y2Layout->addWidget(y2Label_); y2Layout->addWidget(autoRangeY2_);
    auto* y2Range = new QHBoxLayout; y2Range->addWidget(new QLabel("Minimum")); y2Range->addWidget(y2Min_); y2Range->addWidget(new QLabel("Maximum")); y2Range->addWidget(y2Max_); y2Layout->addLayout(y2Range);
    auto* y2TickRow = new QHBoxLayout; y2TickRow->addWidget(new QLabel("Axis ticks")); y2TickRow->addWidget(y2Ticks_); y2Layout->addLayout(y2TickRow);
    auto* y2GapRow = new QHBoxLayout; y2GapRow->addWidget(new QLabel("Number gap")); y2GapRow->addWidget(y2NumberGap_); y2GapRow->addWidget(new QLabel("Title gap")); y2GapRow->addWidget(y2TitleGap_); y2Layout->addLayout(y2GapRow);
    y2Group->setEnabled(false); auto* y2AxisWithToggle = new QVBoxLayout; y2AxisWithToggle->addWidget(useSecondYAxis_); y2AxisWithToggle->addWidget(y2Group); auto* y2Wrap = new QWidget; y2Wrap->setLayout(y2AxisWithToggle); axes->addWidget(y2Wrap); layout->addLayout(axes);
    status_ = new QLabel("Choose a database or open a report file."); status_->setObjectName("profilesStatus"); status_->setWordWrap(true); layout->addWidget(status_);
    plot_ = new ProfilePlot; plot_->setObjectName("profilePlot"); layout->addWidget(plot_, 1);
    layout->addWidget(singleRunGroup_);
    connect(source_, &QComboBox::currentIndexChanged, this, [this] {
        openReport_->setEnabled(source_->currentIndex() == 1);
        if (source_->currentIndex() == 0) loadDatabase();
        else if (source_->currentIndex() == 1 && !reportFile_.isEmpty()) loadReport(reportFile_);
        else if (source_->currentIndex() == 2) showLiveProfiles();
    });
    connect(openReport_, &QPushButton::clicked, this, [this] { const auto f = QFileDialog::getOpenFileName(this, "Open profile report", reportFile_, "Reports (*.txt *.csv *.tsv);;All files (*)"); if (!f.isEmpty()) loadReport(f); });
    connect(singleRunExperiment_, &QComboBox::currentTextChanged, this, [this] { refreshSingleRunDatasets(); });
    connect(singleRunDataset_, &QComboBox::currentTextChanged, this, [this] { refreshSingleRunVariables(); });
    connect(singleRunButton_, &QPushButton::clicked, this, [this] { requestSingleRun(); });
    connect(singleRunPause_, &QPushButton::clicked, this, [this] { if (pauseRequested) pauseRequested(); });
    connect(singleRunStop_, &QPushButton::clicked, this, [this] { if (stopRequested) stopRequested(); });
    connect(singleRunSaveInputs_, &QPushButton::clicked, this, [this] { if (saveFittedInputsRequested) saveFittedInputsRequested(); });
    connect(singleRunVariables_, &QTableWidget::cellDoubleClicked, this, [this](int row, int column) {
        if (column == 0 && row >= 0 && singleRunVariables_->item(row, 0) && variableEditRequested)
            variableEditRequested(singleRunVariables_->item(row, 0)->text());
    });
    connect(addSeriesButton_, &QPushButton::clicked, this, [this] { addSeriesRow(); updatePlot(); });
    connect(useSecondYAxis_, &QCheckBox::toggled, y2Group, &QWidget::setEnabled);
    connect(useSecondYAxis_, &QCheckBox::toggled, this, [this](bool enabled) { for (auto* row : seriesRows_) row->yAxis->setEnabled(enabled); updatePlot(); });
    connect(aspectRatio_, &QDoubleSpinBox::valueChanged, this, [this] { updatePlot(); });
    const QList<QWidget*> settingWidgets{autoRangeX_, autoRangeY_, autoRangeY2_, grid_, ticks_, xMin_, xMax_, yMin_, yMax_, y2Min_, y2Max_, xTicks_, yTicks_, y2Ticks_, xNumberGap_, xTitleGap_, yNumberGap_, yTitleGap_, y2NumberGap_, y2TitleGap_};
    for (auto* widget : settingWidgets) if (auto* check = qobject_cast<QCheckBox*>(widget)) connect(check, &QCheckBox::toggled, this, [this] { updatePlot(); }); else if (auto* real = qobject_cast<QDoubleSpinBox*>(widget)) connect(real, &QDoubleSpinBox::valueChanged, this, [this] { updatePlot(); }); else if (auto* integer = qobject_cast<QSpinBox*>(widget)) connect(integer, &QSpinBox::valueChanged, this, [this] { updatePlot(); });
    for (auto* edit : {title_, xLabel_, yLabel_, y2Label_}) connect(edit, &QLineEdit::textChanged, this, [this] { updatePlot(); });
    connect(legendPosition_, &QComboBox::currentIndexChanged, this, [this] { updatePlot(); });
    connect(exportImage_, &QPushButton::clicked, this, [this] { QString name = plot_->displayTitle().trimmed(); if (name.isEmpty()) name = "profile-graph"; name.replace(QRegularExpression("[<>:\"/\\\\|?*\\x00-\\x1f]"), "_"); while (name.endsWith('.') || name.endsWith(' ')) name.chop(1); const auto filename = QFileDialog::getSaveFileName(this, "Export profile graph", name + ".png", "PNG image (*.png)"); if (filename.isEmpty()) return; QImage image(plot_->size(), QImage::Format_ARGB32_Premultiplied); image.fill(Qt::transparent); plot_->render(&image); QString output = filename; if (QFileInfo(output).suffix().isEmpty()) output += ".png"; status_->setText(image.save(output, "PNG") ? "Graph image exported: " + output : "Could not export graph image: " + output); });
    connect(source_, &QComboBox::currentIndexChanged, this, [this] { if (seriesRows_.isEmpty()) addSeriesRow(); });
    loadDatabase();
}

void ProfilesTab::setSingleRunEnabled(bool enabled)
{
    singleRunExperiment_->setEnabled(enabled);
    singleRunDataset_->setEnabled(enabled);
    singleRunSolver_->setEnabled(enabled);
    singleRunVariables_->setEnabled(enabled);
    singleRunButton_->setEnabled(enabled && singleRunDataset_->count() > 0);
}

void ProfilesTab::setSingleRunControlState(bool canPause, bool canStop, bool canSaveInputs)
{
    singleRunPause_->setEnabled(canPause);
    singleRunStop_->setEnabled(canStop);
    singleRunSaveInputs_->setEnabled(canSaveInputs);
}

void ProfilesTab::setSingleRunStatus(const QString& message)
{
    singleRunStatus_->setText(message);
}

void ProfilesTab::setSingleRunContext(const QString& filename, const model_controls::Snapshot& controls)
{
    const auto previousExperiment = singleRunExperiment_->currentText();
    const auto previousDataset = singleRunDataset_->currentText();
    const bool unchanged = singleRunDatabase_ == filename && singleRunControls_ == controls;
    singleRunDatabase_ = filename;
    singleRunControls_ = controls;
    if (unchanged) return;

    QStringList experiments;
    for (const auto& row : controls.rows)
        if (row.name == "experiment_name") experiments = row.value.value_or(QString{}).split(' ', Qt::SkipEmptyParts);
    experiments.removeDuplicates();
    singleRunExperiment_->blockSignals(true);
    singleRunExperiment_->clear(); singleRunExperiment_->addItems(experiments);
    int experimentIndex = singleRunExperiment_->findText(previousExperiment);
    if (experimentIndex < 0 && !experiments.isEmpty()) experimentIndex = 0;
    singleRunExperiment_->setCurrentIndex(experimentIndex);
    singleRunExperiment_->blockSignals(false);

    bool solver = true;
    for (const auto& row : controls.rows)
        if (row.name == "run_solver") solver = row.value.value_or("true") == "true";
    singleRunSolver_->setChecked(solver);
    refreshSingleRunDatasets(previousDataset);
}

void ProfilesTab::refreshSingleRunDatasets(const QString& preferredDataset)
{
    const auto wanted = preferredDataset.isEmpty() ? singleRunDataset_->currentText() : preferredDataset;
    const auto datasets = model_controls::profileDatasets(singleRunDatabase_, singleRunExperiment_->currentText());
    singleRunDataset_->blockSignals(true);
    singleRunDataset_->clear(); singleRunDataset_->addItems(datasets);
    int index = singleRunDataset_->findText(wanted);
    if (index < 0 && !datasets.isEmpty()) index = 0;
    singleRunDataset_->setCurrentIndex(index);
    singleRunDataset_->blockSignals(false);
    singleRunStatus_->setText(datasets.isEmpty() ? "This experiment has no profile datasets." : QString{});
    singleRunButton_->setEnabled(singleRunExperiment_->isEnabled() && !datasets.isEmpty());
    refreshSingleRunVariables();
}

void ProfilesTab::refreshSingleRunVariables()
{
    singleRunVariables_->setRowCount(0);
    if (singleRunDatabase_.isEmpty() || singleRunExperiment_->currentText().isEmpty() || singleRunDataset_->currentText().isEmpty()) return;
    for (const auto& name : model_controls::applicableProfileVariables(singleRunDatabase_, singleRunControls_, singleRunExperiment_->currentText(), singleRunDataset_->currentText())) {
        const int row = singleRunVariables_->rowCount(); singleRunVariables_->insertRow(row);
        auto* nameItem = new QTableWidgetItem(name); nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable); singleRunVariables_->setItem(row, 0, nameItem);
        auto* fitItem = new QTableWidgetItem; fitItem->setFlags((fitItem->flags() | Qt::ItemIsUserCheckable) & ~Qt::ItemIsEditable); fitItem->setCheckState(Qt::Checked); singleRunVariables_->setItem(row, 1, fitItem);
        auto* startCell = new QWidget(singleRunVariables_);
        auto* startLayout = new QHBoxLayout(startCell); startLayout->setContentsMargins(0, 0, 0, 0); startLayout->setSpacing(4);
        auto* startValue = new QLineEdit(startCell); startValue->setObjectName("singleRunStartValue_" + name);
        auto* validator = new QDoubleValidator(-1e300, 1e300, 15, startValue); validator->setNotation(QDoubleValidator::ScientificNotation); startValue->setValidator(validator);
        const auto initialValue = model_controls::profileVariableInitialValue(singleRunDatabase_, singleRunControls_, singleRunExperiment_->currentText(), singleRunDataset_->currentText(), name);
        if (initialValue) startValue->setText(QString::number(*initialValue, 'g', 15));
        startValue->setProperty("singleRunInitialValue", initialValue ? QVariant::fromValue(*initialValue) : QVariant{});
        auto* previousValue = new QPushButton(startCell); previousValue->setObjectName("singleRunPreviousValue_" + name); previousValue->hide();
        startLayout->addWidget(startValue, 1); startLayout->addWidget(previousValue);
        singleRunVariables_->setCellWidget(row, 2, startCell);
        connect(previousValue, &QPushButton::clicked, this, [startValue, previousValue] {
            bool ok = false;
            const double current = startValue->text().toDouble(&ok);
            if (!ok) return;
            const double alternate = previousValue->property("alternateStartingValue").toDouble();
            startValue->setText(QString::number(alternate, 'g', 15));
            previousValue->setProperty("alternateStartingValue", current);
            previousValue->setText("Switch to: " + QString::number(current, 'g', 8));
            previousValue->setToolTip("Switch between the previous start and the saved fitted value.");
        });
    }
}

void ProfilesTab::acceptSavedFittedInputs()
{
    for (int row = 0; row < singleRunVariables_->rowCount(); ++row) {
        auto* cell = singleRunVariables_->cellWidget(row, 2);
        auto* startValue = cell ? cell->findChild<QLineEdit*>() : nullptr;
        auto* previousValue = cell ? cell->findChild<QPushButton*>() : nullptr;
        if (!startValue || !previousValue) continue;
        const auto variable = singleRunVariables_->item(row, 0)->text();
        const auto fitted = model_controls::profileVariableInitialValue(singleRunDatabase_, singleRunControls_,
            singleRunExperiment_->currentText(), singleRunDataset_->currentText(), variable);
        if (!fitted) continue;
        bool ok = false;
        const double oldValue = startValue->text().toDouble(&ok);
        const QString fittedText = QString::number(*fitted, 'g', 15);
        if (ok && oldValue != *fitted) {
            previousValue->setProperty("alternateStartingValue", oldValue);
            previousValue->setText("Switch to: " + QString::number(oldValue, 'g', 8));
            previousValue->setToolTip("Switch between the previous start and the saved fitted value.");
            previousValue->show();
        } else {
            previousValue->hide();
        }
        startValue->setText(fittedText);
        startValue->setProperty("singleRunInitialValue", *fitted);
    }
    singleRunStatus_->setText("Saved fitted values are now the starting values. Previous starts are available beside changed values.");
}

void ProfilesTab::requestSingleRun()
{
    if (singleRunExperiment_->currentText().isEmpty() || singleRunDataset_->currentText().isEmpty()) return;
    SingleProfileRunRequest request;
    request.experiment = singleRunExperiment_->currentText(); request.dataset = singleRunDataset_->currentText();
    request.runSolver = singleRunSolver_->isChecked();
    for (int row = 0; row < singleRunVariables_->rowCount(); ++row) {
        if (singleRunVariables_->item(row, 1)->checkState() != Qt::Checked) continue;
        const auto variable = singleRunVariables_->item(row, 0)->text();
        auto* cell = singleRunVariables_->cellWidget(row, 2);
        auto* startValue = cell ? cell->findChild<QLineEdit*>() : nullptr;
        bool ok = false; const double value = startValue ? startValue->text().toDouble(&ok) : 0.0;
        if (!ok || !std::isfinite(value)) {
            QMessageBox::warning(this, "Invalid starting value", "Enter a finite starting value for " + variable + "."); return;
        }
        request.variables.push_back(variable);
        const QVariant original = startValue->property("singleRunInitialValue");
        if (!original.isValid() || value != original.toDouble()) request.initialValues[variable.toStdString()] = value;
    }
    if (singleRunRequested) singleRunRequested(request);
}

void ProfilesTab::addSeriesRow(const QJsonObject& value)
{
    auto* row = new SeriesRow;
    row->dataset = new QComboBox; row->profile = new QComboBox; row->yAxis = new QComboBox; row->yAxis->addItems({"Left", "Right"}); row->yAxis->setCurrentIndex(value.value("yAxis").toInt(0)); row->yAxis->setEnabled(useSecondYAxis_->isChecked()); row->lineColorButton = new QPushButton("Line"); row->fillColorButton = new QPushButton("Fill");
    static const QColor defaults[] = {QColor("#2563eb"), QColor("#16a34a"), QColor("#9333ea"), QColor("#ea580c"), QColor("#0891b2"), QColor("#db2777"), QColor("#65a30d"), QColor("#4f46e5")};
    const QColor defaultLine = defaults[seriesRows_.size() % (sizeof(defaults) / sizeof(defaults[0]))];
    QColor defaultFill = defaultLine; defaultFill.setAlpha(72);
    row->lineColor = QColor(value.value("lineColor").toString(defaultLine.name(QColor::HexArgb))); row->fillColor = QColor(value.value("fillColor").toString(defaultFill.name(QColor::HexArgb)));
    row->fill = new QCheckBox; row->fill->setChecked(value.value("fill").toBool(false)); row->legend = new QLineEdit; row->legend->setPlaceholderText("Profile name"); row->lastProfile = value.value("profile").toString();
    auto style = [](QPushButton* button, const QColor& color) { button->setStyleSheet(QString("background-color:%1").arg(color.name(QColor::HexArgb))); };
    style(row->lineColorButton, row->lineColor); style(row->fillColorButton, row->fillColor);
    const int index = seriesTable_->rowCount(); seriesTable_->insertRow(index); seriesTable_->setCellWidget(index, 0, row->dataset); seriesTable_->setCellWidget(index, 1, row->profile); seriesTable_->setCellWidget(index, 2, row->yAxis); seriesTable_->setCellWidget(index, 3, row->lineColorButton); seriesTable_->setCellWidget(index, 4, row->fillColorButton); seriesTable_->setCellWidget(index, 5, row->fill); seriesTable_->setCellWidget(index, 6, row->legend);
    auto* remove = new QPushButton("Remove"); seriesTable_->setCellWidget(index, 7, remove); seriesRows_.push_back(row);
    refreshSeriesRows();
    if (!value.value("datasetId").toString().isEmpty()) { int ds = row->dataset->findData(value.value("datasetId").toString()); if (ds >= 0) row->dataset->setCurrentIndex(ds); }
    if (!row->lastProfile.isEmpty()) { int p = row->profile->findText(row->lastProfile); if (p >= 0) row->profile->setCurrentIndex(p); }
    row->legend->setText(value.value("legend").toString());
    connect(row->dataset, &QComboBox::currentIndexChanged, this, [this, row] { refreshSeriesRows(); updatePlot(); });
    connect(row->yAxis, &QComboBox::currentIndexChanged, this, [this] { updatePlot(); });
    connect(row->profile, &QComboBox::currentTextChanged, this, [this, row](const QString& text) { const QString oldProfile = row->lastProfile; if (row->legend->text().isEmpty() || row->legend->text() == oldProfile) row->legend->setText(text); row->lastProfile = text; updatePlot(); });
    connect(row->legend, &QLineEdit::textChanged, this, [this] { updatePlot(); }); connect(row->fill, &QCheckBox::toggled, this, [this] { updatePlot(); });
    connect(row->lineColorButton, &QPushButton::clicked, this, [this, row] { chooseLineColor(row); }); connect(row->fillColorButton, &QPushButton::clicked, this, [this, row] { chooseFillColor(row); });
    connect(remove, &QPushButton::clicked, this, [this, row] { removeSeriesRow(row); });
}

void ProfilesTab::removeSeriesRow(SeriesRow* row)
{
    const int index = seriesRows_.indexOf(row); if (index < 0) return; seriesTable_->removeRow(index); seriesRows_.removeAt(index); delete row->dataset; delete row->profile; delete row->yAxis; delete row->lineColorButton; delete row->fillColorButton; delete row->fill; delete row->legend; delete row; updatePlot();
}

void ProfilesTab::refreshSeriesRows()
{
    for (auto* row : seriesRows_) {
        const QString dsId = row->dataset->currentData().toString(); const QString profile = row->profile->currentText();
        row->dataset->blockSignals(true); row->dataset->clear(); for (int i = 0; i < datasetIds_.size(); ++i) row->dataset->addItem(datasetLabels_.value(i), datasetIds_[i]);
        int ds = row->dataset->findData(dsId);
        // Presets saved before version browsing stored only experiment and weight.
        // Dataset versions are ordered newest first within that pair.
        if (ds < 0) {
            const auto legacyKey = dsId.split(QChar(0x1f));
            if (legacyKey.size() == 2) {
                const QString prefix = dsId + QChar(0x1f);
                for (int i = 0; i < datasetIds_.size(); ++i) {
                    if (datasetIds_[i].startsWith(prefix)) { ds = i; break; }
                }
            }
        }
        row->dataset->setCurrentIndex(ds >= 0 ? ds : (datasetIds_.isEmpty() ? -1 : 0)); row->dataset->blockSignals(false);
        row->profile->blockSignals(true); row->profile->clear(); row->profile->addItems(profiles_); int p = row->profile->findText(profile.isEmpty() ? row->lastProfile : profile); row->profile->setCurrentIndex(p >= 0 ? p : (profiles_.isEmpty() ? -1 : 0)); row->profile->blockSignals(false);
        row->lastProfile = row->profile->currentText();
    }
}

void ProfilesTab::chooseLineColor(SeriesRow* row) { auto color = QColorDialog::getColor(row->lineColor, this, "Profile line color"); if (color.isValid()) { row->lineColor = color; row->lineColorButton->setStyleSheet(QString("background-color:%1").arg(color.name(QColor::HexArgb))); updatePlot(); } }
void ProfilesTab::chooseFillColor(SeriesRow* row) { auto color = QColorDialog::getColor(row->fillColor, this, "Profile fill color", QColorDialog::ShowAlphaChannel); if (color.isValid()) { row->fillColor = color; row->fillColorButton->setStyleSheet(QString("background-color:%1").arg(color.name(QColor::HexArgb))); updatePlot(); } }

QJsonObject ProfilesTab::settings() const
{
    QJsonArray rows; for (auto* row : seriesRows_) rows.append(QJsonObject{{"datasetId", row->dataset->currentData().toString()}, {"profile", row->profile->currentText()}, {"yAxis", row->yAxis->currentIndex()}, {"lineColor", row->lineColor.name(QColor::HexArgb)}, {"fillColor", row->fillColor.name(QColor::HexArgb)}, {"fill", row->fill->isChecked()}, {"legend", row->legend->text()}});
    return {{"series", rows}, {"title", title_->text()}, {"aspectRatio", aspectRatio_->value()}, {"legendPosition", legendPosition_->currentIndex()}, {"xLabel", xLabel_->text()}, {"yLabel", yLabel_->text()}, {"useSecondYAxis", useSecondYAxis_->isChecked()}, {"y2Label", y2Label_->text()}, {"autoRangeY2", autoRangeY2_->isChecked()}, {"y2Min", y2Min_->value()}, {"y2Max", y2Max_->value()}, {"y2Ticks", y2Ticks_->value()}, {"y2NumberGap", y2NumberGap_->value()}, {"y2TitleGap", y2TitleGap_->value()}, {"autoRange", autoRangeX_->isChecked() && autoRangeY_->isChecked()}, {"autoRangeX", autoRangeX_->isChecked()}, {"autoRangeY", autoRangeY_->isChecked()}, {"xMin", xMin_->value()}, {"xMax", xMax_->value()}, {"yMin", yMin_->value()}, {"yMax", yMax_->value()}, {"xTicks", xTicks_->value()}, {"yTicks", yTicks_->value()}, {"xNumberGap", xNumberGap_->value()}, {"xTitleGap", xTitleGap_->value()}, {"yNumberGap", yNumberGap_->value()}, {"yTitleGap", yTitleGap_->value()}, {"grid", grid_->isChecked()}, {"ticks", ticks_->isChecked()}};
}

void ProfilesTab::setSettings(const QJsonObject& value)
{
    title_->setText(value.value("title").toString()); aspectRatio_->setValue(value.value("aspectRatio").toDouble(1.6)); legendPosition_->setCurrentIndex(value.value("legendPosition").toInt(1)); xLabel_->setText(value.value("xLabel").toString()); yLabel_->setText(value.value("yLabel").toString()); useSecondYAxis_->setChecked(value.value("useSecondYAxis").toBool(false)); y2Label_->setText(value.value("y2Label").toString()); autoRangeY2_->setChecked(value.value("autoRangeY2").toBool(true)); y2Min_->setValue(value.value("y2Min").toDouble(0)); y2Max_->setValue(value.value("y2Max").toDouble(1)); y2Ticks_->setValue(value.value("y2Ticks").toInt(5)); y2NumberGap_->setValue(value.value("y2NumberGap").toDouble(5)); y2TitleGap_->setValue(value.value("y2TitleGap").toDouble(8));
    const bool legacy = value.value("autoRange").toBool(true); autoRangeX_->setChecked(value.value("autoRangeX").toBool(legacy)); autoRangeY_->setChecked(value.value("autoRangeY").toBool(legacy));
    xMin_->setValue(value.value("xMin").toDouble(0)); xMax_->setValue(value.value("xMax").toDouble(1)); yMin_->setValue(value.value("yMin").toDouble(0)); yMax_->setValue(value.value("yMax").toDouble(1)); xTicks_->setValue(value.value("xTicks").toInt(5)); yTicks_->setValue(value.value("yTicks").toInt(5));
    xNumberGap_->setValue(value.value("xNumberGap").toDouble(5)); xTitleGap_->setValue(value.value("xTitleGap").toDouble(6)); yNumberGap_->setValue(value.value("yNumberGap").toDouble(5)); yTitleGap_->setValue(value.value("yTitleGap").toDouble(8)); grid_->setChecked(value.value("grid").toBool(true)); ticks_->setChecked(value.value("ticks").toBool(true));
    while (!seriesRows_.isEmpty()) removeSeriesRow(seriesRows_.back());
    const auto savedRows = value.value("series").toArray(); for (const auto& item : savedRows) addSeriesRow(item.toObject());
    if (savedRows.isEmpty()) {
        QStringList old; for (const auto& item : value.value("profiles").toArray()) old.push_back(item.toString()); if (old.isEmpty() && !value.value("profile").toString().isEmpty()) old.push_back(value.value("profile").toString());
        if (old.isEmpty()) addSeriesRow(); else for (const auto& profile : old) addSeriesRow(QJsonObject{{"profile", profile}, {"datasetId", datasetIds_.value(value.value("block").toInt(0))}, {"legend", profile}, {"lineColor", value.value("lineColor").toString("#ff2563eb")}, {"fillColor", value.value("fillColor").toString("#482563eb")}, {"fill", value.value("fill").toBool(false)}});
    }
    updatePlot();
}

void ProfilesTab::setDatabase(const QString& filename)
{
    if (database_ != filename) {
        liveProfiles_.clear();
        if (source_->currentIndex() == 2) source_->setCurrentIndex(0);
    }
    database_ = filename;
    if (source_->currentIndex() == 0) loadDatabase();
}

void ProfilesTab::beginLiveRun()
{
    liveProfiles_.clear();
    source_->setCurrentIndex(2);
    showLiveProfiles();
}

void ProfilesTab::setLiveProfiles(std::vector<tsensor_workflow::live_profile_snapshot> profiles)
{
    liveProfiles_ = std::move(profiles);
    if (source_->currentIndex() == 2) showLiveProfiles();
}

void ProfilesTab::showLiveProfiles()
{
    datasetLabels_.clear(); datasetIds_.clear();
    profiles_ = {"Experimental", "Numeric model"};
    for (const auto& profile : liveProfiles_) {
        const QString name = QString::fromStdString(profile.experiment_name);
        const QString wt = QString::fromStdString(profile.weight_percent);
        const QString id = name + QChar(0x1f) + wt;
        if (datasetIds_.contains(id)) continue;
        const QString profileName = wt.isEmpty() ? name : name + " · " + wt;
        datasetLabels_ << profileName + " · live";
        datasetIds_ << id;
    }
    refreshSeriesRows();
    if (seriesRows_.isEmpty() && !datasetIds_.isEmpty()) addSeriesRow();
    status_->setText(datasetIds_.isEmpty() ? "Waiting for the first model evaluation…" : "Live preview of the latest completed evaluation; save profiles after the run to keep them.");
    updatePlot();
}

void ProfilesTab::loadDatabase()
{
    datasetLabels_.clear(); datasetIds_.clear(); profiles_.clear();
    if (database_.isEmpty()) { status_->setText("Choose a database file."); refreshSeriesRows(); updatePlot(); return; }
    sqlite3* raw = nullptr; const auto encoded = database_.toUtf8(); if (sqlite3_open_v2(encoded.constData(), &raw, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) { const QString error = raw ? QString::fromUtf8(sqlite3_errmsg(raw)) : "Cannot open database"; if (raw) sqlite3_close(raw); status_->setText(error); return; }
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close); sqlite3_stmt* statement = nullptr;
    const char* query = "SELECT EXPERIMENT_NAME, WT_PERCENT, RUN_AT_UTC_US FROM profile_records ORDER BY EXPERIMENT_NAME, WT_PERCENT, RUN_AT_UTC_US DESC";
    if (sqlite3_prepare_v2(db.get(), query, -1, &statement, nullptr) != SQLITE_OK) { status_->setText("No saved model profiles are available in this database."); return; }
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> stmt(statement, sqlite3_finalize);
    while (sqlite3_step(stmt.get()) == SQLITE_ROW) {
        const auto* rawName = sqlite3_column_text(stmt.get(), 0);
        const auto* rawWt = sqlite3_column_text(stmt.get(), 1);
        const QString name = rawName ? QString::fromUtf8(reinterpret_cast<const char*>(rawName)) : QString{};
        const QString wt = rawWt ? QString::fromUtf8(reinterpret_cast<const char*>(rawWt)) : QString{};
        const sqlite3_int64 runAtUtcUs = sqlite3_column_int64(stmt.get(), 2);
        const auto savedAt = QDateTime::fromMSecsSinceEpoch(runAtUtcUs / 1000, QTimeZone::UTC);
        const QString fractionalMicroseconds = QString::number(std::abs(runAtUtcUs % 1000)).rightJustified(3, '0');
        const QString version = savedAt.toString("yyyy-MM-dd HH:mm:ss.zzz") + fractionalMicroseconds + " UTC";
        const QString profileName = wt.isEmpty() ? name : name + " · " + wt;
        datasetLabels_ << (profileName + " · " + version);
        datasetIds_ << name + QChar(0x1f) + wt + QChar(0x1f) + QString::number(runAtUtcUs);
    }
    profiles_ = {"Experimental", "Numeric model", "Experimental derivative", "Numeric derivative", "Experimental difference", "Numeric difference", "Free dye", "Bound dye", "Total dye", "Unbound beads", "Bound beads", "Total beads", "Error"};
    refreshSeriesRows(); if (seriesRows_.isEmpty() && !datasetIds_.isEmpty()) addSeriesRow(); status_->setText(datasetIds_.isEmpty() ? "No saved model profiles are available in this database." : "Saved profiles from the selected database."); updatePlot();
}

void ProfilesTab::loadReport(const QString& filename)
{
    reportFile_ = filename;
    try { QFile file(filename); if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error(file.errorString().toStdString()); const auto rows = report_format::parse(QString::fromUtf8(file.readAll())); datasetLabels_.clear(); datasetIds_.clear(); profiles_.clear(); int blocks = 0;
        for (const auto& row : rows) { if (row.isEmpty()) continue; if (row[0] == "sec") { datasetLabels_ << QString("Profile block %1 · %2 sec").arg(++blocks).arg(row.size() > 1 ? row[1] : QString{}); datasetIds_ << QString::number(blocks - 1); } else if (row[0] != "res_time" && row.size() >= 10 && !profiles_.contains(row[8])) profiles_ << row[8]; }
        refreshSeriesRows(); if (seriesRows_.isEmpty() && !datasetIds_.isEmpty()) addSeriesRow(); status_->setText("Loaded report: " + filename); source_->setCurrentIndex(1); updatePlot();
    } catch (const std::exception& error) { status_->setText("Cannot load profile report: " + QString::fromUtf8(error.what())); }
}

void ProfilesTab::updatePlot()
{
    QString autoTitle; QStringList names; for (auto* row : seriesRows_) { const QString name = row->legend->text().trimmed().isEmpty() ? row->profile->currentText() : row->legend->text().trimmed(); if (!name.isEmpty()) names << name; }
    if (!datasetLabels_.isEmpty()) autoTitle = datasetLabels_.value(seriesRows_.isEmpty() ? 0 : datasetLabels_.indexOf(seriesRows_.front()->dataset->currentText())) + " · " + names.join(", ");
    QString yAuto = names.size() == 1 ? names.front() : (names.isEmpty() ? QString{} : "Profiles");
    plot_->setAppearance(grid_->isChecked(), ticks_->isChecked(), autoRangeX_->isChecked(), autoRangeY_->isChecked(), useSecondYAxis_->isChecked(), autoRangeY2_->isChecked(), xMin_->value(), xMax_->value(), yMin_->value(), yMax_->value(), y2Min_->value(), y2Max_->value(), xTicks_->value(), yTicks_->value(), y2Ticks_->value(), xLabel_->text().isEmpty() ? "X (µm)" : xLabel_->text(), yLabel_->text().isEmpty() ? yAuto : yLabel_->text(), y2Label_->text().isEmpty() ? "Right Y axis" : y2Label_->text(), title_->text().trimmed().isEmpty() ? autoTitle : title_->text().trimmed(), aspectRatio_->value(), legendPosition_->currentIndex(), qRound(xNumberGap_->value()), qRound(xTitleGap_->value()), qRound(yNumberGap_->value()), qRound(yTitleGap_->value()), qRound(y2NumberGap_->value()), qRound(y2TitleGap_->value()));
    QVector<ProfileSeriesData> series;
    if (source_->currentIndex() == 1) {
        try { QFile file(reportFile_); if (!file.open(QIODevice::ReadOnly)) { plot_->setData(autoTitle, "X (µm)", yAuto, {}); return; } const auto rows = report_format::parse(QString::fromUtf8(file.readAll()));
            struct ReportBlock { QVector<double> x; QHash<QString, QVector<double>> y; }; QVector<ReportBlock> blocks; ReportBlock* block = nullptr;
            QVector<double> currentAxis;
            for (const auto& row : rows) { if (row.isEmpty()) continue; if (row[0] == "res_time") { currentAxis.clear(); for (int i = 9; i < row.size(); ++i) { bool ok = false; double v = row[i].toDouble(&ok); currentAxis << (ok ? v : std::numeric_limits<double>::quiet_NaN()); } } else if (row[0] == "sec") { blocks.push_back({}); block = &blocks.back(); block->x = currentAxis; } else if (block && row.size() >= 10) { auto& values = block->y[row[8]]; for (int i = 9; i < row.size(); ++i) { bool ok = false; double v = row[i].toDouble(&ok); values << (ok ? v : std::numeric_limits<double>::quiet_NaN()); } } }
            for (auto* row : seriesRows_) { bool ok = false; int b = row->dataset->currentData().toInt(&ok); if (!ok || b < 0 || b >= blocks.size()) continue; const QString profile = row->profile->currentText(); auto values = blocks[b].y.value(profile); int n = std::min(blocks[b].x.size(), values.size()); if (n < 2) continue; ProfileSeriesData data; data.name = row->legend->text().trimmed().isEmpty() ? profile : row->legend->text().trimmed(); data.lineColor = row->lineColor; data.fillColor = row->fillColor; data.filled = row->fill->isChecked(); data.yAxis = row->yAxis->currentIndex(); data.x = blocks[b].x; data.y = values; data.x.resize(n); data.y.resize(n); series << data; }
        } catch (...) {}
    } else if (source_->currentIndex() == 2) {
        for (auto* row : seriesRows_) {
            const auto key = row->dataset->currentData().toString().split(QChar(0x1f));
            const int col = profiles_.indexOf(row->profile->currentText());
            if (key.size() != 2 || (col != 0 && col != 1)) continue;
            const auto found = std::find_if(liveProfiles_.begin(), liveProfiles_.end(), [&](const auto& profile) {
                return QString::fromStdString(profile.experiment_name) == key[0] &&
                       QString::fromStdString(profile.weight_percent) == key[1];
            });
            if (found == liveProfiles_.end()) continue;
            const auto& values = col == 0 ? found->experimental : found->numeric;
            const int count = static_cast<int>(std::min(found->x.size(), values.size()));
            if (count < 2) continue;
            ProfileSeriesData data;
            const QString profile = row->profile->currentText();
            data.name = row->legend->text().trimmed().isEmpty() ? profile : row->legend->text().trimmed();
            data.lineColor = row->lineColor; data.fillColor = row->fillColor; data.filled = row->fill->isChecked(); data.yAxis = row->yAxis->currentIndex();
            data.x.reserve(count); data.y.reserve(count);
            for (int i = 0; i < count; ++i) { data.x << found->x[static_cast<std::size_t>(i)]; data.y << values[static_cast<std::size_t>(i)]; }
            series << data;
        }
    } else {
        static const char* columns[] = {"Experimental", "Numeric", "Experimental_Derivative", "Numeric_Derivative", "Experimental_Difference", "Numeric_Difference", "Free_Dye", "Bound_Dye", "Total_Dye", "Unbound_Beads", "Bound_Beads", "Total_Beads", "Error"};
        sqlite3* db = nullptr; const auto encoded = database_.toUtf8(); if (sqlite3_open_v2(encoded.constData(), &db, SQLITE_OPEN_READONLY, nullptr) == SQLITE_OK) {
            for (auto* row : seriesRows_) { const auto key = row->dataset->currentData().toString().split(QChar(0x1f)); const int col = profiles_.indexOf(row->profile->currentText()); if (key.size() != 3 || col < 0 || col >= int(std::size(columns))) continue;
                const QByteArray query = QByteArray("SELECT X,") + columns[col] + " FROM profile_record_samples WHERE EXPERIMENT_NAME=? AND WT_PERCENT=? AND RUN_AT_UTC_US=? ORDER BY X"; sqlite3_stmt* stmt = nullptr; if (sqlite3_prepare_v2(db, query.constData(), -1, &stmt, nullptr) != SQLITE_OK) continue;
                const auto e = key[0].toUtf8(), w = key[1].toUtf8(), timestamp = key[2].toUtf8(); sqlite3_bind_text(stmt, 1, e.constData(), -1, SQLITE_TRANSIENT); sqlite3_bind_text(stmt, 2, w.constData(), -1, SQLITE_TRANSIENT); sqlite3_bind_int64(stmt, 3, timestamp.toLongLong());
                ProfileSeriesData data; const QString profile = row->profile->currentText(); data.name = row->legend->text().trimmed().isEmpty() ? profile : row->legend->text().trimmed(); data.lineColor = row->lineColor; data.fillColor = row->fillColor; data.filled = row->fill->isChecked(); data.yAxis = row->yAxis->currentIndex();
                while (sqlite3_step(stmt) == SQLITE_ROW) { if (sqlite3_column_type(stmt, 1) == SQLITE_NULL) continue; const auto* txt = sqlite3_column_text(stmt, 1); bool valid = false; const double y = txt ? QString::fromUtf8(reinterpret_cast<const char*>(txt)).trimmed().toDouble(&valid) : 0; const double x = sqlite3_column_double(stmt, 0); if (valid && std::isfinite(x) && std::isfinite(y)) { data.x << x; data.y << y; } }
                sqlite3_finalize(stmt); if (data.x.size() >= 2) series << data;
            } sqlite3_close(db);
        }
    }
    plot_->setData(title_->text().trimmed().isEmpty() ? autoTitle : title_->text().trimmed(), xLabel_->text().isEmpty() ? "X (µm)" : xLabel_->text(), yLabel_->text().isEmpty() ? yAuto : yLabel_->text(), series);
    if (!seriesRows_.isEmpty() && !datasetIds_.isEmpty() && series.isEmpty()) status_->setText("No numeric samples are available for the selected profiles.");
}
