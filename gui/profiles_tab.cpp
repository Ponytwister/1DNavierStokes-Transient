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
#include <QLineEdit>
#include <QSpinBox>
#include <QFileInfo>
#include <QImage>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <sqlite3.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

ProfilePlot::ProfilePlot(QWidget* parent) : QWidget(parent) { setMinimumSize(480, 320); }

void ProfilePlot::setData(const QString& title, const QString& xLabel, const QString& yLabel,
                          const QVector<double>& x, const QVector<double>& y)
{
    title_ = title; xLabel_ = xLabel; yLabel_ = yLabel; x_ = x; y_ = y; update();
}

void ProfilePlot::setAppearance(const QColor& line, const QColor& fill, bool filled, bool grid,
                                bool ticks, bool autoRange, double xmin, double xmax, double ymin, double ymax,
                                int xTickCount, int yTickCount, const QString& xLabel, const QString& yLabel)
{
    lineColor_ = line; fillColor_ = fill; filled_ = filled; grid_ = grid; ticks_ = ticks;
    autoRange_ = autoRange; xmin_ = xmin; xmax_ = xmax; ymin_ = ymin; ymax_ = ymax;
    xTickCount_ = xTickCount; yTickCount_ = yTickCount;
    if (!xLabel.isEmpty()) xLabel_ = xLabel;
    if (!yLabel.isEmpty()) yLabel_ = yLabel;
    update();
}

void ProfilePlot::paintEvent(QPaintEvent*)
{
    QPainter painter(this); painter.fillRect(rect(), palette().base()); painter.setRenderHint(QPainter::Antialiasing);
    const QRectF area(76, 42, width() - 108, height() - 104);
    painter.setPen(palette().text().color()); painter.drawText(QRectF(8, 8, width() - 16, 24), Qt::AlignCenter, title_);
    if (x_.size() < 2 || x_.size() != y_.size() || area.width() <= 0 || area.height() <= 0) {
        painter.drawText(area, Qt::AlignCenter, "Select profile data to plot"); return;
    }
    const auto [xmin, xmax] = std::minmax_element(x_.cbegin(), x_.cend());
    const auto [ymin, ymax] = std::minmax_element(y_.cbegin(), y_.cend());
    double lowX = autoRange_ ? *xmin : xmin_, highX = autoRange_ ? *xmax : xmax_;
    double lowY = autoRange_ ? *ymin : ymin_, highY = autoRange_ ? *ymax : ymax_;
    double dx = highX - lowX, dy = highY - lowY;
    if (!std::isfinite(dx) || !std::isfinite(dy)) return;
    if (dx == 0) dx = 1; if (dy == 0) dy = std::max(1.0, std::abs(lowY) * 0.1);
    const double x0 = autoRange_ ? lowX - dx * .02 : lowX, x1 = autoRange_ ? highX + dx * .02 : highX;
    const double y0 = autoRange_ ? lowY - dy * .08 : lowY, y1 = autoRange_ ? highY + dy * .08 : highY;
    if (x1 <= x0 || y1 <= y0) return;
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
        painter.drawText(QRectF(8, y - 10, 62, 20), Qt::AlignRight | Qt::AlignVCenter, QString::number(value, 'g', 4));
    }
    for (int i = 0; i <= xTickCount_; ++i) {
        const qreal x = area.left() + area.width() * i / xTickCount_;
        const double value = x0 + (x1 - x0) * i / xTickCount_;
        if (ticks_) painter.drawLine(QPointF(x, area.bottom()), QPointF(x, area.bottom() + 4));
        painter.drawText(QRectF(x - 35, area.bottom() + 5, 70, 16), Qt::AlignHCenter | Qt::AlignTop, QString::number(value, 'g', 4));
    }
    painter.drawText(QRectF(area.left(), area.bottom() + 22, area.width(), 20), Qt::AlignCenter, xLabel_);
    painter.save(); painter.translate(18, area.center().y()); painter.rotate(-90);
    painter.drawText(QRectF(-area.height()/2, -10, area.height(), 20), Qt::AlignCenter, yLabel_); painter.restore();
    QPainterPath path;
    for (qsizetype i = 0; i < x_.size(); ++i) {
        if (!std::isfinite(x_[i]) || !std::isfinite(y_[i])) continue;
        const QPointF point(area.left() + (x_[i] - x0) / (x1 - x0) * area.width(),
                            area.bottom() - (y_[i] - y0) / (y1 - y0) * area.height());
        if (path.elementCount() == 0) path.moveTo(point); else path.lineTo(point);
    }
    painter.save(); painter.setClipRect(area);
    if (filled_) {
        QPainterPath fillPath = path;
        const double baseline = std::clamp(0.0, y0, y1);
        const qreal baseY = area.bottom() - (baseline - y0) / (y1 - y0) * area.height();
        if (path.elementCount() > 0) {
            const auto first = path.elementAt(0), last = path.elementAt(path.elementCount() - 1);
            fillPath.lineTo(QPointF(last.x, baseY)); fillPath.lineTo(QPointF(first.x, baseY)); fillPath.closeSubpath();
        }
        painter.fillPath(fillPath, fillColor_);
    }
    painter.setPen(QPen(lineColor_.isValid() ? lineColor_ : QColor(37, 99, 235), 2)); painter.drawPath(path);
    painter.restore();
}

ProfilesTab::ProfilesTab(QWidget* parent) : QWidget(parent)
{
    auto* layout = new QVBoxLayout(this); auto* controls = new QHBoxLayout;
    source_ = new QComboBox; source_->setObjectName("profilesSource"); source_->addItems({"Database", "Report file"});
    block_ = new QComboBox; block_->setObjectName("profilesBlock");
    profile_ = new QComboBox; profile_->setObjectName("profilesType");
    openReport_ = new QPushButton("Open report..."); openReport_->setObjectName("openProfileReportButton");
    exportImage_ = new QPushButton("Export image..."); exportImage_->setObjectName("exportProfileImageButton");
    controls->addWidget(new QLabel("Source:")); controls->addWidget(source_); controls->addWidget(new QLabel("Experiment:"));
    controls->addWidget(block_, 1); controls->addWidget(new QLabel("Profile:")); controls->addWidget(profile_, 1); controls->addWidget(openReport_); controls->addWidget(exportImage_);
    layout->addLayout(controls);
    auto* options = new QHBoxLayout;
    auto spin = [](const char* name, double value) {
        auto* box = new QDoubleSpinBox; box->setObjectName(name); box->setDecimals(6); box->setRange(-1e12, 1e12); box->setSingleStep(0.1); box->setValue(value); box->setKeyboardTracking(false); return box;
    };
    xLabel_ = new QLineEdit; xLabel_->setObjectName("profileXLabel"); xLabel_->setPlaceholderText("X axis label (automatic)");
    yLabel_ = new QLineEdit; yLabel_->setObjectName("profileYLabel"); yLabel_->setPlaceholderText("Y axis label (automatic)");
    options->addWidget(xLabel_, 1); options->addWidget(yLabel_, 1);
    autoRange_ = new QCheckBox("Auto range"); autoRange_->setObjectName("profileAutoRange"); autoRange_->setChecked(true);
    grid_ = new QCheckBox("Grid"); grid_->setObjectName("profileGrid"); grid_->setChecked(true);
    ticks_ = new QCheckBox("Tick marks"); ticks_->setObjectName("profileTicks"); ticks_->setChecked(true);
    fill_ = new QCheckBox("Fill under curve"); fill_->setObjectName("profileFill");
    options->addWidget(autoRange_); options->addWidget(grid_); options->addWidget(ticks_); options->addWidget(fill_);
    lineColorButton_ = new QPushButton("Line color"); lineColorButton_->setObjectName("profileLineColor");
    fillColorButton_ = new QPushButton("Fill color"); fillColorButton_->setObjectName("profileFillColor");
    options->addWidget(lineColorButton_); options->addWidget(fillColorButton_);
    layout->addLayout(options);
    auto* ranges = new QHBoxLayout;
    xMin_ = spin("profileXMin", 0); xMax_ = spin("profileXMax", 1); yMin_ = spin("profileYMin", 0); yMax_ = spin("profileYMax", 1);
    xTicks_ = new QSpinBox; xTicks_->setObjectName("profileXTickCount"); xTicks_->setRange(2, 12); xTicks_->setValue(5);
    yTicks_ = new QSpinBox; yTicks_->setObjectName("profileYTickCount"); yTicks_->setRange(2, 12); yTicks_->setValue(5);
    ranges->addWidget(new QLabel("X min")); ranges->addWidget(xMin_); ranges->addWidget(new QLabel("X max")); ranges->addWidget(xMax_);
    ranges->addWidget(new QLabel("Y min")); ranges->addWidget(yMin_); ranges->addWidget(new QLabel("Y max")); ranges->addWidget(yMax_);
    ranges->addWidget(new QLabel("X ticks")); ranges->addWidget(xTicks_); ranges->addWidget(new QLabel("Y ticks")); ranges->addWidget(yTicks_);
    layout->addLayout(ranges);
    status_ = new QLabel("Choose a database or open a report file."); status_->setObjectName("profilesStatus"); status_->setWordWrap(true); layout->addWidget(status_);
    plot_ = new ProfilePlot; plot_->setObjectName("profilePlot"); layout->addWidget(plot_, 1);
    connect(source_, &QComboBox::currentIndexChanged, this, [this] { openReport_->setEnabled(source_->currentIndex() == 1); if (source_->currentIndex() == 0) loadDatabase(); else if (!reportFile_.isEmpty()) loadReport(reportFile_); });
    connect(openReport_, &QPushButton::clicked, this, [this] { const auto f = QFileDialog::getOpenFileName(this, "Open profile report", reportFile_, "Reports (*.txt *.csv *.tsv);;All files (*)"); if (!f.isEmpty()) loadReport(f); });
    connect(block_, &QComboBox::currentIndexChanged, this, [this] { updatePlot(); });
    connect(profile_, &QComboBox::currentIndexChanged, this, [this] { updatePlot(); });
    const QList<QWidget*> settingWidgets{autoRange_, grid_, ticks_, fill_, xMin_, xMax_, yMin_, yMax_, xTicks_, yTicks_};
    for (auto* widget : settingWidgets)
        if (auto* check = qobject_cast<QCheckBox*>(widget)) connect(check, &QCheckBox::toggled, this, [this] { updatePlot(); });
        else if (auto* real = qobject_cast<QDoubleSpinBox*>(widget)) connect(real, &QDoubleSpinBox::valueChanged, this, [this] { updatePlot(); });
        else if (auto* integer = qobject_cast<QSpinBox*>(widget)) connect(integer, &QSpinBox::valueChanged, this, [this] { updatePlot(); });
    connect(xLabel_, &QLineEdit::textChanged, this, [this] { updatePlot(); });
    connect(yLabel_, &QLineEdit::textChanged, this, [this] { updatePlot(); });
    connect(lineColorButton_, &QPushButton::clicked, this, [this] { chooseLineColor(); });
    connect(fillColorButton_, &QPushButton::clicked, this, [this] { chooseFillColor(); });
    connect(exportImage_, &QPushButton::clicked, this, [this] {
        const auto filename = QFileDialog::getSaveFileName(this, "Export profile graph", {}, "PNG image (*.png)");
        if (filename.isEmpty()) return;
        QImage image(plot_->size(), QImage::Format_ARGB32_Premultiplied); image.fill(Qt::transparent); plot_->render(&image);
        QString output = filename; if (QFileInfo(output).suffix().isEmpty()) output += ".png";
        if (!image.save(output, "PNG")) status_->setText("Could not export graph image: " + output);
        else status_->setText("Graph image exported: " + output);
    });
    const auto setColorStyle = [](QPushButton* button, const QColor& color) { button->setStyleSheet(QString("background-color: %1").arg(color.name(QColor::HexArgb))); };
    setColorStyle(lineColorButton_, lineColor_); setColorStyle(fillColorButton_, fillColor_);
}

void ProfilesTab::chooseLineColor() {
    const auto color = QColorDialog::getColor(lineColor_, this, "Profile line color");
    if (color.isValid()) { lineColor_ = color; lineColorButton_->setStyleSheet(QString("background-color: %1").arg(color.name())); updatePlot(); }
}

void ProfilesTab::chooseFillColor() {
    auto color = QColorDialog::getColor(fillColor_, this, "Profile fill color", QColorDialog::ShowAlphaChannel);
    if (color.isValid()) { fillColor_ = color; fillColorButton_->setStyleSheet(QString("background-color: %1").arg(color.name(QColor::HexArgb))); updatePlot(); }
}

QJsonObject ProfilesTab::settings() const {
    return {{"block", block_->currentIndex()}, {"profile", profile_->currentText()},
        {"xLabel", xLabel_->text()}, {"yLabel", yLabel_->text()}, {"autoRange", autoRange_->isChecked()},
        {"xMin", xMin_->value()}, {"xMax", xMax_->value()}, {"yMin", yMin_->value()}, {"yMax", yMax_->value()},
        {"xTicks", xTicks_->value()}, {"yTicks", yTicks_->value()}, {"grid", grid_->isChecked()}, {"ticks", ticks_->isChecked()},
        {"fill", fill_->isChecked()}, {"lineColor", lineColor_.name(QColor::HexArgb)}, {"fillColor", fillColor_.name(QColor::HexArgb)}};
}

void ProfilesTab::setSettings(const QJsonObject& value) {
    xLabel_->setText(value.value("xLabel").toString()); yLabel_->setText(value.value("yLabel").toString());
    autoRange_->setChecked(value.value("autoRange").toBool(true));
    xMin_->setValue(value.value("xMin").toDouble(0)); xMax_->setValue(value.value("xMax").toDouble(1));
    yMin_->setValue(value.value("yMin").toDouble(0)); yMax_->setValue(value.value("yMax").toDouble(1));
    xTicks_->setValue(value.value("xTicks").toInt(5)); yTicks_->setValue(value.value("yTicks").toInt(5));
    grid_->setChecked(value.value("grid").toBool(true)); ticks_->setChecked(value.value("ticks").toBool(true)); fill_->setChecked(value.value("fill").toBool(false));
    const QColor line(value.value("lineColor").toString("#ff2563eb")), fill(value.value("fillColor").toString("#482563eb"));
    if (line.isValid()) lineColor_ = line; if (fill.isValid()) fillColor_ = fill;
    lineColorButton_->setStyleSheet(QString("background-color: %1").arg(lineColor_.name(QColor::HexArgb)));
    fillColorButton_->setStyleSheet(QString("background-color: %1").arg(fillColor_.name(QColor::HexArgb)));
    const int block = value.value("block").toInt(-1);
    if (block >= 0 && block < block_->count()) block_->setCurrentIndex(block);
    const int profile = profile_->findText(value.value("profile").toString());
    if (profile >= 0) profile_->setCurrentIndex(profile);
    updatePlot();
}

void ProfilesTab::setDatabase(const QString& filename) { database_ = filename; if (source_->currentIndex() == 0) loadDatabase(); }

void ProfilesTab::loadDatabase()
{
    block_->clear(); profile_->clear();
    if (database_.isEmpty()) { status_->setText("Choose a database file."); return; }
    sqlite3* raw = nullptr;
    const auto encoded = database_.toUtf8();
    if (sqlite3_open_v2(encoded.constData(), &raw, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
        const QString error = raw ? QString::fromUtf8(sqlite3_errmsg(raw)) : "Cannot open database"; if (raw) sqlite3_close(raw); status_->setText(error); return;
    }
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    sqlite3_stmt* statement = nullptr;
    const char* query = "SELECT DISTINCT EXPERIMENT_NAME, WT_PERCENT FROM profile_record_samples ORDER BY EXPERIMENT_NAME, WT_PERCENT";
    if (sqlite3_prepare_v2(db.get(), query, -1, &statement, nullptr) != SQLITE_OK) { status_->setText("No saved model profiles are available in this database."); return; }
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> stmt(statement, sqlite3_finalize);
    while (sqlite3_step(stmt.get()) == SQLITE_ROW) {
        const QString name = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(stmt.get(), 0)));
        const QString wt = QString::fromUtf8(reinterpret_cast<const char*>(sqlite3_column_text(stmt.get(), 1)));
        block_->addItem(wt.isEmpty() ? name : name + " · " + wt, name + QChar(0x1f) + wt);
    }
    if (block_->count() == 0) { status_->setText("No saved model profiles are available in this database."); return; }
    const QStringList labels = {"Experimental", "Numeric model", "Experimental derivative", "Numeric derivative", "Experimental difference", "Numeric difference", "Free dye", "Bound dye", "Total dye", "Unbound beads", "Bound beads", "Total beads", "Error"};
    for (const auto& label : labels) profile_->addItem(label);
    status_->setText("Saved profiles from the selected database."); updatePlot();
}

void ProfilesTab::loadReport(const QString& filename)
{
    reportFile_ = filename;
    try {
        QFile file(filename); if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error(file.errorString().toStdString());
        const auto rows = report_format::parse(QString::fromUtf8(file.readAll()));
        blocks_.clear(); profiles_.clear();
        for (const auto& row : rows) {
            if (row.isEmpty()) continue;
            if (row[0] == "sec")
                blocks_.push_back(QString("Profile block %1 · %2 sec").arg(blocks_.size()+1).arg(row.size()>1?row[1]:QString{}));
            else if (row[0] != "res_time" && row.size() >= 10 && !profiles_.contains(row[8]))
                profiles_.push_back(row[8]);
        }
        block_->clear(); for (const auto& label : blocks_) block_->addItem(label);
        profile_->clear(); for (const auto& label : profiles_) profile_->addItem(label);
        status_->setText("Loaded report: " + filename); source_->setCurrentIndex(1); updatePlot();
    } catch (const std::exception& error) { status_->setText("Cannot load profile report: " + QString::fromUtf8(error.what())); }
}

void ProfilesTab::updatePlot()
{
    const auto setAppearance = [this](const QString& xLabel, const QString& yLabel) {
        plot_->setAppearance(lineColor_, fillColor_, fill_->isChecked(), grid_->isChecked(), ticks_->isChecked(), autoRange_->isChecked(),
            xMin_->value(), xMax_->value(), yMin_->value(), yMax_->value(), xTicks_->value(), yTicks_->value(),
            xLabel_->text().isEmpty() ? xLabel : xLabel_->text(), yLabel_->text().isEmpty() ? yLabel : yLabel_->text());
    };
    if (block_->count() == 0 || profile_->count() == 0) return;
    const QString type = profile_->currentText();
    if (source_->currentIndex() == 1) {
        try {
            QFile file(reportFile_); if (!file.open(QIODevice::ReadOnly)) return;
            const auto rows = report_format::parse(QString::fromUtf8(file.readAll()));
            int b = -1;
            QVector<double> axis, selectedAxis, values;
            bool found = false;
            for (const auto& row : rows) {
                if (row.isEmpty()) continue;
                if (row[0] == "res_time") {
                    axis.clear();
                    for (int i = 9; i < row.size(); ++i) {
                        bool ok = false;
                        const double value = row[i].toDouble(&ok);
                        axis.push_back(ok ? value : std::numeric_limits<double>::quiet_NaN());
                    }
                } else if (row[0] == "sec") {
                    ++b;
                    if (b == block_->currentIndex()) { selectedAxis = axis; values.clear(); }
                } else if (b == block_->currentIndex() && row.size() >= 10 && row[8] == type) {
                    values.clear();
                    for (int i = 9; i < row.size(); ++i) {
                        bool ok = false;
                        const double value = row[i].toDouble(&ok);
                        values.push_back(ok ? value : std::numeric_limits<double>::quiet_NaN());
                    }
                    found = true;
                }
            }
            const int count = std::min(selectedAxis.size(), values.size());
            selectedAxis.resize(count); values.resize(count);
            plot_->setData(block_->currentText() + " · " + type, "X (µm)", type, selectedAxis, values);
            setAppearance("X (µm)", type);
            if (!found) status_->setText("The selected profile is not present in this report block.");
        } catch (...) {}
        return;
    }
    const auto key = block_->currentData().toString().split(QChar(0x1f)); if(key.size()!=2) return;
    static const char* columns[] = {"Experimental","Numeric","Experimental_Derivative","Numeric_Derivative","Experimental_Difference","Numeric_Difference","Free_Dye","Bound_Dye","Total_Dye","Unbound_Beads","Bound_Beads","Total_Beads","Error"};
    sqlite3* db=nullptr; const auto encoded=database_.toUtf8();
    if(sqlite3_open_v2(encoded.constData(),&db,SQLITE_OPEN_READONLY,nullptr)!=SQLITE_OK){if(db)sqlite3_close(db);return;}
    const QByteArray query = QByteArray("SELECT X,") + columns[profile_->currentIndex()] + " FROM profile_record_samples WHERE EXPERIMENT_NAME=? AND WT_PERCENT=? AND RUN_AT_UTC_US=(SELECT MAX(RUN_AT_UTC_US) FROM profile_record_samples WHERE EXPERIMENT_NAME=? AND WT_PERCENT=?) ORDER BY X";
    sqlite3_stmt* st=nullptr; if(sqlite3_prepare_v2(db,query.constData(),-1,&st,nullptr)!=SQLITE_OK){sqlite3_close(db);return;}
    sqlite3_bind_text(st,1,key[0].toUtf8().constData(),-1,SQLITE_TRANSIENT); sqlite3_bind_text(st,2,key[1].toUtf8().constData(),-1,SQLITE_TRANSIENT); sqlite3_bind_text(st,3,key[0].toUtf8().constData(),-1,SQLITE_TRANSIENT); sqlite3_bind_text(st,4,key[1].toUtf8().constData(),-1,SQLITE_TRANSIENT);
    QVector<double> x,y; while(sqlite3_step(st)==SQLITE_ROW){if(sqlite3_column_type(st,1)==SQLITE_NULL)continue;x.push_back(sqlite3_column_double(st,0));y.push_back(sqlite3_column_double(st,1));}
    sqlite3_finalize(st); sqlite3_close(db); plot_->setData(block_->currentText()+" · "+type,"X (µm)",type,x,y);
    setAppearance("X (µm)", type);
}
