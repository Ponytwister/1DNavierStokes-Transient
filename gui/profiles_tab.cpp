#include "profiles_tab.h"

#include "report_tab.h"
#include <QComboBox>
#include <QFileDialog>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QFile>
#include <QPainterPath>
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
    double dx = *xmax - *xmin, dy = *ymax - *ymin;
    if (!std::isfinite(dx) || !std::isfinite(dy)) return;
    if (dx == 0) dx = 1; if (dy == 0) dy = std::max(1.0, std::abs(*ymin) * 0.1);
    const double x0 = *xmin - dx * .02, x1 = *xmax + dx * .02;
    const double y0 = *ymin - dy * .08, y1 = *ymax + dy * .08;
    painter.setPen(QPen(palette().mid().color(), 1)); painter.drawRect(area);
    painter.setPen(palette().text().color());
    painter.drawText(QRectF(8, area.top(), 62, 20), Qt::AlignRight, QString::number(y1, 'g', 4));
    painter.drawText(QRectF(8, area.bottom() - 20, 62, 20), Qt::AlignRight, QString::number(y0, 'g', 4));
    painter.drawText(QRectF(area.left(), area.bottom() + 5, area.width(), 20), Qt::AlignCenter, xLabel_);
    painter.save(); painter.translate(18, area.center().y()); painter.rotate(-90);
    painter.drawText(QRectF(-area.height()/2, -10, area.height(), 20), Qt::AlignCenter, yLabel_); painter.restore();
    QPainterPath path;
    for (qsizetype i = 0; i < x_.size(); ++i) {
        if (!std::isfinite(x_[i]) || !std::isfinite(y_[i])) continue;
        const QPointF point(area.left() + (x_[i] - x0) / (x1 - x0) * area.width(),
                            area.bottom() - (y_[i] - y0) / (y1 - y0) * area.height());
        if (path.elementCount() == 0) path.moveTo(point); else path.lineTo(point);
    }
    painter.setPen(QPen(QColor(37, 99, 235), 2)); painter.drawPath(path);
}

ProfilesTab::ProfilesTab(QWidget* parent) : QWidget(parent)
{
    auto* layout = new QVBoxLayout(this); auto* controls = new QHBoxLayout;
    source_ = new QComboBox; source_->setObjectName("profilesSource"); source_->addItems({"Database", "Report file"});
    block_ = new QComboBox; block_->setObjectName("profilesBlock");
    profile_ = new QComboBox; profile_->setObjectName("profilesType");
    openReport_ = new QPushButton("Open report..."); openReport_->setObjectName("openProfileReportButton");
    controls->addWidget(new QLabel("Source:")); controls->addWidget(source_); controls->addWidget(new QLabel("Experiment:"));
    controls->addWidget(block_, 1); controls->addWidget(new QLabel("Profile:")); controls->addWidget(profile_, 1); controls->addWidget(openReport_);
    layout->addLayout(controls);
    status_ = new QLabel("Choose a database or open a report file."); status_->setObjectName("profilesStatus"); status_->setWordWrap(true); layout->addWidget(status_);
    plot_ = new ProfilePlot; plot_->setObjectName("profilePlot"); layout->addWidget(plot_, 1);
    connect(source_, &QComboBox::currentIndexChanged, this, [this] { openReport_->setEnabled(source_->currentIndex() == 1); if (source_->currentIndex() == 0) loadDatabase(); else if (!reportFile_.isEmpty()) loadReport(reportFile_); });
    connect(openReport_, &QPushButton::clicked, this, [this] { const auto f = QFileDialog::getOpenFileName(this, "Open profile report", reportFile_, "Reports (*.txt *.csv *.tsv);;All files (*)"); if (!f.isEmpty()) loadReport(f); });
    connect(block_, &QComboBox::currentIndexChanged, this, [this] { updatePlot(); });
    connect(profile_, &QComboBox::currentIndexChanged, this, [this] { updatePlot(); });
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
}
