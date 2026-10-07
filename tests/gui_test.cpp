#include <main_window.h>
#include <report_tab.h>
#include <QClipboard>
#include <QTableView>
#include <controls_dialog.h>
#include <model_controls.h>
#include <setup_file.h>
#include <QCheckBox>
#include <QComboBox>
#include <QListWidget>
#include <QToolButton>
#include <QMenu>
#include <QAction>
#include <QMenuBar>
#include <gtest/gtest.h>
#include <QApplication>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QFile>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QPixmap>
#include <QPlainTextEdit>
#include <QTableWidget>
#include <QHeaderView>
#include <QFormLayout>
#include <quantity_units.h>
#include <experiment_selections.h>
#include <QTabWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <thread>

namespace {
TEST(Report, PreservesUnitsNamesPrecisionAndUnequalProfileLengths)
{
    const QString header = "res_time  bind_ratio(p1)  forward_reaction_rate_1  equalibrium_constant_1  dye_conc.  bead_conc.  bead_surface_area  D-A  profile_type  0  10  \r\n";
    const QString units = "sec  mol bead  rate  keq  mg/ml  wt%  mol  deriv  Channel_Width_(um)->  1  \r\n";
    const QString profile = "1e-09  2  3  4  5  6  7  -  run name_Experimental_Profile  1.23456e-12    \r\n";
    auto rows = report_format::parse(header + units + profile + "\r\n" + header);
    ASSERT_EQ(rows.size(), 5);
    EXPECT_EQ(rows[0].size(), 11);
    EXPECT_EQ(rows[1][1], "mol bead");
    EXPECT_EQ(rows[2].size(), 10);
    EXPECT_EQ(rows[2][8], "run name_Experimental_Profile");
    EXPECT_TRUE(rows[3].isEmpty());
    EXPECT_EQ(rows[0], rows[4]);
    EXPECT_TRUE(report_format::tsv(rows).contains("1e-09\t2\t3\t4\t5\t6\t7\t-\trun name_Experimental_Profile\t1.23456e-12\r\n\r\n"));
    EXPECT_THROW(report_format::parse("unrelated text"), std::runtime_error);
    EXPECT_THROW(report_format::parse(header + "1  2\n"), std::runtime_error);
}

TEST(Report, CopiesImportedReportAndClearsStaleDataOnFailure)
{
    QTemporaryDir directory;
    const auto filename = directory.filePath("report.txt");
    QFile file(filename);
    ASSERT_TRUE(file.open(QIODevice::WriteOnly));
    file.write("res_time  b  k  eq  d  c  s  D-A  profile_type  0  \n"); file.close();
    ReportTab tab;
    ASSERT_TRUE(tab.loadFile(filename));
    auto* copy = tab.findChild<QPushButton*>("copyReportButton");
    ASSERT_NE(copy, nullptr); copy->click();
    EXPECT_EQ(QApplication::clipboard()->text(), "res_time\tb\tk\teq\td\tc\ts\tD-A\tprofile_type\t0\r\n");
    EXPECT_FALSE(tab.loadFile(directory.filePath("missing.txt")));
    EXPECT_FALSE(copy->isEnabled());
    EXPECT_EQ(tab.findChild<QTableView*>("reportTable")->model(), nullptr);
}

TEST(Report, OptionalExternalSample)
{
    const auto filename = qEnvironmentVariable("NAVIER_TEST_REPORT");
    if (filename.isEmpty()) GTEST_SKIP() << "Set NAVIER_TEST_REPORT to check an external report without copying it into the repository.";
    QFile file(filename);
    ASSERT_TRUE(file.open(QIODevice::ReadOnly));
    const auto source = QString::fromUtf8(file.readAll());
    const auto rows = report_format::parse(source);
    EXPECT_EQ(report_format::parse(report_format::tsv(rows)), rows);
    ReportTab tab;
    ASSERT_TRUE(tab.loadFile(filename));
    auto* table = tab.findChild<QTableView*>("reportTable");
    ASSERT_NE(table->model(), nullptr);
    EXPECT_EQ(table->model()->rowCount(), rows.size());
    tab.findChild<QPushButton*>("copyReportButton")->click();
    EXPECT_EQ(QApplication::clipboard()->text(), report_format::tsv(rows));
}

TEST(Report, SelectsBlocksTypesAndColumnsWithTheirOwnAxes)
{
    using namespace report_format;
    const QString header = "res_time  b  k  eq  d  c  s  D-A  profile_type  0  10\n";
    const QString units = "sec  b  k  eq  d  c  s  deriv  Channel_Width_(um)->  1\n";
    const auto rows = parse(header + units +
        "1  2  3  4  5  6  7  -  first_Experimental_Profile  0.12345\n\n" + units +
        "2  2  3  4  5  6  7  -  second_Experimental_Profile  9.87654\n"
        "2  2  3  4  5  6  7  -  second_Numeric_Model_Profile  8.76543\n\n" +
        "res_time  b  k  eq  d  c  s  D-A  profile_type  0  20\n" + units +
        "3  2  3  4  5  6  7  -  third_Experimental_Profile  7.65432\n\n");
    Options options;
    options.blocks = {1, 2}; options.types = {"Experimental_Profile"}; options.columns = {0, 8};
    auto selected = select(rows, options);
    ASSERT_EQ(selected.size(), 8);
    EXPECT_EQ(selected[0], QStringList({"res_time", "profile_type", "0", "10"}));
    EXPECT_EQ(selected[1], QStringList({"sec", "Channel_Width_(um)->", "1"}));
    EXPECT_EQ(selected[2], QStringList({"2", "second_Experimental_Profile", "9.87654"}));
    EXPECT_EQ(selected[4], QStringList({"res_time", "profile_type", "0", "20"}));
    EXPECT_EQ(selected[6][1], "third_Experimental_Profile");
    options.headers = options.units = options.blankRows = false;
    options.numberFormat = 'f'; options.decimals = 2;
    selected = select(rows, options);
    ASSERT_EQ(selected.size(), 2);
    EXPECT_EQ(selected[0], QStringList({"2.00", "second_Experimental_Profile", "9.88"}));
    options.blocks.clear();
    EXPECT_TRUE(select(rows, options).isEmpty());
    EXPECT_EQ(rows[6].back(), "8.76543"); // Formatting never mutates source data.
}

TEST(Report, ExportAndClipboardUseSelectionsAndPreview)
{
    QTemporaryDir directory;
    const auto filename = directory.filePath("source.txt");
    const QByteArray source = "res_time  b  k  eq  d  c  s  D-A  profile_type  0\n"
        "sec  b  k  eq  d  c  s  deriv  Channel_Width_(um)->  1\n"
        "1  2  3  4  5  6  7  -  sample_Experimental_Profile  1.23456\n"
        "1  2  3  4  5  6  7  -  sample_Numeric_Model_Profile  9.87654\n\n";
    QFile file(filename); ASSERT_TRUE(file.open(QIODevice::WriteOnly));
    ASSERT_EQ(file.write(source), source.size()); file.close();
    ReportTab tab; ASSERT_TRUE(tab.loadFile(filename));
    auto* types = tab.findChild<QListWidget*>("reportTypes");
    ASSERT_EQ(types->count(), 2); types->item(0)->setCheckState(Qt::Unchecked);
    auto* columns = tab.findChild<QListWidget*>("reportColumns");
    for (int i = 0; i < 8; ++i) columns->item(i)->setCheckState(Qt::Unchecked);
    for (const auto* name : {"reportHeaders", "reportUnits", "reportBlankRows"})
        tab.findChild<QCheckBox*>(name)->setChecked(false);
    const QString expected = "sample_Numeric_Model_Profile\t9.87654\r\n";
    auto* model = tab.findChild<QTableView*>("reportTable")->model();
    ASSERT_EQ(model->rowCount(), 1); ASSERT_EQ(model->columnCount(), 2);
    EXPECT_EQ(model->data(model->index(0, 1)).toString(), "9.87654");
    tab.findChild<QPushButton*>("copyReportButton")->click();
    EXPECT_EQ(QApplication::clipboard()->text(), expected);
    const auto output = directory.filePath("selected.tsv");
    ASSERT_TRUE(tab.saveFile(output));
    QFile exported(output); ASSERT_TRUE(exported.open(QIODevice::ReadOnly));
    EXPECT_EQ(exported.readAll(), expected.toUtf8());
    ASSERT_TRUE(file.open(QIODevice::ReadOnly)); EXPECT_EQ(file.readAll(), source);
    types->item(1)->setCheckState(Qt::Unchecked);
    EXPECT_FALSE(tab.findChild<QPushButton*>("saveReportButton")->isEnabled());
    EXPECT_FALSE(tab.saveFile(directory.filePath("empty.tsv")));
    EXPECT_FALSE(QFile::exists(directory.filePath("empty.tsv")));
}

TEST(Report, SplitMetadataUsesAllSixValuesBeforeFilteringAndKeepsBlocksIndependent)
{
    using namespace report_format;
    const QString header = "res_time  b  k  eq  d  c  s  D-A  profile_type  0  10\n";
    const QString units = "sec  b  k  eq  d  c  s  deriv  Channel_Width_(um)->  1\n";
    const QString prefix = "1  2  3  4  5  6  7  ";
    const QString block = units + prefix + "11.125  Experimental_Derivative  101\n" +
        prefix + "22.25  Numeric_Derivative  102\n" +
        prefix + "33.5  experiment_Experimental_Profile  103\n" +
        prefix + "44.75  experiment_Numeric_Model_Profile  104\n" +
        prefix + "55.875  Experimental_Difference  105\n" +
        prefix + "66.126  Numeric_Difference  106\n" +
        prefix + "-  Free_Dye  107  108\n\n";
    const auto rows = parse(header + block + header + units +
        prefix + "99  Experimental_Derivative  201\n" + prefix + "-  Free_Dye  207\n\n");
    Options options; options.blocks = {0, 1}; options.types = {"Free_Dye"}; options.columns = {7, 8};
    options.metadataLayout = MetadataLayout::separateBlock;
    const auto selected = select(rows, options);
    ASSERT_EQ(selected.size(), 8);
    EXPECT_EQ(selected[0].mid(0, 6), metricNames());
    EXPECT_EQ(selected[1].mid(0, 6), QStringList({"-", "-", "-", "-", "-", "-"}));
    EXPECT_EQ(selected[2], QStringList({"11.125", "22.25", "33.5", "44.75", "55.875", "66.126", "Free_Dye", "107", "108"}));
    EXPECT_EQ(selected[6], QStringList({"99", "-", "-", "-", "-", "-", "Free_Dye", "207"}));
    options.metrics = {1, 5}; options.numberFormat = 'f'; options.decimals = 2;
    EXPECT_EQ(select(rows, options)[2], QStringList({"22.25", "66.13", "Free_Dye", "107.00", "108.00"}));
    options.metadataLayout = MetadataLayout::stacked; options.numberFormat = 0;
    EXPECT_EQ(select(rows, options)[2], QStringList({"-", "Free_Dye", "107", "108"}));
    EXPECT_EQ(rows[2][7], "11.125");
}

TEST(Report, SplitMetadataPreviewClipboardAndFileAgreeForImportedAndGeneratedSources)
{
    QTemporaryDir directory;
    const QString source = "res_time  b  k  eq  d  c  s  D-A  profile_type  0\n"
        "sec  b  k  eq  d  c  s  deriv  Channel_Width_(um)->  1\n"
        "1  2  3  4  5  6  7  12.5  Experimental_Derivative  1\n"
        "1  2  3  4  5  6  7  25  example_Numeric_Model_Profile  2\n"
        "1  2  3  4  5  6  7  -  Free_Dye  3\n\n";
    const auto inputPath = directory.filePath("input.txt");
    QFile input(inputPath); ASSERT_TRUE(input.open(QIODevice::WriteOnly));
    input.write(source.toUtf8()); input.close();
    for (bool generated : {false, true}) {
        ReportTab tab;
        ASSERT_TRUE(generated ? tab.loadGenerated(source, "test run", "report") : tab.loadFile(inputPath));
        auto* layout = tab.findChild<QComboBox*>("reportMetadataLayout");
        layout->setCurrentIndex(1);
        auto* columns = tab.findChild<QListWidget*>("reportColumns");
        ASSERT_EQ(columns->count(), 14);
        for (int i = 0; i < columns->count(); ++i) {
            const auto field = columns->item(i)->data(Qt::UserRole).toInt();
            columns->item(i)->setCheckState(field == 9 || field == 12 || field == 8 ? Qt::Checked : Qt::Unchecked);
        }
        auto* types = tab.findChild<QListWidget*>("reportTypes");
        for (int i = 0; i < types->count(); ++i)
            types->item(i)->setCheckState(types->item(i)->text() == "Free_Dye" ? Qt::Checked : Qt::Unchecked);
        for (const auto* name : {"reportUnits", "reportBlankRows"}) tab.findChild<QCheckBox*>(name)->setChecked(false);
        const QString expected = "exp_DA\tmodel_integral\tprofile_type\t0\r\n12.5\t25\tFree_Dye\t3\r\n";
        tab.findChild<QPushButton*>("copyReportButton")->click();
        EXPECT_EQ(QApplication::clipboard()->text(), expected);
        auto* model = tab.findChild<QTableView*>("reportTable")->model();
        ASSERT_EQ(model->columnCount(), 4); ASSERT_EQ(model->rowCount(), 2);
        EXPECT_EQ(model->data(model->index(1, 1)).toString(), "25");
        const auto outputPath = directory.filePath(generated ? "generated.tsv" : "imported.tsv");
        ASSERT_TRUE(tab.saveFile(outputPath));
        QFile output(outputPath); ASSERT_TRUE(output.open(QIODevice::ReadOnly));
        EXPECT_EQ(output.readAll(), expected.toUtf8());
        layout->setCurrentIndex(0); layout->setCurrentIndex(1);
        tab.findChild<QPushButton*>("copyReportButton")->click();
        EXPECT_EQ(QApplication::clipboard()->text(), expected);
        if (generated) {
            ASSERT_TRUE(tab.loadGenerated(source, "same run", "report"));
            tab.findChild<QPushButton*>("copyReportButton")->click();
            EXPECT_EQ(QApplication::clipboard()->text(), expected);
        }
    }
}

template<class T> T* widget(QWidget& window, const char* name) {
    auto* result = window.findChild<T*>(name);
    if (!result) throw std::runtime_error(name);
    return result;
}
bool until(const std::function<bool()>& ready) {
    QElapsedTimer clock; clock.start();
    while (!ready() && clock.elapsed() < 15000) {
        QApplication::processEvents(); std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return ready();
}
struct Inputs {
    QTemporaryDir directory;
    QString database;
    Inputs() {
        if (!directory.isValid()) throw std::runtime_error("Temporary directory failed");
        database = directory.filePath("inputs.db");
        QFile fixture(QString::fromUtf8(TSENSOR_FIXTURE_DIR) + "/workflow.sql");
        if (!fixture.open(QIODevice::ReadOnly)) throw std::runtime_error("Fixture missing");
        execute(fixture.readAll());
        execute("UPDATE model_controls SET value='3' WHERE criterion='debug_level'");
    }
    double execute(const QByteArray& query) {
        sqlite3* raw = nullptr;
        if (sqlite3_open(database.toUtf8().constData(), &raw) != SQLITE_OK) throw std::runtime_error("Open failed");
        std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
        double value = 0;
        auto callback = [](void* ptr, int, char** values, char**) {
            if (values[0]) *static_cast<double*>(ptr) = std::stod(values[0]);
            return 0;
        };
        if (sqlite3_exec(db.get(), query.constData(), callback, &value, nullptr) != SQLITE_OK)
            throw std::runtime_error(sqlite3_errmsg(db.get()));
        return value;
    }
    void choose(MainWindow& window) {
        widget<QLineEdit>(window, "databasePath")->setText(database);
        widget<QLineEdit>(window, "outputPath")->setText(directory.filePath("reports"));
    }
};

model_controls::Row& control(std::vector<model_controls::Row>& rows, const QString& name) {
    for (auto& row : rows) if (row.name == name) return row;
    throw std::runtime_error("Control missing");
}

TEST(Gui, FileMenuContainsSetupAndPathControls)
{
    MainWindow window;
    auto* file = window.menuBar()->actions().front()->menu();
    ASSERT_NE(file, nullptr);
    EXPECT_TRUE(file->actions().contains(widget<QAction>(window, "openSetupAction")));
    EXPECT_TRUE(file->actions().contains(widget<QAction>(window, "saveSetupAction")));
    EXPECT_TRUE(widget<QAction>(window, "openSetupAction")->isEnabled());
    EXPECT_FALSE(widget<QAction>(window, "saveSetupAction")->isEnabled());
    EXPECT_EQ(window.centralWidget()->findChild<QLineEdit*>("databasePath"), nullptr);
    EXPECT_EQ(window.centralWidget()->findChild<QLineEdit*>("outputPath"), nullptr);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    EXPECT_EQ(tabs->tabText(1), "Model controls");
    EXPECT_FALSE(tabs->isTabEnabled(1));
    Inputs input; input.choose(window);
    EXPECT_TRUE(tabs->isTabEnabled(1));
    EXPECT_TRUE(widget<QAction>(window, "saveSetupAction")->isEnabled());
}

TEST(Gui, GeneratesSelectedReportFromRunBeforeAnyFileExport)
{
    Inputs input;
    MainWindow window; input.choose(window);
    auto* generate = widget<QPushButton>(window, "generateReportButton");
    EXPECT_FALSE(generate->isEnabled());
    widget<QComboBox>(window, "reportNumberFormat")->setCurrentIndex(2);
    widget<QLineEdit>(window, "outputPath")->clear();
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(generate->isEnabled());
    const auto solutions = input.execute("SELECT count(*) FROM solutions");
    generate->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_FALSE(QFile::exists(input.directory.filePath("reports")));
    EXPECT_EQ(input.execute("SELECT count(*) FROM solutions"), solutions);
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    auto* types = widget<QListWidget>(window, "reportTypes");
    QStringList available;
    for (int i = 0; i < types->count(); ++i) {
        available.push_back(types->item(i)->text());
        types->item(i)->setCheckState(types->item(i)->text() == "Bound_Beads_(wt%)" ? Qt::Checked : Qt::Unchecked);
    }
    EXPECT_TRUE(available.contains("Bound_Beads_(wt%)"));
    EXPECT_TRUE(available.contains("Analytical_Zero_(umol)"));
    EXPECT_TRUE(available.contains("Species:FITC_(umol)"));
    // Regeneration preserves report choices and still does not write a raw report.
    generate->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    types = widget<QListWidget>(window, "reportTypes");
    for (int i = 0; i < types->count(); ++i)
        EXPECT_EQ(types->item(i)->checkState() == Qt::Checked, types->item(i)->text() == "Bound_Beads_(wt%)");
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    ASSERT_EQ(tabs->tabText(tabs->currentIndex()), "Report");
    auto* report = static_cast<ReportTab*>(tabs->currentWidget());
    const auto filename = input.directory.filePath("first-report.txt");
    ASSERT_TRUE(report->saveFile(filename));
    QFile output(filename); ASSERT_TRUE(output.open(QIODevice::ReadOnly));
    const auto contents = output.readAll();
    EXPECT_TRUE(contents.contains("Bound_Beads_(wt%)"));
    EXPECT_TRUE(contents.contains("0.000000e+00"));
    EXPECT_FALSE(contents.contains("Analytical_Zero"));
    EXPECT_FALSE(contents.contains("Species:"));
    widget<QLineEdit>(window, "databasePath")->clear();
    EXPECT_FALSE(generate->isEnabled());
    EXPECT_FALSE(widget<QPushButton>(window, "saveReportButton")->isEnabled());
}

TEST(Gui, ExperimentsTabReadsAllRowsWithoutModelWritesAndClearsStaleData)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME) VALUES('second')");
    MainWindow window; input.choose(window);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    EXPECT_EQ(tabs->tabText(2), "Experiments");
    tabs->setCurrentIndex(2);
    auto* table = widget<QTableWidget>(window, "experimentsTable");
    ASSERT_EQ(table->rowCount(), 2);
    EXPECT_EQ(table->item(0, 0)->text(), "second");
    EXPECT_EQ(table->item(0, 1)->text(), "NULL");
    EXPECT_EQ(table->item(1, 0)->text(), "uniform");
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    widget<QLineEdit>(window, "databasePath")->setText(input.directory.filePath("missing.db"));
    EXPECT_EQ(table->rowCount(), 0);
    EXPECT_FALSE(QFile::exists(input.directory.filePath("missing.db")));
    EXPECT_TRUE(widget<QLabel>(window, "experimentsStatus")->text().contains("Cannot load"));
    input.choose(window);
    EXPECT_EQ(table->rowCount(), 2);
}

TEST(Gui, RawProfilesBrowseRefreshAndClearWithoutWrites)
{
    Inputs input;
    MainWindow window; input.choose(window);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    EXPECT_EQ(tabs->tabText(6), "Raw profiles");
    tabs->setCurrentIndex(6);
    auto* table = widget<QTableWidget>(window, "raw_profileTable");
    ASSERT_EQ(table->rowCount(), 4);
    EXPECT_EQ(table->horizontalHeaderItem(0)->text(), "NAME");
    EXPECT_EQ(table->item(0, 0)->text(), "uniform");
    EXPECT_EQ(table->editTriggers(), QAbstractItemView::NoEditTriggers);
    EXPECT_TRUE(widget<QPushButton>(window, "raw_profileAddButton")->isEnabled());
    input.execute("INSERT INTO raw_profile(NAME) VALUES('second')");
    widget<QPushButton>(window, "raw_profileRefreshButton")->click();
    ASSERT_EQ(table->rowCount(), 5);
    EXPECT_EQ(table->item(0, 1)->toolTip(), "SQL NULL (no value)");
    EXPECT_EQ(input.execute("SELECT count(*) FROM solutions"), 1);
    widget<QLineEdit>(window, "databasePath")->setText(input.directory.filePath("missing.db"));
    EXPECT_EQ(table->rowCount(), 0);
    EXPECT_FALSE(QFile::exists(input.directory.filePath("missing.db")));
    EXPECT_TRUE(widget<QLabel>(window, "raw_profileStatus")->text().contains("Cannot load"));
    input.choose(window);
    EXPECT_EQ(table->rowCount(), 5);
}

TEST(Gui, RawProfilesAddValidateAndModifyCompositeIdentity)
{
    Inputs input;
    input.execute("CREATE TABLE edited_profiles(NAME TEXT NOT NULL, WT_PERCENT DOUBLE NOT NULL, "
                  "CHANNEL_LEFT_EDGE INT, CHANNEL_RIGHT_EDGE INT, INTENSITY_ARRAY TEXT, ENTRANCE_CONC ANY, "
                  "INLET_COND_ID INTEGER, OMIT ANY, INDEPENDENT_PARAMETERS_TO_SOLVE_FOR TEXT, "
                  "LEFT_EDGE NUMERIC, WIDTH NUMERIC, EXTRA_DATA BLOB, PRIMARY KEY(NAME,WT_PERCENT)); "
                  "INSERT INTO edited_profiles(NAME,WT_PERCENT,CHANNEL_LEFT_EDGE,CHANNEL_RIGHT_EDGE,INTENSITY_ARRAY,INLET_COND_ID,LEFT_EDGE,WIDTH,OMIT) "
                  "SELECT NAME,WT_PERCENT,CHANNEL_LEFT_EDGE,CHANNEL_RIGHT_EDGE,INTENSITY_ARRAY,INLET_COND_ID,LEFT_EDGE,WIDTH,OMIT FROM raw_profile; "
                  "DROP TABLE raw_profile; ALTER TABLE edited_profiles RENAME TO raw_profile; "
                  "INSERT INTO experiments(NAME) VALUES('second')");
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs"); tabs->setCurrentIndex(6);
    auto* table = widget<QTableWidget>(window, "raw_profileTable");
    auto* add = widget<QPushButton>(window, "raw_profileAddButton");
    auto* modify = widget<QPushButton>(window, "raw_profileModifyButton");
    auto fill = [](QWidget& editor) {
        auto* name = widget<QComboBox>(editor, "NAME");
        EXPECT_FALSE(name->isEditable()); EXPECT_EQ(name->count(), 2);
        name->setCurrentText("uniform");
        widget<QLineEdit>(editor, "WT_PERCENT")->setText("5");
        widget<QLineEdit>(editor, "CHANNEL_LEFT_EDGE")->setText("0");
        widget<QLineEdit>(editor, "CHANNEL_RIGHT_EDGE")->setText("3");
        widget<QPlainTextEdit>(editor, "INTENSITY_ARRAY")->setPlainText("1  2 3");
        widget<QLineEdit>(editor, "INLET_COND_ID")->setText("1");
    };
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        fill(*editor);
        EXPECT_EQ(widget<QComboBox>(*editor, "OMIT")->currentText(), "NULL");
        auto* save = widget<QPushButton>(*editor, "saveRawProfileButton");
        auto* error = widget<QLabel>(*editor, "rawProfileEditorStatus");
        auto* left = widget<QLineEdit>(*editor, "CHANNEL_LEFT_EDGE");
        left->setText("-1"); save->click(); EXPECT_TRUE(error->text().contains("non-negative"));
        left->setText("0.5"); save->click(); EXPECT_TRUE(error->text().contains("integer"));
        left->setText("0");
        auto* right = widget<QLineEdit>(*editor, "CHANNEL_RIGHT_EDGE");
        right->setText("4"); save->click(); EXPECT_TRUE(error->text().contains("no greater"));
        right->setText("3");
        auto* intensity = widget<QPlainTextEdit>(*editor, "INTENSITY_ARRAY");
        intensity->setPlainText("1 nan 3"); save->click(); EXPECT_TRUE(error->text().contains("finite"));
        intensity->setPlainText("1  2 3");
        auto* name = widget<QComboBox>(*editor, "NAME");
        name->setCurrentIndex(-1); save->click(); EXPECT_TRUE(error->text().contains("NAME"));
        name->setCurrentText("uniform");
        auto* choices = widget<QListWidget>(*editor, "rawProfileParameterChoices");
        EXPECT_EQ(choices->count(), model_controls::solvableParameters(input.database).size());
        choices->findItems("left_edge", Qt::MatchExactly).front()->setCheckState(Qt::Checked);
        choices->findItems("width", Qt::MatchExactly).front()->setCheckState(Qt::Checked);
        save->click();
        EXPECT_FALSE(editor->isVisible()) << error->text().toStdString();
        if (editor->isVisible()) widget<QPushButton>(*editor, "cancelRawProfileButton")->click();
    });
    add->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile"), 5);
    EXPECT_EQ(input.execute("SELECT OMIT IS NULL FROM raw_profile WHERE WT_PERCENT='5.0'"), 1);
    EXPECT_EQ(input.execute("SELECT INTENSITY_ARRAY=('1'||char(9)||'2'||char(9)||'3') FROM raw_profile WHERE WT_PERCENT='5.0'"), 1);
    EXPECT_EQ(input.execute("SELECT INDEPENDENT_PARAMETERS_TO_SOLVE_FOR='left_edge width' FROM raw_profile WHERE WT_PERCENT='5.0'"), 1);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr); fill(*editor);
        widget<QPushButton>(*editor, "saveRawProfileButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "rawProfileEditorStatus")->text().contains("already exists"));
        widget<QPushButton>(*editor, "cancelRawProfileButton")->click();
    });
    add->click(); EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile"), 5);
    input.execute("UPDATE raw_profile SET EXTRA_DATA=x'0011' WHERE WT_PERCENT='5.0'");
    widget<QPushButton>(window, "raw_profileRefreshButton")->click();
    for (int row = 0; row < table->rowCount(); ++row)
        if (table->item(row, 1)->text().toDouble() == 5) { table->setCurrentCell(row, 0); table->selectRow(row); }
    ASSERT_TRUE(modify->isEnabled());
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QComboBox>(*editor, "NAME")->setCurrentText("second");
        widget<QLineEdit>(*editor, "WT_PERCENT")->setText("6");
        widget<QLineEdit>(*editor, "CHANNEL_RIGHT_EDGE")->setText("2"); // Below sample count is allowed.
        widget<QComboBox>(*editor, "OMIT")->setCurrentText("true");
        auto* choices = widget<QListWidget>(*editor, "rawProfileParameterChoices");
        for (int i = 0; i < choices->count(); ++i) choices->item(i)->setCheckState(Qt::Unchecked);
        widget<QPushButton>(*editor, "saveRawProfileButton")->click();
        EXPECT_FALSE(editor->isVisible()) << widget<QLabel>(*editor, "rawProfileEditorStatus")->text().toStdString();
        if (editor->isVisible()) widget<QPushButton>(*editor, "cancelRawProfileButton")->click();
    });
    modify->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile WHERE NAME='uniform'"), 4);
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile WHERE NAME='second' AND WT_PERCENT='6.0' AND CHANNEL_RIGHT_EDGE=2 AND OMIT='true' AND INDEPENDENT_PARAMETERS_TO_SOLVE_FOR IS NULL AND EXTRA_DATA=x'0011'"), 1);
    EXPECT_EQ(input.execute("SELECT count(*) FROM solutions"), 1);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
}

TEST(Gui, RawProfilesCancelConflictsAndMissingExperimentDoNotOverwrite)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(6);
    auto* table = widget<QTableWidget>(window, "raw_profileTable");
    table->setCurrentCell(0, 0); table->selectRow(0);
    const auto weight = table->item(0, 1)->text();
    auto* modify = widget<QPushButton>(window, "raw_profileModifyButton");
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "CHANNEL_RIGHT_EDGE")->setText("8");
        widget<QPushButton>(*editor, "cancelRawProfileButton")->click();
    });
    modify->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile WHERE CHANNEL_RIGHT_EDGE=9"), 4);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        input.execute(("UPDATE raw_profile SET OMIT='true' WHERE WT_PERCENT='" + weight + "'").toUtf8());
        widget<QLineEdit>(*editor, "CHANNEL_RIGHT_EDGE")->setText("8");
        widget<QPushButton>(*editor, "saveRawProfileButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "rawProfileEditorStatus")->text().contains("changed in the database"));
        widget<QPushButton>(*editor, "cancelRawProfileButton")->click();
    });
    modify->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile WHERE CHANNEL_RIGHT_EDGE=9"), 4);
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile WHERE OMIT='true'"), 1);
    widget<QPushButton>(window, "raw_profileRefreshButton")->click();
    table->setCurrentCell(0, 0); table->selectRow(0);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        input.execute("DELETE FROM experiments");
        widget<QLineEdit>(*editor, "CHANNEL_RIGHT_EDGE")->setText("8");
        widget<QPushButton>(*editor, "saveRawProfileButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "rawProfileEditorStatus")->text().contains("experiment is missing"));
        widget<QPushButton>(*editor, "cancelRawProfileButton")->click();
    });
    modify->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM raw_profile WHERE CHANNEL_RIGHT_EDGE=9"), 4);
}

TEST(Gui, VariableFilterCombinesSelectionsAndExcludesOmittedProfiles)
{
    Inputs input;
    input.execute("ALTER TABLE experiments ADD COLUMN PARAMETERS_TO_SOLVE_FOR TEXT;"
                  "ALTER TABLE raw_profile ADD COLUMN INDEPENDENT_PARAMETERS_TO_SOLVE_FOR TEXT;"
                  "UPDATE experiments SET PARAMETERS_TO_SOLVE_FOR=' p1  keq1 ';"
                  "INSERT INTO experiments(NAME,PARAMETERS_TO_SOLVE_FOR) VALUES('uniform-extra','kon1');"
                  "UPDATE raw_profile SET INDEPENDENT_PARAMETERS_TO_SOLVE_FOR='width' WHERE INLET_COND_ID=1;"
                  "UPDATE raw_profile SET INDEPENDENT_PARAMETERS_TO_SOLVE_FOR='left_edge',OMIT=' TrUe ' WHERE INLET_COND_ID=2;"
                  "UPDATE raw_profile SET INDEPENDENT_PARAMETERS_TO_SOLVE_FOR='QE1',OMIT=NULL WHERE INLET_COND_ID=3;"
                  "INSERT INTO raw_profile(NAME,INDEPENDENT_PARAMETERS_TO_SOLVE_FOR) VALUES('uniform-extra','kon1');"
                  "INSERT INTO alglib_input VALUES('p1',1,0,2,1),('width',1,0,2,1),('left_edge',1,0,2,1),"
                  "('QE1',1,0,2,1),('kon1',1,0,2,1)");
    MainWindow window; input.choose(window);
    auto* table = widget<QTableWidget>(window, "alglib_inputTable");
    auto* toggle = widget<QCheckBox>(window, "alglib_inputSelectedOnly");
    auto visible = [&] {
        QStringList names;
        for (int row = 0; row < table->rowCount(); ++row)
            if (!table->isRowHidden(row)) names << table->item(row, 0)->text();
        return names;
    };
    EXPECT_TRUE(toggle->isChecked());
    EXPECT_EQ(visible(), (QStringList{"QE1", "keq1", "p1", "width"}));
    toggle->setChecked(false);
    EXPECT_EQ(visible().size(), 6);
    toggle->setChecked(true);
    input.execute("UPDATE raw_profile SET OMIT='false' WHERE INLET_COND_ID=2");
    widget<QPushButton>(window, "alglib_inputRefreshButton")->click();
    EXPECT_EQ(visible(), (QStringList{"QE1", "keq1", "left_edge", "p1", "width"}));
    auto controls = model_controls::load(input.database);
    control(controls.rows, "universal_solve_for").value = std::nullopt;
    control(controls.rows, "experiment_name").value = "uniform-extra";
    EXPECT_EQ(model_controls::selectedVariables(input.database, controls), QStringList{"kon1"});
    control(controls.rows, "experiment_name").value = "";
    EXPECT_TRUE(model_controls::selectedVariables(input.database, controls).isEmpty());
    input.execute("DROP TABLE raw_profile");
    widget<QPushButton>(window, "alglib_inputRefreshButton")->click();
    EXPECT_TRUE(visible().isEmpty());
    EXPECT_TRUE(widget<QLabel>(window, "alglib_inputFilterStatus")->text().contains("Cannot load"));
    toggle->setChecked(false);
    EXPECT_EQ(visible().size(), 6);
}

TEST(Gui, ExperimentFiltersFollowAppliedSelectionAndIndependentToggles)
{
    Inputs input;
    // Acceptance now starts a run, so the selected experiment needs complete inputs.
    input.execute("INSERT INTO experiments SELECT 'second',LOW_REF_LEFT,LOW_REF_RIGHT,HIGH_REF_LEFT,HIGH_REF_RIGHT,"
                  "DEFAULT_NORMALIZATION,SPECIES,REACTIONS,ENTRANCE_FLOWRATE,EDGES,SPECIE_INLET_CONC_UNITS,"
                  "SPECIE_MODEL_CONC_UNITS,WIDTH FROM experiments WHERE NAME='uniform';"
                  "INSERT INTO experiments(NAME) VALUES('uniform-extra');"
                  "INSERT INTO raw_profile SELECT 'second',WT_PERCENT,CHANNEL_LEFT_EDGE,CHANNEL_RIGHT_EDGE,"
                  "INTENSITY_ARRAY,INLET_COND_ID,LEFT_EDGE,WIDTH,OMIT FROM raw_profile WHERE NAME='uniform' LIMIT 1;"
                  "INSERT INTO raw_profile(NAME) VALUES('uniform-extra'),(NULL)");
    const auto original = model_controls::load(input.database);
    MainWindow window; input.choose(window);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    auto visibleNames = [](QTableWidget* table) {
        QStringList names;
        for (int row = 0; row < table->rowCount(); ++row)
            if (!table->isRowHidden(row)) names << table->item(row, 0)->text();
        return names;
    };
    auto* experiments = widget<QTableWidget>(window, "experimentsTable");
    auto* profiles = widget<QTableWidget>(window, "raw_profileTable");
    auto* experimentFilter = widget<QCheckBox>(window, "experimentsSelectedOnly");
    auto* profileFilter = widget<QCheckBox>(window, "raw_profileSelectedOnly");
    EXPECT_TRUE(experimentFilter->isChecked()); EXPECT_TRUE(profileFilter->isChecked());
    tabs->setCurrentIndex(2);
    EXPECT_EQ(visibleNames(experiments), QStringList{"uniform"});
    tabs->setCurrentIndex(6);
    EXPECT_EQ(visibleNames(profiles), (QStringList{"uniform", "uniform", "uniform", "uniform"}));
    profileFilter->setChecked(false);
    EXPECT_EQ(visibleNames(profiles).size(), 7);
    EXPECT_EQ(visibleNames(experiments), QStringList{"uniform"});
    profileFilter->setChecked(true);

    tabs->setCurrentIndex(1);
    ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "useControlsButton")->isEnabled(); }));
    auto* choices = widget<QListWidget>(window, "experimentChoices");
    for (int i = 0; i < choices->count(); ++i)
        choices->item(i)->setCheckState(choices->item(i)->text() == "second" ? Qt::Checked : Qt::Unchecked);
    QTimer::singleShot(0, [] {
        auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
        ASSERT_NE(box, nullptr); box->button(QMessageBox::Yes)->click();
    });
    widget<QPushButton>(window, "runButton")->click();
    ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "runButton")->isEnabled(); }));
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_EQ(visibleNames(experiments), QStringList{"second"});
    EXPECT_EQ(visibleNames(profiles), QStringList{"second"});
    tabs->setCurrentIndex(6);
    widget<QPushButton>(window, "raw_profileRefreshButton")->click();
    EXPECT_EQ(visibleNames(profiles), QStringList{"second"});
    tabs->setCurrentIndex(2);
    widget<QPushButton>(window, "experimentsRefreshButton")->click();
    EXPECT_EQ(visibleNames(experiments), QStringList{"second"});
    experimentFilter->setChecked(false);
    EXPECT_EQ(visibleNames(experiments).size(), 3);
    EXPECT_EQ(visibleNames(profiles), QStringList{"second"});
    experimentFilter->setChecked(true);

    Inputs other;
    other.execute("UPDATE model_controls SET value='' WHERE criterion='experiment_name'");
    other.choose(window);
    EXPECT_TRUE(visibleNames(experiments).isEmpty());
    tabs->setCurrentIndex(6);
    EXPECT_TRUE(visibleNames(profiles).isEmpty());
    other.execute("DROP TABLE model_controls");
    widget<QPushButton>(window, "raw_profileRefreshButton")->click();
    EXPECT_TRUE(widget<QLabel>(window, "raw_profileFilterStatus")->text().contains("Cannot load"));
    profileFilter->setChecked(false);
    EXPECT_EQ(visibleNames(profiles).size(), 4);
}

TEST(Gui, ReferenceFiltersUseSelectedExperimentListsAndRefresh)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME,REACTIONS,SPECIES) VALUES"
                  "('second',' ExtraReaction  FITC_40nm_1 ',' ExtraSpecies FITC '),"
                  "('unused','UnusedReaction','UnusedSpecies'); "
                  "INSERT INTO reactions(REACTION_NAME) VALUES('ExtraReaction'),('UnusedReaction'),('FITC_40nm_1-extra'); "
                  "INSERT INTO species(SPECIES_NAME) VALUES('ExtraSpecies'),('UnusedSpecies'),('FITC-extra')");
    const auto original = model_controls::load(input.database);
    MainWindow window; input.choose(window);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    auto visibleNames = [](QTableWidget* table) {
        QStringList names;
        for (int row = 0; row < table->rowCount(); ++row)
            if (!table->isRowHidden(row)) names << table->item(row, 0)->text();
        return names;
    };
    auto* reactions = widget<QTableWidget>(window, "reactionsTable");
    auto* species = widget<QTableWidget>(window, "speciesTable");
    auto* reactionFilter = widget<QCheckBox>(window, "reactionsSelectedOnly");
    auto* speciesFilter = widget<QCheckBox>(window, "speciesSelectedOnly");
    EXPECT_TRUE(reactionFilter->isChecked()); EXPECT_TRUE(speciesFilter->isChecked());
    tabs->setCurrentIndex(3);
    EXPECT_EQ(visibleNames(reactions), QStringList{"FITC_40nm_1"});
    tabs->setCurrentIndex(4);
    EXPECT_EQ(visibleNames(species), (QStringList{"40nm_Bound_Dye_1", "FITC", "PS_40nm"}));
    reactionFilter->setChecked(false);
    EXPECT_EQ(visibleNames(reactions).size(), 4);
    EXPECT_EQ(visibleNames(species).size(), 3);
    reactions->setCurrentCell(0, 0); reactions->selectRow(0);
    EXPECT_TRUE(widget<QPushButton>(window, "reactionsModifyButton")->isEnabled());
    reactionFilter->setChecked(true);
    EXPECT_FALSE(widget<QPushButton>(window, "reactionsModifyButton")->isEnabled());

    tabs->setCurrentIndex(1);
    ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "useControlsButton")->isEnabled(); }));
    auto* choices = widget<QListWidget>(window, "experimentChoices");
    choices->findItems("second", Qt::MatchExactly).front()->setCheckState(Qt::Checked);
    QTimer::singleShot(0, [] {
        auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
        ASSERT_NE(box, nullptr); box->button(QMessageBox::Yes)->click();
    });
    widget<QPushButton>(window, "runButton")->click();
    ASSERT_TRUE(until([&] { return widget<QPushButton>(window, "runButton")->isEnabled(); }));
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_EQ(visibleNames(reactions), (QStringList{"ExtraReaction", "FITC_40nm_1"}));
    EXPECT_EQ(visibleNames(species), (QStringList{"40nm_Bound_Dye_1", "ExtraSpecies", "FITC", "PS_40nm"}));
    input.execute("UPDATE experiments SET REACTIONS=NULL,SPECIES='' WHERE NAME='second'");
    tabs->setCurrentIndex(3);
    widget<QPushButton>(window, "reactionsRefreshButton")->click();
    EXPECT_EQ(visibleNames(reactions), QStringList{"FITC_40nm_1"});
    tabs->setCurrentIndex(4);
    widget<QPushButton>(window, "speciesRefreshButton")->click();
    EXPECT_EQ(visibleNames(species).size(), 3);
    speciesFilter->setChecked(false);
    EXPECT_EQ(visibleNames(species).size(), 6);
    EXPECT_EQ(visibleNames(reactions).size(), 1);
    speciesFilter->setChecked(true);
    EXPECT_EQ(input.execute("SELECT count(*) FROM solutions"), 1);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);

    Inputs other;
    other.execute("UPDATE model_controls SET value='missing' WHERE criterion='experiment_name'");
    other.choose(window);
    EXPECT_TRUE(visibleNames(species).isEmpty());
    tabs->setCurrentIndex(3);
    EXPECT_TRUE(visibleNames(reactions).isEmpty());
    other.execute("DROP TABLE experiments");
    widget<QPushButton>(window, "reactionsRefreshButton")->click();
    EXPECT_TRUE(widget<QLabel>(window, "reactionsFilterStatus")->text().contains("Cannot load"));
    reactionFilter->setChecked(false);
    EXPECT_EQ(visibleNames(reactions).size(), 1);
    tabs->setCurrentIndex(4);
    widget<QPushButton>(window, "speciesRefreshButton")->click();
    EXPECT_TRUE(widget<QLabel>(window, "speciesFilterStatus")->text().contains("Cannot load"));
    speciesFilter->setChecked(false);
    EXPECT_EQ(visibleNames(species).size(), 3);
}

TEST(Gui, AllDatabaseTabsArePreloadedAndRetainedDuringNavigation)
{
    Inputs input;
    MainWindow window; input.choose(window);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    ASSERT_EQ(tabs->currentIndex(), 0);
    const char* tables[] = {"experimentsTable", "reactionsTable", "speciesTable", "alglib_inputTable", "raw_profileTable"};
    // Remove the disposable database after loading: even first visits must use memory.
    ASSERT_TRUE(QFile::remove(input.database));
    for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i < 5; ++i) {
            auto* table = widget<QTableWidget>(window, tables[i]);
            ASSERT_GT(table->rowCount(), 0);
            auto* first = table->item(0, 0);
            table->selectRow(0);
            tabs->setCurrentIndex(i + 2);
            EXPECT_EQ(table->item(0, 0), first);
            EXPECT_EQ(table->currentRow(), 0);
            ASSERT_GT(table->selectedItems().size(), 0);
        }
    }
    widget<QLineEdit>(window, "databasePath")->clear();
    for (const auto* name : tables) EXPECT_EQ(widget<QTableWidget>(window, name)->rowCount(), 0);
}

TEST(Gui, ReferenceTabsBrowseRefreshAndSwitchDatabasesWithoutWrites)
{
    Inputs input;
    input.execute("INSERT INTO reactions(REACTION_NAME) VALUES('!Unassigned')");
    input.execute("INSERT INTO species(SPECIES_NAME) VALUES('!Unassigned')");
    MainWindow window; input.choose(window);
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    struct Spec { int index; const char* title; const char* table; const char* refresh; const char* firstColumn; int rows; };
    for (const auto& spec : {
        Spec{3, "Reactions", "reactionsTable", "reactionsRefreshButton", "REACTION_NAME", 2},
        Spec{4, "Species", "speciesTable", "speciesRefreshButton", "SPECIES_NAME", 4}}) {
        EXPECT_EQ(tabs->tabText(spec.index), spec.title);
        tabs->setCurrentIndex(spec.index);
        auto* table = widget<QTableWidget>(window, spec.table);
        ASSERT_EQ(table->rowCount(), spec.rows);
        EXPECT_EQ(table->horizontalHeaderItem(0)->text(), spec.firstColumn);
        EXPECT_EQ(table->item(0, 0)->text(), "!Unassigned");
        EXPECT_EQ(table->item(0, 1)->text(), "NULL");
        EXPECT_EQ(table->item(0, 1)->toolTip(), "SQL NULL (no value)");
        EXPECT_EQ(table->editTriggers(), QAbstractItemView::NoEditTriggers);
        for (int row = 0; row < table->rowCount(); ++row)
            for (int col = 0; col < table->columnCount(); ++col)
                EXPECT_FALSE(table->item(row, col)->flags() & Qt::ItemIsEditable);
        // Refresh preserves model controls and causes no model initialization.
        widget<QPushButton>(window, spec.refresh)->click();
        EXPECT_EQ(table->rowCount(), spec.rows);
    }
    input.execute("UPDATE reactions SET Ks='2 3' WHERE REACTION_NAME='!Unassigned'");
    tabs->setCurrentIndex(3);
    auto* reactions = widget<QTableWidget>(window, "reactionsTable");
    EXPECT_EQ(reactions->item(0, 3)->text(), "NULL"); // Navigation retains the cached snapshot.
    input.execute("UPDATE reactions SET Ks='4 5' WHERE REACTION_NAME='!Unassigned'");
    widget<QPushButton>(window, "reactionsRefreshButton")->click();
    EXPECT_EQ(reactions->item(0, 3)->text(), "4 5");
    EXPECT_EQ(input.execute("SELECT count(*) FROM solutions"), 1);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    EXPECT_DOUBLE_EQ(input.execute("SELECT DIFFUSION_RATE FROM species WHERE SPECIES_NAME='FITC'"), 4.9e-10);
    Inputs other;
    other.choose(window);
    ASSERT_EQ(reactions->rowCount(), 1);
    EXPECT_EQ(reactions->item(0, 0)->text(), "FITC_40nm_1");
    auto* species = widget<QTableWidget>(window, "speciesTable");
    EXPECT_EQ(species->rowCount(), 3); // Hidden tabs are preloaded for the new database.
    tabs->setCurrentIndex(4);
    EXPECT_EQ(species->rowCount(), 3);
}

TEST(Gui, ReferenceTabsClearMissingInputsAndHandleEmptyTables)
{
    Inputs input;
    MainWindow window;
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    tabs->setCurrentIndex(3);
    EXPECT_FALSE(widget<QPushButton>(window, "reactionsRefreshButton")->isEnabled());
    input.choose(window);
    auto* reactions = widget<QTableWidget>(window, "reactionsTable");
    ASSERT_EQ(reactions->rowCount(), 1);
    input.execute("DROP TABLE reactions");
    widget<QPushButton>(window, "reactionsRefreshButton")->click();
    EXPECT_EQ(reactions->rowCount(), 0);
    EXPECT_EQ(reactions->columnCount(), 0);
    EXPECT_TRUE(widget<QLabel>(window, "reactionsStatus")->text().contains("Cannot load reactions"));
    tabs->setCurrentIndex(4);
    auto* species = widget<QTableWidget>(window, "speciesTable");
    ASSERT_EQ(species->rowCount(), 3);
    input.execute("DELETE FROM species");
    widget<QPushButton>(window, "speciesRefreshButton")->click();
    EXPECT_EQ(species->rowCount(), 0);
    EXPECT_EQ(species->columnCount(), 7);
    const auto missing = input.directory.filePath("missing-reference.db");
    widget<QLineEdit>(window, "databasePath")->setText(missing);
    EXPECT_EQ(species->rowCount(), 0);
    EXPECT_EQ(species->columnCount(), 0);
    EXPECT_FALSE(QFile::exists(missing));
    EXPECT_TRUE(widget<QLabel>(window, "speciesStatus")->text().contains("Cannot load species"));
}

TEST(Gui, AlglibBrowseAddValidateAndModify)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs"); tabs->setCurrentIndex(5);
    EXPECT_EQ(tabs->tabText(5), "Variables");
    auto* table = widget<QTableWidget>(window, "alglib_inputTable");
    ASSERT_EQ(table->rowCount(), 1); ASSERT_EQ(table->columnCount(), 5);
    EXPECT_EQ(table->item(0, 0)->text(), "keq1");
    EXPECT_EQ(table->horizontalHeaderItem(1)->text(), "INITIAL VALUE");
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "VARIABLE")->setText("new_variable");
        widget<QLineEdit>(*editor, "INITIAL VALUE")->setText("1");
        widget<QLineEdit>(*editor, "LOWER BOUND")->setText("0");
        widget<QLineEdit>(*editor, "UPPER BOUND")->setText("2");
        auto* scale = widget<QLineEdit>(*editor, "SCALE");
        auto* save = widget<QPushButton>(*editor, "saveReferenceButton");
        scale->setText("nan"); save->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("finite"));
        scale->setText("0"); save->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("nonzero"));
        scale->setText("1");
        widget<QLineEdit>(*editor, "LOWER BOUND")->setText("3"); save->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("must not exceed"));
        widget<QLineEdit>(*editor, "LOWER BOUND")->setText("0");
        widget<QLineEdit>(*editor, "INITIAL VALUE")->setText("3"); save->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("within the bounds"));
        widget<QLineEdit>(*editor, "INITIAL VALUE")->setText("1"); save->click();
    });
    widget<QPushButton>(window, "alglib_inputAddButton")->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM alglib_input"), 2);
    table->setCurrentCell(0, 0); table->selectRow(0);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        EXPECT_TRUE(widget<QLineEdit>(*editor, "VARIABLE")->isReadOnly());
        widget<QLineEdit>(*editor, "INITIAL VALUE")->setText("1.25");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
    });
    widget<QPushButton>(window, "alglib_inputModifyButton")->click();
    EXPECT_DOUBLE_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input WHERE VARIABLE='keq1'"), 1.25);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    input.execute("DROP TABLE alglib_input");
    widget<QPushButton>(window, "alglib_inputRefreshButton")->click();
    EXPECT_EQ(table->rowCount(), 0);
    EXPECT_FALSE(widget<QPushButton>(window, "alglib_inputAddButton")->isEnabled());
    EXPECT_TRUE(widget<QLabel>(window, "alglib_inputStatus")->text().contains("Cannot load alglib_input"));
}

TEST(Gui, SpeciesAddModifyCancelAndConcurrentChange)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs"); tabs->setCurrentIndex(4);
    auto* add = widget<QPushButton>(window, "speciesAddButton");
    auto* modify = widget<QPushButton>(window, "speciesModifyButton");
    auto* table = widget<QTableWidget>(window, "speciesTable");
    EXPECT_TRUE(add->isEnabled()); EXPECT_FALSE(modify->isEnabled());
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "SPECIES_NAME")->setText("NewSpecies");
        widget<QLineEdit>(*editor, "SPECIES_TYPE")->setText("molecule");
        widget<QCheckBox>(*editor, "QENull")->setChecked(false);
        widget<QLineEdit>(*editor, "QE")->setText("nan");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(editor->isVisible());
        widget<QLineEdit>(*editor, "QE")->setText("2.5");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
    });
    add->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM species WHERE SPECIES_NAME='NewSpecies'"), 1);
    EXPECT_DOUBLE_EQ(input.execute("SELECT QE FROM species WHERE SPECIES_NAME='NewSpecies'"), 2.5);
    EXPECT_EQ(input.execute("SELECT DIFFUSION_RATE IS NULL FROM species WHERE SPECIES_NAME='NewSpecies'"), 1);
    widget<QCheckBox>(window, "speciesSelectedOnly")->setChecked(false);
    auto selectNew = [&] {
        for (int row = 0; row < table->rowCount(); ++row)
            if (table->item(row, 0)->text() == "NewSpecies") { table->setCurrentCell(row, 0); table->selectRow(row); return; }
        ADD_FAILURE() << "New species missing";
    };
    selectNew(); EXPECT_TRUE(modify->isEnabled());
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        EXPECT_TRUE(widget<QLineEdit>(*editor, "SPECIES_NAME")->isReadOnly());
        widget<QLineEdit>(*editor, "QE")->setText("4");
        widget<QPushButton>(*editor, "cancelReferenceButton")->click();
    });
    modify->click();
    EXPECT_DOUBLE_EQ(input.execute("SELECT QE FROM species WHERE SPECIES_NAME='NewSpecies'"), 2.5);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "QE")->setText("4");
        input.execute("UPDATE species SET QE=3 WHERE SPECIES_NAME='NewSpecies'");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("changed in the database"));
        widget<QPushButton>(*editor, "cancelReferenceButton")->click();
    });
    modify->click();
    EXPECT_DOUBLE_EQ(input.execute("SELECT QE FROM species WHERE SPECIES_NAME='NewSpecies'"), 3);
    widget<QPushButton>(window, "speciesRefreshButton")->click(); selectNew();
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "QE")->setText("4");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
    });
    modify->click();
    EXPECT_DOUBLE_EQ(input.execute("SELECT QE FROM species WHERE SPECIES_NAME='NewSpecies'"), 4);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
}

TEST(Gui, ReactionsAddValidateDuplicatesAndModify)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(3);
    auto* add = widget<QPushButton>(window, "reactionsAddButton");
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "REACTION_NAME")->setText("NewReaction");
        widget<QLineEdit>(*editor, "SPECIES")->setText("FITC Missing");
        widget<QLineEdit>(*editor, "COEFFICIENTS")->setText("-1 1");
        widget<QLineEdit>(*editor, "Ks")->setText("0 1 2");
        widget<QLineEdit>(*editor, "EXPONENTS")->setText("1 1");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("incorrect number"));
        widget<QLineEdit>(*editor, "Ks")->setText("0 1");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("Unknown or ambiguous species"));
        widget<QLineEdit>(*editor, "SPECIES")->setText("FITC PS_40nm");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
    });
    add->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM reactions WHERE REACTION_NAME='NewReaction'"), 1);
    widget<QCheckBox>(window, "reactionsSelectedOnly")->setChecked(false);
    auto* table = widget<QTableWidget>(window, "reactionsTable");
    for (int row = 0; row < table->rowCount(); ++row)
        if (table->item(row, 0)->text() == "NewReaction") { table->setCurrentCell(row, 0); table->selectRow(row); }
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        EXPECT_TRUE(widget<QLineEdit>(*editor, "REACTION_NAME")->isReadOnly());
        widget<QLineEdit>(*editor, "Ks")->setText("2 3");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
    });
    widget<QPushButton>(window, "reactionsModifyButton")->click();
    EXPECT_EQ(input.execute("SELECT Ks='2 3' FROM reactions WHERE REACTION_NAME='NewReaction'"), 1);
    QTimer::singleShot(0, &window, [&] {
        auto* editor = QApplication::activeModalWidget(); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "REACTION_NAME")->setText("NewReaction");
        widget<QLineEdit>(*editor, "SPECIES")->setText("FITC");
        widget<QLineEdit>(*editor, "COEFFICIENTS")->setText("-1");
        widget<QLineEdit>(*editor, "Ks")->setText("0 1");
        widget<QLineEdit>(*editor, "EXPONENTS")->setText("1");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("already exists"));
        widget<QPushButton>(*editor, "cancelReferenceButton")->click();
    });
    add->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM reactions WHERE REACTION_NAME='NewReaction'"), 1);
    EXPECT_EQ(input.execute("SELECT count(*) FROM solutions"), 1);
}

TEST(Gui, ExperimentEditorUsesOrderedSelectionsAndConvertsDisplayUnits)
{
    Inputs input;
    QFile migration(QString::fromUtf8(TSENSOR_FIXTURE_DIR) + "/../../migrations/001_channel_dimensions.sql");
    ASSERT_TRUE(migration.open(QIODevice::ReadOnly)); input.execute(migration.readAll());
    input.execute("ALTER TABLE experiments ADD COLUMN PARAMETERS_TO_SOLVE_FOR TEXT");
    MainWindow window; input.choose(window); window.show();
    widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(2);
    auto* table = widget<QTableWidget>(window, "experimentsTable"); table->selectRow(0);
    EXPECT_EQ(table->horizontalHeaderItem(table->horizontalHeader()->logicalIndex(1))->text(), "Parameters to solve for");
    QTimer::singleShot(0, [&] {
        auto* editor = window.findChild<QDialog*>("experimentsEditor"); ASSERT_NE(editor, nullptr);
        EXPECT_TRUE(editor->findChildren<QCheckBox*>().isEmpty());
        EXPECT_NE(editor->findChild<QToolButton*>("SPECIES"), nullptr);
        EXPECT_NE(editor->findChild<QToolButton*>("REACTIONS"), nullptr);
        auto* inlet = widget<QComboBox>(*editor, "SPECIE_INLET_CONC_UNITS:FITC");
        EXPECT_EQ(inlet->currentData().toString(), "mg/ml");
        EXPECT_EQ(inlet->findData("um2/ul"), -1); // Molecules cannot use particle surface-area units.
        EXPECT_GE(widget<QComboBox>(*editor, "SPECIE_MODEL_CONC_UNITS:PS_40nm")->findData("um2/ul"), 0);
        auto* species = widget<QListWidget>(*editor, "experimentSpeciesChoices");
        for (int i = 0; i < species->count(); ++i) if (species->item(i)->text() == "FITC") {
            species->item(i)->setCheckState(Qt::Unchecked); species->item(i)->setCheckState(Qt::Checked);
        }
        EXPECT_EQ(inlet->currentData().toString(), "mg/ml");
        auto* form = qobject_cast<QFormLayout*>(widget<QLineEdit>(*editor, "NAME")->parentWidget()->layout());
        ASSERT_NE(form, nullptr);
        for (int i = 0; i < form->rowCount(); ++i) {
            auto* label = qobject_cast<QLabel*>(form->itemAt(i, QFormLayout::LabelRole)->widget());
            ASSERT_NE(label, nullptr);
            EXPECT_EQ(label->text(), table->horizontalHeaderItem(table->horizontalHeader()->logicalIndex(i))->text());
        }
        auto* widthUnits = widget<QComboBox>(*editor, "CHANNEL_WIDTHUnits");
        widthUnits->setCurrentText("mm");
        EXPECT_NEAR(widget<QLineEdit>(*editor, "CHANNEL_WIDTH")->text().toDouble(), .5, 1e-12);
        widget<QLineEdit>(*editor, "CHANNEL_WIDTH")->setText("2");
        widget<QComboBox>(*editor, "CHANNEL_HEIGHTUnits")->setCurrentIndex(3); // micrometers
        widget<QLineEdit>(*editor, "CHANNEL_HEIGHT")->setText("80");
        widget<QComboBox>(*editor, "CHANNEL_LENGTHUnits")->setCurrentText("cm");
        widget<QLineEdit>(*editor, "CHANNEL_LENGTH")->setText("5");
        widget<QComboBox>(*editor, "ENTRANCE_FLOWRATEUnits")->setCurrentIndex(5); // microliters/min
        EXPECT_NEAR(widget<QLineEdit>(*editor, "ENTRANCE_FLOWRATE")->text().toDouble(), 30, 1e-10);
        widget<QLineEdit>(*editor, "ENTRANCE_FLOWRATE")->setText("60 120");
        widget<QLineEdit>(*editor, "PARAMETERS_TO_SOLVE_FOR")->setText("keq1");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        if (editor->isVisible()) { ADD_FAILURE() << widget<QLabel>(*editor, "referenceEditorStatus")->text().toStdString(); editor->reject(); }
    });
    widget<QPushButton>(window, "experimentsModifyButton")->click();
    EXPECT_NEAR(input.execute("SELECT CHANNEL_WIDTH FROM experiments"), .002, 1e-12);
    EXPECT_NEAR(input.execute("SELECT CHANNEL_HEIGHT FROM experiments"), 80e-6, 1e-12);
    EXPECT_NEAR(input.execute("SELECT CHANNEL_LENGTH FROM experiments"), .05, 1e-12);
    EXPECT_EQ(input.execute("SELECT count(*) FROM experiments WHERE PARAMETERS_TO_SOLVE_FOR='keq1' AND SPECIES='FITC PS_40nm 40nm_Bound_Dye_1' AND SPECIE_INLET_CONC_UNITS='mg/ml wt% mg/ml'"), 1);
    EXPECT_EQ(input.execute("SELECT count(*) FROM experiments WHERE abs(CAST(substr(ENTRANCE_FLOWRATE,1,instr(ENTRANCE_FLOWRATE,' ')-1) AS REAL)-1e-9)<1e-20 AND abs(CAST(substr(ENTRANCE_FLOWRATE,instr(ENTRANCE_FLOWRATE,' ')+1) AS REAL)-2e-9)<1e-20"), 1);
}

TEST(Gui, ExperimentReactionMismatchOffersCancelAddAndRemoveWithoutPartialWrites)
{
    Inputs input;
    QFile migration(QString::fromUtf8(TSENSOR_FIXTURE_DIR) + "/../../migrations/001_channel_dimensions.sql");
    ASSERT_TRUE(migration.open(QIODevice::ReadOnly)); input.execute(migration.readAll());
    MainWindow window; input.choose(window); window.show();
    widget<QTabWidget>(window, "mainTabs")->setCurrentIndex(2);
    auto* table = widget<QTableWidget>(window, "experimentsTable"); table->selectRow(0);
    QTimer::singleShot(0, [&] {
        auto* editor = window.findChild<QDialog*>("experimentsEditor"); ASSERT_NE(editor, nullptr);
        auto* choices = widget<QListWidget>(*editor, "experimentSpeciesChoices");
        QListWidgetItem* removed = nullptr;
        for (int i = 0; i < choices->count(); ++i) if (choices->item(i)->text() == "PS_40nm") removed = choices->item(i);
        ASSERT_NE(removed, nullptr); removed->setCheckState(Qt::Unchecked);
        auto answer = [&](const QString& text) {
            QTimer::singleShot(0, [&, text] {
                auto* warning = editor->findChild<QMessageBox*>("reactionSpeciesWarning"); ASSERT_NE(warning, nullptr);
                for (auto* button : warning->buttons()) if (button->text() == text) { button->click(); return; }
                ADD_FAILURE() << "Missing warning action"; warning->reject();
            });
            widget<QPushButton>(*editor, "saveReferenceButton")->click();
            EXPECT_TRUE(editor->isVisible());
            EXPECT_EQ(input.execute("SELECT count(*) FROM experiments WHERE SPECIES='FITC PS_40nm 40nm_Bound_Dye_1' AND REACTIONS='FITC_40nm_1'"), 1);
        };
        answer("Cancel"); EXPECT_EQ(removed->checkState(), Qt::Unchecked);
        answer("Add species"); EXPECT_EQ(removed->checkState(), Qt::Checked);
        EXPECT_EQ(widget<QComboBox>(*editor, "SPECIE_INLET_CONC_UNITS:PS_40nm")->currentData().toString(), "wt%");
        removed->setCheckState(Qt::Unchecked); answer("Remove reactions");
        EXPECT_EQ(widget<QListWidget>(*editor, "experimentReactionChoices")->item(0)->checkState(), Qt::Unchecked);
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        if (editor->isVisible()) { ADD_FAILURE(); editor->reject(); }
    });
    widget<QPushButton>(window, "experimentsModifyButton")->click();
    EXPECT_EQ(input.execute("SELECT count(*) FROM experiments WHERE SPECIES='FITC 40nm_Bound_Dye_1' AND REACTIONS IS NULL AND SPECIE_INLET_CONC_UNITS='mg/ml mg/ml' AND SPECIE_MODEL_CONC_UNITS='umol umol'"), 1);
}

TEST(Gui, QuantityConversionsRejectInvalidAndOverflowingValues)
{
    EXPECT_EQ(convertQuantityList("60 120", 1e-9 / 60).split(' ').size(), 2);
    EXPECT_NEAR(convertQuantityList("60", 1e-9 / 60).toDouble(), 1e-9, 1e-22);
    EXPECT_THROW(convertQuantityList("nan", 1), std::invalid_argument);
    EXPECT_THROW(convertQuantityList("1e308", 1000), std::invalid_argument);
    EXPECT_THROW(convertQuantityList("1e-300", 1e-30), std::invalid_argument);
    QLineEdit line("0.00012345678901234567");
    std::unique_ptr<QComboBox> combo(quantityUnits("CHANNEL_WIDTH", &line));
    for (int repeat = 0; repeat < 5; ++repeat)
        for (int i = 0; i < combo->count(); ++i) combo->setCurrentIndex(i);
    EXPECT_EQ(quantityStoredText(&line, combo.get()), "0.00012345678901234567");
}

TEST(Gui, ExperimentSelectionsRequireUnitsForNewSpeciesAndRecheckDatabaseReferences)
{
    Inputs input; QWidget parent;
    ExperimentSelections selections(input.database, {}, &parent);
    auto* species = widget<QListWidget>(parent, "experimentSpeciesChoices");
    for (int i = 0; i < species->count(); ++i)
        if (species->item(i)->text() == "FITC") species->item(i)->setCheckState(Qt::Checked);
    EXPECT_THROW(selections.value("SPECIE_INLET_CONC_UNITS"), std::runtime_error);
    widget<QComboBox>(parent, "SPECIE_INLET_CONC_UNITS:FITC")->setCurrentText("mg/ml");
    widget<QComboBox>(parent, "SPECIE_MODEL_CONC_UNITS:FITC")->setCurrentText("umol");
    EXPECT_EQ(selections.value("SPECIE_INLET_CONC_UNITS"), "mg/ml");
    EXPECT_EQ(selections.value("SPECIE_MODEL_CONC_UNITS"), "umol");
    input.execute("UPDATE species SET SPECIES_TYPE='particle' WHERE SPECIES_NAME='FITC'");
    sqlite3* raw = nullptr;
    ASSERT_EQ(sqlite3_open(input.database.toUtf8().constData(), &raw), SQLITE_OK);
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    EXPECT_THROW(selections.validateReferences(raw), std::runtime_error);
}

TEST(Gui, ChannelDimensionsValidateSaveConflictRollbackAndDiscard)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME) VALUES('second')");
    QFile migration(QString::fromUtf8(TSENSOR_FIXTURE_DIR) + "/../../migrations/001_channel_dimensions.sql");
    ASSERT_TRUE(migration.open(QIODevice::ReadOnly)); input.execute(migration.readAll());
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs"); tabs->setCurrentIndex(2);
    auto* table = widget<QTableWidget>(window, "experimentsTable");
    const QStringList expectedHeaders{"Name", "Reactions", "Species", "Specie inlet conc units",
        "Specie model conc units", "Edges", "Width", "Normalization", "Low ref left", "Low ref right",
        "High ref left", "High ref right", "Channel width (m)", "Channel height (m)",
        "Channel length (m)", "Entrance flowrate"};
    for (int i = 0; i < expectedHeaders.size(); ++i)
        EXPECT_EQ(table->horizontalHeaderItem(table->horizontalHeader()->logicalIndex(i))->text(), expectedHeaders[i]);
    EXPECT_EQ(window.findChild<QPushButton*>("saveDimensionsButton"), nullptr);
    auto* modify = widget<QPushButton>(window, "experimentsModifyButton");
    EXPECT_FALSE(modify->isEnabled());
    widget<QCheckBox>(window, "experimentsSelectedOnly")->setChecked(false);
    table->selectRow(1);
    EXPECT_FALSE(table->item(1, 0)->flags() & Qt::ItemIsEditable);
    QTimer::singleShot(0, [&] {
        auto* editor = window.findChild<QDialog*>("experimentsEditor"); ASSERT_NE(editor, nullptr);
        EXPECT_TRUE(widget<QLineEdit>(*editor, "NAME")->isReadOnly());
        auto* width = widget<QLineEdit>(*editor, "CHANNEL_WIDTH");
        for (const auto invalid : {"nan", "0", "-1"}) {
            width->setText(invalid);
            widget<QPushButton>(*editor, "saveReferenceButton")->click();
            EXPECT_FALSE(widget<QLabel>(*editor, "referenceEditorStatus")->text().isEmpty());
            EXPECT_TRUE(editor->isVisible());
        }
        width->setText(".001");
        input.execute("UPDATE experiments SET CHANNEL_WIDTH=.002 WHERE NAME='uniform'");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("changed in the database"));
        widget<QPushButton>(*editor, "cancelReferenceButton")->click();
    });
    modify->click();
    EXPECT_DOUBLE_EQ(input.execute("SELECT CHANNEL_WIDTH FROM experiments WHERE NAME='uniform'"), .002);
    EXPECT_DOUBLE_EQ(input.execute("SELECT CHANNEL_WIDTH FROM experiments WHERE NAME='second'"), 5e-4);
    widget<QPushButton>(window, "experimentsRefreshButton")->click(); table->selectRow(1);
    QTimer::singleShot(0, [&] {
        auto* editor = window.findChild<QDialog*>("experimentsEditor"); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "CHANNEL_WIDTH")->setText(".001");
        widget<QLineEdit>(*editor, "CHANNEL_HEIGHT")->setText(".00008");
        widget<QLineEdit>(*editor, "CHANNEL_LENGTH")->setText(".05");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        if (editor->isVisible()) { ADD_FAILURE(); editor->reject(); }
    });
    modify->click();
    QTimer::singleShot(0, [&] {
        auto* editor = window.findChild<QDialog*>("experimentsEditor"); ASSERT_NE(editor, nullptr);
        widget<QLineEdit>(*editor, "NAME")->setText("uniform");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        EXPECT_TRUE(widget<QLabel>(*editor, "referenceEditorStatus")->text().contains("already exists"));
        widget<QLineEdit>(*editor, "NAME")->setText("new-experiment");
        widget<QPushButton>(*editor, "saveReferenceButton")->click();
        if (editor->isVisible()) { ADD_FAILURE(); editor->reject(); }
    });
    widget<QPushButton>(window, "experimentsAddButton")->click();
    EXPECT_DOUBLE_EQ(input.execute("SELECT CHANNEL_WIDTH FROM experiments WHERE NAME='new-experiment'"), 5e-4);
    EXPECT_DOUBLE_EQ(input.execute("SELECT CHANNEL_WIDTH FROM experiments WHERE NAME='uniform'"), .001);
    EXPECT_DOUBLE_EQ(input.execute("SELECT CHANNEL_HEIGHT FROM experiments WHERE NAME='uniform'"), .00008);
    EXPECT_DOUBLE_EQ(input.execute("SELECT CHANNEL_LENGTH FROM experiments WHERE NAME='uniform'"), .05);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    tsensor_workflow::run_session session(std::filesystem::path(input.database.toStdWString()));
    session.load_inputs();
    EXPECT_DOUBLE_EQ(session.parameters().experiment_runs.front().W, .001);
    EXPECT_DOUBLE_EQ(session.parameters().experiment_runs.front().H, .00008);
    EXPECT_DOUBLE_EQ(session.parameters().experiment_runs.front().L, .05);
    if (const auto capture = qEnvironmentVariable("NAVIER_EXPERIMENTS_TAB_CAPTURE"); !capture.isEmpty()) {
        QApplication::processEvents(); EXPECT_TRUE(window.grab().save(capture));
    }
}

TEST(SetupFile, RoundTripAndTransactionalRestore)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('future_control',NULL)");
    const auto original = model_controls::load(input.database);
    const auto filename = input.directory.filePath(QString::fromUtf8("saved setup ü.navier.json"));
    setup_file::Setup setup{input.database, input.directory.filePath("output with spaces"), original.rows};
    setup_file::save(filename, setup);
    const auto loaded = setup_file::load(filename);
    EXPECT_EQ(loaded.database, setup.database);
    EXPECT_EQ(loaded.outputDirectory, setup.outputDirectory);
    EXPECT_EQ(loaded.controls, original.rows);
    input.execute("UPDATE model_controls SET value='7' WHERE criterion='max_iterations'");
    const auto changed = model_controls::load(input.database);
    setup_file::restore(loaded, changed);
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    EXPECT_THROW(setup_file::restore(loaded, changed), std::runtime_error);
    auto incompatible = loaded;
    incompatible.controls.pop_back();
    EXPECT_THROW(setup_file::restore(incompatible, original), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_THROW(setup_file::save(input.database, setup), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    setup.controls.push_back(setup.controls.front());
    EXPECT_THROW(setup_file::save(filename, setup), std::runtime_error);
    EXPECT_EQ(setup_file::load(filename).controls, original.rows);
}

TEST(SetupFile, ValidatesFileAndResolvesRelativePaths)
{
    QTemporaryDir directory;
    const auto filename = directory.filePath("setup.json");
    auto write = [&](const QByteArray& bytes) {
        QFile file(filename); ASSERT_TRUE(file.open(QIODevice::WriteOnly)); file.write(bytes);
    };
    write(R"({"format":"navier-setup","version":1,"database":"inputs.db","outputDirectory":"reports","controls":[{"name":"max_iterations","value":"3"}]})");
    const auto loaded = setup_file::load(filename);
    EXPECT_EQ(loaded.database, directory.filePath("inputs.db"));
    EXPECT_EQ(loaded.outputDirectory, directory.filePath("reports"));
    for (const auto& bytes : {
        QByteArray("{"),
        QByteArray(R"({"format":"navier-setup","version":2})"),
        QByteArray(R"({"format":"navier-setup","version":1,"database":"a","outputDirectory":"b","controls":[{"name":"max_iterations","value":"-1"}]})"),
        QByteArray(R"({"format":"navier-setup","version":1,"database":"a","outputDirectory":"b","controls":[{"name":"max_iterations"}]})"),
        QByteArray(R"({"format":"navier-setup","version":1,"database":"a","outputDirectory":"b","controls":[{"name":"x","value":null},{"name":"x","value":""}]})")}) {
        write(bytes);
        EXPECT_THROW(setup_file::load(filename), std::runtime_error);
    }
    EXPECT_THROW(setup_file::load(directory.filePath("missing.json")), std::runtime_error);
    EXPECT_THROW(setup_file::save(directory.filePath("missing/setup.json"), loaded), std::runtime_error);
}

TEST(ModelControls, SupportsRealColumnNamesNullsAndPreservesOtherData)
{
    Inputs input;
    input.execute("CREATE TABLE controls_real(Parameter TEXT NOT NULL, Setting ANY);"
                  "INSERT INTO controls_real SELECT * FROM model_controls;"
                  "DROP TABLE model_controls; ALTER TABLE controls_real RENAME TO model_controls;"
                  "INSERT INTO model_controls VALUES('run_solver',NULL),('save_normalized_profiles','false')");
    auto original = model_controls::load(input.database);
    EXPECT_EQ(original.nameColumn, "Parameter");
    EXPECT_FALSE(control(original.rows, "run_solver").value.has_value());
    auto edited = original.rows;
    control(edited, "run_solver").value = "false";
    control(edited, "max_iterations").value = "12";
    model_controls::save(input.database, original, edited);
    auto loaded = model_controls::load(input.database);
    EXPECT_EQ(loaded.rows, edited);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    control(edited, "universal_solve_for").value.reset();
    model_controls::save(input.database, loaded, edited);
    loaded = model_controls::load(input.database);
    EXPECT_FALSE(control(loaded.rows, "universal_solve_for").value.has_value());
}

TEST(ModelControls, ValidationConflictAndTransactionRollback)
{
    Inputs input;
    auto original = model_controls::load(input.database);
    auto edited = original.rows;
    control(edited, "width resolution (X)").value = "0";
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    EXPECT_THROW(model_controls::validate({"convergence_epsx", "nan"}), std::runtime_error);
    EXPECT_THROW(model_controls::validate({"max_iterations", "3abc"}), std::runtime_error);
    EXPECT_THROW(model_controls::validate({"run_solver", "yes"}), std::runtime_error);
    EXPECT_NO_THROW(model_controls::validate({"convergence_epsx", "1e-12"}));
    edited = original.rows;
    control(edited, "debug_level").value = "4";
    control(edited, "max_iterations").value = "12";
    input.execute("CREATE TRIGGER reject_control BEFORE UPDATE ON model_controls "
                  "WHEN NEW.criterion='max_iterations' BEGIN SELECT RAISE(ABORT,'test failure'); END");
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original); // Earlier debug_level update rolled back.
    input.execute("DROP TRIGGER reject_control; UPDATE model_controls SET value='5' WHERE criterion='debug_level'");
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='debug_level'"), 5);
    EXPECT_EQ(input.execute("SELECT value FROM model_controls WHERE criterion='max_iterations'"), 3);
    EXPECT_THROW(model_controls::load(input.directory.filePath("missing.db")), std::runtime_error);
    EXPECT_FALSE(QFile::exists(input.directory.filePath("missing.db")));
}

TEST(ModelControls, UnknownRowsAndAmbiguousNamesCannotBeEdited)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('future_control','keep me')");
    auto original = model_controls::load(input.database);
    auto edited = original.rows;
    control(edited, "future_control").value = "changed";
    EXPECT_THROW(model_controls::save(input.database, original, edited), std::runtime_error);
    EXPECT_EQ(model_controls::load(input.database), original);
    input.execute("CREATE TABLE duplicate_controls(Parameter TEXT, Setting ANY);"
                  "INSERT INTO duplicate_controls VALUES('debug_level',1),('debug_level',2);"
                  "DROP TABLE model_controls; ALTER TABLE duplicate_controls RENAME TO model_controls");
    EXPECT_THROW(model_controls::load(input.database), std::runtime_error);
}

TEST(ModelControls, RequiredValuesAndStrictNumericBounds)
{
    for (const auto* value : {"", "-1", "7", "1.5", "nan"})
        EXPECT_THROW(model_controls::validate({"debug_level", value}), std::runtime_error);
    for (const auto* value : {"0", "6"})
        EXPECT_NO_THROW(model_controls::validate({"debug_level", value}));
    for (const auto* value : {"", "0", "-1e-9", "1e-3", "1", "nan", "inf"})
        EXPECT_THROW(model_controls::validate({"convergence_epsx", value}), std::runtime_error);
    EXPECT_NO_THROW(model_controls::validate({"convergence_epsx", "9.99e-4"}));
    for (const auto* name : {"debug_level", "experiment_name", "run_solver", "convergence_epsx"})
        EXPECT_THROW(model_controls::validate({name, std::nullopt}), std::runtime_error);
    EXPECT_NO_THROW(model_controls::validate({"universal_solve_for", std::nullopt}));
    EXPECT_NO_THROW(model_controls::validate({"universal_solve_for", ""}));
}

TEST(Gui, InvalidUnchangedFieldsBlockBothActions)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value=NULL WHERE criterion='experiment_name'");
    ControlsDialog dialog(input.database); dialog.show();
    auto* use = widget<QPushButton>(dialog, "useControlsButton");
    auto* save = widget<QPushButton>(dialog, "saveControlsButton");
    ASSERT_TRUE(until([&] { return use->isEnabled(); }));
    const auto original = model_controls::load(input.database);
    for (auto* button : {use, save}) {
        button->click();
        EXPECT_TRUE(dialog.isVisible());
        EXPECT_TRUE(widget<QLabel>(dialog, "controlsStatus")->text().contains("experiment_name"));
        EXPECT_EQ(model_controls::load(input.database), original);
    }
    widget<QListWidget>(dialog, "experimentChoices")->item(0)->setCheckState(Qt::Checked);
    auto* parameters = widget<QListWidget>(dialog, "parameterChoices");
    for (int i = 0; i < parameters->count(); ++i) parameters->item(i)->setCheckState(Qt::Unchecked);
    use->click();
    EXPECT_EQ(dialog.result(), QDialog::Accepted);
    auto values = dialog.values().rows;
    EXPECT_FALSE(control(values, "universal_solve_for").value);
}

TEST(Gui, ExperimentChecklistSupportsMultipleSelectionsAndExplicitPersistence)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME) VALUES('second'),('third')");
    const auto original = model_controls::load(input.database);
    auto choose = [](ControlsDialog& dialog, const QString& name, Qt::CheckState state) {
        auto* list = widget<QListWidget>(dialog, "experimentChoices");
        const auto items = list->findItems(name, Qt::MatchExactly);
        ASSERT_EQ(items.size(), 1);
        items.front()->setCheckState(state);
    };
    {
        ControlsDialog dialog(input.database); dialog.show();
        ASSERT_TRUE(until([&] { return widget<QPushButton>(dialog, "useControlsButton")->isEnabled(); }));
        EXPECT_EQ(widget<QListWidget>(dialog, "experimentChoices")->count(), 3);
        EXPECT_EQ(dialog.findChild<QLineEdit*>("experiment_name"), nullptr);
        choose(dialog, "second", Qt::Checked);
        dialog.reject();
        EXPECT_EQ(model_controls::load(input.database), original);
    }
    ControlsDialog dialog(input.database); dialog.show();
    ASSERT_TRUE(until([&] { return widget<QPushButton>(dialog, "useControlsButton")->isEnabled(); }));
    choose(dialog, "second", Qt::Checked);
    EXPECT_TRUE(widget<QToolButton>(dialog, "experiment_name")->text().contains("second"));
    if (const auto capture = qEnvironmentVariable("NAVIER_EXPERIMENTS_CAPTURE"); !capture.isEmpty()) {
        auto* picker = widget<QToolButton>(dialog, "experiment_name");
        picker->menu()->popup(picker->mapToGlobal(QPoint(0, picker->height())));
        QApplication::processEvents();
        EXPECT_TRUE(picker->menu()->grab().save(capture));
        picker->menu()->hide();
    }
    widget<QPushButton>(dialog, "useControlsButton")->click();
    ASSERT_EQ(dialog.result(), QDialog::Accepted);
    auto memory = dialog.values();
    EXPECT_EQ(control(memory.rows, "experiment_name").value, "uniform second");
    EXPECT_EQ(model_controls::load(input.database), original);
    ControlsDialog reopened(input.database, nullptr, memory); reopened.show();
    auto* save = widget<QPushButton>(reopened, "saveControlsButton");
    ASSERT_TRUE(until([&] { return save->isEnabled(); }));
    choose(reopened, "uniform", Qt::Unchecked);
    choose(reopened, "second", Qt::Unchecked);
    for (auto* button : {save, widget<QPushButton>(reopened, "useControlsButton")}) {
        button->click(); EXPECT_TRUE(reopened.isVisible());
        EXPECT_TRUE(widget<QLabel>(reopened, "controlsStatus")->text().contains("experiment_name"));
    }
    choose(reopened, "second", Qt::Checked);
    choose(reopened, "third", Qt::Checked);
    save->click();
    ASSERT_TRUE(until([&] { return reopened.result() == QDialog::Accepted; }));
    auto stored = model_controls::load(input.database);
    EXPECT_EQ(control(stored.rows, "experiment_name").value, "second third");
}

TEST(Gui, MissingSavedExperimentMustBeDeselected)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='missing uniform' WHERE criterion='experiment_name'");
    ControlsDialog dialog(input.database); dialog.show();
    auto* use = widget<QPushButton>(dialog, "useControlsButton");
    ASSERT_TRUE(until([&] { return use->isEnabled(); }));
    use->click(); EXPECT_TRUE(dialog.isVisible());
    EXPECT_TRUE(widget<QLabel>(dialog, "controlsStatus")->text().contains("missing"));
    widget<QListWidget>(dialog, "experimentChoices")->findItems("missing", Qt::MatchExactly).front()->setCheckState(Qt::Unchecked);
    use->click(); EXPECT_EQ(dialog.result(), QDialog::Accepted);
}

TEST(ModelControls, ParameterChoicesDependOnAnyExperimentsReactionCount)
{
    Inputs input;
    const QStringList base{"p1", "kon1", "keq1", "left_edge", "width", "QE1"};
    EXPECT_EQ(model_controls::solvableParameters(input.database), base);
    input.execute("INSERT INTO experiments(NAME,REACTIONS) VALUES('other',NULL)");
    EXPECT_EQ(model_controls::solvableParameters(input.database), base);
    input.execute("UPDATE experiments SET REACTIONS='  first  ' WHERE NAME='other'");
    EXPECT_EQ(model_controls::solvableParameters(input.database), base);
    input.execute("UPDATE experiments SET REACTIONS='first second' WHERE NAME='other'");
    auto expanded = base; expanded.append({"p2", "kon2", "keq2", "QE2"});
    EXPECT_EQ(model_controls::solvableParameters(input.database), expanded);
}

TEST(Gui, ParameterChecklistSupportsMemoryDefaultsAndEmptySelection)
{
    Inputs input;
    input.execute("INSERT INTO experiments(NAME,REACTIONS) VALUES('other','first second')");
    const auto original = model_controls::load(input.database);
    ControlsDialog dialog(input.database); dialog.show();
    auto* use = widget<QPushButton>(dialog, "useControlsButton");
    ASSERT_TRUE(until([&] { return use->isEnabled(); }));
    EXPECT_EQ(dialog.findChild<QLineEdit*>("universal_solve_for"), nullptr);
    auto* choices = widget<QListWidget>(dialog, "parameterChoices");
    ASSERT_EQ(choices->count(), 10);
    choices->findItems("kon2", Qt::MatchExactly).front()->setCheckState(Qt::Checked);
    use->click(); ASSERT_EQ(dialog.result(), QDialog::Accepted);
    auto memory = dialog.values();
    EXPECT_EQ(control(memory.rows, "universal_solve_for").value, "keq1 kon2");
    EXPECT_EQ(model_controls::load(input.database), original);
    ControlsDialog reopened(input.database, nullptr, memory); reopened.show();
    auto* save = widget<QPushButton>(reopened, "saveControlsButton");
    ASSERT_TRUE(until([&] { return save->isEnabled(); }));
    auto* storedChoices = widget<QListWidget>(reopened, "parameterChoices");
    EXPECT_EQ(storedChoices->findItems("kon2", Qt::MatchExactly).front()->checkState(), Qt::Checked);
    save->click(); ASSERT_TRUE(until([&] { return reopened.result() == QDialog::Accepted; }));
    auto stored = model_controls::load(input.database);
    EXPECT_EQ(control(stored.rows, "universal_solve_for").value, "keq1 kon2");
    ControlsDialog empty(input.database); empty.show();
    auto* emptySave = widget<QPushButton>(empty, "saveControlsButton");
    ASSERT_TRUE(until([&] { return emptySave->isEnabled(); }));
    auto* emptyChoices = widget<QListWidget>(empty, "parameterChoices");
    for (int i = 0; i < emptyChoices->count(); ++i) emptyChoices->item(i)->setCheckState(Qt::Unchecked);
    emptySave->click(); ASSERT_TRUE(until([&] { return empty.result() == QDialog::Accepted; }));
    stored = model_controls::load(input.database);
    EXPECT_FALSE(control(stored.rows, "universal_solve_for").value);
}

TEST(Gui, PresetButtonCreatesUpdatesAndLoadsSetupWithoutDatabaseWrites)
{
    Inputs input;
    const auto defaults = model_controls::load(input.database);
    const auto filename = input.directory.filePath("preset.navier.json");
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    auto chooseFile = [&](const QString& file) {
        QTimer::singleShot(0, [file] {
            auto* dialog = qobject_cast<QFileDialog*>(QApplication::activeModalWidget());
            ASSERT_NE(dialog, nullptr);
            dialog->selectFile(file);
            QMetaObject::invokeMethod(dialog, "accept", Qt::DirectConnection);
        });
    };
    tabs->setCurrentIndex(1);
    auto* editor = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
    ASSERT_TRUE(until([&] { return editor->ready(); }));
    auto* button = widget<QPushButton>(*editor, "useControlsButton");
    EXPECT_EQ(button->text(), "Save preset...");
    widget<QLineEdit>(*editor, "max_iterations")->setText("9");
    chooseFile(filename); button->click();
    auto saved = setup_file::load(filename);
    EXPECT_EQ(control(saved.controls, "max_iterations").value, std::optional<QString>("9"));
    EXPECT_EQ(button->text(), "Update preset");
    widget<QLineEdit>(*editor, "max_iterations")->setText("11");
    button->click();
    saved = setup_file::load(filename);
    EXPECT_EQ(control(saved.controls, "max_iterations").value, std::optional<QString>("11"));
    chooseFile(filename); widget<QAction>(window, "openSetupAction")->trigger();
    QApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    tabs->setCurrentIndex(1);
    editor = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
    ASSERT_TRUE(until([&] { return editor->ready(); }));
    EXPECT_EQ(widget<QLineEdit>(*editor, "max_iterations")->text(), "11");
    EXPECT_EQ(widget<QPushButton>(*editor, "useControlsButton")->text(), "Update preset");
    widget<QLineEdit>(*editor, "max_iterations")->setText("13");
    widget<QPushButton>(*editor, "useControlsButton")->click();
    saved = setup_file::load(filename);
    EXPECT_EQ(control(saved.controls, "max_iterations").value, std::optional<QString>("13"));
    EXPECT_EQ(model_controls::load(input.database), defaults);
}

TEST(Gui, ModelControlsDoNotLockNavigationOrClose)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    auto* run = widget<QPushButton>(window, "runButton");
    tabs->setCurrentIndex(1);
    auto* editor = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
    ASSERT_NE(editor, nullptr);
    ASSERT_TRUE(until([&] { return editor->ready(); }));
    EXPECT_TRUE(run->isEnabled());
    EXPECT_TRUE(tabs->isTabEnabled(2));
    EXPECT_TRUE(widget<QAction>(window, "openSetupAction")->isEnabled());
    widget<QLineEdit>(*editor, "max_iterations")->setText("9");
    tabs->setCurrentIndex(2);
    EXPECT_EQ(tabs->currentIndex(), 2);
    EXPECT_TRUE(run->isEnabled());
    tabs->setCurrentIndex(1);
    EXPECT_EQ(widget<QLineEdit>(*editor, "max_iterations")->text(), "9");
    EXPECT_TRUE(window.close());
}

TEST(Gui, UntouchedControlsPreserveOptionalValuesAndFormatting)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='' WHERE criterion='universal_solve_for';"
                  "UPDATE model_controls SET value=' uniform  ' WHERE criterion='experiment_name'");
    const auto original = model_controls::load(input.database);
    ControlsDialog editor(input.database);
    ASSERT_TRUE(until([&] { return editor.ready(); }));
    EXPECT_EQ(editor.draft(false).rows, original.rows);
}

TEST(Gui, RunPromptsOnlyForChangedControlsAndPreservesDefaults)
{
    Inputs input;
    input.execute("INSERT INTO model_controls VALUES('run_solver','false')");
    const auto defaults = model_controls::load(input.database);
    MainWindow window; input.choose(window); window.show();
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    auto* run = widget<QPushButton>(window, "runButton");
    tabs->setCurrentIndex(1);
    auto* editor = dynamic_cast<ControlsDialog*>(window.findChild<QDialog*>("modelControlsEditor"));
    ASSERT_TRUE(until([&] { return editor->ready(); }));
    // Merely visiting the tab must not require acceptance.
    run->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    auto* iterations = widget<QLineEdit>(*editor, "max_iterations");
    iterations->setText("9");
    auto answer = [](QMessageBox::StandardButton button) {
        QTimer::singleShot(0, [button] {
            auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
            ASSERT_NE(box, nullptr);
            EXPECT_EQ(box->windowTitle(), "Accept model controls");
            box->button(button)->click();
        });
    };
    answer(QMessageBox::Cancel); run->click();
    EXPECT_TRUE(run->isEnabled());
    EXPECT_EQ(iterations->text(), "9");
    iterations->setText("bad");
    answer(QMessageBox::Yes); run->click();
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("max_iterations"));
    iterations->setText("9");
    answer(QMessageBox::Yes); run->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Optimizer iterations: Not run"));
    run->click(); // Accepted edits do not prompt again.
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_EQ(model_controls::load(input.database), defaults);
}

TEST(ModelControls, WorkerUsesSnapshotBeforeLoadingExperimentInputs)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='false' WHERE criterion='save_model_profiles'");
    const auto defaults = model_controls::load(input.database);
    auto rows = defaults.rows;
    control(rows, "max_iterations").value = "19";
    control(rows, "width resolution (X)").value = "7";
    control(rows, "save_model_profiles").value.reset();
    tsensor_workflow::control_values values;
    for (const auto& row : rows)
        values.emplace_back(row.name.toStdString(), row.value ? std::optional<std::string>(row.value->toStdString()) : std::nullopt);
    values.emplace_back("run_solver", "false");
    tsensor_workflow::background_runner runner;
    runner.start(std::filesystem::path(input.database.toStdWString()), {}, values);
    runner.wait();
    auto result = runner.take_result();
    if (result.failure) std::rethrow_exception(result.failure);
    ASSERT_NE(result.session, nullptr);
    EXPECT_EQ(result.session->parameters().max_iterations, 19);
    EXPECT_EQ(result.session->parameters().X, 7);
    EXPECT_TRUE(result.session->parameters().save_model_profiles);
    EXPECT_FALSE(result.result->optimizer_ran);
    EXPECT_EQ(model_controls::load(input.database), defaults);
}

TEST(Gui, RunExportAndExplicitSaves)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    auto* exportButton = widget<QPushButton>(window, "exportButton");
    bool sawLiveParameters = false;
    auto* liveTable = widget<QTableWidget>(window, "parameterTable");
    QObject::connect(liveTable, &QTableWidget::itemChanged, &window, [&](QTableWidgetItem* item) {
        // The final result uses a different summary. Observe the queued
        // optimizer report even when this small solve finishes between polls.
        if (item->column() == 3 && !item->text().isEmpty() &&
            widget<QLabel>(window, "resultSummary")->text().startsWith("Completed model evaluations:")) {
            sawLiveParameters = true;
            EXPECT_FALSE(run->isEnabled());
            EXPECT_NEAR(item->text().toDouble(), 1, 1e-10);
            EXPECT_DOUBLE_EQ(liveTable->item(item->row(), 2)->text().toDouble(), 1);
        }
    });
    EXPECT_FALSE(exportButton->isEnabled());
    run->click();
    EXPECT_FALSE(widget<QLineEdit>(window, "databasePath")->isEnabled());
    EXPECT_FALSE(run->isEnabled());
    EXPECT_FALSE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(3));
    EXPECT_FALSE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(4));
    EXPECT_FALSE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(5));
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(3));
    EXPECT_TRUE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(4));
    EXPECT_TRUE(widget<QTabWidget>(window, "mainTabs")->isTabEnabled(5));
    ASSERT_TRUE(exportButton->isEnabled());
    EXPECT_TRUE(sawLiveParameters);
    EXPECT_EQ(widget<QTableWidget>(window, "parameterTable")->rowCount(), 1);
    auto* parameters = widget<QTableWidget>(window, "parameterTable");
    ASSERT_EQ(parameters->columnCount(), 4);
    EXPECT_EQ(parameters->horizontalHeaderItem(2)->text(), "Initial value");
    EXPECT_DOUBLE_EQ(parameters->item(0, 2)->text().toDouble(), 1);
    EXPECT_NEAR(parameters->item(0, 3)->text().toDouble(), 1, 1e-10);
    EXPECT_FALSE(widget<QPlainTextEdit>(window, "progressLog")->toPlainText().contains("Initial parameter "));
    EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Termination code:"));
    EXPECT_TRUE(widget<QLabel>(window, "resultSummary")->text().contains("Sum of squared residuals:"));
    if (const auto capture = qEnvironmentVariable("NAVIER_GUI_CAPTURE"); !capture.isEmpty()) {
        QApplication::processEvents();
        EXPECT_TRUE(window.grab().save(capture));
    }
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    exportButton->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("reports/uniform,.txt")));
    auto* tabs = widget<QTabWidget>(window, "mainTabs");
    EXPECT_EQ(tabs->tabText(tabs->currentIndex()), "Report");
    ASSERT_NE(widget<QTableView>(window, "reportTable")->model(), nullptr);
    EXPECT_TRUE(widget<QPushButton>(window, "copyReportButton")->isEnabled());
    widget<QPushButton>(window, "profilesButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 36);
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    QTimer::singleShot(0, [] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget()))
            box->button(QMessageBox::Yes)->click();
    });
    widget<QPushButton>(window, "inputsButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_EQ(input.execute("SELECT [INITIAL VALUE] FROM alglib_input"), 1);
    EXPECT_EQ(input.execute("SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    // Changing databases must not leave actions attached to an old session.
    widget<QLineEdit>(window, "databasePath")->setText("another.db");
    EXPECT_FALSE(exportButton->isEnabled());
    EXPECT_EQ(widget<QTableWidget>(window, "parameterTable")->rowCount(), 0);
}

TEST(Gui, InitialValuesRemainVisibleWhenSolveFails)
{
    Inputs input;
    // Invalid optimizer settings fail after inputs (and initial values) are loaded.
    input.execute("UPDATE model_controls SET value='-1' WHERE criterion='max_iterations'");
    MainWindow window; input.choose(window);
    auto* run = widget<QPushButton>(window, "runButton");
    run->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    auto* table = widget<QTableWidget>(window, "parameterTable");
    ASSERT_EQ(table->rowCount(), 1);
    EXPECT_EQ(table->item(0, 1)->text(), "keq1");
    EXPECT_DOUBLE_EQ(table->item(0, 2)->text().toDouble(), 1);
    EXPECT_TRUE(table->item(0, 3)->text().isEmpty());
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("negative MaxIts"));
}

TEST(Gui, CancelAndCloseDuringCalculation)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='100000' WHERE criterion='length/time resolution (Z)'");
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    bool heartbeat = false;
    run->click();
    QTimer::singleShot(0, &window, [&] {
        heartbeat = true;
        EXPECT_FALSE(run->isEnabled());
        widget<QPushButton>(window, "cancelButton")->click();
        // Programmatic changes can refresh controls even while inputs are locked.
        widget<QLineEdit>(window, "outputPath")->setText(input.directory.filePath("cancelled"));
        EXPECT_FALSE(widget<QPushButton>(window, "cancelButton")->isEnabled());
    });
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(heartbeat);
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("cancelled"));
    EXPECT_FALSE(widget<QPushButton>(window, "exportButton")->isEnabled());
    run->click();
    EXPECT_FALSE(window.close()); // Close deferred while the worker stops.
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
}

TEST(Gui, FailureCanBeCorrectedAndRetried)
{
    Inputs input;
    input.execute("UPDATE alglib_input SET SCALE=0");
    MainWindow window; input.choose(window);
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_FALSE(widget<QPushButton>(window, "inputsButton")->isEnabled());
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("zero passed to SCALE"));
    input.execute("UPDATE alglib_input SET SCALE=1");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(widget<QPushButton>(window, "exportButton")->isEnabled());
    widget<QLineEdit>(window, "outputPath")->setText(input.database); // File, not directory.
    widget<QPushButton>(window, "exportButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("export", Qt::CaseInsensitive));
    widget<QLineEdit>(window, "outputPath")->setText(input.directory.filePath("retry"));
    widget<QPushButton>(window, "exportButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("retry/uniform,.txt")));
}
TEST(Gui, ExplicitProfileSaveIgnoresLegacyFlagAndCloseDuringExport)
{
    Inputs input;
    input.execute("UPDATE model_controls SET value='false' WHERE criterion='save_model_profiles'");
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    ASSERT_TRUE(widget<QPushButton>(window, "profilesButton")->isEnabled());
    EXPECT_EQ(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    widget<QPushButton>(window, "profilesButton")->click();
    ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    EXPECT_GT(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    EXPECT_EQ(input.execute("SELECT value='false' FROM model_controls WHERE criterion='save_model_profiles'"), 1);
    ASSERT_TRUE(widget<QPushButton>(window, "exportButton")->isEnabled());
    widget<QPushButton>(window, "exportButton")->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("reports/uniform,.txt")));
    EXPECT_GT(input.execute("SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
}

TEST(Gui, FailedExportAbortsCloseAndAllowsRetry)
{
    Inputs input;
    MainWindow window; input.choose(window); window.show();
    auto* run = widget<QPushButton>(window, "runButton");
    auto* exportButton = widget<QPushButton>(window, "exportButton");
    run->click(); ASSERT_TRUE(until([&] { return run->isEnabled(); }));
    widget<QLineEdit>(window, "outputPath")->setText(input.database);
    exportButton->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return run->isEnabled() || !window.isVisible(); }));
    ASSERT_TRUE(window.isVisible());
    EXPECT_TRUE(widget<QLabel>(window, "runStatus")->text().contains("export", Qt::CaseInsensitive));
    ASSERT_TRUE(exportButton->isEnabled());
    widget<QLineEdit>(window, "outputPath")->setText(input.directory.filePath("retry"));
    exportButton->click();
    EXPECT_FALSE(window.close());
    ASSERT_TRUE(until([&] { return !window.isVisible(); }));
    EXPECT_TRUE(QFile::exists(input.directory.filePath("retry/uniform,.txt")));
}
} // namespace

int main(int argc, char** argv) {
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
    QApplication application(argc, argv);
    QApplication::setQuitOnLastWindowClosed(false);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
