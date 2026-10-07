#include <background_runner.h>
#include <gtest/gtest.h>
#include <atomic>
#include <sstream>
#include <type_traits>

namespace {
using namespace tsensor_workflow;

void sql(sqlite3* db, const std::string& text)
{
    if (sqlite3_exec(db, text.c_str(), nullptr, nullptr, nullptr) != SQLITE_OK)
        throw std::runtime_error(sqlite3_errmsg(db));
}

double scalar(sqlite3* db, const char* query)
{
    sqlite3_stmt* raw = nullptr;
    const int rc = sqlite3_prepare_v2(db, query, -1, &raw, nullptr);
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(raw, sqlite3_finalize);
    if (rc != SQLITE_OK || sqlite3_step(raw) != SQLITE_ROW)
        throw std::runtime_error(sqlite3_errmsg(db));
    return sqlite3_column_double(raw, 0);
}

struct fixture {
    std::filesystem::path root, database;
    fixture() {
        // Atomic directory creation avoids collisions across parallel test processes.
        for (int i = 0; ; ++i) {
            root = std::filesystem::temp_directory_path() / ("navier-workflow-" + std::to_string(i));
            if (std::filesystem::create_directory(root)) break;
        }
        database = root / "inputs.db";
        sqlite3* raw = nullptr;
        const auto name = database.u8string();
        const int rc = sqlite3_open(reinterpret_cast<const char*>(name.c_str()), &raw);
        std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
        if (rc != SQLITE_OK) throw std::runtime_error("Cannot create test database");
        std::ifstream input(std::filesystem::path(TSENSOR_FIXTURE_DIR) / "workflow.sql");
        if (!input) throw std::runtime_error("Missing workflow fixture");
        sql(db.get(), std::string(std::istreambuf_iterator<char>(input), {}));
    }
    ~fixture() { std::error_code error; std::filesystem::remove_all(root, error); }
};

void near(double actual, double expected)
{
    ASSERT_TRUE(std::isfinite(actual));
    EXPECT_NEAR(actual, expected, 1e-12 + 1e-10 * std::abs(expected));
}

void migrate_channels(sqlite3* db) {
    std::ifstream file(std::filesystem::path(TSENSOR_FIXTURE_DIR) / "../../migrations/001_channel_dimensions.sql");
    if (!file) throw std::runtime_error("Missing channel migration");
    sql(db, std::string(std::istreambuf_iterator<char>(file), {}));
}

TEST(Workflow, NoSelectedParametersEvaluatesAndGeneratesReport) {
    for (const char* selection : {"NULL", "''"}) {
        for (const char* optimize : {"true", "false"}) {
            SCOPED_TRACE(std::string(selection) + " run_solver=" + optimize);
            fixture files;
            int optimizer_events = 0;
            run_session session(files.database, [&](const progress_event& event) {
                if (event.kind == event_kind::optimizer_progress) ++optimizer_events;
            });
            sql(session.database(), std::string("UPDATE model_controls SET value=") + selection +
                " WHERE criterion='universal_solve_for';"
                "INSERT INTO model_controls VALUES('run_solver','" + optimize + "');"
                "DROP TABLE alglib_input;"); // No optimizer inputs are needed.
            session.load_inputs();
            ASSERT_TRUE(session.parameters().solvables.empty());
            const auto result = session.run();
            EXPECT_EQ(session.state(), session_state::completed);
            EXPECT_FALSE(result.optimizer_ran);
            EXPECT_FALSE(result.optimizer_iterations.has_value());
            EXPECT_FALSE(result.termination_type.has_value());
            EXPECT_TRUE(result.parameters.empty());
            EXPECT_EQ(result.residual_evaluations, 1);
            EXPECT_EQ(optimizer_events, 0);
            ASSERT_TRUE(result.sum_squared_residuals.has_value());
            near(*result.sum_squared_residuals, 0);
            // Uniform dye, no reaction and no-flux walls preserve concentration.
            for (const auto& experiment : session.parameters().experiments) {
                ASSERT_FALSE(experiment.model_profile.empty());
                for (double value : experiment.model_profile) near(value, 20);
                for (int x = 0; x < experiment.window_size; ++x)
                    near(experiment.channel_position[x], x * experiment.run->W * 1e6 / 8);
            }
            EXPECT_NE(session.generate_report().find("Species:FITC"), std::string::npos);
            EXPECT_GT(std::filesystem::file_size(session.export_results(files.root / "reports")), 0u);
            session.save_model_profiles();
            near(scalar(session.database(), "SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 36);
            near(scalar(session.database(), "SELECT max(abs(Numeric-0.002)) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
            near(scalar(session.database(), "SELECT count(*) FROM parameter_solutions WHERE SOLVE_SETTING_ID<>99"), 0);
        }
    }
}

TEST(Workflow, ReactionKsAliasesResolveThroughSolveForVariables) {
    fixture files;
    run_session session(files.database);
    sql(session.database(),
        "UPDATE model_controls SET value='association_rate' WHERE criterion='universal_solve_for';"
        "INSERT INTO alglib_input VALUES('association_rate',0.75,0.5,2,1);"
        "UPDATE reactions SET Ks='#association_rate 1' WHERE REACTION_NAME='FITC_40nm_1'");
    session.load_inputs();

    auto& p = session.parameters();
    ASSERT_EQ(p.solvables.size(), 1u);
    ASSERT_EQ(p.initial_values_alglib.size(), 1u);
    near(p.initial_values_alglib.front(), 0.75);
    auto& reaction = p.experiment_runs.front().reactions.front();
    EXPECT_EQ(reaction.k_alias[0], "association_rate");
    EXPECT_EQ(p.experiment_runs.front().reaction_variables.at("association_rate").source, &p.solvables.front());
    near(reaction.k[0].value(), 0.75);
}

TEST(Workflow, ReactionCoefficientAndExponentAliasesSupportSigns) {
    fixture files;
    run_session session(files.database);
    sql(session.database(),
        "UPDATE model_controls SET value='p1' WHERE criterion='universal_solve_for';"
        "INSERT INTO alglib_input VALUES('p1',0.75,0.5,2,1);"
        "UPDATE reactions SET COEFFICIENTS='#p1 -1 -#p1', EXPONENTS='#p1 2 -#p1' "
        "WHERE REACTION_NAME='FITC_40nm_1'");
    session.load_inputs();

    auto& p = session.parameters();
    auto& run = p.experiment_runs.front();
    auto& reaction = run.reactions.front();
    auto* fitc = &run.species.at(specie_index(&run, "FITC"));
    auto* bound = &run.species.at(specie_index(&run, "40nm_Bound_Dye_1"));
    ASSERT_EQ(p.solvables.size(), 1u);
    ASSERT_EQ(p.initial_values_alglib.size(), 1u);
    near(reaction.coef.at(fitc).value(), 0.75);
    near(reaction.coef.at(bound).value(), -0.75);
    near(reaction.exp.at(fitc).value(), 0.75);
    near(reaction.exp.at(bound).value(), -0.75);

    p.row_count = 0; // Exercise the solver's parameter-update step without a transient solve.
    alglib::real_1d_array control_parameters, residuals;
    control_parameters.setlength(1);
    control_parameters[0] = 0.625;
    alglib_solver(control_parameters, residuals, &p);
    near(reaction.coef.at(fitc).value(), 0.625);
    near(reaction.coef.at(bound).value(), -0.625);
    near(reaction.exp.at(fitc).value(), 0.625);
    near(reaction.exp.at(bound).value(), -0.625);
}

TEST(Workflow, ReactionKsAliasesRequireVariablesTableEntry) {
    fixture files;
    run_session session(files.database);
    sql(session.database(),
        "UPDATE model_controls SET value='association_rate' WHERE criterion='universal_solve_for';"
        "UPDATE reactions SET Ks='#association_rate 1' WHERE REACTION_NAME='FITC_40nm_1'");

    EXPECT_THROW(session.load_inputs(), std::runtime_error);
}

TEST(Workflow, ReactionKsAliasesMustBeSelectedForSolving) {
    fixture files;
    run_session session(files.database);
    sql(session.database(),
        "UPDATE reactions SET Ks='#association_rate 1' WHERE REACTION_NAME='FITC_40nm_1'");

    EXPECT_THROW(session.load_inputs(), std::runtime_error);
}

TEST(Workflow, RawProfileEntranceConcentrationOverridesLegacyInletId) {
    fixture files;
    run_session session(files.database);
    sql(session.database(), "UPDATE raw_profile SET ENTRANCE_CONC=5, ENTRANCE_CONC_UNITS='umol' WHERE WT_PERCENT='1.0'");
    migrate_channels(session.database());
    session.load_inputs();
    auto& parameters = session.parameters();
    const auto profile = std::find_if(parameters.experiments.begin(), parameters.experiments.end(),
        [](const experiment_struct& item) { return item.second_name == "1.0"; });
    ASSERT_NE(profile, parameters.experiments.end());
    auto* fitc = &profile->run->species.at(profile->run->FITC);
    ASSERT_FALSE(profile->entrances.empty());
    ASSERT_TRUE(profile->entrances.front().CONC.contains(fitc));
    // 5 umol of a 100 g/mol molecule is 0.0005 mg/ml; this replaces ID 1's 0.002 mg/ml.
    EXPECT_NEAR(profile->entrances.front().CONC.at(fitc), 0.0005,
                1e-12 + 1e-10 * 0.0005);
}

TEST(Workflow, OmittedProfileParametersNeverReachGuiEventsOrSolver) {
    for (bool omit_first : {false, true}) {
        SCOPED_TRACE(omit_first);
        fixture files;
        std::vector<progress_event> events;
        run_session session(files.database, [&](const progress_event& event) {
            events.push_back(event);
        });
        auto* db = session.database();
        // Exercise either column order: omission must be known before registration.
        sql(db, "ALTER TABLE raw_profile RENAME TO original_profiles");
        sql(db, std::string("CREATE TABLE raw_profile AS SELECT ") +
            (omit_first ? "OMIT," : "") +
            "NAME,WT_PERCENT,CHANNEL_LEFT_EDGE,CHANNEL_RIGHT_EDGE,INTENSITY_ARRAY,"
            "INLET_COND_ID,'width' AS INDEPENDENT_PARAMETERS_TO_SOLVE_FOR,LEFT_EDGE,WIDTH" +
            (omit_first ? "" : ",OMIT") + " FROM original_profiles");
        sql(db, "UPDATE raw_profile SET OMIT='TrUe' WHERE WT_PERCENT='4.0';"
                "INSERT INTO alglib_input VALUES('width',8,7,9,1);"
                "INSERT INTO model_controls VALUES('run_solver','true')");
        session.load_inputs();
        auto& p = session.parameters();
        ASSERT_EQ(p.solvables.size(), 4u); // Universal keq1 and three active widths.
        EXPECT_EQ(p.initial_values_alglib.size(), 4u);
        EXPECT_EQ(p.scale.size(), 4u);
        EXPECT_EQ(p.low_bound.size(), 4u);
        EXPECT_EQ(p.up_bound.size(), 4u);
        for (const auto& parameter : p.solvables)
            EXPECT_NE(parameter.source_name, "uniform_4.0");
        ASSERT_TRUE(p.experiments.back().omit);
        EXPECT_EQ(p.experiments.back().width.source, &p.experiments.back().run->width);
        const auto result = session.run();
        EXPECT_GT(result.residual_evaluations, 0);
        ASSERT_EQ(result.parameters.size(), 4u);
        bool initialized = false, progressed = false;
        for (const auto& event : events) {
            if (event.kind != event_kind::parameters_initialized &&
                event.kind != event_kind::optimizer_progress) continue;
            initialized |= event.kind == event_kind::parameters_initialized;
            progressed |= event.kind == event_kind::optimizer_progress;
            EXPECT_EQ(event.parameters.size(), 4u);
            for (const auto& parameter : event.parameters)
                EXPECT_NE(parameter.source, "uniform_4.0");
        }
        EXPECT_TRUE(initialized);
        EXPECT_TRUE(progressed);
    }
}

TEST(ScatterControl, DatabaseAndInMemoryChoicesReachResiduals) {
    // Uniform dye and beads, zero reaction rate, and no-flux walls give a
    // constant analytical solution. At the Gaussian center the correction
    // reduces to 1 - amplitude + center * slope, independently of the solver.
    for (bool in_memory : {false, true}) {
        for (double diameter : {20.0, 40.0}) {
            const double center = diameter < 30 ? 0.0258645848310996 : 0.0711222856018783;
            const double amplitude = diameter < 30 ? 0.238826108843563 : 0.448191328804794;
            const double slope = diameter < 30 ? 1.40044738887211 : 2.14776259822044;
            fixture files;
            // Repeated fresh sessions also check that on does not leak into off.
            for (const std::string choice : {"none", "NS_ND", "none"}) {
                SCOPED_TRACE(choice + " diameter=" + std::to_string(diameter) +
                             " in_memory=" + std::to_string(in_memory));
                run_session session(files.database);
                auto* db = session.database();
                sql(db, "UPDATE species SET PARTICLE_DIAMETER=" + std::to_string(diameter) +
                        " WHERE SPECIES_NAME='PS_40nm'; DELETE FROM inlet_conditions WHERE SPECIES_NAME='PS_40nm'");
                std::ostringstream beads;
                beads.precision(17);
                beads << "INSERT INTO inlet_conditions SELECT INLET_COND_ID," << center
                      << ",1,'PS_40nm' FROM inlet_conditions WHERE SPECIES_NAME='FITC'";
                sql(db, beads.str());
                sql(db, "UPDATE model_controls SET value='" +
                        (in_memory ? (choice == "none" ? std::string("NS_ND") : std::string("none")) : choice) +
                        "' WHERE criterion='scatter_correction_type'");
                if (in_memory) {
                    control_values controls{{"experiment_name", "uniform"},
                        {"debug_level", "7"}, {"width resolution (X)", "8"},
                        {"length/time resolution (Z)", "4"}, {"run_solver", "true"},
                        {"universal_solve_for", "keq1"}, {"max_iterations", "3"},
                        {"scatter_correction_type", choice}};
                    session.load_inputs(controls);
                } else {
                    sql(db, "INSERT OR REPLACE INTO model_controls VALUES('run_solver','true')");
                    session.load_inputs();
                }
                EXPECT_EQ(session.parameters().scatter_correction_type, choice);
                const auto result = session.run();
                EXPECT_GT(result.residual_evaluations, 0);
                const double factor = choice == "none" ? 1 : 1 - amplitude + center * slope;
                // ALGLIB squares the signed profile discrepancies exactly once.
                ASSERT_TRUE(result.sum_squared_residuals.has_value());
                near(*result.sum_squared_residuals,
                     session.parameters().total_window_size * std::pow(factor - 1, 2));
                for (const auto& exp : session.parameters().experiments) {
                    for (double value : exp.model_profile) near(value, 20);
                    for (double value : exp.error) near(value, (factor - 1) * (factor - 1));
                }
            }
        }
    }
}

TEST(Workflow, SignedResidualsGiveLeastSquaresObjectiveAndPreserveReportedErrors)
{
    fixture files;
    std::vector<progress_event> reports;
    run_session session(files.database, [&](const progress_event& event) {
        if (event.kind == event_kind::optimizer_progress) reports.push_back(event);
    });
    sql(session.database(), "INSERT INTO model_controls VALUES('run_solver','true')");
    session.load_inputs();
    auto& p = session.parameters();
    ASSERT_EQ(p.experiments.size(), 4u);
    // Uniform equilibrium gives model/dye_conc = 1 independently of keq1.
    // Choose normalized observations to exercise both signs, zero, and omission.
    const double observations[] = {0.5, 3.0, 1.0, 5.0};
    const double expected[] = {0.5, -2.0, 0.0, 0.0};
    double expected_sum = 0;
    for (std::size_t row = 0; row < p.experiments.size(); ++row) {
        auto& exp = p.experiments[row];
        std::fill(exp.raw_experimental_profile.begin(), exp.raw_experimental_profile.end(), observations[row]);
        exp.omit = row == 3;
        expected_sum += exp.window_size * expected[row] * expected[row];
    }
    alglib::real_1d_array controls, residuals;
    controls.setcontent(p.initial_values_alglib.size(), p.initial_values_alglib.data());
    residuals.setlength(p.total_window_size);
    alglib_solver(controls, residuals, &p);
    for (std::size_t row = 0; row < p.experiments.size(); ++row) {
        const auto& exp = p.experiments[row];
        for (int i = 0; i < exp.window_size; ++i) {
            near(residuals[exp.window_start + i], expected[row]);
            near(exp.error[i], expected[row] * expected[row]);
        }
    }
    const auto result = session.run();
    ASSERT_TRUE(result.sum_squared_residuals.has_value());
    near(*result.sum_squared_residuals, expected_sum);
    ASSERT_FALSE(reports.empty());
    for (const auto& report : reports) {
        ASSERT_TRUE(report.sum_squared_residuals.has_value());
        near(*report.sum_squared_residuals, expected_sum);
    }
}

static_assert(std::is_const_v<decltype(parameters_t::W)>);
static_assert(std::is_const_v<decltype(experiment_run_struct::W)>);
static_assert(std::is_const_v<decltype(experiment_run_struct::H)>);
static_assert(std::is_const_v<decltype(experiment_run_struct::L)>);

TEST(ChannelDimensions, DefaultsAndInvalidConstructorValues) {
    experiment_run_struct initialized(channel_dimensions(.001, .00008, .05));
    near(initialized.dye_conc, 0); near(initialized.dt, 0);
    EXPECT_EQ(initialized.number_of_species, 0);
    parameters_t p;
    near(p.W, 5e-4); near(p.H, 4e-5); near(p.L, .025);
    parameters_t custom(channel_dimensions(.001, .00008, .05));
    near(custom.W, .001); near(custom.H, .00008); near(custom.L, .05);
    EXPECT_THROW(channel_dimensions(0, 1, 1), std::invalid_argument);
    EXPECT_THROW(channel_dimensions(1, -1, 1), std::invalid_argument);
    EXPECT_THROW(channel_dimensions(1, 1, std::numeric_limits<double>::infinity()), std::invalid_argument);
}

TEST(ChannelDimensions, MixedExperimentsUseIndependentGeometryAndUniformAnalyticalSolution) {
    fixture files;
    run_session session(files.database);
    auto* db = session.database();
    sql(db, "CREATE TEMP TABLE second AS SELECT * FROM experiments; UPDATE second SET NAME='wide'; "
            "INSERT INTO experiments SELECT * FROM second; "
            "CREATE TEMP TABLE profiles AS SELECT * FROM raw_profile; UPDATE profiles SET NAME='wide'; "
            "INSERT INTO raw_profile SELECT * FROM profiles; "
            "UPDATE model_controls SET value='uniform wide' WHERE criterion='experiment_name'");
    migrate_channels(db);
    sql(db, "UPDATE experiments SET CHANNEL_WIDTH=.001, CHANNEL_HEIGHT=.00008, CHANNEL_LENGTH=.05 WHERE NAME='wide'");
    session.load_inputs();
    const auto& runs = session.parameters().experiment_runs;
    ASSERT_EQ(runs.size(), 2u);
    // Volume/flow = 1 s and 8 s; four axial steps imply dt = .25 s and 2 s.
    near(runs[0].dt, .25); near(runs[1].dt, 2);
    near(runs[0].species.at(0).r, 4.9e-10 * .25 / (5e-4 * 5e-4) * 64);
    near(runs[1].species.at(0).r, 4.9e-10 * 2 / (.001 * .001) * 64);
    session.run();
    for (const auto& exp : session.parameters().experiments) {
        near(exp.scale_factor, exp.run->name == "wide" ? 125 : 62.5);
        // Constant inlet + no reaction + no-flux walls stays uniform for either geometry.
        for (double value : exp.model_profile) near(value, 20);
        for (double value : exp.error) near(value, 0);
    }
    const auto report = session.export_results(files.root);
    std::ifstream input(report); std::string line;
    std::vector<double> last_positions;
    while (std::getline(input, line)) {
        if (line.rfind("res_time", 0) != 0) continue;
        std::istringstream row(line); std::string label;
        for (int i = 0; i < 9; ++i) row >> label;
        double position = 0, last = 0;
        while (row >> position) last = position;
        last_positions.push_back(last);
    }
    ASSERT_EQ(last_positions.size(), 2u);
    near(last_positions[0], 437.5); near(last_positions[1], 875);
    session.save_model_profiles();
    near(scalar(db, "SELECT MAX(X) FROM model_profile JOIN solutions USING(SOLUTION_ID) WHERE EXPERIMENT_NAME='wide'"), 1000);
}

TEST(ChannelDimensions, RejectsPartialAndInvalidSchemasBeforeCreatingSolutionRecords) {
    fixture files;
    {
        run_session session(files.database);
        sql(session.database(), "ALTER TABLE experiments ADD COLUMN CHANNEL_WIDTH REAL");
        EXPECT_THROW(session.load_inputs(), workflow_error);
        near(scalar(session.database(), "SELECT count(*) FROM solutions"), 1);
    }
    {
        run_session session(files.database);
        sql(session.database(), "ALTER TABLE experiments ADD COLUMN CHANNEL_HEIGHT REAL; ALTER TABLE experiments ADD COLUMN CHANNEL_LENGTH REAL");
        EXPECT_THROW(session.load_inputs(), workflow_error);
        near(scalar(session.database(), "SELECT count(*) FROM solutions"), 1);
    }
}

void verify_solution(run_session& session, const run_result& result)
{
    EXPECT_EQ(session.state(), session_state::completed);
    EXPECT_TRUE(result.optimizer_ran);
    ASSERT_TRUE(result.termination_type);
    EXPECT_GT(*result.termination_type, 0);
    EXPECT_GT(result.residual_evaluations, 0);
    ASSERT_EQ(result.parameters.size(), 1u);
    EXPECT_EQ(result.parameters[0].name, "keq1");
    // keq1 is intentionally unidentifiable when kon=0; only workflow is tested.
    near(result.parameters[0].value, 1);
    auto& p = session.parameters();
    ASSERT_EQ(p.experiments.size(), 4u);
    for (const auto& exp : p.experiments) {
        ASSERT_EQ(exp.species_out.size(), 3u);
        for (std::size_t species = 0; species < 3; ++species) {
            ASSERT_EQ(exp.species_out[species].size(), 8u);
            for (double value : exp.species_out[species]) near(value, species == 0 ? 20 : 0);
        }
        for (double value : exp.model_profile) near(value, 20);
        for (double value : exp.experimental_profile) near(value, 1);
        for (double value : exp.error) near(value, 0);
    }
}

void verify_persistence(run_session& session, const std::filesystem::path& directory)
{
    auto* db = session.database();
    near(scalar(db, "SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    near(scalar(db, "SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    auto path = session.export_results(directory);
    EXPECT_EQ(path.parent_path(), directory);
    std::ifstream report(path);
    std::string line;
    int dye_rows = 0;
    while (std::getline(report, line)) {
        std::istringstream row(line);
        std::string field;
        for (int i = 0; i < 9 && row >> field; ++i) {
            if (i == 8 && field == "Free_Dye") {
                double value;
                int count = 0;
                while (row >> value) { near(value, .002); ++count; }
                EXPECT_EQ(count, 8);
                ++dye_rows;
            }
        }
    }
    EXPECT_EQ(dye_rows, 4);
    for (int repeat = 0; repeat < 2; ++repeat) {
        session.save_model_profiles();
        near(scalar(db, "SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 36);
        near(scalar(db, "SELECT count(*) FROM parameter_solutions WHERE SOLVE_SETTING_ID<>99"), 1);
        near(scalar(db, "SELECT max(abs(Free_Dye-0.002)) FROM model_profile WHERE SOLUTION_ID<>99 AND Free_Dye<>''"), 0);
        near(scalar(db, "SELECT max(abs(Numeric-0.002)) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
        near(scalar(db, "SELECT max(abs(Error)) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
    }
    near(scalar(db, "SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    session.save_fitted_parameters();
    near(scalar(db, "SELECT [INITIAL VALUE] FROM alglib_input"), session.parameters().initial_values_alglib[0]);
    near(scalar(db, "SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99"), 42);
    near(scalar(db, "SELECT VALUE FROM parameter_solutions WHERE SOLUTION_ID=99"), 42);
    near(scalar(db, "SELECT MODEL_DA FROM solutions WHERE SOLUTION_ID=99"), 4);
}

TEST(Workflow, InitialParameterSnapshotPrecedesSolveAndOwnsStartingValues)
{
    for (bool optimize : {false, true}) {
        fixture files;
        std::vector<progress_event> events;
        run_session session(files.database, [&](const progress_event& event) { events.push_back(event); });
        sql(session.database(), std::string("INSERT INTO model_controls VALUES('run_solver','") +
            (optimize ? "true" : "false") + "')");
        session.load_inputs();
        auto initial = std::find_if(events.begin(), events.end(), [](const progress_event& event) {
            return event.kind == event_kind::parameters_initialized;
        });
        ASSERT_NE(initial, events.end());
        ASSERT_EQ(initial->parameters.size(), 1u);
        EXPECT_EQ(initial->parameters[0].source, "global");
        EXPECT_EQ(initial->parameters[0].name, "keq1");
        // The fixture's explicit reaction Ks="0 1" initializes the linked keq1
        // to 1, taking precedence over the Variables table's 0.75 default.
        near(initial->parameters[0].initial_value, 1);
        EXPECT_FALSE(std::any_of(events.begin(), events.end(), [](const progress_event& event) {
            return event.message.find("Initial parameter ") != std::string::npos;
        }));
        const auto result = session.run();
        ASSERT_EQ(result.parameters.size(), 1u);
        EXPECT_EQ(result.optimizer_ran, optimize);
        if (!optimize) EXPECT_EQ(result.residual_evaluations, 1);
        for (const auto& experiment : session.parameters().experiments)
            for (double value : experiment.model_profile) near(value, 20);
        near(result.parameters[0].initial_value, 1);
        near(result.parameters[0].value, 1);
        ASSERT_TRUE(result.sum_squared_residuals.has_value());
        near(*result.sum_squared_residuals, 0); // Uniform equilibrium matches exactly.
        std::size_t reports = 0;
        for (const auto& event : events) {
            if (event.kind != event_kind::optimizer_progress) continue;
            ++reports;
            ASSERT_TRUE(event.sum_squared_residuals.has_value());
            near(*event.sum_squared_residuals, 0);
            ASSERT_EQ(event.parameters.size(), 1u);
            EXPECT_EQ(event.parameters[0].source, "global");
            EXPECT_EQ(event.parameters[0].name, "keq1");
            near(event.parameters[0].initial_value, 1);
            near(event.parameters[0].value, 1);
            EXPECT_GT(event.evaluations.value_or(0), 0);
        }
        EXPECT_EQ(reports > 0, optimize);
        session.parameters().initial_values_alglib[0] = 2;
        near(result.parameters[0].initial_value, 1); // Owned, not a live reference.
        for (const auto& event : events)
            if (event.kind == event_kind::optimizer_progress)
                near(event.parameters[0].initial_value, 1);
    }
}

TEST(Workflow, SynchronousLoadSolveExportSaveAndReload)
{
    fixture files;
    {
        run_session session(files.database);
        session.load_inputs();
        const auto result = session.run();
        verify_solution(session, result);
        verify_persistence(session, files.root / "reports");
    }
    run_session again(files.database);
    again.load_inputs();
    verify_solution(again, again.run());
    near(scalar(again.database(), "SELECT count(*) FROM solve_settings"), 2);
    near(scalar(again.database(), "SELECT count(*) FROM solutions"), 5);
}

TEST(Workflow, GeneratedReportSeparatesAnalyticalZeroBoundBeadsAndSpecies)
{
    fixture files;
    run_session session(files.database);
    EXPECT_THROW(session.generate_report(), workflow_error);
    session.load_inputs();
    EXPECT_THROW(session.generate_report(), workflow_error);
    session.run();
    auto& p = session.parameters();
    p.row_count = 1;
    auto& exp = p.experiments.front();
    auto& run = *exp.run;
    // Independent export fixture: 12 bound dye / 3 dye per bead = 4 mg/ml;
    // density 2 g/ml means 1 wt% = 20 mg/ml, so bound beads = 0.2 wt%.
    run.solution_density = 2;
    auto& bead = run.species.at(run.PS_beads);
    bead.model_units = "mg/ml"; bead.input_units = "wt%";
    run.reactions.at(run.FITC_Bead_1).coef.at(&run.species.at(run.FITC)).value() = -3;
    exp.species_out.at(run.Bound_Dye_1).assign(p.X, 12);
    exp.species_out.at(run.PS_beads).assign(p.X, 2);
    for (int i = 0; i < exp.window_size; ++i) exp.analytical_zero[i] = .125 * i;
    exp.channel_position[0] = 1e-310;
    const auto changes = sqlite3_total_changes(session.database());
    const auto generated = session.generate_report();
    EXPECT_EQ(sqlite3_total_changes(session.database()), changes);
    EXPECT_EQ(exp.channel_position[0], 1e-310);
    EXPECT_FALSE(std::filesystem::exists(files.root / "reports"));
    auto profile = [](const std::string& report, const std::string& name) {
        std::istringstream input(report);
        std::string line;
        while (std::getline(input, line)) {
            std::istringstream row(line); std::string label;
            for (int i = 0; i < 9; ++i) row >> label;
            if (label != name) continue;
            std::vector<double> values; double value;
            while (row >> value) values.push_back(value);
            return values;
        }
        return std::vector<double>{};
    };
    const auto bound = profile(generated, "Bound_Beads_(wt%)");
    ASSERT_EQ(bound.size(), p.X);
    for (double value : bound) near(value, .2);
    const auto analytic = profile(generated, "Analytical_Zero_(umol)");
    ASSERT_EQ(analytic.size(), exp.window_size);
    for (int i = 0; i < exp.window_size; ++i) near(analytic[i], .125 * i);
    const auto species = profile(generated, "Species:PS_40nm_(mg/ml)");
    ASSERT_EQ(species.size(), p.X);
    for (double value : species) near(value, 2);
    const auto total = profile(generated, "Total_Beads_(wt%)");
    ASSERT_EQ(total.size(), p.X);
    for (double value : total) near(value, .3);
    // Existing CLI/export contract stays unchanged, including the historical label.
    const auto legacyPath = session.export_results(files.root / "legacy");
    std::ifstream legacyFile(legacyPath);
    const std::string legacy(std::istreambuf_iterator<char>(legacyFile), {});
    EXPECT_EQ(profile(legacy, "Bound_Beads_(wt%)"), analytic);
    EXPECT_EQ(legacy.find("Species:"), std::string::npos);
    EXPECT_EQ(legacy.find("Analytical_Zero"), std::string::npos);
}

TEST(Workflow, BackgroundSuccessOutlivesRunner)
{
    fixture files;
    background_result outcome;
    {
        background_runner runner;
        runner.start(files.database);
        runner.wait();
        ASSERT_EQ(runner.status(), background_state::completed);
        outcome = runner.take_result();
    }
    ASSERT_FALSE(outcome.failure);
    ASSERT_TRUE(outcome.session);
    ASSERT_TRUE(outcome.result);
    verify_solution(*outcome.session, *outcome.result);
    verify_persistence(*outcome.session, files.root / "background reports");
}

TEST(Workflow, CancelDuringRealSolveThenRunAgain)
{
    fixture files;
    background_runner runner;
    std::atomic<bool> accepted = false;
    runner.start(files.database, [&](const progress_event& event) {
        if (event.kind == event_kind::evaluation && event.evaluations == 1)
            accepted = runner.request_cancel();
    });
    runner.wait();
    EXPECT_TRUE(accepted);
    EXPECT_EQ(runner.status(), background_state::cancelled);
    auto cancelled = runner.take_result();
    EXPECT_FALSE(cancelled.session);
    EXPECT_FALSE(cancelled.result);
    ASSERT_TRUE(cancelled.failure);
    try { std::rethrow_exception(cancelled.failure); }
    catch (const workflow_error& error) { EXPECT_EQ(error.code, error_code::cancelled); }
    {
        run_session inspect(files.database);
        near(scalar(inspect.database(), "SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
        near(scalar(inspect.database(), "SELECT [INITIAL VALUE] FROM alglib_input"), .75);
    }
    runner.start(files.database);
    runner.wait();
    ASSERT_EQ(runner.status(), background_state::completed);
    auto completed = runner.take_result();
    ASSERT_TRUE(completed.session);
    verify_solution(*completed.session, *completed.result);
}

TEST(Workflow, InvalidInputsFailWithoutResultsAndAllowFreshRun)
{
    fixture files;
    {
        run_session edit(files.database);
        sql(edit.database(), "UPDATE alglib_input SET SCALE=0");
    }
    background_runner runner;
    runner.start(files.database);
    runner.wait();
    EXPECT_EQ(runner.status(), background_state::failed);
    auto failed = runner.take_result();
    EXPECT_FALSE(failed.session);
    ASSERT_TRUE(failed.failure);
    try { std::rethrow_exception(failed.failure); }
    catch (const workflow_error& error) {
        EXPECT_EQ(error.action, operation::load_inputs);
        EXPECT_EQ(error.code, error_code::invalid_input);
    }
    {
        run_session edit(files.database);
        near(scalar(edit.database(), "SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99"), 0);
        near(scalar(edit.database(), "SELECT [INITIAL VALUE] FROM alglib_input"), .75);
        sql(edit.database(), "UPDATE alglib_input SET SCALE=1");
    }
    runner.start(files.database);
    runner.wait();
    ASSERT_EQ(runner.status(), background_state::completed);
    auto completed = runner.take_result();
    ASSERT_TRUE(completed.session);
    verify_solution(*completed.session, *completed.result);
}
} // namespace
