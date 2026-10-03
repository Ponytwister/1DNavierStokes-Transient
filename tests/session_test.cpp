#include <workflow.h>
#include <gtest/gtest.h>

#include <type_traits>

namespace {
using tsensor_workflow::run_session;
using tsensor_workflow::session_state;
using tsensor_workflow::error_code;
using tsensor_workflow::operation;

template<class Action>
void expect_workflow_error(Action action, error_code code, operation expected_operation)
{
    try {
        action();
        FAIL() << "Expected a workflow error";
    } catch (const tsensor_workflow::workflow_error& error) {
        EXPECT_EQ(error.code, code);
        EXPECT_EQ(error.action, expected_operation);
    }
}

static_assert(!std::is_copy_constructible_v<run_session>);
static_assert(!std::is_move_constructible_v<run_session>);
static_assert(!std::is_copy_constructible_v<parameters_t>);
static_assert(!std::is_move_constructible_v<parameters_t>);

void execute(sqlite3* db, const char* sql)
{
    const int rc = sqlite3_exec(db, sql, nullptr, nullptr, nullptr);
    if (rc != SQLITE_OK) { throw std::runtime_error(sqlite3_errmsg(db)); }
}

void create_inputs(run_session& session, int cells, double initial)
{
    execute(session.database(), R"(
        CREATE TABLE model_controls (criterion TEXT PRIMARY KEY, value TEXT);
        CREATE TABLE alglib_input (
            VARIABLE TEXT, "INITIAL VALUE" REAL, "LOWER BOUND" REAL,
            "UPPER BOUND" REAL, SCALE REAL);
    )");
    const auto sql = "INSERT INTO model_controls VALUES ('width resolution (X)', '"
        + std::to_string(cells) + "'), ('universal_solve_for', 'width');"
        + "INSERT INTO alglib_input VALUES ('width', " + std::to_string(initial)
        + ", 0.1, 10, 1);";
    execute(session.database(), sql.c_str());
    session.parameters().debug_level = 7;
    read_model_parameters_from_db(session.parameters(), session.database());
    read_alglib_values_from_db(session.parameters(), session.database());
}

TEST(RunSession, IndependentInputsBuffersAndParameterLinks)
{
    run_session first(":memory:");
    create_inputs(first, 4, 2.0);
    auto& p = first.parameters();
    ASSERT_EQ(p.solvables.size(), 1u);
    p.solvables.front().value() = 2.0;
    p.experiment_runs.emplace_back();
    p.experiment_runs.front().width.source = &p.solvables.front();
    p.experiments.resize(1);
    p.experiments.front().run = &p.experiment_runs.front();
    p.experiments.front().width.source = &p.experiment_runs.front().width;
    p.experiments.front().species_out.assign(2, std::vector<double>(4, 0.0));
    auto* linked_source = &p.solvables.front();

    {
        run_session second(":memory:");
        create_inputs(second, 9, 5.0);
        EXPECT_EQ(second.parameters().X, 9);
        EXPECT_EQ(second.parameters().initial_values_alglib, std::vector<double>({5.0}));
        EXPECT_EQ(p.X, 4);
        EXPECT_EQ(p.initial_values_alglib, std::vector<double>({2.0}));
        p.experiments.front().width.value() = 3.0;
        EXPECT_DOUBLE_EQ(p.solvables.front().value(), 3.0);
        EXPECT_DOUBLE_EQ(second.parameters().initial_values_alglib.front(), 5.0);
    }
    EXPECT_EQ(p.experiment_runs.front().width.source, linked_source);
    EXPECT_DOUBLE_EQ(p.experiments.front().width.value(), 3.0);
    EXPECT_EQ(p.experiments.front().species_out[1], std::vector<double>(4, 0.0));
    run_session next(":memory:");
    EXPECT_TRUE(next.parameters().solvables.empty());
    EXPECT_TRUE(next.parameters().initial_values_alglib.empty());
    EXPECT_TRUE(next.parameters().state.empty());
    EXPECT_EQ(next.parameters().iterations, 0);
    EXPECT_EQ(next.parameters().SOLVE_SETTING_ID, 0);
    EXPECT_FALSE(next.parameters().SOLUTION_ID_RECURSIVE_CALL);
}

TEST(RunSession, InletLoaderOwnsZeroInitializedConcentrationBuffers)
{
    run_session session(":memory:");
    auto& p = session.parameters();
    p.debug_level = 7;
    p.X = 4;
    // The existing loader prints diagnostic details for row 3.
    p.row_count = 4;
    p.experiment_runs.emplace_back();
    auto& run = p.experiment_runs.front();
    run.number_of_species = 1;
    run.FITC = 0;
    run.dye_conc_mgml = 0;
    run.species.resize(1);
    run.species.front().name = "FITC";
    run.species.front().type = 1;
    run.species.front().molecular_weight = 100;
    run.species.front().input_units = "mg/ml";
    run.species.front().model_units = "umol";
    p.experiments.resize(4);
    for (auto& experiment : p.experiments) {
        experiment.run = &run;
        experiment.INLET_COND_ID = 1;
        experiment.entrances.resize(1);
    }
    execute(session.database(), R"(
        CREATE TABLE inlet_conditions (
            INLET_COND_ID INTEGER, SPECIE_CONC REAL,
            ENTRANCE_NUMBER INTEGER, SPECIES_NAME TEXT);
        INSERT INTO inlet_conditions VALUES (1, 2.0, 1, 'FITC');
    )");
    read_inlet_cond_from_db(p, session.database());
    for (const auto& experiment : p.experiments) {
        ASSERT_EQ(experiment.species_out.size(), 1u);
        EXPECT_EQ(experiment.species_out.front(), std::vector<double>(4, 0.0));
        EXPECT_DOUBLE_EQ(experiment.entrances.front().CONC.at(&run.species.front()), 2.0);
    }
}

TEST(RunSession, FailedLoadCannotBeReusedOrSaved)
{
    run_session session(":memory:");
    session.parameters().debug_level = 7;
    expect_workflow_error([&] { session.run(); }, error_code::invalid_state, operation::solve);
    expect_workflow_error([&] { session.export_results("unused"); }, error_code::invalid_state, operation::export_results);
    expect_workflow_error([&] { session.save_fitted_parameters(); }, error_code::invalid_state, operation::save_fitted_parameters);
    expect_workflow_error([&] { session.load_inputs(); }, error_code::database, operation::load_inputs);
    EXPECT_EQ(session.state(), session_state::failed);
    expect_workflow_error([&] { session.load_inputs(); }, error_code::invalid_state, operation::load_inputs);
    expect_workflow_error([&] { session.run(); }, error_code::invalid_state, operation::solve);
    expect_workflow_error([&] { session.save_model_profiles(); }, error_code::invalid_state, operation::save_model_profiles);
    EXPECT_EQ(sqlite3_next_stmt(session.database(), nullptr), nullptr);
}

TEST(RunSession, CallbackExceptionFinalizesStatementBeforeRethrowing)
{
    run_session session(":memory:");
    session.parameters().debug_level = 7;
    execute(session.database(), R"(
        CREATE TABLE model_controls (criterion TEXT, value TEXT);
        INSERT INTO model_controls VALUES ('width resolution (X)', 'invalid integer');
    )");
    expect_workflow_error([&] { session.load_inputs(); }, error_code::invalid_input, operation::load_inputs);
    EXPECT_EQ(session.state(), session_state::failed);
    EXPECT_EQ(sqlite3_next_stmt(session.database(), nullptr), nullptr);
    // SQLite remains usable: the failure did not unwind through its C stack.
    EXPECT_NO_THROW(execute(session.database(), "DROP TABLE model_controls;"));
}

TEST(RunSession, RepeatedSuccessAndPartialFailureReleaseDatabaseMemory)
{
    ASSERT_EQ(sqlite3_initialize(), SQLITE_OK);
    const auto baseline = sqlite3_memory_used();
    for (int iteration = 0; iteration < 12; ++iteration) {
        {
            run_session session(":memory:");
            ASSERT_GT(sqlite3_memory_used(), baseline);
            create_inputs(session, 4, 2.0);
            // Fail in the SQLite callback after all optimizer buffers allocate.
            execute(session.database(), "UPDATE alglib_input SET SCALE = 0;");
            EXPECT_THROW(read_alglib_values_from_db(session.parameters(), session.database()),
                         std::runtime_error);
            EXPECT_EQ(session.parameters().scale.size(), 1u);
            EXPECT_EQ(sqlite3_next_stmt(session.database(), nullptr), nullptr);
        }
        EXPECT_EQ(sqlite3_memory_used(), baseline);
        {
            run_session session(":memory:");
            session.parameters().debug_level = 7;
            // Exercise allocated SQLite error strings as well as callbacks.
            EXPECT_THROW(session.load_inputs(), std::runtime_error);
        }
        EXPECT_EQ(sqlite3_memory_used(), baseline);
    }
}

TEST(RunSession, ProfileSamplingReportsTheExperimentBeforeReadingPastItsEnd)
{
    run_session session(":memory:");
    auto& p = session.parameters();
    p.debug_level = 7;
    p.row_count = 0; // Isolate profile preparation; no transient workers.
    p.experiment_runs.resize(1);
    p.experiment_runs.front().name = "sample-run";
    p.experiments.resize(1);
    auto& exp = p.experiments.front();
    exp.run = &p.experiment_runs.front();
    exp.second_name = "sample-profile";
    exp.window_size = 5;
    exp.raw_experimental_profile = {10, 20, 30, 40, 50};
    exp.experimental_profile.resize(5);
    exp.channel_position.resize(5);
    exp.scale_factor = 1;
    exp.left_edge.value() = 1.2;
    exp.width.value() = 3.2;
    alglib::real_1d_array controls, residuals;
    // ceil(3 + 1.2) == 5, despite left_edge + width < sample count.
    try {
        alglib_solver(controls, residuals, &p);
        FAIL() << "Expected profile-domain error";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("sample-run:sample-profile"), std::string::npos);
        EXPECT_NE(std::string(error.what()).find("available samples=5"), std::string::npos);
    }
    EXPECT_EQ(p.iterations, 0);
    // Last valid sample is 4; the existing tail repeats that sample.
    exp.width.value() = 3;
    EXPECT_NO_THROW(alglib_solver(controls, residuals, &p));
    EXPECT_EQ(exp.experimental_profile, (std::vector<double>{30, 40, 50, 50, 50}));
    exp.left_edge.value() = 1;
    exp.width.value() = 4;
    EXPECT_NO_THROW(alglib_solver(controls, residuals, &p));
    EXPECT_EQ(exp.experimental_profile, (std::vector<double>{20, 30, 40, 50, 50}));
}

TEST(RunSession, SolverWorkersReturnExceptionsToTheirCaller)
{
    run_session session(":memory:");
    auto& p = session.parameters();
    p.debug_level = 7;
    p.row_count = 2;
    // Each worker fails at experiments.at(row), before touching model inputs.
    alglib::real_1d_array controls, residuals;
    EXPECT_THROW(alglib_solver(controls, residuals, &p), std::out_of_range);
    EXPECT_EQ(p.iterations, 0);
    EXPECT_TRUE(p.state.empty());
    EXPECT_THROW(alglib_solver(controls, residuals, nullptr), std::invalid_argument);
}
} // namespace
