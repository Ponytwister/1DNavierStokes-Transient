#include <workflow.h>
#include <gtest/gtest.h>

#include <type_traits>

namespace {
using tsensor_workflow::run_session;
using tsensor_workflow::session_state;

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
    EXPECT_THROW(session.run(), std::logic_error);
    EXPECT_THROW(session.export_results("unused"), std::logic_error);
    EXPECT_THROW(session.save_fitted_parameters(), std::logic_error);
    EXPECT_THROW(session.load_inputs(), std::runtime_error); // Missing schema.
    EXPECT_EQ(session.state(), session_state::failed);
    EXPECT_THROW(session.load_inputs(), std::logic_error);
    EXPECT_THROW(session.run(), std::logic_error);
    EXPECT_THROW(session.save_model_profiles(), std::logic_error);
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
    EXPECT_THROW(session.load_inputs(), std::invalid_argument);
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
