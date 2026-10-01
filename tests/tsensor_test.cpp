#include <tsensor.h>
#include <gtest/gtest.h>
#include <iterator>
#include <memory>
#include <new>
#include <sstream>

// The application owns this global; tests supply their own instance.
parameters_t p;

namespace {
constexpr double absolute_tolerance = 1e-12;
constexpr double relative_tolerance = 1e-10;

void expect_numeric(double actual, double expected)
{
    ASSERT_TRUE(std::isfinite(actual));
    EXPECT_NEAR(actual, expected,
                absolute_tolerance + relative_tolerance * std::abs(expected));
}

std::string fixture_text(const char* name)
{
    std::ifstream input(std::string(TSENSOR_FIXTURE_DIR) + "/" + name);
    if (!input) {
        throw std::runtime_error(std::string("Cannot open fixture: ") + name);
    }
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

TEST(Solvable, LinkedWritesReachSourceAndUnlinkKeepsValue)
{
    solvable root(2.0);
    solvable middle(&root);
    solvable leaf(&middle);
    leaf.value() = 7.0;
    EXPECT_DOUBLE_EQ(root.value(), 7.0);
    leaf.unlink();
    root.value() = 9.0;
    EXPECT_FALSE(leaf.is_linked());
    EXPECT_DOUBLE_EQ(leaf.value(), 7.0);
}

TEST(Solvable, CyclesAreRejected)
{
    solvable first(1.0);
    solvable second(&first);
    first.source = &second;
    EXPECT_THROW(first.value(), std::runtime_error);
    EXPECT_THROW(first.init(), std::runtime_error);
}

TEST(Interpolation, LinearProfileIncludesEndpointsAndReversedCoordinates)
{
    // Independent analytic profile y = 3x + 2.
    for (double x : {0.0, 0.25, 0.5, 1.0}) {
        expect_numeric(lin_interpolate(x, 0.0, 2.0, 1.0, 5.0), 3.0 * x + 2.0);
        expect_numeric(lin_interpolate(x, 1.0, 5.0, 0.0, 2.0), 3.0 * x + 2.0);
    }
}

TEST(Interpolation, ExtrapolationIsRejected)
{
    EXPECT_THROW(lin_interpolate(-0.1, 0, 2, 1, 5), std::runtime_error);
    EXPECT_THROW(lin_interpolate(1.1, 0, 2, 1, 5), std::runtime_error);
}

TEST(UnitConversion, MolecularMassAndRoundTrip)
{
    experiment_run_struct run{};
    run.species.resize(1);
    run.species[0].type = 1;
    run.species[0].molecular_weight = 100.0;
    const double forward = unit_conversion(&run, 0, "mg/ml", "umol");
    expect_numeric(forward, 10000.0);
    expect_numeric(forward * unit_conversion(&run, 0, "umol", "mg/ml"), 1.0);
    EXPECT_THROW(unit_conversion(&run, 0, "unsupported", "mg/ml"), std::runtime_error);
}

class ModelInputs : public ::testing::Test {
protected:
    void SetUp() override
    {
        // parameters_t has const members; reconstruct to isolate each test.
        p.~parameters_t();
        new (&p) parameters_t{};
        p.debug_level = 7;
    }
};

TEST_F(ModelInputs, InletProfileMatchesAnalyticFixture)
{
    std::istringstream input(fixture_text("inlet_profile.txt"));
    int cells;
    double molecular_weight, flow_a, flow_b, concentration_a, concentration_b;
    ASSERT_TRUE(static_cast<bool>(input >> cells >> molecular_weight >> flow_a >> flow_b
                                       >> concentration_a >> concentration_b));
    ASSERT_GT(cells, 0);
    p.X = cells;
    experiment_run_struct run{};
    run.number_of_species = 1;
    run.total_flowrate = flow_a + flow_b;
    run.species.resize(1);
    auto& species = run.species[0];
    species.type = 1;
    species.molecular_weight = molecular_weight;
    species.input_units = "mg/ml";
    species.model_units = "umol";
    experiment_struct experiment{};
    experiment.run = &run;
    experiment.entrances = {{{{&species, concentration_a}}, flow_a},
                           {{{&species, concentration_b}}, flow_b}};
    std::vector<double> solution(2 * cells, -1.0);
    set_inlet_conc(&experiment, solution.data());
    for (int i = 0; i < cells; ++i) {
        double expected;
        ASSERT_TRUE(static_cast<bool>(input >> expected));
        SCOPED_TRACE(i);
        expect_numeric(solution[i], expected);
        expect_numeric(solution[cells + i], expected);
    }
    std::string extra;
    EXPECT_FALSE(static_cast<bool>(input >> extra));
}

TEST_F(ModelInputs, LoadsControlsFromDisposableDatabase)
{
    sqlite3* raw = nullptr;
    const int result = sqlite3_open(":memory:", &raw);
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
    ASSERT_EQ(result, SQLITE_OK);
    const auto sql = fixture_text("model_controls.sql");
    ASSERT_EQ(sqlite3_exec(db.get(), sql.c_str(), nullptr, nullptr, nullptr), SQLITE_OK)
        << sqlite3_errmsg(db.get());
    read_model_parameters_from_db(db.get());
    EXPECT_EQ(p.X, 4);
    EXPECT_EQ(p.Z, 2);
    EXPECT_TRUE(p.disable_reactions);
    EXPECT_FALSE(p.save_model_profiles);
    EXPECT_EQ(p.max_iterations, 3);
    expect_numeric(p.convergence_epsx, 1e-9);
}
} // namespace
