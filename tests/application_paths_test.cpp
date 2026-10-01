#include <application_options.h>
#include <workflow.h>
#include <gtest/gtest.h>

#include <chrono>
#include <iterator>

namespace {
namespace fs = std::filesystem;

TEST(ApplicationOptions, PreservesLegacyDefaults)
{
    const char* args[] = {"Navier"};
    const auto options = tsensor_workflow::parse_options(1, args);
    EXPECT_EQ(options.database, fs::path("../../navier.db"));
    EXPECT_EQ(options.output_directory, fs::path(".."));
    EXPECT_FALSE(options.help);
}

TEST(ApplicationOptions, AcceptsExplicitPathsWithSpaces)
{
    const char* args[] = {"Navier", "--output-dir", "results with spaces",
                          "--database", "data with spaces/experiment.db"};
    const auto options = tsensor_workflow::parse_options(5, args);
    EXPECT_EQ(options.database, fs::path("data with spaces/experiment.db"));
    EXPECT_EQ(options.output_directory, fs::path("results with spaces"));
}

TEST(ApplicationOptions, RejectsUnknownAndMissingArguments)
{
    const char* unknown[] = {"Navier", "--databse", "data.db"};
    EXPECT_THROW(tsensor_workflow::parse_options(3, unknown), std::invalid_argument);
    const char* missing[] = {"Navier", "--database"};
    EXPECT_THROW(tsensor_workflow::parse_options(2, missing), std::invalid_argument);
    const char* empty[] = {"Navier", "--output-dir", ""};
    EXPECT_THROW(tsensor_workflow::parse_options(3, empty), std::invalid_argument);
    const char* next_option[] = {"Navier", "--database", "--output-dir", "results"};
    EXPECT_THROW(tsensor_workflow::parse_options(4, next_option), std::invalid_argument);
}

class ApplicationPaths : public ::testing::Test {
protected:
    parameters_t p{};
    fs::path root;
    fs::path original_directory;
    bool owns_directory = false;

    void SetUp() override
    {
        original_directory = fs::current_path();
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        root = original_directory / ("path-test-" + std::to_string(stamp));
        owns_directory = fs::create_directory(root);
        ASSERT_TRUE(owns_directory);
        p.X = 1;
        // Header-only export: no solver or experimental data is needed.
        p.row_count = 0;
        p.experiment_runs.emplace_back();
        p.experiment_runs.back().name = "sample run";
    }

    void TearDown() override
    {
        fs::current_path(original_directory);
        std::error_code error;
        if (owns_directory) {
            fs::remove_all(root, error);
        }
    }
};

TEST_F(ApplicationPaths, ExportsToAbsoluteDirectoryFromAnotherWorkingDirectory)
{
    fs::create_directory(root / "launch");
    fs::current_path(root / "launch");
    const auto destination = root / "results with spaces" / "nested";
    const auto result = tsensor_workflow::export_results(p, destination);
    EXPECT_EQ(result, destination / "sample run,.txt");
    ASSERT_TRUE(fs::is_regular_file(result));
    std::ifstream input(result);
    const std::string text{std::istreambuf_iterator<char>(input), {}};
    EXPECT_EQ(text.find("res_time  bind_ratio(p1)"), 0u);
    EXPECT_FALSE(fs::exists(root / "sample run,.txt"));
    EXPECT_EQ(fs::current_path(), root / "launch");
}

TEST_F(ApplicationPaths, RelativeOutputUsesWorkingDirectory)
{
    fs::current_path(root);
    tsensor_workflow::export_results(p, "relative results");
    EXPECT_TRUE(fs::is_regular_file(root / "relative results" / "sample run,.txt"));
}

TEST_F(ApplicationPaths, ReportsInvalidDirectoryAndFileOpenFailures)
{
    const auto blocked = root / "a file";
    { std::ofstream file(blocked); file << "keep"; }
    EXPECT_THROW(tsensor_workflow::export_results(p, blocked / "results"), fs::filesystem_error);
    EXPECT_THROW(tsensor_workflow::export_results(p, {}), std::invalid_argument);
    // A directory at the exact filename causes a deterministic file-open error.
    fs::create_directory(root / "sample run,.txt");
    EXPECT_THROW(tsensor_workflow::export_results(p, root), std::runtime_error);
}

TEST_F(ApplicationPaths, RejectsExperimentNamesThatEscapeOutputDirectory)
{
    p.experiment_runs.back().name = "../escaped";
    EXPECT_THROW(tsensor_workflow::export_results(p, root / "results"), std::invalid_argument);
    EXPECT_FALSE(fs::exists(root / "escaped,.txt"));
}

} // namespace
