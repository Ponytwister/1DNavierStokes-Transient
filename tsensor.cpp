#include <workflow.h>
#include <application_options.h>
parameters_t p;

int
main(int argc, char* argv[])
{
    sqlite3* db = nullptr; // Owned by this entry point.
    try
    {
        auto options = tsensor_workflow::parse_options(argc, argv);
        if (options.help) {
            std::cout << "Usage: Navier [--database PATH] [--output-dir DIRECTORY]\n"
                      << "       Navier --help\n"
                      << "Defaults: --database ../../navier.db --output-dir ..\n"
                      << "Relative paths use the launch working directory.\n";
            return 0;
        }
        options.database = std::filesystem::absolute(options.database).lexically_normal();
        options.output_directory = std::filesystem::absolute(options.output_directory).lexically_normal();
        // SQLite expects UTF-8 filenames; filesystem paths retain native encoding.
        const auto database_utf8 = options.database.u8string();
        const int rc = sqlite3_open_v2(reinterpret_cast<const char*>(database_utf8.c_str()),
                                      &db, SQLITE_OPEN_READWRITE, nullptr);
        if (rc != SQLITE_OK) {
            std::cerr << "Cannot open database " << options.database << ": "
                      << sqlite3_errmsg(db) << '\n';
            sqlite3_close(db);
            return 1;
        }
        std::cout << "Database opened: " << options.database << '\n';
        // Detect invalid directory paths before loading or running the model.
        std::filesystem::create_directories(options.output_directory);
        std::cout << "Output directory: " << options.output_directory << '\n';

        tsensor_workflow::load_inputs(db);
        tsensor_workflow::run();
        tsensor_workflow::export_results(options.output_directory);
        tsensor_workflow::save_model_profiles(db);

        std::string accept_alglib_values;
        std::cout << "Save parameter solutions to db as initial values? (y/n) ";
        std::cin >> accept_alglib_values;
        if (accept_alglib_values == "y") {
            tsensor_workflow::save_fitted_parameters(db);
        }

        tsensor_workflow::release_run_resources();
    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        sqlite3_close(db);
        return 1;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        sqlite3_close(db);
        return 1;
    }
    
    
    sqlite3_close(db);
    std::cout << "Database closed." << std::endl;
    return 0;
}
