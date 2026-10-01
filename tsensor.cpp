#include <workflow.h>
#include <application_options.h>

int
main(int argc, char* argv[])
{
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
        tsensor_workflow::run_session session(options.database);
        std::cout << "Database opened: " << options.database << '\n';
        // Detect invalid directory paths before loading or running the model.
        std::filesystem::create_directories(options.output_directory);
        std::cout << "Output directory: " << options.output_directory << '\n';

        session.load_inputs();
        session.run();
        session.export_results(options.output_directory);
        session.save_model_profiles();

        std::string accept_alglib_values;
        std::cout << "Save parameter solutions to db as initial values? (y/n) ";
        std::cin >> accept_alglib_values;
        if (accept_alglib_values == "y") {
            session.save_fitted_parameters();
        }

    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        return 1;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    
    
    std::cout << "Database closed." << std::endl;
    return 0;
}
