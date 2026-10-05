#include <workflow.h>
#include <application_options.h>
#include <cstdio>
#include <functional>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {
bool terminal_output() {
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

class ProgressPrinter {
public:
    ~ProgressPrinter() { finish_line(); }
    void operator()(const tsensor_workflow::progress_event& event) {
        using tsensor_workflow::event_kind;
        if (event.kind == event_kind::parameters_initialized) return;
        if (event.kind == event_kind::message && event.detail_level < 3) return;
        if (event.kind == event_kind::evaluation) {
            if (interactive_) {
                std::cout << '\r' << "Model evaluations completed: " << event.evaluations.value() << std::flush;
                pending_line_ = true;
            }
            return;
        }
        finish_line();
        if (event.kind == event_kind::optimizer_progress) {
            std::cout << "Sum of squared residuals: " << event.sum_squared_residuals.value() << '\n';
            for (const auto& parameter : event.parameters)
                std::cout << '(' << parameter.source << ':' << parameter.name << ':' << parameter.value << ") ";
            std::cout << '\n';
            return;
        }
        // Failure exceptions are printed once by main's error handler.
        if (event.kind == event_kind::failed) return;
        std::cout << event.message << '\n';
    }
    void final_count(std::size_t count) {
        if (!interactive_) std::cout << "Model evaluations completed: " << count << '\n';
    }
private:
    void finish_line() {
        if (pending_line_) { std::cout << '\n'; pending_line_ = false; }
    }
    bool interactive_ = terminal_output();
    bool pending_line_ = false;
};
} // namespace

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
        ProgressPrinter print_progress;
        tsensor_workflow::run_session session(options.database, std::ref(print_progress));
        std::cout << "Database opened: " << options.database << '\n';
        // Detect invalid directory paths before loading or running the model.
        std::filesystem::create_directories(options.output_directory);
        std::cout << "Output directory: " << options.output_directory << '\n';

        session.load_inputs();
        const auto result = session.run();
        print_progress.final_count(result.residual_evaluations);
        if (result.sum_squared_residuals)
            std::cout << "Sum of squared residuals: " << *result.sum_squared_residuals << '\n';
        if (result.optimizer_ran) {
            std::cout << "Optimizer iterations: " << result.optimizer_iterations.value()
                      << "; termination code: " << result.termination_type.value() << '\n';
        }
        for (const auto& parameter : result.parameters) {
            std::cout << '(' << parameter.source << ':' << parameter.name << ':' << parameter.value << ") ";
        }
        std::cout << '\n';
        session.export_results(options.output_directory);
        session.save_model_profiles();

        std::string accept_alglib_values;
        std::cout << "Save parameter solutions to db as initial values? (y/n) ";
        std::cin >> accept_alglib_values;
        if (accept_alglib_values == "y") {
            session.save_fitted_parameters();
        }
        if (session.progress_failure()) {
            std::cerr << "Progress display failed; the model operations continued.\n";
        }
    }
    catch (const tsensor_workflow::workflow_error& error)
    {
        std::cerr << tsensor_workflow::operation_name(error.action) << " failed: " << error.what() << '\n';
        return 1;
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
