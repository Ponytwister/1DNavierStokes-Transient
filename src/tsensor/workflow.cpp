#include <workflow.h>
#include <type_traits>

namespace tsensor_workflow {

namespace {
[[noreturn]] void report_failure(parameters_t& p, const workflow_error& error)
{
    publish_event(p, {event_kind::failed, error.action, error.what(), 3, {},
                      error.code, error.sqlite_code});
    throw error;
}

template<class Action>
auto perform(parameters_t& p, operation action, Action&& work)
{
    p.active_operation = action;
    publish_event(p, {event_kind::started, action, operation_name(action)});
    try {
        check_cancellation(p);
        if constexpr (std::is_void_v<std::invoke_result_t<Action>>) {
            work();
            publish_event(p, {event_kind::completed, action,
                              std::string(operation_name(action)) + " completed"});
        } else {
            auto result = work();
            publish_event(p, {event_kind::completed, action,
                              std::string(operation_name(action)) + " completed"});
            return result;
        }
    } catch (const workflow_error& error) {
        report_failure(p, error);
    } catch (const alglib::ap_error& error) {
        report_failure(p, {action, error_code::solver, error.msg});
    } catch (const std::filesystem::filesystem_error& error) {
        report_failure(p, {action, error_code::io, error.what()});
    } catch (const std::ios_base::failure& error) {
        report_failure(p, {action, error_code::io, error.what()});
    } catch (const std::invalid_argument& error) {
        report_failure(p, {action, error_code::invalid_input, error.what()});
    } catch (const std::out_of_range& error) {
        report_failure(p, {action, error_code::invalid_input, error.what()});
    } catch (const std::exception& error) {
        const auto code = action == operation::load_inputs ? error_code::invalid_input
                        : action == operation::solve ? error_code::solver : error_code::internal;
        report_failure(p, {action, code, error.what()});
    } catch (...) {
        report_failure(p, {action, error_code::internal, "Unknown model failure"});
    }
}
} // namespace

run_session::run_session(const std::filesystem::path& database_path, progress_callback progress,
                         std::stop_token cancellation)
{
    parameters_.cancellation = cancellation;
    parameters_.progress = std::move(progress);
    perform(parameters_, operation::open_database, [&] {
        if (database_path.empty()) {
            throw std::invalid_argument("Database path must not be empty");
        }
        const auto utf8 = database_path.u8string();
        sqlite3* raw = nullptr;
        const int rc = sqlite3_open_v2(reinterpret_cast<const char*>(utf8.c_str()),
                                      &raw, SQLITE_OPEN_READWRITE, nullptr);
        database_.reset(raw);
        if (rc != SQLITE_OK) {
            throw workflow_error(operation::open_database, error_code::database,
                "Cannot open database " + database_path.string() + ": " + sqlite3_errmsg(raw), rc);
        }
    });
}

void run_session::require_state(session_state expected, operation action) const
{
    if (state_ != expected) {
        throw workflow_error(action, error_code::invalid_state,
            "Invalid session operation order; use a fresh session for a new run");
    }
}

void run_session::load_inputs(const std::optional<control_values>& controls)
{
    perform(parameters_, operation::load_inputs, [&] {
        require_state(session_state::empty, operation::load_inputs);
        state_ = session_state::failed;
        tsensor_workflow::load_inputs(parameters_, database_.get(), controls);
        state_ = session_state::loaded;
    });
}

run_result run_session::run()
{
    return perform(parameters_, operation::solve, [&] {
        require_state(session_state::loaded, operation::solve);
        state_ = session_state::failed;
        auto result = tsensor_workflow::run(parameters_);
        state_ = session_state::completed;
        return result;
    });
}

std::filesystem::path run_session::export_results(const std::filesystem::path& directory)
{
    return perform(parameters_, operation::export_results, [&] {
        require_state(session_state::completed, operation::export_results);
        return tsensor_workflow::export_results(parameters_, directory);
    });
}

void run_session::save_model_profiles()
{
    perform(parameters_, operation::save_model_profiles, [&] {
        require_state(session_state::completed, operation::save_model_profiles);
        tsensor_workflow::save_model_profiles(parameters_, database_.get());
    });
}

void run_session::save_fitted_parameters()
{
    perform(parameters_, operation::save_fitted_parameters, [&] {
        require_state(session_state::completed, operation::save_fitted_parameters);
        tsensor_workflow::save_fitted_parameters(parameters_, database_.get());
    });
}

void load_inputs(parameters_t& p, sqlite3* db, const std::optional<control_values>& controls)
{
    if (controls) {
        for (const auto& [name, value] : *controls) {
            auto key = name;
            auto setting = value;
            char* row[] = {key.data(), setting ? setting->data() : nullptr};
            read_model_parameters_db_callback(&p, 2, row, nullptr);
        }
    } else read_model_parameters_from_db(p, db);
    initialize_channel_dimensions(p, db);
    lines_from_profile_text(p, db);
    read_exp_parameters_from_db(p, db);
    get_solve_settings_ID_from_db(p, db);
    read_specie_and_reaction_values_from_db(p, db);
    read_raw_profiles_from_db(p, db);
    read_inlet_cond_from_db(p, db);
    get_SOLUTION_IDs_from_db(p, db);
    //get_solvable_initial_values_from_db(p, db);
    read_alglib_values_from_db(p, db);
    progress_event initialized{event_kind::parameters_initialized, operation::load_inputs, {}};
    for (std::size_t i = 0; i < p.solvables.size(); ++i) {
        const auto& parameter = p.solvables[i];
        initialized.parameters.push_back({parameter.source_name, parameter.name,
                                         p.initial_values_alglib[i], p.initial_values_alglib[i]});
    }
    publish_event(p, std::move(initialized));
    check_cancellation(p);
    normalize_profile(p);
    check_cancellation(p);
}

run_result run(parameters_t& p)
{
    check_cancellation(p);
    /*
    double epsx = 1e-11;
    finishes if |v|<=EpsX is fulfilled
    |.| means Euclidian norm
    v - scaled step vector
    v[i]=dx[i]/s[i]
    dx - step vector
    dx=X(k+1)-X(k)
    s - scaling coefficients set by MinLMSetScale()
    Recommended values: 1E-9 ... 1E-12.
    */
    run_result result;
    const auto initial_values = p.initial_values_alglib;
    double DiffStep = 0.0001;
    alglib::real_1d_array control_parameters;
    control_parameters.setcontent(p.solvables.size(), p.initial_values_alglib.data());
    alglib::real_1d_array s;
    s.setcontent(p.solvables.size(), p.scale.data());
    alglib::real_1d_array bndl;
    bndl.setcontent(p.solvables.size(), p.low_bound.data());
    alglib::real_1d_array bndu;
    bndu.setcontent(p.solvables.size(), p.up_bound.data());
    alglib::minlmstate state;
    alglib::minlmreport rep;
    alglib::minlmcreatev(p.solvables.size(), p.total_window_size, control_parameters, DiffStep, state);
    alglib::minlmsetbc(state, bndl, bndu);
    alglib::minlmsetcond(state, p.convergence_epsx, p.max_iterations);
    alglib::minlmsetscale(state, s);
    alglib::minlmsetnonmonotonicsteps(state, 2);
    if (p.run_solver) {
        struct progress_context {
            parameters_t& parameters;
            run_result& result;
            const std::vector<double>& initial_values;
        } context{p, result, initial_values};
        const auto evaluate = [](const alglib::real_1d_array& x, alglib::real_1d_array& fi, void* ptr) {
            alglib_solver(x, fi, &static_cast<progress_context*>(ptr)->parameters);
        };
        const auto report = [](const alglib::real_1d_array& x, double objective, void* ptr) {
            auto& context = *static_cast<progress_context*>(ptr);
            auto& p = context.parameters;
            check_cancellation(p);
            context.result.sum_squared_residuals = objective;
            progress_event event{event_kind::optimizer_progress, operation::solve, {}};
            event.evaluations = p.iterations;
            event.sum_squared_residuals = objective;
            for (std::size_t i = 0; i < p.solvables.size(); ++i) {
                const auto& parameter = p.solvables[i];
                event.parameters.push_back({parameter.source_name, parameter.name, x[i], context.initial_values[i]});
            }
            publish_event(p, std::move(event));
        };
        // Reports include the initial point and internal optimizer steps, not
        // finite-difference trial evaluations. The objective is sum(fi[i]^2).
        alglib::minlmsetxrep(state, true);
        alglib::minlmoptimize(state, +evaluate, +report, &context);
        alglib::minlmresults(state, control_parameters, rep);
        result.optimizer_ran = true;
        result.optimizer_iterations = rep.iterationscount;
        result.termination_type = rep.terminationtype;

        for (int index = 0; index < p.solvables.size(); index++) {
            p.initial_values_alglib[index] = control_parameters[index];
        }

    } else {
        alglib::real_1d_array residuals;
        alglib_solver(control_parameters, residuals, &p);
        double objective = 0;
        for (alglib::ae_int_t i = 0; i < residuals.length(); ++i)
            objective += residuals[i] * residuals[i];
        result.sum_squared_residuals = objective;
    }
    check_cancellation(p);
    result.residual_evaluations = p.iterations;
    for (std::size_t i = 0; i < p.solvables.size(); ++i) {
        const auto& parameter = p.solvables[i];
        result.parameters.push_back({parameter.source_name, parameter.name, p.initial_values_alglib[i], initial_values[i]});
    }
    return result;
}

std::filesystem::path export_results(parameters_t& p, const std::filesystem::path& output_directory)
{
    if (output_directory.empty()) {
        throw std::invalid_argument("Output directory must not be empty");
    }
    p.output_file_name.clear();
    for (int run = 0; run < p.experiment_runs.size(); run++) {
        experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
        p.output_file_name.append(run_ptr->name + ",");
    }
    // Experiment names must not redirect the export outside the chosen directory.
    if (p.output_file_name.find_first_of("/\\:") != std::string::npos) {
        throw std::invalid_argument("Experiment names must not contain path separators or colons");
    }
    std::filesystem::create_directories(output_directory);
    const auto output_path = output_directory / (p.output_file_name + ".txt");
    try {
        save_excel_output(p, output_path);
    } catch (const std::ios_base::failure& error) {
        throw workflow_error(p.active_operation, error_code::io,
                             "Cannot export results to " + output_path.string() + ": " + error.what());
    }
    return output_path;
}

void save_model_profiles(parameters_t& p, sqlite3* db)
{
    if (p.save_model_profiles) {
        write_model_profile_to_db(p, db);
    }
}

void save_fitted_parameters(parameters_t& p, sqlite3* db)
{
    write_alglib_values_to_db(p, db);
}

} // namespace tsensor_workflow
