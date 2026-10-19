#include <workflow.h>
#include <type_traits>
#include <sqlite3.h>
#include <sstream>
#include <map>

namespace tsensor_workflow {

namespace {
void record_run(sqlite3* db, parameters_t& p, const control_values& supplied,
                const run_result& result, const char* status)
{
    auto exec = [db](const char* sql) {
        char* error = nullptr;
        const int rc = sqlite3_exec(db, sql, nullptr, nullptr, &error);
        std::string message = error ? error : sqlite3_errmsg(db);
        sqlite3_free(error);
        if (rc != SQLITE_OK) throw std::runtime_error("Cannot save run history: " + message);
    };
    exec("CREATE TABLE IF NOT EXISTS run_history (ID INTEGER PRIMARY KEY, COMPLETED_AT_UTC TEXT NOT NULL, EXPERIMENTS TEXT, STATUS TEXT NOT NULL, RESIDUAL_EVALUATIONS INTEGER, ITERATIONS INTEGER, RESIDUAL_SUM REAL);");
    exec("CREATE TABLE IF NOT EXISTS run_history_values (RUN_ID INTEGER NOT NULL REFERENCES run_history(ID) ON DELETE CASCADE, CATEGORY TEXT NOT NULL, ITEM TEXT NOT NULL, FIELD TEXT NOT NULL, VALUE TEXT, PRIMARY KEY(RUN_ID,CATEGORY,ITEM,FIELD));");
    sqlite3_stmt* versionRaw = nullptr;
    if (sqlite3_prepare_v2(db, "PRAGMA user_version", -1, &versionRaw, nullptr) == SQLITE_OK) {
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> version(versionRaw, sqlite3_finalize);
        if (sqlite3_step(versionRaw) == SQLITE_ROW && sqlite3_column_int(versionRaw, 0) < 6)
            exec("PRAGMA user_version = 6");
    }
    exec("BEGIN IMMEDIATE");
    try {
        sqlite3_stmt* raw = nullptr;
        const char* insert = "INSERT INTO run_history(COMPLETED_AT_UTC,EXPERIMENTS,STATUS,RESIDUAL_EVALUATIONS,ITERATIONS,RESIDUAL_SUM) VALUES(datetime('now'),?1,?2,?3,?4,?5)";
        if (sqlite3_prepare_v2(db, insert, -1, &raw, nullptr) != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(raw, sqlite3_finalize);
        std::string experiments;
        for (const auto& run : p.experiment_runs) { if (!experiments.empty()) experiments += ", "; experiments += run.name; }
        sqlite3_bind_text(raw,1,experiments.c_str(),-1,SQLITE_TRANSIENT);
        sqlite3_bind_text(raw,2,status,-1,SQLITE_STATIC);
        sqlite3_bind_int(raw,3,result.residual_evaluations);
        if (result.optimizer_iterations) sqlite3_bind_int64(raw,4,*result.optimizer_iterations); else sqlite3_bind_null(raw,4);
        if (result.sum_squared_residuals) sqlite3_bind_double(raw,5,*result.sum_squared_residuals); else sqlite3_bind_null(raw,5);
        if (sqlite3_step(raw) != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
        const auto id = sqlite3_last_insert_rowid(db);
        sqlite3_finalize(statement.release());
        raw = nullptr;
        if (sqlite3_prepare_v2(db,"INSERT INTO run_history_values VALUES(?1,?2,?3,?4,?5)",-1,&raw,nullptr)!=SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> value(raw, sqlite3_finalize);
        auto add = [&](const char* category, const std::string& item, const char* field, const std::optional<std::string>& v) {
            sqlite3_reset(raw); sqlite3_clear_bindings(raw); sqlite3_bind_int64(raw,1,id);
            sqlite3_bind_text(raw,2,category,-1,SQLITE_STATIC); sqlite3_bind_text(raw,3,item.c_str(),-1,SQLITE_TRANSIENT);
            sqlite3_bind_text(raw,4,field,-1,SQLITE_STATIC);
            if (v) sqlite3_bind_text(raw,5,v->c_str(),-1,SQLITE_TRANSIENT); else sqlite3_bind_null(raw,5);
            if (sqlite3_step(raw)!=SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
        };
        auto number = [](double v) { std::ostringstream out; out.precision(17); out << v; return out.str(); };
        std::map<std::string, std::optional<std::string>> controls;
        sqlite3_stmt* controlRaw = nullptr;
        if (sqlite3_prepare_v2(db, "SELECT * FROM model_controls", -1, &controlRaw, nullptr) == SQLITE_OK) {
            std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> controlStatement(controlRaw, sqlite3_finalize);
            if (sqlite3_column_count(controlRaw) >= 2) {
                while (sqlite3_step(controlRaw) == SQLITE_ROW) {
                    const auto* name = reinterpret_cast<const char*>(sqlite3_column_text(controlRaw, 0));
                    const auto* value = reinterpret_cast<const char*>(sqlite3_column_text(controlRaw, 1));
                    if (name) controls[name] = sqlite3_column_type(controlRaw, 1) == SQLITE_NULL
                        ? std::nullopt : std::optional<std::string>(value ? value : "");
                }
            }
        }
        for (const auto& [name, value] : supplied) controls[name] = value;
        for (const auto& [name, val] : controls) add("control",name,"value",val);
        for (auto& run : p.experiment_runs) {
            for (auto& s : run.species) {
                add("species",run.name + "/" + s.name,"QE",number(s.QE.value()));
                add("species",run.name + "/" + s.name,"DIFFUSION_RATE",number(s.diffusion_rate));
                add("species",run.name + "/" + s.name,"PARTICLE_DIAMETER",number(s.diameter));
                add("species",run.name + "/" + s.name,"PARTICLE_DENSITY",number(s.particle_density));
                add("species",run.name + "/" + s.name,"MOLECULAR_WEIGHT",number(s.molecular_weight));
                add("species",run.name + "/" + s.name,"TYPE",std::to_string(s.type));
                add("species",run.name + "/" + s.name,"QE_ALIAS",s.qe_alias);
                add("species",run.name + "/" + s.name,"DIFFUSION_ALIAS",s.diffusion_alias);
                add("species",run.name + "/" + s.name,"DIAMETER_ALIAS",s.diameter_alias);
                add("species",run.name + "/" + s.name,"INPUT_UNITS",s.input_units);
                add("species",run.name + "/" + s.name,"MODEL_UNITS",s.model_units);
            }
            for (auto& reaction : run.reactions) {
                add("reaction",run.name + "/" + reaction.name,"K1",number(reaction.k[0].value()));
                add("reaction",run.name + "/" + reaction.name,"K2",number(reaction.k[1].value()));
                add("reaction",run.name + "/" + reaction.name,"K1_ALIAS",reaction.k_alias[0]);
                add("reaction",run.name + "/" + reaction.name,"K2_ALIAS",reaction.k_alias[1]);
                for (auto* s : reaction.specie_vect) {
                    const auto item = run.name + "/" + reaction.name + "/" + s->name;
                    if (reaction.coef.contains(s)) {
                        add("reaction",item,"COEFFICIENT",number(reaction.coef.at(s).value()));
                        const auto alias = reaction.coef_alias.find(s);
                        if (alias != reaction.coef_alias.end()) add("reaction",item,"COEFFICIENT_ALIAS",alias->second.name);
                    }
                    if (reaction.exp.contains(s)) {
                        add("reaction",item,"EXPONENT",number(reaction.exp.at(s).value()));
                        const auto alias = reaction.exp_alias.find(s);
                        if (alias != reaction.exp_alias.end()) add("reaction",item,"EXPONENT_ALIAS",alias->second.name);
                    }
                }
            }
        }
        exec("COMMIT");
    } catch (...) { sqlite3_exec(db,"ROLLBACK",nullptr,nullptr,nullptr); throw; }
}

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
                         std::stop_token cancellation, std::stop_token pause_request,
                         bool capture_live_profiles,
                         std::optional<single_profile_selection> single_profile)
{
    parameters_.cancellation = cancellation;
    parameters_.pause_request = pause_request;
    parameters_.progress = std::move(progress);
    parameters_.capture_live_profiles = capture_live_profiles;
    if (single_profile) {
        parameters_.single_profile_experiment = std::move(single_profile->experiment_name);
        parameters_.single_profile_dataset = std::move(single_profile->dataset_name);
        parameters_.single_profile_variables = std::move(single_profile->variables);
        parameters_.single_profile_initial_values = std::move(single_profile->initial_values);
    }
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
        loaded_controls_ = controls.value_or(control_values{});
        tsensor_workflow::load_inputs(parameters_, database_.get(), controls);
        state_ = session_state::loaded;
    });
}

run_result run_session::run()
{
    return perform(parameters_, operation::solve, [&] {
        require_state(session_state::loaded, operation::solve);
        state_ = session_state::failed;
        run_result result;
        try {
            result = tsensor_workflow::run(parameters_);
        } catch (...) {
            const auto failure = std::current_exception();
            try {
                run_result unsuccessful;
                record_run(database_.get(), parameters_, loaded_controls_, unsuccessful,
                    parameters_.cancellation.stop_requested() ? "cancelled" : "failed");
            } catch (...) {}
            std::rethrow_exception(failure);
        }
        record_run(database_.get(), parameters_, loaded_controls_, result, result.paused ? "paused" : "completed");
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

std::filesystem::path run_session::export_results_as(const std::filesystem::path& output_file)
{
    return perform(parameters_, operation::export_results, [&] {
        require_state(session_state::completed, operation::export_results);
        return tsensor_workflow::export_results_as(parameters_, output_file);
    });
}

void run_session::save_model_profiles()
{
    perform(parameters_, operation::save_model_profiles, [&] {
        require_state(session_state::completed, operation::save_model_profiles);
        tsensor_workflow::save_model_profiles(parameters_, database_.get());
    });
}

std::string run_session::generate_report()
{
    return perform(parameters_, operation::export_results, [&] {
        require_state(session_state::completed, operation::export_results);
        return generate_excel_report(parameters_);
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
    prepare_parameter_scope_selection(p, db);
    initialize_channel_dimensions(p, db);
    lines_from_profile_text(p, db);
    if (p.single_profile_experiment && p.experiments.size() != 1)
        throw std::runtime_error("The selected experiment dataset is no longer available.");
    read_exp_parameters_from_db(p, db);
    read_specie_and_reaction_values_from_db(p, db);
    read_raw_profiles_from_db(p, db);
    initialize_profile_entrance_concentrations(p, db);
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
    if (!p.solvables.empty())
        control_parameters.setcontent(p.solvables.size(), p.initial_values_alglib.data());
    if (p.run_solver && !p.solvables.empty()) {
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
        struct progress_context {
            parameters_t& parameters;
            run_result& result;
            const std::vector<double>& initial_values;
            alglib::minlmstate* state;
        } context{p, result, initial_values, &state};
        const auto evaluate = [](const alglib::real_1d_array& x, alglib::real_1d_array& fi, void* ptr) {
            alglib_solver(x, fi, &static_cast<progress_context*>(ptr)->parameters);
        };
        const auto report = [](const alglib::real_1d_array& x, double objective, void* ptr) {
            auto& context = *static_cast<progress_context*>(ptr);
            auto& p = context.parameters;
            check_cancellation(p);
            if (!context.result.original_sum_squared_residuals)
                context.result.original_sum_squared_residuals = objective;
            context.result.sum_squared_residuals = objective;
            progress_event event{event_kind::optimizer_progress, operation::solve, {}};
            event.evaluations = p.iterations;
            event.sum_squared_residuals = objective;
            for (std::size_t i = 0; i < p.solvables.size(); ++i) {
                const auto& parameter = p.solvables[i];
                event.parameters.push_back({parameter.source_name, parameter.name, x[i], context.initial_values[i]});
            }
            publish_event(p, std::move(event));
            // ALGLIB permits a termination request from another thread. Apply it
            // at a reporting point so the returned vector and profile are a
            // complete optimizer iterate, not a finite-difference trial.
            if (p.pause_request.stop_requested()) alglib::minlmrequesttermination(*context.state);
        };
        // Reports include the initial point and internal optimizer steps, not
        // finite-difference trial evaluations. The objective is sum(fi[i]^2).
        alglib::minlmsetxrep(state, true);
        alglib::minlmoptimize(state, +evaluate, +report, &context);
        alglib::minlmresults(state, control_parameters, rep);
        result.optimizer_ran = true;
        result.optimizer_iterations = rep.iterationscount;
        result.termination_type = rep.terminationtype;
        result.paused = p.pause_request.stop_requested() && rep.terminationtype == 8;

        for (int index = 0; index < p.solvables.size(); index++) {
            p.initial_values_alglib[index] = control_parameters[index];
        }

    } else {
        alglib::real_1d_array residuals;
        residuals.setlength(p.total_window_size);
        alglib_solver(control_parameters, residuals, &p);
        double objective = 0;
        for (alglib::ae_int_t i = 0; i < residuals.length(); ++i)
            objective += residuals[i] * residuals[i];
        result.original_sum_squared_residuals = objective;
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
    return export_results_as(p, output_path);
}

std::filesystem::path export_results_as(parameters_t& p, const std::filesystem::path& output_file)
{
    if (output_file.empty() || output_file.filename().empty()) {
        throw std::invalid_argument("Output file must not be empty");
    }
    const auto output_directory = output_file.parent_path();
    if (!output_directory.empty())
        std::filesystem::create_directories(output_directory);
    try {
        save_excel_output(p, output_file);
    } catch (const std::ios_base::failure& error) {
        throw workflow_error(p.active_operation, error_code::io,
                             "Cannot export results to " + output_file.string() + ": " + error.what());
    }
    return output_file;
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
