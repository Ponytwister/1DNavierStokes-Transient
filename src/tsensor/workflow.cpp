#include <workflow.h>

namespace tsensor_workflow {

run_session::run_session(const std::filesystem::path& database_path)
{
    if (database_path.empty()) {
        throw std::invalid_argument("Database path must not be empty");
    }
    const auto utf8 = database_path.u8string();
    sqlite3* raw = nullptr;
    const int rc = sqlite3_open_v2(reinterpret_cast<const char*>(utf8.c_str()),
                                  &raw, SQLITE_OPEN_READWRITE, nullptr);
    database_.reset(raw); // Also owns the handle returned on an open failure.
    if (rc != SQLITE_OK) {
        throw std::runtime_error("Cannot open database " + database_path.string()
                                 + ": " + sqlite3_errmsg(raw));
    }
}

void run_session::require_state(session_state expected) const
{
    if (state_ != expected) {
        throw std::logic_error("Invalid session operation order; use a fresh session for a new run");
    }
}

void run_session::load_inputs()
{
    require_state(session_state::empty);
    state_ = session_state::failed;
    tsensor_workflow::load_inputs(parameters_, database_.get());
    state_ = session_state::loaded;
}

void run_session::run()
{
    require_state(session_state::loaded);
    state_ = session_state::failed;
    tsensor_workflow::run(parameters_);
    state_ = session_state::completed;
}

std::filesystem::path run_session::export_results(const std::filesystem::path& directory)
{
    require_state(session_state::completed);
    return tsensor_workflow::export_results(parameters_, directory);
}

void run_session::save_model_profiles()
{
    require_state(session_state::completed);
    tsensor_workflow::save_model_profiles(parameters_, database_.get());
}

void run_session::save_fitted_parameters()
{
    require_state(session_state::completed);
    tsensor_workflow::save_fitted_parameters(parameters_, database_.get());
}

void load_inputs(parameters_t& p, sqlite3* db)
{
    read_model_parameters_from_db(p, db);
    lines_from_profile_text(p, db);
    read_exp_parameters_from_db(p, db);
    get_solve_settings_ID_from_db(p, db);
    read_specie_and_reaction_values_from_db(p, db);
    read_raw_profiles_from_db(p, db);
    read_inlet_cond_from_db(p, db);
    get_SOLUTION_IDs_from_db(p, db);
    //get_solvable_initial_values_from_db(p, db);
    read_alglib_values_from_db(p, db);
    normalize_profile(p);
}

void run(parameters_t& p)
{
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
        std::cout << "minlmoptimize: iteration ";
        alglib::minlmoptimize(state, alglib_solver, nullptr, &p);   // Optimize
        std::cout << p.iterations - 1 << ". Done" << std::endl;
        alglib::minlmresults(state, control_parameters, rep);

        for (int index = 0; index < p.solvables.size(); index++) {
            p.initial_values_alglib[index] = control_parameters[index];
        }

        for (int i = 0; i < p.solvables.size(); i++) {
            auto& s = p.solvables.at(i);
            std::cout << "(" << s.source_name << ":" << s.name << ":" << s.value() << ")";
            if (i != p.solvables.size() - 1) {
                std::cout << " ";
            }
        }
        std::cout << std::endl;
    } else {
        alglib::real_1d_array residuals;
        alglib_solver(control_parameters, residuals, &p);
    };
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
        throw std::runtime_error("Cannot export results to " + output_path.string()
                                 + ": " + error.what());
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
