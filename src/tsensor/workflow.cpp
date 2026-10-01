#include <workflow.h>

namespace tsensor_workflow {

void load_inputs(sqlite3* db)
{
    read_model_parameters_from_db(db);
    lines_from_profile_text(db);
    read_exp_parameters_from_db(db);
    get_solve_settings_ID_from_db(db);
    read_specie_and_reaction_values_from_db(db);
    read_raw_profiles_from_db(db);
    read_inlet_cond_from_db(db);
    get_SOLUTION_IDs_from_db(db);
    //get_solvable_initial_values_from_db(db);
    read_alglib_values_from_db(db);
    normalize_profile();
}

void run()
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
    control_parameters.setcontent(p.solvables.size(), p.initial_values_alglib);
    alglib::real_1d_array s;
    s.setcontent(p.solvables.size(), p.scale);
    alglib::real_1d_array bndl;
    bndl.setcontent(p.solvables.size(), p.low_bound);
    alglib::real_1d_array bndu;
    bndu.setcontent(p.solvables.size(), p.up_bound);
    alglib::minlmstate state;
    alglib::minlmreport rep;
    alglib::minlmcreatev(p.solvables.size(), p.total_window_size, control_parameters, DiffStep, state);
    alglib::minlmsetbc(state, bndl, bndu);
    alglib::minlmsetcond(state, p.convergence_epsx, p.max_iterations);
    alglib::minlmsetscale(state, s);
    alglib::minlmsetnonmonotonicsteps(state, 2);
    if (p.run_solver) {
        std::cout << "minlmoptimize: iteration ";
        alglib::minlmoptimize(state, alglib_solver);   // Optimize
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
        void *ptr = nullptr;
        alglib_solver(control_parameters, residuals, ptr);
    };
}

std::filesystem::path export_results(const std::filesystem::path& output_directory)
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
        save_excel_output(output_path);
    } catch (const std::ios_base::failure& error) {
        throw std::runtime_error("Cannot export results to " + output_path.string()
                                 + ": " + error.what());
    }
    return output_path;
}

void save_model_profiles(sqlite3* db)
{
    if (p.save_model_profiles) {
        write_model_profile_to_db(db);
    }
}

void save_fitted_parameters(sqlite3* db)
{
    write_alglib_values_to_db(db);
}

void release_run_resources()
{
    // cleanup for potential loop
    p.output_file_name.clear();
    p.scatter_correction_type.clear();

    p.experiment_runs.clear();
    p.experiment_runs.shrink_to_fit();
    p.solve_for.clear();
    p.solve_for.shrink_to_fit();
    p.global_solve_for.clear();
    p.global_solve_for.shrink_to_fit();
    p.experiments.clear();
    p.experiments.shrink_to_fit();

    delete[] p.initial_values_alglib;
    delete[] p.low_bound;
    delete[] p.up_bound;
    delete[] p.scale;
}

} // namespace tsensor_workflow
