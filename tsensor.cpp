#include <tsensor.h>
parameters_t p;

int
main()
{
    std::cout << "Openning Database navier.db: ";
    sqlite3 *db; // Database connection handle
    int rc = sqlite3_open("../../navier.db", &db);
    char* zErrMsg = 0; // For error messages
    if (rc != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
    } else {
        fprintf(stdout, "Database opened successfully\n");
    }

    read_model_parameters_from_db(db);
    lines_from_profile_text(db);
    read_exp_parameters_from_db(db);
    get_solve_settings_ID_from_db(db);
    read_specie_and_reaction_values_from_db(db);
    read_raw_profiles_from_db(db);
    read_inlet_cond_from_db(db);
    get_SOLUTION_IDs_from_db(db);
    read_alglib_values_from_db(db);
    normalize_profile();
    
    try
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
        double DiffStep = 0.000001;
        alglib::real_1d_array control_parameters;
        control_parameters.setcontent(p.number_of_variables, p.initial_values_alglib);
        alglib::real_1d_array s;
        s.setcontent(p.number_of_variables, p.scale);
        alglib::real_1d_array bndl;
        bndl.setcontent(p.number_of_variables, p.low_bound);
        alglib::real_1d_array bndu;
        bndu.setcontent(p.number_of_variables, p.up_bound);
        alglib::minlmstate state;
        alglib::minlmreport rep;
        alglib::minlmcreatev(p.number_of_variables, p.total_window_size, control_parameters, DiffStep, state);
        alglib::minlmsetbc(state, bndl, bndu);
        alglib::minlmsetcond(state, p.convergence_epsx, p.max_iterations);
        alglib::minlmsetscale(state, s);
        alglib::minlmsetnonmonotonicsteps(state, 2);
        if (p.run_solver) {
            std::cout << "minlmoptimize: iteration ";
            alglib::minlmoptimize(state, alglib_solver);   // Optimize
            std::cout << p.iterations - 1 << ". Done" << std::endl;
            alglib::minlmresults(state, control_parameters, rep);
            for (int index = 0; index < p.number_of_variables; index++) {
                p.initial_values_alglib[index] = control_parameters[index];
            }
            std::cout << "global(";
            for (int index = 0; index < p.global_solve_for.size(); index++) {
                std::cout << p.solve_for[index] << ":" << control_parameters[index];
                if (index != p.global_solve_for.size() - 1) {
                    std::cout << " ";
                }
            }
            std::cout << ")";
            int displacement = p.global_solve_for.size();
            for (int run = 0; run < p.experiment_runs.size(); run++) {
                experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
                for (int index = 0; index < p.global_solve_for.size(); index++) {
                    variable_location(p.solve_for[index], run_ptr) = control_parameters[index];
                }
                std::cout << " " + run_ptr->name + "(";
                for (int i = 0; i < run_ptr->solve_for.size(); i++) {
                    std::cout << run_ptr->solve_for.at(i) << ":" << control_parameters[displacement + i];
                    variable_location(run_ptr->solve_for.at(i), run_ptr) = control_parameters[displacement + i];
                    if (i != run_ptr->solve_for.size() - 1) {
                        std::cout << " ";
                    }
                }
                displacement += run_ptr->solve_for.size();
                std::cout << ")";
            }
            std::cout << std::endl;
        };

        p.output_file_name.clear();
        for (int run = 0; run < p.experiment_runs.size(); run++) {
            experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
            p.output_file_name.append(run_ptr->name + ",");
        }
        save_excel_output("../" + p.output_file_name + ".txt");
        //save_csv_out("../csv_out/" + p.output_file_name + ".csv");

        if (p.save_model_profiles) {
            write_model_profile_to_db(db);
        }
        
        std::string accept_alglib_values;
        std::cout << "Save parameter solutions to db as initial values? (y/n) ";
        std::cin >> accept_alglib_values;
        if (accept_alglib_values == "y") {
            write_alglib_values_to_db(db);
        }

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
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        return 1;
    }
    
    
    sqlite3_close(db);
    std::cout << "Database closed." << std::endl;
    return 0;
}