#include <tsensor.h>
#include <sstream>
#include <bits/stdc++.h>
#include <algorithm>
#include <exception>

namespace {
// Own SQLite's allocated error text, including when a caller throws.
class sqlite_error_message {
    char* message_ = nullptr;
public:
    sqlite_error_message() = default;
    sqlite_error_message(const sqlite_error_message&) = delete;
    sqlite_error_message& operator=(const sqlite_error_message&) = delete;
    ~sqlite_error_message() { sqlite3_free(message_); }
    char** out() {
        sqlite3_free(message_);
        message_ = nullptr;
        return &message_;
    }
};

using row_callback = int (*)(void*, int, char**, char**);

std::string sql_diagnostic(const char* sql)
{
    constexpr std::size_t maximum_length = 2048;
    const std::string statement = sql ? sql : "";
    if (statement.size() <= maximum_length) return statement;
    return statement.substr(0, maximum_length) + "\n... [" +
        std::to_string(statement.size() - maximum_length) + " SQL characters omitted]";
}

struct value_alias_token { double value = 0.0d; std::string alias; double sign = 1.0d; bool has_value = false; };
value_alias_token parse_value_alias(const std::string& token, const std::string& description)
{
    value_alias_token result;
    std::string value = token;
    if (!value.empty() && value.front() == '-') { result.sign = -1.0d; value.erase(value.begin()); }
    const auto marker = value.find('#');
    if (marker != std::string::npos) {
        result.alias = value.substr(marker + 1);
        if (result.alias.empty() || result.alias.find_first_of(" \t\r\n") != std::string::npos)
            throw std::runtime_error(description + " alias must include a variable name with no spaces");
        value.resize(marker);
    }
    if (!value.empty()) {
        std::size_t consumed = 0;
        result.value = result.sign * std::stod(value, &consumed);
        if (consumed != value.size()) throw std::runtime_error(description + " must be a number optionally followed by #variable_name");
        result.has_value = true;
    }
    return result;
}
bool alias_is_selected(const parameters_t& p, const experiment_run_struct& run, const std::string& alias)
{
    return std::find(p.global_solve_for.begin(), p.global_solve_for.end(), alias) != p.global_solve_for.end() ||
           std::find(run.solve_for.begin(), run.solve_for.end(), alias) != run.solve_for.end();
}

std::string value_with_alias(double value, const std::string& alias, double sign = 1.0d)
{
    if (alias.empty()) return double_to_string(value);
    const std::string numeric = double_to_string(sign < 0.0d ? std::abs(value) : value);
    return (sign < 0.0d ? "-" : "") + numeric + "#" + alias;
}

void execute_sql(parameters_t& p, sqlite3* db, const char* sql,
                row_callback callback, char** error)
{
    check_cancellation(p);
    struct callback_context {
        parameters_t& parameters;
        row_callback callback;
        std::exception_ptr failure;
    } context{p, callback, {}};
    // Never unwind C++ exceptions through sqlite3_exec. Returning nonzero lets
    // SQLite finalize its statement before the original exception is rethrown.
    const auto invoke = [](void* data, int count, char** values, char** columns) noexcept -> int {
        auto& context = *static_cast<callback_context*>(data);
        try {
            check_cancellation(context.parameters);
            return context.callback(&context.parameters, count, values, columns);
        } catch (...) {
            context.failure = std::current_exception();
            return 1;
        }
    };
    const int rc = sqlite3_exec(db, sql, callback ? +invoke : nullptr, &context, error);
    if (context.failure) { std::rethrow_exception(context.failure); }
    if (rc != SQLITE_OK) {
        throw tsensor_workflow::workflow_error(p.active_operation,
            tsensor_workflow::error_code::database,
            std::string("Error in executing SQL: ") + (error && *error ? *error : sqlite3_errmsg(db))
                + "\n" + sql_diagnostic(sql), rc);
    }
}
} // namespace

static void write_excel_report(parameters_t& p, std::ostream& fout, bool extended)
{
    auto write_header = [&](double width) {
        fout << "res_time"                  << "  "; //1
        fout << "bind_ratio(p1)"            << "  "; //2
        fout << "forward_reaction_rate_1"   << "  "; //3
        fout << "equalibrium_constant_1"    << "  "; //4
        fout << "dye_conc."                 << "  "; //5
        fout << "bead_conc."                << "  "; //6
        fout << "bead_surface_area"         << "  "; //7
        fout << "D-A"                       << "  "; //8
        fout << "profile_type"              << "  "; //9
        const double model_size_ratio = width * 1.0e6 / p.X;
        for (int x = 0; x < p.X; x++) {
            double channel_position = (double)x * model_size_ratio;
            if (std::abs(channel_position) < 1e-307) {
                channel_position = 0;
            }
            fout << channel_position << "  ";
        }
        fout << '\n';


    };
    std::optional<double> header_width;
    if (!extended && p.row_count == 0) write_header(p.experiment_runs.empty() ? p.W : p.experiment_runs.front().W);

    double dye_bead_ratio = 2000.0d;
    double total_dye = 0.0d;
    experiment_struct* exp_ptr;
    experiment_run_struct* run_ptr;
    specie_struct* FITC_ptr;
    specie_struct* PS_beads_ptr;
    specie_struct* Bound_Dye_1_ptr;
    reaction_struct* FITC_Bead_1_ptr;
    
    std::vector<std::string> specie_name_vect = {"Free_Dye", "Bound_Dye", "Total_Dye", "Unbound_Beads_(wt%)", "Bound_Beads_(wt%)", "Total_Beads_(wt%)", "Experimental_Derivative", "Numeric_Derivative", "Experimental_Profile", "Numeric_Model_Profile", "Experimental_Difference", "Numeric_Difference", "Error"};
    
    for (int row = 0; row < p.row_count; row++) {
        check_cancellation(p);
        exp_ptr = &p.experiments.at(row);
        run_ptr = exp_ptr->run;
        if (!extended && (!header_width || *header_width != run_ptr->W)) {
            write_header(run_ptr->W);
            header_width = run_ptr->W;
        }
        ptrdiff_t FITC = run_ptr->FITC;
        FITC_ptr = &run_ptr->species.at(FITC);
        ptrdiff_t PS_beads = run_ptr->PS_beads;
        PS_beads_ptr = &run_ptr->species.at(PS_beads);
        ptrdiff_t Bound_Dye_1 = run_ptr->Bound_Dye_1;
        Bound_Dye_1_ptr = &run_ptr->species.at(Bound_Dye_1);
        ptrdiff_t FITC_Bead_1 = run_ptr->FITC_Bead_1;
        FITC_Bead_1_ptr = &run_ptr->reactions.at(FITC_Bead_1);
        auto profile_names = specie_name_vect;
        if (extended) {
            // The legacy Bound_Beads row contains analytical_zero, not beads.
            // Keep its sampled axis and give the true model-grid beads their own row.
            profile_names[4] = "Analytical_Zero_(" + FITC_ptr->model_units + ")";
            profile_names[3] = "Unbound_Beads_(" + PS_beads_ptr->input_units + ")";
            profile_names[5] = "Total_Beads_(" + PS_beads_ptr->input_units + ")";
            profile_names.push_back("Bound_Beads_(" + PS_beads_ptr->input_units + ")");
            for (const auto& species : run_ptr->species)
                profile_names.push_back("Species:" + species.name + "_(" + species.model_units + ")");
        }
        const auto& model = exp_ptr->model_profile;
        const auto& experiment = exp_ptr->experimental_profile;
        const auto& out = exp_ptr->species_out;
        if (!extended) {
        fout << "sec"                       << "  "; //1
        fout << "(" + FITC_ptr->model_units + "_PSbead)/(" + PS_beads_ptr->model_units << "_FITC)  "; //2
        fout << "forward_reaction_rate_1"   << "  "; //3
        fout << "keq_1"                     << "  "; //4
        fout << FITC_ptr->model_units       << "  "; //5
        fout << PS_beads_ptr->input_units   << "  "; //6
        fout << PS_beads_ptr->model_units   << "  "; //7
        fout << "deriv"                     << "  "; //8
        fout << "Channel_Width_(um)->"      << "  "; //9
        for (int x = 0; x < exp_ptr->window_size; x++) {
            const double position = exp_ptr->channel_position.at(x);
            fout << (std::abs(position) < 1e-307 ? 0 : position) << "  ";
        }
        fout << '\n';
        }
        double beads_sa = 0;
        for (auto specie_key : exp_ptr->entrances.at(0).CONC) {
            if (specie_key.first == PS_beads_ptr) {
                beads_sa = specie_key.second * unit_conversion(run_ptr, PS_beads, PS_beads_ptr->input_units, PS_beads_ptr->model_units);
                break;
            }
        }

        const double FITC_unit = unit_conversion(run_ptr, FITC, FITC_ptr->model_units, FITC_ptr->input_units);
        const double bead_unit = unit_conversion(run_ptr, PS_beads, PS_beads_ptr->model_units, PS_beads_ptr->input_units);
        const double bead_coefficient = std::abs(FITC_Bead_1_ptr->coef.at(FITC_ptr).value());
        const double p1 = variable_location("p1", run_ptr).value();
        const double kon1 = variable_location("kon1", run_ptr).value();
        const double keq1 = variable_location("keq1", run_ptr).value();

        auto write_profile_axis = [&](const std::vector<double>& axis, const std::string& axis_name) {
            fout << "res_time  bind_ratio(p1)  forward_reaction_rate_1  equalibrium_constant_1  dye_conc.  bead_conc.  bead_surface_area  D-A  profile_type  ";
            for (const auto position : axis) fout << (std::abs(position) < 1e-307 ? 0 : position) << "  ";
            fout << '\n';
            fout << "sec  (" + FITC_ptr->model_units + "_PSbead)/(" + PS_beads_ptr->model_units
                 << "_FITC)  forward_reaction_rate_1  keq_1  " << FITC_ptr->model_units << "  "
                 << PS_beads_ptr->input_units << "  " << PS_beads_ptr->model_units
                 << "  deriv  " << axis_name << "  ";
            for (const auto position : axis) fout << (std::abs(position) < 1e-307 ? 0 : position) << "  ";
            fout << '\n';
        };
        auto sampled_axis = [&]() {
            std::vector<double> axis;
            axis.reserve(exp_ptr->channel_position.size());
            for (const auto position : exp_ptr->channel_position) axis.push_back(position);
            return axis;
        };
        auto model_axis = [&]() {
            std::vector<double> axis;
            axis.reserve(p.X);
            const double scale = run_ptr->W * 1.0e6 / p.X;
            for (int x = 0; x < p.X; ++x) axis.push_back(x * scale);
            return axis;
        };

        int width;
        if (extended) {
            // The Report tab uses these source and rescaled views; legacy
            // exports retain their historical profile set.
            const auto source_count = exp_ptr->loaded_experimental_profile.size();
            if (source_count > 0) {
                std::vector<double> source_axis;
                source_axis.reserve(source_count);
                for (std::size_t x = 0; x < source_count; ++x) source_axis.push_back(static_cast<double>(x));
                auto write_source_profile = [&](const std::string& type, const std::vector<double>& values) {
                    write_profile_axis(source_axis, "Raw_Sample_Index");
                    fout << p.Z * run_ptr->dt << "  " << p1 << "  " << kon1 << "  " << keq1 << "  "
                         << run_ptr->dye_conc << "  " << exp_ptr->second_name << "  " << beads_sa
                         << "  -  " << run_ptr->name << "_" << exp_ptr->second_name << "wt%_" << type << "  ";
                    for (std::size_t x = 0; x < source_count; ++x) fout << values.at(x) << "    ";
                    fout << '\n';
                };
                write_source_profile("Raw_Profile", exp_ptr->loaded_experimental_profile);
                write_source_profile("Normalized_Profile", exp_ptr->normalized_experimental_profile);
            }

            std::vector<double> rescaled_raw_profile;
            rescaled_raw_profile.reserve(exp_ptr->window_size);
            int last_raw_sample = 0;
            for (int x = 0; x < exp_ptr->window_size; ++x) {
                if (x < exp_ptr->width.value()) {
                    const auto padded_index = static_cast<int>(std::ceil(x + exp_ptr->left_edge.value()));
                    if (!exp_ptr->loaded_experimental_profile.empty()) {
                        last_raw_sample = std::clamp(padded_index - p.exp_left_padding, 0,
                            static_cast<int>(exp_ptr->loaded_experimental_profile.size()) - 1);
                    }
                }
                rescaled_raw_profile.push_back(exp_ptr->loaded_experimental_profile.empty() ? 0.0
                    : exp_ptr->loaded_experimental_profile.at(static_cast<std::size_t>(last_raw_sample)));
            }
            write_profile_axis(sampled_axis(), "Channel_Width_(um)->");
            fout << p.Z * run_ptr->dt << "  " << p1 << "  " << kon1 << "  " << keq1 << "  "
                 << run_ptr->dye_conc << "  " << exp_ptr->second_name << "  " << beads_sa
                 << "  -  " << run_ptr->name << "_" << exp_ptr->second_name
                 << "wt%_Rescaled_Raw_Profile  ";
            for (const auto value : rescaled_raw_profile) fout << value << "    ";
            fout << '\n';
        }
        for (int j = 0; j < profile_names.size(); j++) {
            if (extended) {
                const bool uses_sample_axis = (j == 4 || (j >= 6 && j <= 12));
                write_profile_axis(uses_sample_axis ? sampled_axis() : model_axis(),
                    uses_sample_axis ? "Channel_Width_(um)->" : "Model_Grid_(um)->");
            }
            fout << p.Z * run_ptr->dt                   << "  "; //1
            fout << p1                                       << "  "; //2
            fout << kon1                                     << "  "; //3
            fout << keq1                                     << "  "; //4
            fout << run_ptr->dye_conc                   << "  "; //5
            fout << exp_ptr->second_name                << "  "; //6
            fout << beads_sa                            << "  "; //7
            if (j < 6 || j > 11) {
                fout << "-"                             << "  "; //8
            } else if (j == 6) {
                fout << exp_ptr->exp_DA                 << "  "; //8
            } else if (j == 7) {
                fout << exp_ptr->model_DA               << "  "; //8
            } else if (j == 8) {
                if (exp_ptr->zero_row_ptr != NULL) {
                    fout << exp_ptr->exp_integral       << "  "; //8
                } else {
                    fout << "-"                         << "  "; //8
                }
            } else if (j == 9) {
                if (exp_ptr->zero_row_ptr != NULL) {
                    fout << exp_ptr->model_integral     << "  "; //8
                } else {
                    fout << "-"                         << "  "; //8
                }
            } else if (j == 10) {
                fout << exp_ptr->analytic_exp_integral  << "  "; //8
            } else if (j == 11) {
                fout << exp_ptr->analytic_model_integral<< "  "; //8
            }
            
            if (j == 8 || j == 9) {
                fout << run_ptr->name << "_" << exp_ptr->second_name << "wt%_"; //9
            }
            fout << profile_names.at(j) << "  "; //9
            
            if ((j < 6 && j != 4) || j >= 13) {
                width = p.X;
            } else {
                width = exp_ptr->window_size;
            }
            for (int x = 0; x < width; x++) {
                switch(j) {
                    case 0: // Free_Dye
                        assert (x < p.X);
                        fout << out[FITC][x] * FITC_unit;
                        break;
                    case 1: // Bound_Dye
                        fout << (out[Bound_Dye_1][x]) * FITC_unit;
                        break;
                    case 2: // Total_Dye
                        fout << (out[FITC][x] + out[Bound_Dye_1][x]) * FITC_unit;
                        break;
                    case 3: // Unbound_Beads_(wt%)
                        fout << out[PS_beads][x] * bead_unit;
                        break;
                    case 4: // Analytical zero; legacy files retain the historical bead label.
                        fout << exp_ptr->analytical_zero.at(x);
                        break;   
                    case 5: // Total_Beads_(wt%)
                        fout << (out[PS_beads][x] + out[Bound_Dye_1][x] / bead_coefficient) * bead_unit;
                        break;
                    case 6: // Experimental_Derivative
                        fout << exp_ptr->experimental_derivative.at(x);
                        break;
                    case 7: // Numeric_Derivative
                        fout << exp_ptr->numeric_derivative.at(x);
                        break;    
                    case 8: // Experimental_Profile
                        fout << experiment.at(x) * run_ptr->dye_conc_mgml;
                        break;
                    case 9: // Numeric_Model_Profile
                        fout << model.at(x) * FITC_unit;
                        break;
                    case 10: // Experimental_Difference
                        fout << exp_ptr->experimental_difference.at(x);
                        break;    
                    case 11: // Numeric_Difference
                        fout << exp_ptr->numeric_difference.at(x) * FITC_unit;
                        break;
                    case 12: // Error
                        fout << exp_ptr->error.at(x);
                        break;
                    case 13: { // Bound beads: original stoichiometric conversion.
                        if (!std::isfinite(bead_coefficient) || bead_coefficient == 0)
                            throw std::invalid_argument("Cannot report bound beads with a zero or nonfinite dye coefficient");
                        fout << out.at(Bound_Dye_1).at(x) / bead_coefficient * bead_unit;
                        break;
                    }
                    default: // Additional species profiles stay in declared model units.
                        fout << out.at(j - 14).at(x);
                        break;
                }
                fout << "    ";
            }
            fout << '\n';
        }
        fout << '\n';
    }
}

void save_excel_output(parameters_t& p, const std::filesystem::path& file_name)
{
    add_report(p, 3, "Exporting results to " + file_name.string());
    std::ofstream fout;
    fout.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    fout.open(file_name, std::ofstream::out | std::ofstream::trunc);
    write_excel_report(p, fout, false);
    fout.close();
}

std::string generate_excel_report(parameters_t& p)
{
    // Preserve round-trip double precision for subsequent user-selected formatting.
    // This is an in-memory snapshot: no report file or database writes occur.
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    write_excel_report(p, output, true);
    return output.str();
}

void 
delete_values_from_db(parameters_t& p, sqlite3* db, std::string table, std::string where_conditions) //reading data using callback functions
{
    int row = 0;
    sqlite_error_message errMsg;
    std::string sqltext;
    if (where_conditions.empty()) {
        return;
    }
    sqltext.append("DELETE FROM " + table + " WHERE " + where_conditions + ";");
    execute_sql(p, db, sqltext.c_str(), 0, errMsg.out());
}

int 
read_model_parameters_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    if (argv[1] == NULL) {
        return 0;
    }
    std::string criterion;
    std::string value = argv[1];
    removeSpaces(value);
    criterion = argv[0];

    ptrdiff_t run_index;
    experiment_run_struct* run_ptr;
    std::string s; // variable to store token obtained from the original string
    std::stringstream ss(value); // constructing stream from the string

    if (criterion == "width resolution (X)") {
        p.X = std::stoi(value);
    } else if (criterion == "length/time resolution (Z)") {
        p.Z = std::stoi(value);
    } else if (criterion == "experiment_name") {
        while (getline(ss, s, ' ')) { 
            removeSpaces(s);
            p.experiment_runs.emplace_back();
            p.experiment_runs.back().name = s; 
        }
    } else if (criterion == "exp_left_padding") {
        p.exp_left_padding = std::stoi(value);
    } else if (criterion == "exp_right_padding") {
        p.exp_right_padding = std::stoi(value);
    } else if (criterion == "disable_reactions") {
        std::transform(value.begin(), value.end(), value.begin(), ::tolower);
        if (value == "true") {
            p.disable_reactions = true;
        } else {p.disable_reactions = false;}
    } else if (criterion == "run_solver") {
        if (value == "true" || value == "false") {
            p.run_solver = value == "true";
        }
    } else if (criterion == "scatter_correction_type") {
        p.scatter_correction_type = value;
    } else if (criterion == "save_model_profiles") {
        if (value == "true") {
            p.save_model_profiles = true;
        } else {
            p.save_model_profiles = false;
        }
    } else if (criterion == "max_iterations") {
        p.max_iterations = std::stoi(value);
    } else if (criterion == "convergence_epsx") {
        p.convergence_epsx = std::stod(value);
    } else if (criterion == "disable_reverse_reactions") {
        std::transform(value.begin(), value.end(), value.begin(), ::tolower);
        if (value == "true") {
            p.disable_reverse_reactions = true;
        } else {p.disable_reverse_reactions = false;}
    } else if (criterion == "use_alglib_init_values") {
        p.use_alglib_init_values = false;
        std::transform(value.begin(), value.end(), value.begin(), ::tolower);
        if (value == "true") {
            p.use_alglib_init_values = true;
        } else {p.use_alglib_init_values = false;}
    } else if (criterion == "universal_solve_for") {
        while (getline(ss, s, ' ')) {
            if (p.single_profile_experiment &&
                std::find(p.single_profile_variables.begin(), p.single_profile_variables.end(), s) == p.single_profile_variables.end()) continue;
            if (std::find(p.suppressed_universal_parameters.begin(),
                          p.suppressed_universal_parameters.end(), s) !=
                p.suppressed_universal_parameters.end()) continue;
            if (std::find(p.solve_for.begin(), p.solve_for.end(), s) == p.solve_for.end()) { // to avoid repeats
                p.global_solve_for.push_back(s);
                p.solve_for.push_back(s);
                p.solvables.emplace_back(s, "global");
            }
        }
    } else if (criterion == "debug_level") {
        p.debug_level = std::stoi(value);
    }
    return 0;
}

void 
read_model_parameters_from_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Retrieving model control parameters:");
    sqlite_error_message errMsg;
    std::string sqltext = "SELECT * FROM 'model_controls';";
    const char* sql = sqltext.c_str();
    execute_sql(p, db, sql, read_model_parameters_db_callback, errMsg.out());
    add_finishing_report(p, 3, "Done");

    std::string runs =std::to_string(p.experiment_runs.size()) + " runs (";
    for (int run = 0; run < p.experiment_runs.size(); run++) {
        runs.append(p.experiment_runs.at(run).name);
        if (run != p.experiment_runs.size() - 1) {
            runs.append(" ");
        } else {
            runs.append(")");
        }
    }
    add_finishing_report(p, 3, runs);
}

int
raw_profile_row_count_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    std::string criterion;

    for(int i = 0; i < count; i++) {
        std::string value;
        if (argv[i] != NULL) {value = argv[i];} else {value.clear();}
        removeSpaces(value);
        criterion = columnNames[i];
        experiment_run_struct* run_ptr;
        if (criterion == "NAME") {
            assert(!value.empty());
            std::string run_names;
            for (int run = 0; run < p.experiment_runs.size(); run++) {
                run_names.append("'" + p.experiment_runs.at(run).name + "' ");
                if (p.experiment_runs.at(run).name == value) {
                    run_ptr = &p.experiment_runs.at(run);
                    break;
                } else if (run == p.experiment_runs.size() - 1) {
                    throw std::runtime_error("Could not find run matching NAME:'" + value + "' existing names:" + run_names);
                }
            }
        } else if (criterion == "WT_PERCENT") {
            assert(!value.empty());
            if (p.single_profile_experiment &&
                (run_ptr->name != *p.single_profile_experiment || value != *p.single_profile_dataset)) continue;
            p.experiments.emplace_back();
            p.experiments.back().second_name = value;
            p.experiments.back().run = run_ptr;
        }
    };
    p.row_count = p.experiments.size();
    return 0;
}

void
lines_from_profile_text(parameters_t& p, sqlite3* db)
{
    add_report(p, 3, "Retrieving profile row count");
    sqlite_error_message errMsg;
    std::string sqltext = "SELECT NAME, WT_PERCENT FROM 'raw_profile' WHERE ";
    int number_of_runs = p.experiment_runs.size();
    for (int run = 0; run < number_of_runs; run++) {
        experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
        sqltext.append("(NAME = '" + run_ptr->name + "')");
        if (run != number_of_runs - 1) {
            sqltext.append(" OR ");
        } else {
            sqltext.append(";");
        }
    }
    const char* sql = sqltext.c_str();
    execute_sql(p, db, sql, raw_profile_row_count_db_callback, errMsg.out());
    add_report(p, 3, "Loaded " + std::to_string(p.row_count) + " profile rows");
}

int 
exp_parameters_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    std::string criterion;
    std::string s; // variable to store token obtained from the original string
    ptrdiff_t run_index;
    experiment_run_struct* run_ptr;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (argv[i] == nullptr) continue;
        
        std::string col_val = argv[i];
        assert (!col_val.empty());
        std::stringstream ss(col_val); // constructing stream from the string
        if (criterion == "NAME") { 
            for (ptrdiff_t j = 0; j < p.experiment_runs.size(); j++) {
                if (p.experiment_runs.at(j).name == col_val) {
                    run_index = j;
                    run_ptr = &p.experiment_runs.at(j);
                }
            }
        } else if (criterion == "LOW_REF_LEFT") {
            run_ptr->low_ref_start = std::stoi(col_val);
        } else if (criterion == "LOW_REF_RIGHT") {
            run_ptr->low_ref_end = std::stoi(col_val);
        } else if (criterion == "HIGH_REF_LEFT") {
            run_ptr->high_ref_start = std::stoi(col_val);
        } else if (criterion == "HIGH_REF_RIGHT") {
            run_ptr->high_ref_end = std::stoi(col_val);
        } else if (criterion == "PARAMETERS_TO_SOLVE_FOR") {
            while (getline(ss, s, ' ')) {
                if (p.single_profile_experiment &&
                    std::find(p.single_profile_variables.begin(), p.single_profile_variables.end(), s) == p.single_profile_variables.end()) continue;
                if (std::find(p.suppressed_experiment_parameters[run_ptr->name].begin(),
                              p.suppressed_experiment_parameters[run_ptr->name].end(), s) !=
                        p.suppressed_experiment_parameters[run_ptr->name].end()) continue;
                if (std::find(p.solve_for.begin(), p.solve_for.end(), s) == p.solve_for.end()) { // to avoid repeats in p.solve_for
                    p.solve_for.push_back(s);
                }
                if (std::find(run_ptr->solve_for.begin(), run_ptr->solve_for.end(), s) == run_ptr->solve_for.end()) { // to avoid repeats in run_ptr->solve_for
                    run_ptr->solve_for.push_back(s);
                    p.solvables.emplace_back(s, run_ptr->name);
                }
            }
        
        } else if (criterion == "DEFAULT_NORMALIZATION") {
            run_ptr->normalization_method = col_val;
        } else if (criterion == "SPECIES") {
            while (getline(ss, s, ' ')) { 
                run_ptr->species.emplace_back();
                run_ptr->species.back().name = s;
            }
            run_ptr->number_of_species = run_ptr->species.size();
        } else if (criterion == "REACTIONS") {
            while (getline(ss, s, ' ')) { 
                run_ptr->reactions.emplace_back();
                run_ptr->reactions.back().name = s;
            }
            run_ptr->number_of_reactions = run_ptr->reactions.size();
        } else if (criterion == "ENTRANCE_FLOWRATE") {
            std::vector <std::string> entrance_vector;
            while (getline(ss, s, ' ')) { 
                entrance_vector.emplace_back(s);
                run_ptr->ENTRANCE_FLOWRATE.emplace_back();
                run_ptr->ENTRANCE_FLOWRATE.back() = std::stod(s);
                run_ptr->total_flowrate += run_ptr->ENTRANCE_FLOWRATE.back();
            }

        } else if (criterion == "LEFT_EDGE") {
            std::vector <double> edges;
            while (getline(ss, s, ' ')) { 
                // Store the experiment value in the same unpadded coordinates as raw_profile.LEFT_EDGE.
                edges.emplace_back(std::max(std::stod(s) + p.exp_left_padding, 0.0d));
            }
            run_ptr->left_edge.value() = edges.at(0);
            run_ptr->left_edge.param_init = true;
        } else if (criterion == "SPECIE_INLET_CONC_UNITS") {
            ptrdiff_t specie_index = 0;
            while (getline(ss, s, ' ')) { 
                run_ptr->species.at(specie_index).input_units = s;
                specie_index++;
            }
        } else if (criterion == "SPECIE_MODEL_CONC_UNITS") {
            ptrdiff_t specie_index = 0;
            while (getline(ss, s, ' ')) { 
                run_ptr->species.at(specie_index).model_units = s;
                specie_index++;
            }
        } else if (criterion == "WIDTH") {
            run_ptr->width.value() = std::stod(col_val);
            run_ptr->width.param_init = true;
        }
        
    }
    return 0;
}

void 
read_exp_parameters_from_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Reading experimental parameters from database");
    for (int run = 0; run < p.experiment_runs.size(); run++) {
        experiment_run_struct* run_ptr = &p.experiment_runs.at(run);

        sqlite_error_message errMsg;
        std::string sqltext = "SELECT * FROM 'experiments' WHERE NAME='" + run_ptr->name + "';";
        const char* sql = sqltext.c_str();
        execute_sql(p, db, sql, exp_parameters_db_callback, errMsg.out());
        if (run_ptr->number_of_species == 0) {throw std::runtime_error(run_ptr->name + " number_of_species is zero");}

        double restime      = run_ptr->W * run_ptr->H * run_ptr->L / (run_ptr->total_flowrate); // seconds
        run_ptr->dt         = restime / p.Z; // seconds
    }

    add_report(p, 3, std::to_string(p.solvables.size()) + " solvables for "
                     + std::to_string(p.solve_for.size()) + " variables");
}

int 
specie_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    std::string criterion;
    std::string value;
    ptrdiff_t index;
    experiment_run_struct* run_ptr;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (argv[i] != NULL) {value = argv[i];} else {value.clear();}
        removeSpaces(value);
        specie_struct* specie_ptr;
        if (criterion == "run") {
            if (!value.empty()) {
                run_ptr = &p.experiment_runs.at(std::stoi(value));
            }
        } else if (criterion == "SPECIES_NAME") {
            if (!value.empty()) {
                for (ptrdiff_t j = 0; j < run_ptr->number_of_species; j++) {
                    if (run_ptr->species.at(j).name == value) {
                        index = j;
                        specie_ptr = &run_ptr->species.at(j);
                        break;
                    }
                }
            }
        } else if (criterion == "SPECIES_TYPE") {
            specie_ptr->type = 0;
            if (value == "molecule") {
                specie_ptr->type = 1;
            } else if (value == "particle") {
                specie_ptr->type = 2;
            }
        } else if (criterion == "DIFFUSION_RATE") {
            const auto parsed = value.empty() ? value_alias_token{} : parse_value_alias(value, "Species diffusion");
            specie_ptr->diffusion_alias = parsed.alias;
            specie_ptr->diffusion_rate = parsed.value;
            if (!parsed.alias.empty()) run_ptr->alias_variables.try_emplace(parsed.alias).first->second.param_init |= parsed.has_value;
        } else if (criterion == "QE") {
            const auto parsed = value.empty() ? value_alias_token{1.0d, {}, 1.0d} : parse_value_alias(value, "Species QE");
            specie_ptr->qe_alias = parsed.alias;
            if (!parsed.alias.empty()) {
                auto& alias_variable = run_ptr->alias_variables.try_emplace(parsed.alias).first->second;
                alias_variable.param_init |= parsed.has_value;
                specie_ptr->QE = alias_is_selected(p, *run_ptr, parsed.alias)
                    ? solvable(&alias_variable, parsed.alias) : solvable(parsed.value);
            } else specie_ptr->QE = solvable(parsed.value);
        } else if (criterion == "PARTICLE_DIAMETER") {
            const auto parsed = value.empty() ? value_alias_token{} : parse_value_alias(value, "Species diameter");
            specie_ptr->diameter_alias = parsed.alias;
            specie_ptr->diameter = parsed.value * 1.0e-9d;
            if (!parsed.alias.empty()) run_ptr->alias_variables.try_emplace(parsed.alias).first->second.param_init |= parsed.has_value;
        } else if (criterion == "PARTICLE_DENSITY") {
            if (!value.empty()) {
                specie_ptr->particle_density = std::stod(value);
            }
        } else if (criterion == "MOLECULAR_WEIGHT") {
            if (!value.empty()) {
                specie_ptr->molecular_weight = std::stod(value);
            }
        }
    };
    return 0;
}

int 
reaction_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    std::string criterion;
    std::string value;
    ptrdiff_t react_index;
    specie_struct* specie_ptr;
    experiment_run_struct* run_ptr;
    std::vector<std::ptrdiff_t> specie_vect;
    std::string s; // variable to store token obtained from the original string
    for(int i = 0; i < count; i++) {
        ptrdiff_t iter = 0;
        criterion = columnNames[i];
        if (argv[i] != NULL) {value = argv[i];} else {value.clear();}
        removeSpaces(value);
        std::stringstream ss(value); // constructing stream from the string
        if (criterion == "run") {
            if (!value.empty()) {
                run_ptr = &p.experiment_runs.at(std::stoi(value));
            }
        } else if (criterion == "REACTION_NAME") {
            if (!value.empty()) {
                for (ptrdiff_t j = 0; j < run_ptr->number_of_reactions; j++) {
                    if (run_ptr->reactions.at(j).name == value) {
                        react_index = j;
                    }
                }
            }
        } else if (criterion == "SPECIES") {
            while (getline(ss, s, ' ')) { 
                specie_vect.push_back(specie_index(run_ptr, s));
                for (ptrdiff_t j = 0; j < run_ptr->number_of_species; j++) {
                    if (run_ptr->species.at(j).name == s) {
                        specie_ptr = &run_ptr->species.at(j);
                    }
                }
                run_ptr->reactions.at(react_index).specie_vect.emplace_back(specie_ptr);
            }
        } else if (criterion == "COEFFICIENTS") {
            while (getline(ss, s, ' ')) {
                removeSpaces(s);
                auto& reaction = run_ptr->reactions.at(react_index);
                specie_ptr = reaction.specie_vect.at(iter);
                const auto parsed = parse_value_alias(s, "Reaction coefficient");
                if (!parsed.alias.empty()) {
                    run_ptr->alias_variables.try_emplace(parsed.alias).first->second.param_init |= parsed.has_value;
                    reaction.coef_alias[specie_ptr] = {parsed.alias, parsed.sign};
                } else reaction.coef_alias.erase(specie_ptr);
                auto& coefficient = reaction.coef.try_emplace(specie_ptr).first->second;
                coefficient.value() = parsed.value;
                coefficient.param_init = true;
                coefficient.source_name = run_ptr->name;
                iter++;
            }
        } else if (criterion == "Ks") {
            while (getline(ss, s, ' ')) {
                removeSpaces(s);
                auto& reaction = run_ptr->reactions.at(react_index);
                const auto parsed = parse_value_alias(s, "Reaction Ks");
                reaction.k_alias[iter] = parsed.alias;
                reaction.k[iter].value() = parsed.value;
                reaction.k[iter].param_init = true;
                reaction.k[iter].source_name = run_ptr->name;
                if (!parsed.alias.empty()) run_ptr->alias_variables.try_emplace(parsed.alias).first->second.param_init |= parsed.has_value;
                iter++;
            }
        } else if (criterion == "EXPONENTS") {
            while (getline(ss, s, ' ')) {
                removeSpaces(s);
                auto& reaction = run_ptr->reactions.at(react_index);
                specie_ptr = reaction.specie_vect.at(iter);
                const auto parsed = parse_value_alias(s, "Reaction exponent");
                if (!parsed.alias.empty()) {
                    run_ptr->alias_variables.try_emplace(parsed.alias).first->second.param_init |= parsed.has_value;
                    reaction.exp_alias[specie_ptr] = {parsed.alias, parsed.sign};
                } else reaction.exp_alias.erase(specie_ptr);
                auto& exponent = reaction.exp.try_emplace(specie_ptr).first->second;
                exponent.value() = parsed.value;
                exponent.param_init = true;
                exponent.source_name = run_ptr->name;
                iter++;
            }
        }
    };
    return 0;
}

static void
synchronize_run_aliases(parameters_t& p, experiment_run_struct& run)
{
    for (auto& reaction : run.reactions) {
        for (std::size_t parameter = 0; parameter < std::size(reaction.k_alias); ++parameter) {
            if (!reaction.k_alias[parameter].empty() && run.alias_variables.at(reaction.k_alias[parameter]).source) {
                reaction.k[parameter].value() = run.alias_variables.at(reaction.k_alias[parameter]).value();
            }
        }
        for (const auto& [species, alias] : reaction.coef_alias) {
            auto& variable = run.alias_variables.at(alias.name);
            if (variable.source) reaction.coef.at(species).value() = alias.sign * variable.value();
        }
        for (const auto& [species, alias] : reaction.exp_alias) {
            auto& variable = run.alias_variables.at(alias.name);
            if (variable.source) reaction.exp.at(species).value() = alias.sign * variable.value();
        }
    }
    for (auto& species : run.species) {
        if (!species.diameter_alias.empty() && run.alias_variables.at(species.diameter_alias).source) {
            species.diameter = run.alias_variables.at(species.diameter_alias).value() * 1.0e-9d;
        }
        if (species.type == 2 && (!std::isfinite(species.diameter) || species.diameter <= 0.0d ||
                                  !std::isfinite(species.particle_density) || species.particle_density <= 0.0d)) {
            throw std::runtime_error("Particle diameter or density undefined for species '" + species.name + "'");
        }
        if (!species.diffusion_alias.empty() && run.alias_variables.at(species.diffusion_alias).source) {
            species.diffusion_rate = run.alias_variables.at(species.diffusion_alias).value();
        } else if (species.type == 2) {
            species.diffusion_rate = 1.380649e-23d * run.temperature / (3.0d * M_PI * run.visc * species.diameter);
        }
        if (!std::isfinite(species.diffusion_rate) || species.diffusion_rate < 0.0d) {
            throw std::runtime_error("Species diffusion rate must be finite and non-negative for '" + species.name + "'");
        }
        species.r = species.diffusion_rate * run.dt / (run.W * run.W) * p.X * p.X;
    }
}

static solvable*
find_profile_alias_parameter(parameters_t& p, const experiment_struct& experiment, const std::string& alias)
{
    const std::string profile_source = experiment.run->name + "_" + experiment.second_name;
    for (const std::string& source : {profile_source, experiment.run->name, std::string("global")}) {
        for (auto& parameter : p.solvables) {
            if (parameter.name == alias && parameter.source_name == source) return &parameter;
        }
    }
    return nullptr;
}

static double
profile_alias_value(parameters_t& p, const experiment_struct& experiment, const std::string& alias,
                    double fallback, bool has_fallback, bool use_initial_values)
{
    auto* parameter = find_profile_alias_parameter(p, experiment, alias);
    if (!parameter) {
        if (has_fallback) return fallback;
        throw std::runtime_error("Raw profile ENTRANCE_CONC alias '" + alias + "' has no selected parameter or numeric fallback");
    }
    if (!use_initial_values) return parameter->value();
    for (std::size_t i = 0; i < p.solvables.size(); ++i) {
        if (&p.solvables[i] == parameter) return p.initial_values_alglib.at(i);
    }
    throw std::runtime_error("Raw profile ENTRANCE_CONC alias parameter was not found in the solver inputs");
}

static void
synchronize_profile_entrance_concentrations(parameters_t& p, bool use_initial_values)
{
    for (auto& experiment : p.experiments) {
        if (experiment.omit) continue;
        auto* run = experiment.run;
        for (const auto& override : experiment.entrance_concentrations) {
            if (override.alias.empty()) continue;
            if (override.entrance_number < 1 || override.entrance_number > static_cast<int>(experiment.entrances.size()))
                throw std::runtime_error("ENTRANCE_CONC refers to an inlet that is missing from ENTRANCE_FLOWRATE.");
            const auto species = specie_index(run, override.species_name);
            auto* inlet_species = &run->species.at(species);
            const double concentration = profile_alias_value(p, experiment, override.alias,
                override.concentration, override.has_paired_value, use_initial_values);
            if (!std::isfinite(concentration) || concentration < 0.0)
                throw std::runtime_error("ENTRANCE_CONC must be a finite non-negative number");
            const std::string source_units = override.units.empty() ? inlet_species->model_units : override.units;
            const double input_concentration = concentration *
                unit_conversion(run, species, source_units, inlet_species->input_units);
            if (!std::isfinite(input_concentration))
                throw std::runtime_error("ENTRANCE_CONC conversion produced a non-finite value.");
            experiment.entrances.at(override.entrance_number - 1).CONC[inlet_species] = input_concentration;
        }
    }

    for (auto& run : p.experiment_runs) {
        run.dye_conc_mgml = 0.0;
        for (const auto& experiment : p.experiments) {
            if (experiment.run != &run || run.FITC < 0 || run.FITC >= static_cast<ptrdiff_t>(run.species.size())) continue;
            auto* fitc = &run.species.at(run.FITC);
            for (const auto& entrance : experiment.entrances) {
                const auto concentration = entrance.CONC.find(fitc);
                if (concentration != entrance.CONC.end()) run.dye_conc_mgml = std::max(run.dye_conc_mgml, concentration->second);
            }
        }
        if (run.FITC >= 0 && run.FITC < static_cast<ptrdiff_t>(run.species.size())) {
            run.dye_conc = run.dye_conc_mgml * unit_conversion(&run, run.FITC,
                run.species.at(run.FITC).input_units, run.species.at(run.FITC).model_units);
        }
        if (!std::isfinite(run.dye_conc_mgml) || run.dye_conc_mgml <= 0.0 ||
            !std::isfinite(run.dye_conc) || run.dye_conc <= 0.0) {
            throw std::invalid_argument("Experiment '" + run.name +
                "' has no positive FITC inlet concentration. Check its raw_profile_entrance_concentrations rows.");
        }
    }
}

void 
read_specie_and_reaction_values_from_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Reading specie and reaction values from database:"); // TODO introduce species lookup and unit correction
    add_report(p, 1, "pre_db");
    experiment_run_struct* run_ptr;
    for (int run_index = 0; run_index < p.experiment_runs.size(); run_index++ ) {
        run_ptr = &p.experiment_runs.at(run_index);
        sqlite_error_message errMsg;
        std::string sqltext = "SELECT '" + std::to_string(run_index) + "' AS run, * FROM 'species' WHERE ";
        for (int specie = 0; specie < run_ptr->number_of_species; specie++) {
            sqltext.append("SPECIES_NAME = '" + run_ptr->species.at(specie).name + "'");
            if (specie != run_ptr->number_of_species - 1) {
                sqltext.append(" OR ");
            } else {
                sqltext.append(";");
            }
        }
        execute_sql(p, db, sqltext.c_str(), specie_db_callback, errMsg.out());

        sqltext = "SELECT '" + std::to_string(run_index) + "' AS run, * FROM 'reactions' WHERE ";
        for (int react = 0; react < run_ptr->number_of_reactions; react++) {
            sqltext.append("REACTION_NAME = '" + run_ptr->reactions.at(react).name + "'");
            if (react != run_ptr->number_of_reactions - 1) {
                sqltext.append(" OR ");
            } else {
                sqltext.append(";");
            }
        }
        execute_sql(p, db, sqltext.c_str(), reaction_db_callback, errMsg.out());

        for (const auto& [alias, variable] : run_ptr->alias_variables) {
            if (alias == "left_edge" || alias == "width") {
                throw std::runtime_error("Model alias '" + alias + "' conflicts with a reserved geometry variable");
            }
            const bool globally_solved = std::find(p.global_solve_for.begin(), p.global_solve_for.end(), alias) != p.global_solve_for.end();
            const bool locally_solved = std::find(run_ptr->solve_for.begin(), run_ptr->solve_for.end(), alias) != run_ptr->solve_for.end();
            if (!globally_solved && !locally_solved) {
                if (!variable.param_init) {
                    throw std::runtime_error("Model alias '" + alias + "' must be selected or have a paired numeric value");
                }
                continue;
            }
        }

        for (int item = 0; item < p.solvables.size(); item++) {
            auto& s = p.solvables.at(item);
            if (s.source_name == run_ptr->name || s.source_name == "global") {
                auto& target = variable_location(s.name, run_ptr);
                s.value() = target.value();
                s.param_init = true;
                target.source = &s;
                target.source_name = s.source_name;
                target.name = s.name;
            }
        }
    
        for (int specie = 0; specie < run_ptr->number_of_species; specie++) {
            add_report(p, 0, "specie: " + std::to_string(specie));
            specie_struct* specie_ptr = &run_ptr->species.at(specie);
            if (specie_ptr->type == 2 && !specie_ptr->diameter_alias.empty()) {
                // Diameter and its derived diffusion rate resolve after Variables load.
            } else if (specie_ptr->type == 2 && specie_ptr->diameter != 0.0d && specie_ptr->particle_density != 0.0d) {
                specie_ptr->diffusion_rate = 1.380649e-23d * run_ptr->temperature / (3.0d * M_PI * run_ptr->visc * specie_ptr->diameter);
            } else if (specie_ptr->type == 1) {
                // Skip. pulled dirrectly from db. 
                // TODO approximate value if not given in db read.
            } else {
                throw std::runtime_error("SPECIE:" + std::to_string(specie) + " diffusion_rate undefined (check Type, Diameter, or Density)");
            }
            specie_ptr->r = specie_ptr->diffusion_rate * run_ptr->dt / (run_ptr->W * run_ptr->W) * p.X * p.X;
            pop_report(p, 0);
        }
    }
    pop_report(p, 1);
    add_finishing_report(p, 3, "Done");

    for (int run = 0; run < p.experiment_runs.size(); run++) {
        run_ptr = &p.experiment_runs.at(run);
        run_ptr->FITC = specie_index(run_ptr, "FITC");
        run_ptr->PS_beads = specie_index(run_ptr, "PS_40nm", "PS_20nm");
        run_ptr->Bound_Dye_1 = specie_index(run_ptr, "40nm_Bound_Dye_1", "20nm_Bound_Dye_1");
        //run_ptr->Bound_Dye_2 = specie_index(run_ptr, "40nm_Bound_Dye_2", "20nm_Bound_Dye_2");
        run_ptr->FITC_Bead_1 = reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1");
        //run_ptr->FITC_Bead_2 = reaction_index(run_ptr, "FITC_40nm_2", "FITC_20nm_2");
        for (int react = 0; react < run_ptr->number_of_reactions; react++) {
            add_report(p, 3, "Reaction: " + run_ptr->reactions.at(react).name + ":");
            for (int specie = 0; specie < run_ptr->number_of_species; specie++) {
                add_report(p, 3, std::to_string(run_ptr->reactions.at(react).coef.at(&run_ptr->species.at(specie)).value()) + " " + run_ptr->species.at(specie).name + ",");
            }
            add_finishing_report(p, 3,"");
        }
    }
}

int 
raw_profiles_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    // SELECT * follows schema order; OMIT can come after the fit parameters.
    // Decide before registering any profile-local solver variables. The GUI and
    // ALGLIB both consume the resulting solvables list.
    bool omit = false;
    for (int i = 0; i < count; ++i) {
        if (std::string(columnNames[i]) == "OMIT" && argv[i] != nullptr) {
            std::string value = argv[i];
            removeSpaces(value);
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            omit = value == "true";
        }
    }
    add_report(p, 0, "callback_start");
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    std::string criterion;
    std::string s; // variable to store token obtained from the original string
    std::string value;
    double current_wt_percent;
    experiment_struct* exp_ptr;
    experiment_run_struct* run_ptr;
    specie_struct* specie_ptr;
    int col_index;
    ptrdiff_t row = 0;
    double last_profile_point = 0;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (argv[i] != NULL) {value = argv[i];} else {value.clear();}
        const std::string raw_value = value;
        removeSpaces(value);
        std::stringstream ss(value); // constructing stream from the string
        pop_and_add(p, 0, criterion);
        if (criterion == "NAME") {
            if (p.single_profile_experiment && value != *p.single_profile_experiment) {
                pop_report(p, 0);
                return 0;
            }
            for (ptrdiff_t j = 0; j < p.experiment_runs.size(); j++) {
                if (p.experiment_runs.at(j).name == value) {
                    run_ptr = &p.experiment_runs.at(j);
                }
            }
        } else if (criterion == "WT_PERCENT") {
            if (p.single_profile_dataset && value != *p.single_profile_dataset) {
                pop_report(p, 0);
                return 0;
            }
            if (!value.empty()) {
                for (ptrdiff_t j = 0; j < p.experiments.size(); j++) {
                    if (p.experiments.at(j).second_name == value && run_ptr == p.experiments.at(j).run) {
                        row = j;
                        exp_ptr = &p.experiments.at(j);
                        exp_ptr->omit = p.single_profile_experiment ? false : omit;
                        break;
                    } else if (j == p.experiments.size() - 1) {
                        pop_report(p, 0); // clear criterion
                        return 0; // omit this experiment. 
                    }
                }
                pop_and_add(p, 0, "row:" + std::to_string(row) + " " + exp_ptr->run->name + " " + exp_ptr->second_name);
                add_report(p, 0, criterion); // re-add criterion once row is found
            } else {throw std::runtime_error("No value in secondary name field in raw_profiles");}
        } else if (criterion == "CHANNEL_LEFT_EDGE") {
            if (!value.empty()) {
                exp_ptr->beginning_of_channel = std::stoi(value);
            } else {throw std::runtime_error("No value in CHANNEL_LEFT_EDGE field in raw_profiles");}
        } else if (criterion == "CHANNEL_RIGHT_EDGE") {
            if (!value.empty()) {
                exp_ptr->end_of_channel = std::stoi(value);
            } else {throw std::runtime_error("No value in CHANNEL_RIGHT_EDGE field in raw_profiles");}
        } else if (criterion == "INTENSITY_ARRAY") {
            exp_ptr->window_size = exp_ptr->end_of_channel - exp_ptr->beginning_of_channel + p.exp_left_padding + p.exp_right_padding;
            exp_ptr->window_start = p.total_window_size;
            p.total_window_size += exp_ptr->window_size;
            
            exp_ptr->channel_position.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->raw_experimental_profile.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->experimental_profile.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->model_profile.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->experimental_derivative.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->numeric_derivative.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->error.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->numeric_difference.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->experimental_difference.resize(exp_ptr->window_size, 0.0d);
            exp_ptr->analytical_zero.resize(exp_ptr->window_size, 0.0d);

            int x = 0;
            for (; x < p.exp_left_padding; x++) {
                exp_ptr->raw_experimental_profile.at(x) = 0.0d;
            }

            int j = 0;
            while (getline(ss, s, '	')) { 
                if (j >= (exp_ptr->beginning_of_channel) && j < exp_ptr->end_of_channel) { // Actual data 
                    const auto intensity = std::stod(s);
                    exp_ptr->raw_experimental_profile.at(x) = intensity;
                    exp_ptr->loaded_experimental_profile.push_back(intensity);
                    x++;
                } else if (j >= exp_ptr->end_of_channel || x >= exp_ptr->raw_experimental_profile.size()) {break;}
                j++;
            }

            for (; x < exp_ptr->raw_experimental_profile.size(); x++) {
                exp_ptr->raw_experimental_profile.at(x) = exp_ptr->raw_experimental_profile.at(x - 1);
            }
        } else if (criterion == "INDEPENDENT_PARAMETERS_TO_SOLVE_FOR") {
            add_report(p, 0, "Value=" + value);
            if ((!omit || p.single_profile_experiment) && !value.empty()) {
                while (getline(ss, s, ' ')) {
                    if (p.single_profile_experiment &&
                        std::find(p.single_profile_variables.begin(), p.single_profile_variables.end(), s) == p.single_profile_variables.end()) continue;
                    if (std::find(p.solve_for.begin(), p.solve_for.end(), s) == p.solve_for.end()) { // to avoid repeats in p.solve_for
                        p.solve_for.push_back(s);
                    }
                    p.solvables.emplace_back(solvable(s));
                    p.solvables.back().source_name = run_ptr->name + "_" + exp_ptr->second_name;
                    
                    if (s == "left_edge") {
                        exp_ptr->left_edge.source = &p.solvables.back();
                        exp_ptr->left_edge.value() = run_ptr->left_edge.value();
                        exp_ptr->left_edge.init();
                    } else if (s == "width") {
                        exp_ptr->width.source = &p.solvables.back();
                        exp_ptr->width.value() = run_ptr->width.value();
                        exp_ptr->width.init();         
                    } /*else if (s == "p1") {
                        exp_ptr->p1 = p.solvables.back();
                    } else if (s == "keq1") {
                        exp_ptr->keq1 = p.solvables.back();
                    } else if (s == "QE1") {
                        exp_ptr->QE1 = p.solvables.back();
                    }*/
                }
            }
            pop_report(p, 0);
        } else if (criterion == "LEFT_EDGE") {
            add_report(p, 0, "LEFT_EDGE Value=" + value);
            if (!value.empty() && exp_ptr->left_edge.param_init) {
                exp_ptr->left_edge.value() = std::max(std::stod(value) + p.exp_left_padding, 0.0d);

            }
            pop_report(p, 0);
        } else if (criterion == "WIDTH") {
            add_report(p, 0, "WIDTH Value=" + value);
            if (!value.empty() && exp_ptr->width.param_init) {
                exp_ptr->width.value() = std::stod(value);
            }
            pop_report(p, 0);
        } else if (criterion == "OMIT") {
            // Already applied when the experiment was identified above.
        }
    };
    pop_report(p, 0); // clear final criterion
    pop_report(p, 0); // clear row count
    return 0;
}

void 
read_raw_profiles_from_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Reading experimental profiles from database:");
    for (int i = 0; i < p.experiments.size(); i++) {
        experiment_struct* exp_ptr = &p.experiments.at(i);
        experiment_run_struct* run_ptr = exp_ptr->run;
        exp_ptr->left_edge.source = &run_ptr->left_edge;
        exp_ptr->width.source = &run_ptr->width;
        /*exp_ptr->p1         = solvable(run_ptr->p1, "p1");
        exp_ptr->keq1       = solvable(run_ptr->keq1, "keq1");
        exp_ptr->QE1        = solvable(run_ptr->QE1, "QE1");*/
    }
    add_report(p, 1, "reading_db");
    sqlite_error_message errMsg;
    std::string sqltext = "SELECT * FROM 'raw_profile' WHERE ";
    int number_of_runs = p.experiment_runs.size();
    for (int run = 0; run < number_of_runs; run++) {
        experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
        sqltext.append("(NAME = '" + run_ptr->name + "')");
        if (run != number_of_runs - 1) {
            sqltext.append(" OR ");
        } else {
            sqltext.append(";");
        }
    }
    const char* sql = sqltext.c_str();
    execute_sql(p, db, sql, raw_profiles_db_callback, errMsg.out());

    for (const auto& experiment : p.experiments) {
        if (experiment.omit) continue;
        auto validate_alias = [&](const std::string& alias, bool has_fallback) {
            if (alias.empty()) return;
            const bool selected = alias_is_selected(p, *experiment.run, alias) ||
                find_profile_alias_parameter(p, experiment, alias) != nullptr;
            const bool has_parameter = find_profile_alias_parameter(p, experiment, alias) != nullptr;
            if ((!selected && !has_fallback) ||
                (selected && !has_parameter))
                throw std::runtime_error("Raw profile ENTRANCE_CONC alias '" + alias +
                    "' must be selected for its profile, experiment, or globally, or include a numeric fallback");
        };
        for (const auto& override : experiment.entrance_concentrations)
            validate_alias(override.alias, override.has_paired_value);
    }

    sqlite3_stmt* table_check_raw = nullptr;
    int table_check_rc = sqlite3_prepare_v2(db,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='raw_profile_entrance_concentrations'",
        -1, &table_check_raw, nullptr);
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> table_check(table_check_raw, sqlite3_finalize);
    if (table_check_rc != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
    const int table_step = sqlite3_step(table_check.get());
    if (table_step == SQLITE_ROW) {
        for (auto& experiment : p.experiments) {
            sqlite3_stmt* query_raw = nullptr;
            const int prepared = sqlite3_prepare_v2(db,
                "SELECT ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS FROM raw_profile_entrance_concentrations WHERE NAME IS ? AND WT_PERCENT IS ?",
                -1, &query_raw, nullptr);
            std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> query(query_raw, sqlite3_finalize);
            if (prepared != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
            if (sqlite3_bind_text(query.get(), 1, experiment.run->name.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK ||
                sqlite3_bind_text(query.get(), 2, experiment.second_name.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db));
            int rc;
            while ((rc = sqlite3_step(query.get())) == SQLITE_ROW) {
                const int entrance_number = sqlite3_column_int(query.get(), 0);
                const auto* species = reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 1));
                const auto* units = reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 3));
                const auto* value_text = reinterpret_cast<const char*>(sqlite3_column_text(query.get(), 2));
                if (entrance_number < 1 || !species || !value_text)
                    throw std::runtime_error("Invalid raw profile entrance concentration row.");
                const auto parsed = parse_value_alias(value_text, "Raw profile ENTRANCE_CONC");
                if (!std::isfinite(parsed.value) || parsed.value < 0.0 || parsed.sign < 0.0 ||
                    (parsed.alias.empty() && !parsed.has_value))
                    throw std::runtime_error("Invalid raw profile entrance concentration row.");
                experiment.entrance_concentrations.push_back({entrance_number, species, parsed.value,
                    units ? units : "", parsed.alias, parsed.has_value});
            }
            if (rc != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
        }
    } else if (table_step != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(db));
    pop_and_add(p, 1, "setting experiments.entrances size and flowrates");

    for (int row = 0; row < p.row_count; row++) {
        experiment_struct* exp_ptr = &p.experiments.at(row);
        exp_ptr->entrances.resize(exp_ptr->run->ENTRANCE_FLOWRATE.size());
        for (int entrance = 0; entrance < exp_ptr->entrances.size(); entrance++) {
            exp_ptr->entrances.at(entrance).ENTRANCE_FLOWRATE = exp_ptr->run->ENTRANCE_FLOWRATE.at(entrance);
        }
        for (int zrow = 0; zrow < p.row_count; zrow++) {
            if (exp_ptr->run == p.experiments.at(zrow).run && p.experiments.at(zrow).second_name == "0.0") {
                exp_ptr->zero_row_ptr = &p.experiments.at(zrow);
                break;
            } else if (zrow == p.row_count - 1) {
                exp_ptr->zero_row_ptr = NULL;
            }
        }
    }
    pop_report(p, 1);
    add_finishing_report(p, 3, "Done");
}

void
prepare_parameter_scope_selection(parameters_t& p, sqlite3* db)
{
    struct profile_selection {
        std::string experiment;
        std::string dataset;
        bool omitted = false;
        std::vector<std::string> variables;
    };
    std::map<std::string, std::vector<std::string>> experiment_variables;
    std::vector<profile_selection> profiles;
    const auto column_text = [](sqlite3_stmt* statement, int column) {
        const auto* value = reinterpret_cast<const char*>(sqlite3_column_text(statement, column));
        return value ? std::string(value) : std::string{};
    };
    const auto parse_variables = [](const std::string& text) {
        std::vector<std::string> variables;
        std::stringstream stream(text);
        std::string variable;
        while (stream >> variable) {
            if (std::find(variables.begin(), variables.end(), variable) == variables.end())
                variables.push_back(std::move(variable));
        }
        return variables;
    };
    const auto find_column = [](sqlite3_stmt* statement, const char* name) {
        for (int column = 0; column < sqlite3_column_count(statement); ++column)
            if (std::string(sqlite3_column_name(statement, column)) == name) return column;
        return -1;
    };

    sqlite3_stmt* raw_experiments = nullptr;
    int rc = sqlite3_prepare_v2(db, "SELECT * FROM experiments", -1, &raw_experiments, nullptr);
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> experiments(raw_experiments, sqlite3_finalize);
    if (rc != SQLITE_OK)
        throw std::runtime_error(std::string("Cannot inspect experiment parameter selections: ") + sqlite3_errmsg(db));
    while ((rc = sqlite3_step(experiments.get())) == SQLITE_ROW) {
        const int name_column = find_column(experiments.get(), "NAME");
        const int variables_column = find_column(experiments.get(), "PARAMETERS_TO_SOLVE_FOR");
        if (name_column < 0 || variables_column < 0) continue;
        const auto name = column_text(experiments.get(), name_column);
        if (std::none_of(p.experiment_runs.begin(), p.experiment_runs.end(),
                         [&](const auto& run) { return run.name == name; })) continue;
        experiment_variables[name] = parse_variables(column_text(experiments.get(), variables_column));
    }
    if (rc != SQLITE_DONE)
        throw std::runtime_error(std::string("Cannot read experiment parameter selections: ") + sqlite3_errmsg(db));

    sqlite3_stmt* raw_profiles = nullptr;
    rc = sqlite3_prepare_v2(db, "SELECT * FROM raw_profile", -1, &raw_profiles, nullptr);
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> raw_profile_rows(raw_profiles, sqlite3_finalize);
    if (rc != SQLITE_OK)
        throw std::runtime_error(std::string("Cannot inspect raw-profile parameter selections: ") + sqlite3_errmsg(db));
    while ((rc = sqlite3_step(raw_profile_rows.get())) == SQLITE_ROW) {
        const int name_column = find_column(raw_profile_rows.get(), "NAME");
        const int dataset_column = find_column(raw_profile_rows.get(), "WT_PERCENT");
        if (name_column < 0 || dataset_column < 0) continue;
        profile_selection profile;
        profile.experiment = column_text(raw_profile_rows.get(), name_column);
        profile.dataset = column_text(raw_profile_rows.get(), dataset_column);
        const int omit_column = find_column(raw_profile_rows.get(), "OMIT");
        if (omit_column >= 0) {
            auto omit = column_text(raw_profile_rows.get(), omit_column);
            removeSpaces(omit);
            std::transform(omit.begin(), omit.end(), omit.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            profile.omitted = omit == "true";
        }
        const int variables_column = find_column(raw_profile_rows.get(), "INDEPENDENT_PARAMETERS_TO_SOLVE_FOR");
        if (variables_column >= 0)
            profile.variables = parse_variables(column_text(raw_profile_rows.get(), variables_column));
        if (std::none_of(p.experiment_runs.begin(), p.experiment_runs.end(),
                         [&](const auto& run) { return run.name == profile.experiment; })) continue;
        if (p.single_profile_experiment &&
            (profile.experiment != *p.single_profile_experiment ||
             profile.dataset != *p.single_profile_dataset)) continue;
        if (!p.single_profile_experiment && profile.omitted) continue;
        profiles.push_back(std::move(profile));
    }
    if (rc != SQLITE_DONE)
        throw std::runtime_error(std::string("Cannot read raw-profile parameter selections: ") + sqlite3_errmsg(db));

    std::vector<std::string> covered_everywhere;
    bool first_profile = true;
    for (const auto& profile : profiles) {
        const auto& experiment_scope = experiment_variables[profile.experiment];
        std::vector<std::string> covered = profile.variables;
        for (const auto& variable : experiment_scope)
            if (std::find(covered.begin(), covered.end(), variable) == covered.end())
                covered.push_back(variable);
        if (first_profile) {
            covered_everywhere = std::move(covered);
            first_profile = false;
        } else {
            covered_everywhere.erase(std::remove_if(covered_everywhere.begin(), covered_everywhere.end(),
                [&](const auto& variable) {
                    return std::find(covered.begin(), covered.end(), variable) == covered.end();
                }), covered_everywhere.end());
        }
    }
    for (const auto& [experiment, variables] : experiment_variables) {
        for (const auto& variable : variables) {
            std::size_t profile_count = 0;
            const bool every_profile_selects_it = std::all_of(profiles.begin(), profiles.end(),
                [&](const auto& profile) {
                    if (profile.experiment != experiment) return true;
                    ++profile_count;
                    return std::find(profile.variables.begin(), profile.variables.end(), variable) !=
                           profile.variables.end();
                });
            if (profile_count > 0 && every_profile_selects_it)
                p.suppressed_experiment_parameters[experiment].push_back(variable);
        }
    }
    if (!first_profile) {
        p.suppressed_universal_parameters = std::move(covered_everywhere);
        for (const auto& variable : p.suppressed_universal_parameters) {
            p.global_solve_for.erase(std::remove(p.global_solve_for.begin(), p.global_solve_for.end(), variable),
                                     p.global_solve_for.end());
            p.solve_for.erase(std::remove(p.solve_for.begin(), p.solve_for.end(), variable), p.solve_for.end());
            p.solvables.erase(std::remove_if(p.solvables.begin(), p.solvables.end(), [&](const auto& solvable) {
                return solvable.source_name == "global" && solvable.name == variable;
            }), p.solvables.end());
        }
    }
}


void 
initialize_profile_entrance_concentrations(parameters_t& p, sqlite3* db)
{
    add_report(p, 3, "Initializing profile entrance concentrations:");
    add_report(p, 1, "profile_values:");

    for (int row = 0; row < p.row_count; row++) {
        //initialize the concentration inlet arrays.
        experiment_struct* exp_ptr = &p.experiments.at(row);
        experiment_run_struct* run_ptr = exp_ptr->run;
        for (const auto& override : exp_ptr->entrance_concentrations) {
            if (exp_ptr->entrances.empty()) throw std::runtime_error("Per-species ENTRANCE_CONC requires at least one entrance.");
            if (override.entrance_number < 1 || override.entrance_number > static_cast<int>(exp_ptr->entrances.size()))
                throw std::runtime_error("ENTRANCE_CONC refers to an inlet that is missing from ENTRANCE_FLOWRATE.");
            if (!override.alias.empty()) continue;
            const auto specie = specie_index(run_ptr, override.species_name);
            const auto& inlet_species = run_ptr->species.at(specie);
            const std::string source_units = override.units.empty() ? inlet_species.model_units : override.units;
            const double input_concentration = override.concentration *
                unit_conversion(run_ptr, specie, source_units, inlet_species.input_units);
            if (!std::isfinite(input_concentration))
                throw std::runtime_error("ENTRANCE_CONC conversion produced a non-finite value.");
            exp_ptr->entrances.at(override.entrance_number - 1).CONC[&run_ptr->species.at(specie)] = input_concentration;
        }
        // The inlet and outlet Concentration Arrays.
        exp_ptr->species_out.assign(run_ptr->number_of_species, std::vector<double>(p.X, 0.0));
        pop_and_add(p, 1, "intit row:" + std::to_string(row));
        for (int entrance = 0; entrance < exp_ptr->entrances.size(); entrance++) {
            entrance_struct* entr_ptr = &exp_ptr->entrances.at(entrance);
            for (auto specie : entr_ptr->CONC) {
                if (specie.first->name == "FITC") {
                    run_ptr->dye_conc_mgml = std::max(specie.second, run_ptr->dye_conc_mgml);
                    break;
                }
            }
        }
        run_ptr->dye_conc = run_ptr->dye_conc_mgml * unit_conversion(run_ptr, run_ptr->FITC, run_ptr->species.at(run_ptr->FITC).input_units, run_ptr->species.at(run_ptr->FITC).model_units);
    }
    pop_report(p, 1);
    add_finishing_report(p, 3, "Done");
}

int 
alglib_input_db_callback(void *data, int count, char **argv, char **columnNames)
{
    auto& p = *static_cast<parameters_t*>(data);
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count->is the number of columns
    //columnNames-> array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv->array of pointers to strings obtained as if from [sqlite3_column_text()]
    std::string criterion;
    std::string current_variable;

    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (argv[i] == NULL) {throw std::runtime_error("No " + criterion + "for" + current_variable);}
        if (criterion == "VARIABLE") {
            current_variable = argv[i];
        } else if (criterion == "INITIAL VALUE") {
            p.initial_values_alglib_map.insert({current_variable, std::stod(argv[i])});
        } else if (criterion == "LOWER BOUND") {
            p.low_bound_map.insert({current_variable, std::stod(argv[i])});
        } else if (criterion == "UPPER BOUND") {
            p.up_bound_map.insert({current_variable, std::stod(argv[i])});
        } else if (criterion == "SCALE") {
            if (std::stod(argv[i]) == 0) {
                throw std::runtime_error("zero passed to SCALE from variable:" + current_variable);
            }
            p.scale_map.insert({current_variable, std::stod(argv[i])});
        }
    };
    //}
    return 0;
}

void 
read_alglib_values_from_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Reading ALGLIB inputs from database");
    p.scale.resize(p.solvables.size());
    p.initial_values_alglib.resize(p.solvables.size());
    p.low_bound.resize(p.solvables.size());
    p.up_bound.resize(p.solvables.size());
    if (p.solvables.empty()) {
        for (auto& run : p.experiment_runs) synchronize_run_aliases(p, run);
        synchronize_profile_entrance_concentrations(p, true);
        return;
    }
    
    sqlite_error_message errMsg;
    std::string sqltext = "SELECT * FROM 'alglib_input' WHERE ";
    
    for (int i = 0; i < p.solve_for.size(); i++) {
        sqltext.append("VARIABLE = '" + p.solve_for.at(i) + "' ");
        if (i != p.solve_for.size() - 1) {
            sqltext.append("OR ");
        } else {
            sqltext.append(";");
        }
    }
    
    execute_sql(p, db, sqltext.c_str(), alglib_input_db_callback, errMsg.out());

    for (const auto& run : p.experiment_runs) {
        for (const auto& [alias, variable] : run.alias_variables) {
            if (alias_is_selected(p, run, alias) && !p.initial_values_alglib_map.contains(alias)) {
                throw std::runtime_error("Model alias '" + alias + "' is missing from the Variables table or has no solve-for bounds");
            }
        }
    }

    
    for (int i = 0; i < p.solvables.size(); i++) {
        auto& s = p.solvables.at(i);
        const bool model_alias = std::any_of(p.experiment_runs.begin(), p.experiment_runs.end(), [&](const auto& run) {
            return run.alias_variables.contains(s.name);
        });
        if (model_alias) {
            p.initial_values_alglib[i] = p.initial_values_alglib_map.at(s.name);
            s.value() = p.initial_values_alglib[i];
        } else if ((p.use_alglib_init_values || s.source_name == "global") && !s.param_init) {
            p.initial_values_alglib[i] = p.initial_values_alglib_map.at(s.name);
        } else {
            p.initial_values_alglib[i] = s.value();
        }
        
        p.low_bound[i] = p.low_bound_map.at(s.name);
        p.up_bound[i] = p.up_bound_map.at(s.name); 
        p.scale[i] = p.scale_map.at(s.name);

        if (p.single_profile_experiment) {
            const auto requested = p.single_profile_initial_values.find(s.name);
            if (requested != p.single_profile_initial_values.end()) {
                p.initial_values_alglib[i] = requested->second;
                s.value() = requested->second;
            }
            // Validate profile-derived defaults as well as explicit UI overrides.
            // Otherwise ALGLIB silently projects an out-of-range start to a bound.
            if (!std::isfinite(p.initial_values_alglib[i]) ||
                p.initial_values_alglib[i] < p.low_bound[i] ||
                p.initial_values_alglib[i] > p.up_bound[i])
                throw std::runtime_error("Starting value for '" + s.name + "' must be finite and within its configured bounds");
        }

    }

    for (auto& run : p.experiment_runs) synchronize_run_aliases(p, run);
    synchronize_profile_entrance_concentrations(p, true);
}

void 
write_model_profile_to_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Writing model profiles to database:");
    sqlite_error_message errMsg;
    std::string sqltext;
    struct transaction_guard {
        sqlite3* db;
        bool active = true;
        explicit transaction_guard(sqlite3* connection) : db(connection) {
            if (sqlite3_exec(db, "BEGIN IMMEDIATE", nullptr, nullptr, nullptr) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db));
        }
        void commit() {
            if (sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db));
            active = false;
        }
        ~transaction_guard() { if (active) sqlite3_exec(db, "ROLLBACK", nullptr, nullptr, nullptr); }
    } transaction(db);
    auto quote = [](const std::string& value) {
        std::string escaped;
        escaped.reserve(value.size() + 2);
        escaped.push_back('\'');
        for (char c : value) { escaped.push_back(c); if (c == '\'') escaped.push_back('\''); }
        escaped.push_back('\'');
        return escaped;
    };
    std::string experiments_fitted;
    for (const auto& run : p.experiment_runs) {
        if (!experiments_fitted.empty()) experiments_fitted += " ";
        experiments_fitted += run.name;
    }
    std::string parameters_solved;
    for (const auto& parameter : p.solvables) {
        if (!parameters_solved.empty()) parameters_solved += " ";
        parameters_solved += parameter.name;
    }
    const std::string reactions_enabled = p.disable_reactions ? "all reactions disabled"
        : p.disable_reverse_reactions ? "all reverse reactions disabled" : "all reactions enabled";
    const auto clock_run_at_utc_us = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    sqlite3_stmt* latest_calculation = nullptr;
    int latest_rc = sqlite3_prepare_v2(db,
        "SELECT MAX(RUN_AT_UTC_US) FROM profile_calculations", -1, &latest_calculation, nullptr);
    if (latest_rc != SQLITE_OK)
        throw std::runtime_error(std::string("Error reading latest profile calculation: ") + sqlite3_errmsg(db));
    const int latest_step = sqlite3_step(latest_calculation);
    const sqlite3_int64 previous_run_at_utc_us = latest_step == SQLITE_ROW &&
        sqlite3_column_type(latest_calculation, 0) != SQLITE_NULL
        ? sqlite3_column_int64(latest_calculation, 0) : 0;
    sqlite3_finalize(latest_calculation);
    if (latest_step != SQLITE_ROW)
        throw std::runtime_error(std::string("Error reading latest profile calculation: ") + sqlite3_errmsg(db));
    if (previous_run_at_utc_us == std::numeric_limits<sqlite3_int64>::max())
        throw std::runtime_error("Profile calculation timestamp is out of range.");
    // The viewer selects MAX(RUN_AT_UTC_US) for each profile. Keep IDs strictly
    // increasing even when two saves share a clock tick or the wall clock moves back.
    const auto run_at_utc_us = std::max<sqlite3_int64>(
        static_cast<sqlite3_int64>(clock_run_at_utc_us), previous_run_at_utc_us + 1);
    sqltext = "INSERT INTO profile_calculations(RUN_AT_UTC_US,ALL_EXP_FITTED,PARAMETERS_SOLVED_FOR,REACTIONS_ENABLED,SCATTER_METHOD,X_RESOLUTION,Z_RESOLUTION) VALUES(" +
        std::to_string(run_at_utc_us) + "," +
        quote(experiments_fitted) + "," + quote(parameters_solved) + "," +
        quote(reactions_enabled) + "," + quote(p.scatter_correction_type) + "," +
        std::to_string(p.X) + "," + std::to_string(p.Z) + ");";
    execute_sql(p, db, sqltext.c_str(), 0, errMsg.out());
    for (std::size_t row = 0; row < p.experiments.size(); ++row) {
        const auto& exp = p.experiments[row];
        sqltext = "INSERT INTO profile_records(RUN_AT_UTC_US,EXPERIMENT_NAME,WT_PERCENT,EXP_DA,MODEL_DA,EXP_INTEGRAL,MODEL_INTEGRAL) VALUES(" +
            std::to_string(run_at_utc_us) + "," + quote(exp.run->name) + "," + quote(exp.second_name) + "," +
            quote(double_to_string(exp.exp_DA)) + "," + quote(double_to_string(exp.model_DA)) + "," +
            quote(double_to_string(exp.exp_integral)) + "," + quote(double_to_string(exp.model_integral)) + ");";
        execute_sql(p, db, sqltext.c_str(), 0, errMsg.out());
    }
    for (auto& s : p.solvables) {
        std::string source = "global";
        std::string experiment_name;
        std::string wt_percent;
        for (std::size_t row = 0; row < p.experiments.size(); ++row) {
            const auto& exp = p.experiments[row];
            if (s.source_name == exp.run->name + "_" + exp.second_name) {
                source.clear(); experiment_name = exp.run->name; wt_percent = exp.second_name; break;
            }
            if (s.source_name == exp.run->name) { source = exp.run->name; break; }
        }
        sqltext = "INSERT INTO profile_record_parameters(RUN_AT_UTC_US,EXPERIMENT_NAME,WT_PERCENT,PARAMETER,VALUE,UNITS,SOURCE) VALUES(" +
            std::to_string(run_at_utc_us) + "," + quote(experiment_name) + "," + quote(wt_percent) + "," +
            quote(s.name) + "," + quote(double_to_string(s.value())) + ",NULL," + quote(source) + ");";
        execute_sql(p, db, sqltext.c_str(), 0, errMsg.out());
    }

    sqltext.clear();
    pop_and_add(p, 1, "preparing items into profile_record_samples table");
    struct model_profile_struct {
        std::string RUN_AT_UTC_US;
        std::string EXPERIMENT_NAME;
        std::string WT_PERCENT;
        std::string X;
        std::string Free_Dye;
        std::string Bound_Dye;
        std::string Total_Dye;
        std::string Unbound_Beads;
        std::string Bound_Beads;
        std::string Total_Beads;
        std::string Experimental_Derivative;
        std::string Numeric_Derivative;
        std::string Experimental;
        std::string Numeric;
        std::string Error;
        std::string experimental_difference;
        std::string numeric_difference;
        model_profile_struct(sqlite3_int64 RUN_AT_UTC_US, std::string EXPERIMENT_NAME, std::string WT_PERCENT, double X)
            : RUN_AT_UTC_US(std::to_string(RUN_AT_UTC_US))
            , EXPERIMENT_NAME(std::move(EXPERIMENT_NAME))
            , WT_PERCENT(std::move(WT_PERCENT))
            , X([&] {
                std::ostringstream formatted;
                formatted << std::setprecision(std::numeric_limits<double>::max_digits10) << X;
                return formatted.str();
            }())
        {}
    };
    std::vector<model_profile_struct> model_profiles;
    experiment_run_struct* run_ptr;
    experiment_struct* exp_ptr;
    for (int row = 0; row < p.row_count; row++) {
        add_report(p, 0, "row:" + std::to_string(row));
        add_report(p, 0, "pre-while loop");
        exp_ptr = &p.experiments.at(row);
        run_ptr = exp_ptr->run;
        const double inv_dye_conc = 1.0 / run_ptr->dye_conc;
        ptrdiff_t FITC = run_ptr->FITC;
        specie_struct* FITC_ptr = &run_ptr->species.at(FITC);
        
        ptrdiff_t PS_beads = run_ptr->PS_beads;
        specie_struct* PS_beads_ptr = &run_ptr->species.at(PS_beads);
        ptrdiff_t Bound_Dye_1 = run_ptr->Bound_Dye_1;
        specie_struct* Bound_Dye_1_ptr = &run_ptr->species.at(Bound_Dye_1);
        ptrdiff_t FITC_Bead_1 = run_ptr->FITC_Bead_1;
        reaction_struct* FITC_Bead_1_ptr = &run_ptr->reactions.at(FITC_Bead_1);

        const double FITC_unit = unit_conversion(run_ptr, FITC, FITC_ptr->model_units, FITC_ptr->input_units);
        const double bead_unit = unit_conversion(run_ptr, PS_beads, PS_beads_ptr->model_units, PS_beads_ptr->input_units);
        std::vector<double> model = exp_ptr->model_profile;
        std::vector<double> experiment = exp_ptr->experimental_profile;
        const auto& out = exp_ptr->species_out;

        const double model_scale = run_ptr->W * 1e6 / p.X;
        const double data_scale  = exp_ptr->scale_factor;

        int model_x = 0;
        int data_x  = 0;
        double X_model = 0.0;
        double X_data = exp_ptr->channel_position.at(data_x);
        double last_X = -1;

        assert(data_scale > 0);
        bool model_match = false;
        bool data_match = false;
        pop_and_add(p, 0, "while loop");
        while (X_model < run_ptr->W * 1e6 || X_data < run_ptr->W * 1e6) {
            X_data = exp_ptr->channel_position.at(data_x);
            X_model = model_scale * model_x;
            model_match = X_data >= X_model;
            data_match = X_model >= X_data;
            if ( (model_match && model_x < p.X) || (data_match && data_x < exp_ptr->window_size)) {
                if (std::min(X_data, X_model) != last_X) {
                    model_profiles.emplace_back(run_at_utc_us, run_ptr->name, exp_ptr->second_name, std::min(X_data, X_model));
                    last_X = std::min(X_data, X_model);
                }
            }
            auto& profile = model_profiles.back();
            if (model_match) {
                if (model_x < p.X) {
                    const double free_dye = out[FITC][model_x] * FITC_unit;
                    const double bound_dye = (out[Bound_Dye_1][model_x]) * FITC_unit;
                    profile.Free_Dye                = double_to_string(free_dye);
                    profile.Bound_Dye               = double_to_string(bound_dye);
                    profile.Total_Dye               = double_to_string(free_dye + bound_dye);
                    const double unbound_beads = out[PS_beads][model_x] * bead_unit;
                    const double bound_beads = out[Bound_Dye_1][model_x] / abs(FITC_Bead_1_ptr->coef.at(FITC_ptr).value()) * bead_unit;
                    profile.Unbound_Beads           = double_to_string(unbound_beads);
                    profile.Bound_Beads             = double_to_string(bound_beads);
                    profile.Total_Beads             = double_to_string(unbound_beads + bound_beads);
                }
                model_x++;
            }
            if (data_match) {
                if (data_x < exp_ptr->window_size) {
                    profile.Experimental_Derivative = double_to_string(exp_ptr->experimental_derivative.at(data_x));
                    profile.Numeric_Derivative      = double_to_string(exp_ptr->numeric_derivative.at(data_x));
                    profile.Experimental            = double_to_string(experiment.at(data_x) * run_ptr->dye_conc_mgml);
                    profile.Numeric                 = double_to_string(model.at(data_x) * FITC_unit);
                    profile.Error                   = double_to_string(exp_ptr->error.at(data_x));
                    profile.experimental_difference = double_to_string(exp_ptr->experimental_difference.at(data_x));
                    profile.numeric_difference      = double_to_string(exp_ptr->numeric_difference.at(data_x));
                }
                data_x++;
            }
        }
        pop_report(p, 0);
        pop_report(p, 0);
    }
    pop_and_add(p, 1, "profile_record_samples table");
    if (model_profiles.size() > 0) {
        static constexpr const char* insert_samples =
            "INSERT INTO profile_record_samples (RUN_AT_UTC_US, EXPERIMENT_NAME, WT_PERCENT, X, Free_Dye, Bound_Dye, Total_Dye, Unbound_Beads, Bound_Beads, Total_Beads, Experimental_Derivative, Numeric_Derivative, Experimental, Numeric, Error, Experimental_Difference, Numeric_Difference) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
        sqlite3_stmt* insert_statement = nullptr;
        int rc = sqlite3_prepare_v2(db, insert_samples, -1, &insert_statement, nullptr);
        if (rc != SQLITE_OK) {
            throw tsensor_workflow::workflow_error(p.active_operation,
                tsensor_workflow::error_code::database,
                std::string("Error preparing profile sample insert: ") + sqlite3_errmsg(db), rc);
        }
        struct statement_guard {
            sqlite3_stmt* statement;
            ~statement_guard() { sqlite3_finalize(statement); }
        } insert_guard{insert_statement};
        auto bind_text = [&](int index, const std::string& value) {
            return sqlite3_bind_text(insert_statement, index, value.c_str(),
                static_cast<int>(value.size()), SQLITE_TRANSIENT);
        };
        for (const auto& profile : model_profiles) {
            rc = sqlite3_bind_int64(insert_statement, 1, std::stoll(profile.RUN_AT_UTC_US));
            if (rc == SQLITE_OK) rc = bind_text(2, profile.EXPERIMENT_NAME);
            if (rc == SQLITE_OK) rc = bind_text(3, profile.WT_PERCENT);
            if (rc == SQLITE_OK) rc = bind_text(4, profile.X);
            if (rc == SQLITE_OK) rc = bind_text(5, profile.Free_Dye);
            if (rc == SQLITE_OK) rc = bind_text(6, profile.Bound_Dye);
            if (rc == SQLITE_OK) rc = bind_text(7, profile.Total_Dye);
            if (rc == SQLITE_OK) rc = bind_text(8, profile.Unbound_Beads);
            if (rc == SQLITE_OK) rc = bind_text(9, profile.Bound_Beads);
            if (rc == SQLITE_OK) rc = bind_text(10, profile.Total_Beads);
            if (rc == SQLITE_OK) rc = bind_text(11, profile.Experimental_Derivative);
            if (rc == SQLITE_OK) rc = bind_text(12, profile.Numeric_Derivative);
            if (rc == SQLITE_OK) rc = bind_text(13, profile.Experimental);
            if (rc == SQLITE_OK) rc = bind_text(14, profile.Numeric);
            if (rc == SQLITE_OK) rc = bind_text(15, profile.Error);
            if (rc == SQLITE_OK) rc = bind_text(16, profile.experimental_difference);
            if (rc == SQLITE_OK) rc = bind_text(17, profile.numeric_difference);
            if (rc == SQLITE_OK) rc = sqlite3_step(insert_statement);
            if (rc != SQLITE_DONE) {
                throw tsensor_workflow::workflow_error(p.active_operation,
                    tsensor_workflow::error_code::database,
                    std::string("Error saving profile sample at X=") + profile.X + ": " + sqlite3_errmsg(db) +
                        "\n" + insert_samples, rc);
            }
            rc = sqlite3_reset(insert_statement);
            if (rc != SQLITE_OK) {
                throw tsensor_workflow::workflow_error(p.active_operation,
                    tsensor_workflow::error_code::database,
                    std::string("Error resetting profile sample insert: ") + sqlite3_errmsg(db), rc);
            }
            sqlite3_clear_bindings(insert_statement);
        }
    }

    transaction.commit();
    pop_report(p, 1);
    add_finishing_report(p, 3, "Done");
}

void 
write_alglib_values_to_db(parameters_t& p, sqlite3* db) //reading data using callback functions
{
    add_report(p, 3, "Writing alglib_values to database:");
    add_report(p, 1, "init values and scale");
    sqlite_error_message errMsg;
    std::string sqltext;
    for (int index = 0; index < p.global_solve_for.size(); index++) {
        sqltext.append("UPDATE alglib_input SET 'INITIAL VALUE' = '" + double_to_string(p.initial_values_alglib[index]) + "' WHERE VARIABLE = '" + p.global_solve_for[index] + "'; ");
        if (p.initial_values_alglib[index] != 0) {
            //sqltext.append("UPDATE alglib_input SET 'SCALE' = '" + double_to_string(p.initial_values_alglib[index]) + "' WHERE VARIABLE = '" + p.global_solve_for[index] + "'; ");
        }
    }
    for (int run = 0; run < p.experiment_runs.size(); run++) {
        experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
        pop_and_add(p, 1, "run " + run_ptr->name); // removes init val state
        add_report(p, 0, "QE");
        for (int specie = 0; specie < run_ptr->number_of_species; specie++) { 
            sqltext.append("UPDATE species SET 'QE' = '" + value_with_alias(run_ptr->species.at(specie).QE.value(), run_ptr->species.at(specie).qe_alias) + "' WHERE SPECIES_NAME = '" + run_ptr->species.at(specie).name + "'; ");
        }
        for (int reaction = 0; reaction < run_ptr->number_of_reactions; reaction++) {
            reaction_struct* react_ptr = &run_ptr->reactions.at(reaction);
            pop_and_add(p, 0, "reaction " + react_ptr->name); // removes QE state then previous reaction state
            add_report(p, 0, "ks");
            sqltext.append("UPDATE reactions SET 'Ks' = '" + value_with_alias(react_ptr->k[0].value(), react_ptr->k_alias[0]) + " " + value_with_alias(react_ptr->k[1].value(), react_ptr->k_alias[1]) + "' WHERE REACTION_NAME = '" + react_ptr->name + "'; ");
            pop_and_add(p, 0, "coefs"); // removes ks state
            sqltext.append("UPDATE reactions SET 'COEFFICIENTS' = '");
            for (std::ptrdiff_t specie = 0; specie < react_ptr->specie_vect.size(); specie++) {
                specie_struct* specie_ptr = react_ptr->specie_vect.at(specie);
                const auto alias = react_ptr->coef_alias.find(specie_ptr);
                sqltext.append(value_with_alias(react_ptr->coef.at(specie_ptr).value(),
                    alias == react_ptr->coef_alias.end() ? std::string{} : alias->second.name,
                    alias == react_ptr->coef_alias.end() ? 1.0d : alias->second.sign));
                if (specie != react_ptr->specie_vect.size() - 1) {
                    sqltext.append(" ");
                }
            }
            sqltext.append("' WHERE REACTION_NAME = '" + react_ptr->name + "'; ");
            pop_and_add(p, 0, "exponents"); // removes coef state
            sqltext.append("UPDATE reactions SET 'EXPONENTS' = '");
            for (std::ptrdiff_t specie = 0; specie < react_ptr->specie_vect.size(); specie++) {
                specie_struct* specie_ptr = react_ptr->specie_vect.at(specie);
                const auto alias = react_ptr->exp_alias.find(specie_ptr);
                sqltext.append(value_with_alias(react_ptr->exp.at(specie_ptr).value(),
                    alias == react_ptr->exp_alias.end() ? std::string{} : alias->second.name,
                    alias == react_ptr->exp_alias.end() ? 1.0d : alias->second.sign));
                if (specie != react_ptr->specie_vect.size() - 1) {
                    sqltext.append(" ");
                }
            }
            sqltext.append("' WHERE REACTION_NAME = '" + react_ptr->name + "'; ");
            pop_report(p, 0); // removes exponents state
        }
        pop_and_add(p, 0, "edges"); // removes reaction state
        // Persist the experiment value without synthetic left padding, matching raw_profile.LEFT_EDGE.
        sqltext.append("UPDATE experiments SET 'LEFT_EDGE' = '" + double_to_string(run_ptr->left_edge.value() - p.exp_left_padding) + "' WHERE NAME = '" + run_ptr->name + "'; ");
        sqltext.append("UPDATE experiments SET 'WIDTH' = '" + double_to_string(run_ptr->width.value()) + "' WHERE NAME = '" + run_ptr->name + "'; ");
        pop_report(p, 0); // removes edges state
    }
    add_report(p, 0, "raw profile left_edge and width");
    for (int i = 0; i < p.experiments.size(); i++) {
        experiment_struct* exp_ptr = &p.experiments.at(i);
        if (exp_ptr->omit) continue;
        sqltext.append("UPDATE raw_profile SET 'LEFT_EDGE' = '" + double_to_string(exp_ptr->left_edge.value() - p.exp_left_padding) + "', 'WIDTH' = '" + double_to_string(exp_ptr->width.value()) + "' WHERE NAME = '" + exp_ptr->run->name + "' AND WT_PERCENT = '" + exp_ptr->second_name + "'; ");
    }
    pop_report(p, 0);
    pop_and_add(p, 1, "running sql"); // removes run state
    execute_sql(p, db, sqltext.c_str(), 0, errMsg.out());
    pop_report(p, 1);
    add_finishing_report(p, 3, "Done");
}

void 
normalize_profile(parameters_t& p)
{
    experiment_struct* exp_ptr;
    experiment_run_struct* run_ptr;
    add_report(p, 3, "Normalizing");
    for (int run = 0; run < p.experiment_runs.size(); run++) {
        experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
        add_report(p, 3, "(" + run_ptr->name + ":" + run_ptr->normalization_method + ")");
    }
    add_report(p, 0, std::to_string(0));
    double low_ref, high_ref, peak_ref;
    for (int j = 0; j < p.row_count; j++) {
        pop_and_add(p, 0, std::to_string(j));
        exp_ptr = &p.experiments.at(j);
        run_ptr = exp_ptr->run;
        int exp_size = exp_ptr->raw_experimental_profile.size();
        low_ref = 0.0d;
        high_ref = 0.0d;
        peak_ref = 0.0d;
        for (int i = run_ptr->low_ref_start + p.exp_left_padding; i < run_ptr->high_ref_end + p.exp_left_padding && i < exp_size; i++) {
            if (i <= run_ptr->low_ref_end + p.exp_left_padding) {
                low_ref += exp_ptr->raw_experimental_profile.at(i);
            } else if (i >= run_ptr->high_ref_start + p.exp_left_padding) {
                high_ref += exp_ptr->raw_experimental_profile.at(i);
            }
        }
        for (int i = p.exp_left_padding; i < run_ptr->high_ref_end + p.exp_left_padding && i < exp_size; i++) {
            peak_ref = std::max(exp_ptr->raw_experimental_profile.at(i), peak_ref);
        }
        if (run_ptr->low_ref_end < run_ptr->low_ref_start) {
            low_ref = 0;
        } else {
            low_ref = low_ref / (double)(run_ptr->low_ref_end - run_ptr->low_ref_start + 1);
        }
        
        high_ref = high_ref / (double)(run_ptr->high_ref_end - run_ptr->high_ref_start);
        for (int i = 0; i < exp_size; i++) {
            if (run_ptr->normalization_method == "ridge linear scaling") {
                if (i < run_ptr->low_ref_start + p.exp_left_padding) {
                    exp_ptr->raw_experimental_profile.at(i) = 0;
                } else if (i >= run_ptr->high_ref_end + p.exp_left_padding) {  
                    exp_ptr->raw_experimental_profile.at(i) = 1;
                } else {
                    exp_ptr->raw_experimental_profile.at(i) = (exp_ptr->raw_experimental_profile.at(i) - low_ref) / (high_ref - low_ref);
                }
            } else if (run_ptr->normalization_method == "peak linear scaling") {
                if (i < run_ptr->low_ref_start + p.exp_left_padding) {
                    exp_ptr->raw_experimental_profile.at(i) = 0;
                } else if (i >= exp_size - p.exp_right_padding) {
                    exp_ptr->raw_experimental_profile.at(i) = (high_ref - low_ref) / (peak_ref - low_ref);
                } else {
                    exp_ptr->raw_experimental_profile.at(i) = (exp_ptr->raw_experimental_profile.at(i) - low_ref) / (peak_ref - low_ref);
                }
            }
        }
        const auto left_padding = static_cast<std::size_t>(std::max(0, p.exp_left_padding));
        const auto right_padding = static_cast<std::size_t>(std::max(0, p.exp_right_padding));
        const auto first_sample = std::min(left_padding, exp_ptr->raw_experimental_profile.size());
        const auto last_sample = exp_ptr->raw_experimental_profile.size() >= right_padding
            ? exp_ptr->raw_experimental_profile.size() - right_padding
            : first_sample;
        if (last_sample > first_sample) {
            exp_ptr->normalized_experimental_profile.assign(
                exp_ptr->raw_experimental_profile.begin() + first_sample,
                exp_ptr->raw_experimental_profile.begin() + last_sample);
        } else {
            exp_ptr->normalized_experimental_profile.clear();
        }
    }
    pop_report(p, 0);
    add_finishing_report(p, 3, std::to_string(p.row_count) + " Done");
}

double
scattering_correction(parameters_t& p, double species_1, double species_2, double coef, experiment_run_struct* run_ptr)
{
    if (p.scatter_correction_type == "NS_ND") {
        double bead_wt = (species_1 + species_2 / coef) * unit_conversion(run_ptr, run_ptr->PS_beads, run_ptr->species.at(run_ptr->PS_beads).model_units, "wt%");
        if (run_ptr->species.at(run_ptr->PS_beads).diameter < 30.0e-9d) {
            return -0.238826108843563 * std::exp(-pow(0.0258645848310996 - bead_wt, 2.0d) / (2.0d * pow(0.00418992180425235, 2.0d))) + bead_wt * 1.40044738887211 + 1;
        };
        return -0.448191328804794 * std::exp(-pow(0.0711222856018783 - bead_wt, 2.0d) / (2.0d * pow(0.012623365952763, 2.0d))) + bead_wt * 2.14776259822044 + 1;
    } else {
        return 1.00d;
    }
}

solvable&
variable_location(const std::string& variable_name, experiment_run_struct* run_ptr) {
    if (variable_name == "left_edge") {return run_ptr->left_edge;}
    if (variable_name == "width") {return run_ptr->width;}
    if (run_ptr->alias_variables.contains(variable_name)) return run_ptr->alias_variables.at(variable_name);
    if (variable_name == "FITC_exp") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).exp.at(&run_ptr->species.at(specie_index(run_ptr, "FITC")));}
    if (variable_name == "bead_exp") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).exp.at(&run_ptr->species.at(specie_index(run_ptr, "PS_40nm", "PS_20nm")));}
    if (variable_name == "bound_bead_exp") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).exp.at(&run_ptr->species.at(specie_index(run_ptr, "40nm_Bound_Dye_1", "20nm_Bound_Dye_1")));}
    if (variable_name == "kon1") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).k[0];}
    if (variable_name == "kon2") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_2", "FITC_20nm_2")).k[0];}
    if (variable_name == "keq1") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).k[1];}
    if (variable_name == "keq2") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_2", "FITC_20nm_2")).k[1];}
    if (variable_name == "QE1") {return run_ptr->species.at(specie_index(run_ptr, "40nm_Bound_Dye_1", "20nm_Bound_Dye_1")).QE;}
    if (variable_name == "QE2") {return run_ptr->species.at(specie_index(run_ptr, "40nm_Bound_Dye_2", "20nm_Bound_Dye_2")).QE;}
    if (variable_name == "p1") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).coef.at(&run_ptr->species.at(specie_index(run_ptr, "FITC")));}
    if (variable_name == "p2") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_2", "FITC_20nm_2")).coef.at(&run_ptr->species.at(specie_index(run_ptr, "FITC")));}
    if (variable_name == "ND1") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_1", "FITC_20nm_1")).coef.at(&run_ptr->species.at(specie_index(run_ptr, "40nm_Bound_Dye_1", "20nm_Bound_Dye_1")));}
    if (variable_name == "NDD2") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_2", "FITC_20nm_2")).coef.at(&run_ptr->species.at(specie_index(run_ptr, "40nm_Bound_Dye_2", "20nm_Bound_Dye_2")));}
    if (variable_name == "ND2") {return run_ptr->reactions.at(reaction_index(run_ptr, "FITC_40nm_2", "FITC_20nm_2")).coef.at(&run_ptr->species.at(specie_index(run_ptr, "40nm_Bound_Dye_1", "20nm_Bound_Dye_1")));}
    else {throw std::runtime_error("Unknown variable '" + variable_name + "' called in variable_location");}
}

void
add_report(parameters_t& p, int debug_level, std::string text_to_add)
{
    std::lock_guard<std::recursive_mutex> lock(p.report_mutex);
    if (p.debug_level <= debug_level) {
        p.state.emplace_back(debug_level, text_to_add + " ");
        publish_event(p, {tsensor_workflow::event_kind::message, p.active_operation,
                          std::move(text_to_add), debug_level});
    }
}

void
pop_report(parameters_t& p, int debug_level)
{
    std::lock_guard<std::recursive_mutex> lock(p.report_mutex);
    if (!p.state.empty() && p.state.back().debug_level <= debug_level) {
        p.state.pop_back();
    }
}

void
pop_and_add(parameters_t& p, int debug_level, std::string text_to_add)
{
    std::lock_guard<std::recursive_mutex> lock(p.report_mutex);
    pop_report(p, debug_level);
    add_report(p, debug_level, text_to_add);
}


void
pop_finishing_report(parameters_t& p, int debug_level, std::string text_to_add)
{
    std::lock_guard<std::recursive_mutex> lock(p.report_mutex);
    if (p.debug_level <= debug_level) {
        pop_and_add(p, debug_level, text_to_add);
        p.state.clear();
    }
}

void
add_finishing_report(parameters_t& p, int debug_level, std::string text_to_add)
{
    std::lock_guard<std::recursive_mutex> lock(p.report_mutex);
    if (p.debug_level <= debug_level) {
        add_report(p, debug_level, text_to_add);
        p.state.clear();
    }
}

void removeSpaces(std::string &str)
{
    const auto first = str.find_first_not_of(' ');
    if (first == std::string::npos) {
        str.clear();
        return;
    }
    str.erase(str.find_last_not_of(' ') + 1);
    str.erase(0, first);
}

ptrdiff_t
specie_index(experiment_run_struct* run_ptr, std::string specie_to_find, std::string second_specie_name)
{
    std::string species_list;
    for (ptrdiff_t i = 0; i < run_ptr->species.size(); i++) {
        if (specie_to_find == run_ptr->species.at(i).name || second_specie_name == run_ptr->species.at(i).name) {
            return i; 
        }
        species_list.append(run_ptr->species.at(i).name + ",");
    }
    throw std::runtime_error("Invalid specie name as string in specie_index (" + specie_to_find + " or " + second_specie_name + ") not in (" + species_list + ")");
}

ptrdiff_t
reaction_index(experiment_run_struct* run_ptr, std::string reaction_to_find, std::string second_reaction_name)
{
    for (ptrdiff_t i = 0; i < run_ptr->reactions.size(); i++) {
        if (reaction_to_find == run_ptr->reactions.at(i).name || second_reaction_name == run_ptr->reactions.at(i).name) {
            return i; 
        }
    }
    throw std::runtime_error("Invalid reaction name as string in reaction_index");
}

double
unit_conversion(experiment_run_struct* run_ptr, int specie, std::string current_unit, std::string desired_unit)
{
    if (current_unit == desired_unit) {
        return 1.0d;
    }
    std::vector<std::string> unit = {current_unit, desired_unit};
    double unit_value[2];
    specie_struct* specie_ptr = &run_ptr->species.at(specie);

    if (specie_ptr->type == 2) {
        assert (specie_ptr->diameter != 0.0d); // TODO make runtime error
        assert (specie_ptr->particle_density != 0.0d); // TODO make runtime error
    } else if (specie_ptr->type == 1) {
        assert (specie_ptr->molecular_weight != 0.0d); // TODO make runtime error
    } else {
        throw std::runtime_error("SPECIE:" + specie_ptr->name + " unit_conversion for type " + std::to_string(specie_ptr->type) + " undefined");
    }

    for (int i = 0; i < 2; i++) {
        if (unit[i] == "umol") {
            if (specie_ptr->type == 2) {
                unit_value[i] = 1.00d / 1.0e+6 * (6.022e+23) * ( 4.0d / 3.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 3.0d)) * specie_ptr->particle_density / 1000.0d;
            } else if (specie_ptr->type == 1) {
                unit_value[i] = 1.00d / 1.0e+6 * specie_ptr->molecular_weight;
            } else {
                throw std::runtime_error("SPECIE:" + specie_ptr->name + " unit_conversion for " + unit[i] + " undefined (check specie_type)");
            }
        } else if (unit[i] == "um2/ul") {
            assert (specie_ptr->type == 2); // TODO make runtime error
            //            um2/ul * (1.0e+6ul/l) / (10^12um2/m2) / (SA m2/bead) * (Vol m3/bead) * (dens g/ml) * (10^6 ml/m3)
            unit_value[i] = 1.0d * 1.0e+6 / 1.0e+12 / (4.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 2.0d)) * ( 4.0d / 3.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 3.0d)) * specie_ptr->particle_density * 1.0e+6; 
        } else if (unit[i] == "nm2/ul") {
            assert (specie_ptr->type == 2); // TODO make runtime error
            //            nm2/ul * (1.0e+6ul/l) / (10^18nm2/m2) / (SA m2/bead) * (Vol m3/bead) * (dens g/ml) * (10^6 ml/m3)
            unit_value[i] = 1.0d * 1.0e+6 / 1.0e+18 / (4.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 2.0d)) * ( 4.0d / 3.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 3.0d)) * specie_ptr->particle_density * 1.0e+6; 
        } else if (unit[i] == "mm2/nl") {
            assert (specie_ptr->type == 2); // TODO make runtime error
            //            mm2/nl * (1.0e+9nl/l) / (10^6mm2/m2) / (SA m2/bead) * (Vol m3/bead) * (dens g/ml) * (10^6 ml/m3)
            unit_value[i] = 1.0d * 1.0e+9 / 1.0e+6 / (4.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 2.0d)) * ( 4.0d / 3.0d * M_PI * pow(specie_ptr->diameter / 2.0d, 3.0d)) * specie_ptr->particle_density * 1.0e+6;
        } else if (unit[i] == "mg/ml") {
            unit_value[i] = 1.00d;
        } else if (unit[i] == "wt%") {
            assert (run_ptr->solution_density != 0.0d); // TODO make runtime error
            unit_value[i] = 1 / 100.0d * run_ptr->solution_density * 1000.0d;
        } else if (unit[i] == "g/ml") {
            unit_value[i] = 1000.00d;
        } else {
            throw std::runtime_error("SPECIE:" + specie_ptr->name + " unit_conversion for " + unit[i] + " undefined");
        }
    }
    assert (unit_value[0] != 0.0d); // TODO make runtime error
    assert (unit_value[1] != 0.0d); // TODO make runtime error
    return unit_value[0] / unit_value[1];
}

std::string
double_to_string(double arg)
{
    std::string sign;
    if (std::abs(arg) < 1e-307) {
        arg = 0;
    }
   
    if (arg == 0) {
        return "0";
    } else if (arg < 0) {
        arg *= -1;
        sign = "-";
    }
    int exponent = (int)std::floor(std::log10(arg));
    arg *= std::pow(10 , -exponent);
    return sign + std::to_string(arg) + "e" + std::to_string(exponent);
}

double
lin_interpolate(double x, double x1, double y1, double x2, double y2)
{
    if (x == x1) {
        return y1;
    } else if (x == x2) {
        return y2;
    } else if ((x > x1 && x > x2) || (x < x1 && x < x2)) {
        throw std::runtime_error("lin_interpolate out of range :" + double_to_string(x) + " not between " + double_to_string(x1) + " and " + double_to_string(x2));
    }
    if (x2 == x1) {
        return (y1 + y2) / 2;
    }
    return (y2 * std::abs(x - x1) + y1 * std::abs(x2 - x)) / std::abs(x2 - x1);
}


void
set_inlet_conc(parameters_t& p, experiment_struct* exp_ptr, double* solution)
{
    auto* run_ptr = exp_ptr->run;
    const int species_count = run_ptr->number_of_species;
    const double total_flowrate = run_ptr->total_flowrate;
    if (p.X <= 0 || species_count <= 0 || !std::isfinite(total_flowrate) || total_flowrate <= 0.0) {
        throw std::invalid_argument("Invalid grid or total flowrate while setting inlet concentrations for " +
            run_ptr->name + ":" + exp_ptr->second_name);
    }

    double flowrate_accounted_for = 0.0;
    for (const auto& entrance : exp_ptr->entrances) {
        const int first_cell = static_cast<int>(p.X * flowrate_accounted_for / total_flowrate);
        flowrate_accounted_for += entrance.ENTRANCE_FLOWRATE;
        const int end_cell = static_cast<int>(p.X * flowrate_accounted_for / total_flowrate);

        for (int specie = 0; specie < species_count; specie++) {
            auto* specie_ptr = &run_ptr->species.at(specie);
            double concentration = 0.0;
            const auto concentration_entry = entrance.CONC.find(specie_ptr);
            if (concentration_entry != entrance.CONC.end()) {
                concentration = concentration_entry->second;
            }
            const double model_concentration = concentration * unit_conversion(
                run_ptr, specie, specie_ptr->input_units, specie_ptr->model_units);
            for (int x = first_cell; x < end_cell; x++) {
                const int grid_index = specie * p.X + x;
                solution[grid_index] = model_concentration;
                solution[grid_index + species_count * p.X] = model_concentration;
            }
        }
    }
}

void
model(parameters_t& p, const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, int row)
{

    //add_report(p, 0, "row:" + std::to_string(row) + "_start");
    auto* exp_ptr = &p.experiments.at(row);
    auto* run_ptr = exp_ptr->run;
    const double left_edge = exp_ptr->left_edge.value();
    const double width = exp_ptr->width.value();
    if (!std::isfinite(left_edge) || left_edge < 0.0 || !std::isfinite(width) || width <= 0.0 ||
        exp_ptr->window_size <= 0 || p.X < 2 || run_ptr->number_of_species <= 0 ||
        run_ptr->number_of_reactions <= 0) {
        throw std::invalid_argument("Invalid model dimensions or profile geometry for " +
            run_ptr->name + ":" + exp_ptr->second_name);
    }
    if (static_cast<double>(exp_ptr->window_size) < std::ceil(width)) {
        throw std::invalid_argument("Model profile width exceeds its sample window for " +
            run_ptr->name + ":" + exp_ptr->second_name);
    }
    exp_ptr->scale_factor = run_ptr->W * 1.0e6 / width;
    if (!std::isfinite(exp_ptr->scale_factor) || exp_ptr->scale_factor <= 0.0) {
        throw std::invalid_argument("Invalid model scale factor for " +
            run_ptr->name + ":" + exp_ptr->second_name);
    }
    const int X = p.X;
    const int species_count = run_ptr->number_of_species;
    const int reaction_count = run_ptr->number_of_reactions;
    const std::size_t grid_size = static_cast<std::size_t>(species_count) * X;
    const std::size_t reaction_species_size = static_cast<std::size_t>(reaction_count) * species_count;
    std::vector<double> kon(reaction_count);
    std::vector<double> reverse_kon(reaction_count);
    std::vector<double> reaction_rate(reaction_count);
    std::vector<double> coef(reaction_species_size);
    std::vector<double> reaction_orders(reaction_species_size);
    std::vector<double> specie_rate(species_count);
    std::vector<double> available(species_count);
    std::vector<double> r(species_count);
    std::vector<double> solution_arena(3 * grid_size);
    double* E = solution_arena.data();
    double* solution = E + grid_size;
    double* old_solution = solution + grid_size;
    double* solution_ptr     = solution;
    double* old_solution_ptr = old_solution;
    //The 3 Concentration Arrays.
    auto& species_out = exp_ptr->species_out;
    //pop_and_add(p, 0, "row:" + std::to_string(row) + "_coef");

    
    for (int react = 0; react < reaction_count; react++) {
        reaction_struct* react_ptr = &run_ptr->reactions.at(react);
        kon[react] = run_ptr->dt * react_ptr->k[0].value();
        const double reverse_rate_constant = react_ptr->k[1].value();
        if (!std::isfinite(reverse_rate_constant) || reverse_rate_constant <= 0.0) {
            throw std::invalid_argument("Reaction reverse rate constant must be positive for " +
                run_ptr->name + ":" + react_ptr->name);
        }
        reverse_kon[react] = kon[react] / reverse_rate_constant;
        const std::size_t reaction_offset = static_cast<std::size_t>(react) * species_count;
        for (int specie = 0; specie < species_count; specie++) {
            auto* specie_ptr = &run_ptr->species.at(specie);
            coef[reaction_offset + specie] = react_ptr->coef.at(specie_ptr).value();
            reaction_orders[reaction_offset + specie] = react_ptr->exp.at(specie_ptr).value();
        }
    }
    assert (coef[run_ptr->FITC] == -coef[run_ptr->Bound_Dye_1]);
    //pop_and_add(p, 0, "row:" + std::to_string(row) + "_preinlet");

    set_inlet_conc(p, exp_ptr, solution);

    for (int specie = 0; specie < species_count; specie++) {
        r[specie] = run_ptr->species.at(specie).r;
    }
    //pop_and_add(p, 0, "row:" + std::to_string(row) + "_preZloop");

    for (int z = 0; z < p.Z; z++) {
        check_cancellation(p);
        for (int i = 0; i < X; i++) { 
            for (int specie = 0; specie < species_count; specie++) {
                solution_ptr = &solution[specie * X];
                if (i == 0) { // left edge
                    available[specie] =                                     (1.00d - r[specie]) * solution_ptr[i]          + r[specie] * solution_ptr[i + 1];
                } else if (i == X - 1) { // right edge
                    available[specie] = r[specie] * solution_ptr[i - 1]   + (1.00d - r[specie]) * solution_ptr[i];
                } else { // mid points
                    available[specie] = r[specie] * solution_ptr[i - 1]   + (2.00d - 2.00d * r[specie]) * solution_ptr[i]  + r[specie] * solution_ptr[i + 1];
                }
                available[specie] = std::max(available[specie], 0.0d);
            }
            if (p.disable_reactions) {
                for (int reaction = 0; reaction < reaction_count; reaction++) {
                    reaction_rate[reaction] = 0.0d;
                }
            } else {
                // specie reaction rates
                for (int reaction = 0; reaction < reaction_count; reaction++) {
                    const std::size_t specie_stagger = static_cast<std::size_t>(reaction) * species_count;
                    reaction_rate[reaction] = kon[reaction];
                    double reverse = reverse_kon[reaction];
                    if (p.disable_reverse_reactions) {
                        reverse = 0.0d;
                    }
                    for (int specie = 0; specie < species_count; specie++) {
                        solution_ptr = &solution[specie * X];
                        if (coef[specie_stagger + specie] < 0.0d) {
                            reaction_rate[reaction] *= std::pow(solution_ptr[i], reaction_orders[specie_stagger + specie]); // Forward Reaction
                        } else if (coef[specie_stagger + specie] > 0.0d) {
                            reverse *= std::pow(solution_ptr[i], reaction_orders[specie_stagger + specie]); // Reverse Reaction
                        }
                    }
                    reaction_rate[reaction] -= reverse;
                }
                // limiting reagents
                for (int specie = 0; specie < species_count; specie++) {
                    specie_rate[specie] = 0.0d;
                    for (int reaction = 0; reaction < reaction_count; reaction++) {
                        specie_rate[specie] += coef[specie + reaction * species_count] * reaction_rate[reaction];
                    }
                    int while_loop_iter = 0;
                    while(available[specie] + specie_rate[specie] < 0) { // check for limiting reagent.
                        double total_positive_magnitude = 0.00d;
                        specie_rate[specie] = 0.0d;
                        for (int reaction = 0; reaction < reaction_count; reaction++) {
                            if (reaction_rate[reaction] * coef[specie + reaction * species_count] > 0.0d) {
                                total_positive_magnitude += reaction_rate[reaction] * coef[specie + reaction * species_count];
                            }
                        }
                        for (int reaction = 0; reaction < reaction_count; reaction++) {
                            double reaction_specie_rate = reaction_rate[reaction] * coef[specie + reaction * species_count];
                            if ((reaction_specie_rate >= 0.0d)) {
                                // skip this condition
                            } else if (total_positive_magnitude == 0.00d) {
                                reaction_rate[reaction] = 0.0d;
                            } else if (std::abs(reaction_specie_rate) > total_positive_magnitude && while_loop_iter == 0) {
                                reaction_rate[reaction] = std::copysign(total_positive_magnitude / coef[specie + reaction * species_count], reaction_rate[reaction]);
                            } else {
                                reaction_rate[reaction] = reaction_rate[reaction] * 0.99;
                            }
                            specie_rate[specie] += coef[specie + reaction * species_count] * reaction_rate[reaction];
                            assert (!std::isnan(specie_rate[specie]));
                            assert (!std::isinf(specie_rate[specie]));
                        }
                        assert (while_loop_iter != 1000);
                        while_loop_iter++;
                    }
                }
                // Final rates after limiting 
                for (int specie = 0; specie < species_count; specie++) {
                    specie_rate[specie] = 0.0d;
                    for (int reaction = 0; reaction < reaction_count; reaction++) {
                        specie_rate[specie] += coef[specie + reaction * species_count] * reaction_rate[reaction];
                    }
                }
                //assert (reaction_rate[run_ptr->FITC_Bead_1] >= 0.0d);
                //assert (specie_rate[run_ptr->FITC] < 0.0d);
            }
            for (int specie = 0; specie < species_count; specie++) {
                E[specie * X + i] = available[specie] + specie_rate[specie];
            }
        }
        for (int specie = 0; specie < species_count; specie++) {
            double oneplus_r = 1.0d / (1.0d + r[specie]);
            double specie_total = 0;
            solution_ptr = &solution[specie * X];
            old_solution_ptr = &old_solution[specie * X];
            for (int i = 0; i < X; i++) {
                old_solution_ptr[i] = solution_ptr[i];
                specie_total += old_solution_ptr[i];
            }
            double error = specie_total;
            int while_loop_iter = 0;
            while ((error > p.time_step_convergence * specie_total && while_loop_iter < X * 2) || while_loop_iter < 3) {
                check_cancellation(p);
                error = 0.0d;
                std::swap(old_solution_ptr, solution_ptr);
                for (int i = 0; i < X; i++) {
                    if (i == 0) { // left edge
                        solution_ptr[i] = (E[specie * X + i] + r[specie] * old_solution_ptr[i + 1]) * oneplus_r;
                    } else if (i == X - 1) { // right edge
                        solution_ptr[i] = (E[specie * X + i] + r[specie] * old_solution_ptr[i - 1]) * oneplus_r;
                    } else { // mid points
                        solution_ptr[i] = (E[specie * X + i] + r[specie] * (old_solution_ptr[i - 1] + old_solution_ptr[i + 1])) * oneplus_r / 2;
                    }
                    error += std::abs(solution_ptr[i] - old_solution_ptr[i]);
                }
                while_loop_iter += 1;
            }
            for (int i = 0; i < X; i++) {
                assert (&solution_arena[run_ptr->number_of_species * X] <= &solution_ptr[i]);
                assert (&solution_arena[3 * run_ptr->number_of_species * X] > &solution_ptr[i]);
                solution[specie * X + i] = solution_ptr[i]; // solution_ptr may be pointing to data in old_solution region
            }
            assert (while_loop_iter < X * 1);
            for (int i = 0; i < X; i++) {
                assert (&solution_arena[run_ptr->number_of_species * X] <= &solution_ptr[i]);
                assert (&solution_arena[3 * run_ptr->number_of_species * X] > &solution_ptr[i]);
                assert (&solution_arena[run_ptr->number_of_species * X] <= &old_solution_ptr[i]);
                assert (&solution_arena[3 * run_ptr->number_of_species * X] > &old_solution_ptr[i]);
                assert (!std::isnan(solution_ptr[i]));
                assert (!std::isinf(solution_ptr[i]));

                solution_ptr[i] = std::max(solution_ptr[i], 0.00d);
                if (z == (p.Z - 1)) {
                    species_out[specie][i] = solution_ptr[i];
                }
            }             
        }
    }
    //pop_and_add(p, 0, "row:" + std::to_string(row) + "_postZloop");
    //add_report(p, 0, "row:" + std::to_string(row) + "_postZloop");

    
    //add_report(p, 0, "row:" + std::to_string(row) + "_postunshift");
    for (int i = 0; i < exp_ptr->window_size; i++) {
        exp_ptr->model_profile.at(i) = 0.0d;
        //assert ((species_out[0][i] + species_out[2][i]) / run_ptr->dye_conc < 2.0); // model profile output is to high
    }
    double split = 99.5d;
    int bottom_point, top_point;
    for (int specie = 0; specie < run_ptr->number_of_species; specie++) {
        specie_struct* specie_ptr = &run_ptr->species.at(specie);
        const double QE = specie_ptr->QE.value();
        for (int i = 0; i < exp_ptr->window_size; i++) {
            if (i < width) {
                split = (double)(i) * ((double)X) / width;
                bottom_point = static_cast<int>(floor(split));
                assert (bottom_point <= split);
                top_point = static_cast<int>(ceil(split));
                assert (top_point >= split);
                if (top_point == bottom_point || bottom_point == (X - 1)) {
                    exp_ptr->model_profile.at(i) += species_out[specie][bottom_point] * QE;
                } else {
                    exp_ptr->model_profile.at(i) += lin_interpolate(split, (double)bottom_point, species_out[specie][bottom_point] * QE, (double)top_point, species_out[specie][top_point] * QE);
                }
            } else {
                exp_ptr->model_profile.at(i) += species_out[specie][X - 1] * QE;
            }
            assert (!std::isnan(exp_ptr->model_profile.at(i)));
            assert (!std::isinf(exp_ptr->model_profile.at(i)));
        }
    }
    pop_and_add(p, 0, "row:" + std::to_string(row) + "_residuals");

    for (int i = 0; i < exp_ptr->window_size; i++) {
        double scatter;
        if (i < width) {
            split = (double)(i) * ((double)X) / width;
            bottom_point = static_cast<int>(floor(split));
            top_point = static_cast<int>(ceil(split));
            if (top_point == bottom_point || bottom_point == (X - 1)) {
                scatter = scattering_correction(p, species_out[run_ptr->PS_beads][bottom_point], species_out[run_ptr->Bound_Dye_1][bottom_point], coef[0], run_ptr);
            } else {
                scatter = scattering_correction(p, species_out[run_ptr->PS_beads][bottom_point] * ((double)top_point - split) + species_out[run_ptr->PS_beads][top_point] * (split - (double)bottom_point), species_out[run_ptr->Bound_Dye_1][bottom_point] * ((double)top_point - split) + species_out[run_ptr->Bound_Dye_1][top_point] * (split - (double)bottom_point), coef[0], run_ptr);
            }
        } else {
            scatter = scattering_correction(p, species_out[run_ptr->PS_beads][X - 1], species_out[run_ptr->Bound_Dye_1][X - 1], coef[0], run_ptr);
        }

        assert (i + exp_ptr->window_start < p.total_window_size); // residuals input is out of bounds
        // MinLM minimizes sum(fi^2); supply signed discrepancies, not their squares.
        const double residual = (exp_ptr->model_profile.at(i) * scatter / run_ptr->dye_conc - exp_ptr->experimental_profile.at(i)) * !exp_ptr->omit;
        residuals[i + exp_ptr->window_start] = residual;

        // Preserve squared-error semantics in exported profiles and database rows.
        exp_ptr->error.at(i) = residual * residual;

        if (i == 0 || i == exp_ptr->window_size - 1) {
            exp_ptr->experimental_derivative.at(i) = 0.0;
            exp_ptr->numeric_derivative.at(i) = 0.0;
        } else {
            exp_ptr->experimental_derivative.at(i)  = (exp_ptr->experimental_profile.at(i + 1)  - exp_ptr->experimental_profile.at(i - 1))  / (2 * exp_ptr->scale_factor);
            exp_ptr->numeric_derivative.at(i)       = (exp_ptr->model_profile.at(i + 1)         - exp_ptr->model_profile.at(i - 1))         / (2 * exp_ptr->scale_factor) / run_ptr->dye_conc;
        }
        
        assert (!std::isnan(residuals[i + exp_ptr->window_start]));
        assert (!std::isinf(residuals[i + exp_ptr->window_start]));
    }
    
    {
        pop_and_add(p, 0, "row:" + std::to_string(row) + "_DAintegrals");
        double exp_d = std::numeric_limits<double>::lowest();
        double model_d = std::numeric_limits<double>::lowest();
        double exp_a = std::numeric_limits<double>::max();
        double model_a = std::numeric_limits<double>::max();
        double exp_integral = 0;
        double model_integral = 0;
        double analytic_exp_integral = 0;
        double analytic_model_integral = 0;
        int exp_window_count = 0;
        int model_window_count = 0;
       
        for (int x = 0; x < exp_ptr->window_size; x++) {
            if (exp_ptr->channel_position.at(x) > (double)run_ptr->W * 1e6 * 0.2d && exp_ptr->channel_position.at(x) < (double)run_ptr->W * 1e6 * 0.8d ) { // TODO: make 0.2 & 0.8 parameters. channel width comes from experiment geometry
                exp_d                   = std::max(exp_d, exp_ptr->experimental_derivative.at(x));
                model_d                 = std::max(model_d, exp_ptr->numeric_derivative.at(x));
                exp_a                   = std::min(exp_a, exp_ptr->experimental_derivative.at(x));
                model_a                 = std::min(model_a, exp_ptr->numeric_derivative.at(x));
                exp_window_count++;
                model_window_count++;
                if (exp_ptr->zero_row_ptr != NULL) {
                    auto& zero = exp_ptr->zero_row_ptr;
                    const double x_split = exp_ptr->channel_position.at(x) / zero->scale_factor + left_edge;
                    int bot_zero = 0;
                    int top_zero = zero->channel_position.size() - 1;
                    for (int j = 0; j < zero->channel_position.size(); j++) {
                        if (zero->channel_position.at(j) > exp_ptr->channel_position.at(x)) {
                            bot_zero = std::max(0, j - 1);
                            top_zero = j;
                            break;
                        } else if (zero->channel_position.at(j) == exp_ptr->channel_position.at(x)) {
                            bot_zero = j;
                            top_zero = j;
                            break;
                        }
                    }
                    assert(exp_ptr->channel_position.at(x) >= zero->channel_position.at(bot_zero));
                    assert(exp_ptr->channel_position.at(x) <= zero->channel_position.at(top_zero));
                    const double zero_exp = lin_interpolate(exp_ptr->channel_position.at(x), zero->channel_position.at(bot_zero), zero->experimental_profile.at(bot_zero), zero->channel_position.at(top_zero), zero->experimental_profile.at(top_zero));
                    const double zero_model = lin_interpolate(exp_ptr->channel_position.at(x), zero->channel_position.at(bot_zero), zero->model_profile.at(bot_zero), zero->channel_position.at(top_zero), zero->model_profile.at(top_zero));
                    assert(!isnan(zero_exp));
                    assert(!isnan(zero_model));
                    const double exp_dif                    = std::abs(exp_ptr->experimental_profile.at(x) - zero_exp);
                    exp_ptr->experimental_difference.at(x)  = exp_dif;
                    exp_integral                            += exp_dif;
                    assert(!isnan(exp_integral));
                    const double num_dif                    = std::abs(exp_ptr->model_profile.at(x) - zero_model);
                    exp_ptr->numeric_difference.at(x)       = num_dif;
                    model_integral                          += num_dif;
                    assert(!isnan(model_integral));
                }
                const double analytic_zero              = 1.0d / 2.0d * std::erfc(((250.0d - exp_ptr->channel_position.at(x)) * 1.0e-6d) / std::sqrt(4.0d * (double)p.Z * run_ptr->dt * 0.00000000049d));
                exp_ptr->analytical_zero.at(x)          = analytic_zero * run_ptr->dye_conc;
                const double analytic_exp_dif           = std::abs(exp_ptr->experimental_profile.at(x) - analytic_zero);
                analytic_exp_integral                   += analytic_exp_dif;
                const double analytic_num_dif           = std::abs(exp_ptr->model_profile.at(x) - analytic_zero * run_ptr->dye_conc);
                analytic_model_integral                 += analytic_num_dif;
            }
        }

        exp_ptr->exp_DA = exp_d - exp_a;
        exp_ptr->model_DA = model_d - model_a;
        assert(exp_window_count > 0);
        exp_integral /= (double)p.X * (0.8d - 0.2d) * exp_window_count;
        analytic_exp_integral /= (double)p.X * (0.8d - 0.2d) * exp_window_count;
        //assert(!isnan(exp_integral));
        assert( model_window_count > 0);
        model_integral /= (double)p.X * (0.8d - 0.2d) * model_window_count * run_ptr->dye_conc;
        analytic_model_integral /= (double)p.X * (0.8d - 0.2d) * model_window_count * run_ptr->dye_conc;
        //assert(!isnan(model_integral));
        exp_ptr->exp_integral = exp_integral;
        exp_ptr->model_integral = model_integral;
        exp_ptr->analytic_exp_integral = analytic_exp_integral;
        exp_ptr->analytic_model_integral = analytic_model_integral;
    }
    pop_report(p, 0);
}



void
alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr)
{
    if (!ptr) { throw std::invalid_argument("Missing solver state"); }
    auto& p = *static_cast<parameters_t*>(ptr);
    check_cancellation(p);
    {
        add_report(p, 1, "setting control parmeters");
        for (int i = 0; i < p.solvables.size(); i++) {
            auto& s = p.solvables.at(i);
            s.value() = control_parameters[i];
        }
        for (auto& run : p.experiment_runs) synchronize_run_aliases(p, run);
        synchronize_profile_entrance_concentrations(p, false);
        for (int run = 0; run < p.experiment_runs.size(); run++) {
            experiment_run_struct* run_ptr = &p.experiment_runs.at(run);
            for (int i = 0; i < run_ptr->reactions.size(); i++) {
                const auto& reaction = run_ptr->reactions.at(i);
                if (reaction.coef_alias.empty() &&
                    (reaction.name == "FITC_40nm_1" || reaction.name == "FITC_20nm_1")) {
                    variable_location("ND1", run_ptr).value() = abs(variable_location("p1", run_ptr).value());
                }
                if (reaction.coef_alias.empty() &&
                    (reaction.name == "FITC_40nm_2" || reaction.name == "FITC_20nm_2")) {
                    variable_location("NDD2", run_ptr).value() = abs(variable_location("p2", run_ptr).value()) + abs(variable_location("ND2", run_ptr).value());
                }
            }
        }

        for (int i = 0; i < p.experiments.size(); i++) {
            experiment_struct* exp_ptr = &p.experiments.at(i);
            const double left = exp_ptr->left_edge.value();
            const double width = exp_ptr->width.value();
            // x is integral and x < width: the last sampled x is ceil(width)-1.
            // Check before converting to int or touching any profile buffers.
            const double last_sample = std::ceil(left + std::min(
                static_cast<double>(exp_ptr->window_size), std::ceil(width)) - 1.0);
            if (!std::isfinite(left) || left < 0 || !std::isfinite(width) || width <= 0 ||
                exp_ptr->window_size <= 0 || exp_ptr->raw_experimental_profile.empty() ||
                last_sample >= static_cast<double>(exp_ptr->raw_experimental_profile.size())) {
                throw std::invalid_argument("Profile sampling exceeds available data for " +
                    exp_ptr->run->name + ":" + exp_ptr->second_name +
                    " (left_edge=" + double_to_string(left) + ", width=" + double_to_string(width) +
                    ", last requested index=" + double_to_string(last_sample) +
                    ", available samples=" + std::to_string(exp_ptr->raw_experimental_profile.size()) +
                    "). Adjust left_edge/width fit bounds or provide a larger raw profile window.");
            }
            // Sampling needs the current width's scale on the first evaluation,
            // before model() runs (a model-only run has no later evaluation).
            exp_ptr->scale_factor = exp_ptr->run->W * 1.0e6 / width;
            int last_x = static_cast<int>(ceil(exp_ptr->left_edge.value()));
            for (int x = 0; x < exp_ptr->window_size; x++) {
                int x_unshift = static_cast<int>(ceil(x + exp_ptr->left_edge.value()));
                exp_ptr->channel_position.at(x) = ((double)x_unshift - exp_ptr->left_edge.value()) * exp_ptr->scale_factor;
                if (x == 0) {
                    exp_ptr->channel_position.at(x) = 0;  // distorts the left edge to prevent a crash later. TODO: FIX
                }
                if (x < exp_ptr->width.value()) {
                    exp_ptr->experimental_profile.at(x) = exp_ptr->raw_experimental_profile.at(x_unshift);
                    last_x = x_unshift;
                } else {
                    exp_ptr->experimental_profile.at(x) = exp_ptr->raw_experimental_profile.at(last_x);
                }
                assert (!std::isnan(exp_ptr->experimental_profile.at(x)));
                assert (!std::isinf(exp_ptr->experimental_profile.at(x)));
            }
        }
        
        pop_report(p, 1);
        unsigned long const hardware_threads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::exception_ptr> failures(p.row_count);
        // jthread also joins if creating a later worker throws.
        std::vector<std::jthread> threads(p.row_count);
        // D-A metrics compare each profile with its run's 0% profile. Finish
        // those reference profiles before starting dependent profiles so the
        // first (and only, for run_solver=false) evaluation is deterministic.
        std::vector<int> dependent_rows;
        dependent_rows.reserve(p.row_count);
        for (int row = 0; row < p.row_count; ++row) {
            if (p.experiments.at(row).second_name == "0.0") {
                model(p, control_parameters, residuals, row);
            } else {
                dependent_rows.push_back(row);
            }
        }
        int start_row = 0;
        std::size_t row = 0;
        while (row < dependent_rows.size()) {
            check_cancellation(p);
            start_row = static_cast<int>(row);
            for (int i = 0; i < hardware_threads && row < dependent_rows.size(); i++) {
                const int model_row = dependent_rows[row];
                threads[row] = std::jthread([&p, &control_parameters, &residuals, &failures, model_row] {
                    try {
                        model(p, control_parameters, residuals, model_row);
                    } catch (...) {
                        failures[model_row] = std::current_exception();
                    }
                });
                row++;
            }
            for (int i = 0; i < hardware_threads && static_cast<std::size_t>(start_row + i) < row; i++) {
                threads[start_row + i].join();
            }
            for (std::size_t i = static_cast<std::size_t>(start_row); i < row; ++i) {
                if (failures[dependent_rows[i]]) { std::rethrow_exception(failures[dependent_rows[i]]); }
            }
        }
    }
    p.iterations = p.iterations + 1;
    tsensor_workflow::progress_event event{tsensor_workflow::event_kind::evaluation, p.active_operation,
                      "Residual evaluation completed", 3, p.iterations};
    if (p.capture_live_profiles) {
        for (const auto& experiment : p.experiments) {
            auto* run = experiment.run;
            if (!run || run->FITC < 0 || run->FITC >= static_cast<std::ptrdiff_t>(run->species.size())) continue;
            const auto& fitc = run->species.at(run->FITC);
            const double numeric_unit = unit_conversion(run, run->FITC, fitc.model_units, fitc.input_units);
            const auto sample_count = std::min({experiment.channel_position.size(),
                experiment.experimental_profile.size(), experiment.model_profile.size()});
            tsensor_workflow::live_profile_snapshot snapshot;
            snapshot.experiment_name = run->name;
            snapshot.weight_percent = experiment.second_name;
            snapshot.x.reserve(sample_count);
            snapshot.experimental.reserve(sample_count);
            snapshot.numeric.reserve(sample_count);
            for (std::size_t sample = 0; sample < sample_count; ++sample) {
                snapshot.x.push_back(experiment.channel_position[sample]);
                snapshot.experimental.push_back(experiment.experimental_profile[sample] * run->dye_conc_mgml);
                snapshot.numeric.push_back(experiment.model_profile[sample] * numeric_unit);
            }
            event.live_profiles.push_back(std::move(snapshot));
        }
    }
    publish_event(p, std::move(event));
}
