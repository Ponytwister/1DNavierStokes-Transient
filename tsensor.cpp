#include <workflow.h>
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

    tsensor_workflow::load_inputs(db);

    try
    {
        tsensor_workflow::run();
        tsensor_workflow::export_results();
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
        return 1;
    }
    
    
    sqlite3_close(db);
    std::cout << "Database closed." << std::endl;
    return 0;
}