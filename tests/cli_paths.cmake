# Only help, argument errors, and database-open errors are exercised here.
# Never load model inputs or launch a solver against an existing database.
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef suffix)
set(run_dir "${TEST_ROOT}-${suffix}")
file(MAKE_DIRECTORY "${run_dir}/launch elsewhere")

execute_process(COMMAND "${NAVIER}" --help
    WORKING_DIRECTORY "${run_dir}/launch elsewhere"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Usage: Navier")
    message(FATAL_ERROR "Help failed: ${result}: ${output} ${error}")
endif()

execute_process(COMMAND "${NAVIER}" --database --output-dir results
    WORKING_DIRECTORY "${run_dir}/launch elsewhere"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 1 OR NOT error MATCHES "Missing path after --database")
    message(FATAL_ERROR "Missing argument was not rejected: ${result}: ${output} ${error}")
endif()

set(missing_db "${run_dir}/missing database.db")
execute_process(COMMAND "${NAVIER}" --database "${missing_db}"
    --output-dir "${run_dir}/results with spaces"
    WORKING_DIRECTORY "${run_dir}/launch elsewhere"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 1 OR NOT error MATCHES "Cannot open database")
    message(FATAL_ERROR "Missing database was not rejected: ${result}: ${output} ${error}")
endif()
if(EXISTS "${missing_db}" OR EXISTS "${run_dir}/results with spaces")
    message(FATAL_ERROR "Failed database open created a database or output directory")
endif()
# Only remove the unique disposable directory created above.
file(REMOVE_RECURSE "${run_dir}")
