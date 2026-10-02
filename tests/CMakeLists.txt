include(FetchContent)

# GoogleTest is a build/test dependency, not part of the desktop distribution.
set(INSTALL_GTEST OFF CACHE BOOL "Do not install bundled test libraries" FORCE)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)

add_executable(tsensor_test tsensor_test.cpp application_paths_test.cpp session_test.cpp feedback_test.cpp background_runner_test.cpp workflow_test.cpp)

target_link_libraries(tsensor_test
 PRIVATE
    gtest_main
    tsensor
)

target_compile_definitions(tsensor_test PRIVATE
    TSENSOR_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
include(GoogleTest)
gtest_discover_tests(tsensor_test DISCOVERY_MODE PRE_TEST
    PROPERTIES TIMEOUT 30)

add_test(NAME Application.CommandLinePaths
    COMMAND ${CMAKE_COMMAND} -DNAVIER=$<TARGET_FILE:Navier>
            -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/cli-path-test
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cli_paths.cmake)
set_tests_properties(Application.CommandLinePaths PROPERTIES TIMEOUT 30)

add_executable(workflow_cli_fixture workflow_cli_fixture.cpp)
target_link_libraries(workflow_cli_fixture PRIVATE sqlite3)
add_test(NAME Application.CommandLineWorkflow
    COMMAND ${CMAKE_COMMAND} -DNAVIER=$<TARGET_FILE:Navier>
            -DFIXTURE_TOOL=$<TARGET_FILE:workflow_cli_fixture>
            -DFIXTURE_SQL=${CMAKE_CURRENT_SOURCE_DIR}/fixtures/workflow.sql
            -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/cli-workflow-test
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cli_workflow.cmake)
set_tests_properties(Application.CommandLineWorkflow PROPERTIES TIMEOUT 60)

if(NAVIER_BUILD_GUI)
    add_executable(navier_gui_test gui_test.cpp)
    target_link_libraries(navier_gui_test PRIVATE gtest navier_gui)
    target_compile_definitions(navier_gui_test PRIVATE TSENSOR_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
    gtest_discover_tests(navier_gui_test DISCOVERY_MODE PRE_TEST
        PROPERTIES TIMEOUT 30 ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
endif()
