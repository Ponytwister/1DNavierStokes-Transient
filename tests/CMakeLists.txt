include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)

add_executable(tsensor_test tsensor_test.cpp application_paths_test.cpp session_test.cpp)

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
