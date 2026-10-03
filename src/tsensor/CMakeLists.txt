add_library(tsensor STATIC channel_dimensions.cpp sqlite_interface.cpp workflow.cpp application_options.cpp feedback.cpp background_runner.cpp)
target_link_libraries(tsensor PUBLIC alglib)
target_link_libraries(tsensor PUBLIC sqlite3)

target_include_directories(tsensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
