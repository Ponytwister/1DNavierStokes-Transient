cmake_minimum_required(VERSION 4.2)
project(Navier LANGUAGES C CXX)

enable_testing()
add_executable(${PROJECT_NAME} tsensor.cpp src/tsensor/sqlite_interface.cpp)

add_subdirectory(src)
add_subdirectory(tests)

target_link_libraries(${PROJECT_NAME} PUBLIC tsensor)