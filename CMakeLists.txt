cmake_minimum_required(VERSION 3.30.2)

project(Navier LANGUAGES C CXX)
add_executable(${PROJECT_NAME} Navier.cpp)

add_subdirectory(src)
target_link_libraries(${PROJECT_NAME} PUBLIC alglib)