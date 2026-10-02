cmake_minimum_required(VERSION 4.2)
project(Navier LANGUAGES C CXX)

option(NAVIER_BUILD_GUI "Build the Qt Widgets desktop application" OFF)

include(CTest)
add_executable(${PROJECT_NAME} tsensor.cpp)

add_subdirectory(src)
if(NAVIER_BUILD_GUI)
    add_subdirectory(gui)
endif()
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC tsensor)
