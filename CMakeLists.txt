cmake_minimum_required(VERSION 4.2)
project(Navier LANGUAGES C CXX)

include(CTest)
add_executable(${PROJECT_NAME} tsensor.cpp)

add_subdirectory(src)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC tsensor)
