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

if(NAVIER_BUILD_GUI)
    # Only the desktop application and its runtime files have install rules.
    set(CPACK_GENERATOR ZIP)
    set(CPACK_PACKAGE_NAME NavierGui)
    set(CPACK_PACKAGE_VERSION "0.1.0")
    set(CPACK_PACKAGE_FILE_NAME "NavierGui-0.1.0-windows-x64-${CMAKE_BUILD_TYPE}")
    set(CPACK_PACKAGE_DIRECTORY "${CMAKE_BINARY_DIR}/packages")
    include(CPack)
endif()
