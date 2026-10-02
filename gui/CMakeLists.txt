include(${PROJECT_SOURCE_DIR}/cmake/QtKit.cmake)

if(NOT WIN32 OR NOT MINGW OR NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "The GUI currently requires the pinned Windows x86_64 MinGW kit. See docs/qt-setup.md.")
endif()
if(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
   NOT CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL NAVIER_QT_GNU_VERSION OR
   NOT CMAKE_C_COMPILER_VERSION VERSION_EQUAL NAVIER_QT_GNU_VERSION)
    message(FATAL_ERROR "The GUI requires GNU ${NAVIER_QT_GNU_VERSION} for C and C++. Use the qt-mingw presets.")
endif()
execute_process(COMMAND "${CMAKE_CXX_COMPILER}" -dumpmachine
    OUTPUT_VARIABLE qt_compiler_target OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
if(NOT qt_compiler_target STREQUAL "${NAVIER_QT_ARCHITECTURE}-w64-mingw32")
    message(FATAL_ERROR "Unexpected GUI compiler target: ${qt_compiler_target}")
endif()

find_package(Qt6 ${NAVIER_QT_VERSION} EXACT REQUIRED COMPONENTS ${NAVIER_QT_COMPONENTS})
add_executable(NavierGui WIN32 main.cpp)
target_link_libraries(NavierGui PRIVATE tsensor Qt6::Widgets)
set_target_properties(NavierGui PROPERTIES AUTOMOC ON)

if(BUILD_TESTING)
    add_test(NAME Application.GuiStartup COMMAND NavierGui --smoke-test)
    set_tests_properties(Application.GuiStartup PROPERTIES
        TIMEOUT 15 ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
endif()
