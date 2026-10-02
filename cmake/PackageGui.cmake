# Standard bin/plugins layout lets Qt generate a relocatable qt.conf.
include(GNUInstallDirs)
install(TARGETS NavierGui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
qt_generate_deploy_app_script(TARGET NavierGui OUTPUT_SCRIPT navier_deploy_script
    NO_COMPILER_RUNTIME NO_TRANSLATIONS
    EXCLUDE_PLUGIN_TYPES generic networkinformation tls)
install(SCRIPT ${navier_deploy_script})

# Copy the runtimes from the compiler used to build this application, never PATH.
get_filename_component(navier_compiler_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
foreach(runtime libgcc_s_seh-1.dll libstdc++-6.dll libwinpthread-1.dll)
    install(FILES "${navier_compiler_bin}/${runtime}" DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()
install(DIRECTORY "${navier_compiler_bin}/../licenses/gcc"
    "${navier_compiler_bin}/../licenses/mingw-w64"
    "${navier_compiler_bin}/../licenses/winpthreads" DESTINATION notices/compiler)
install(FILES "${PROJECT_SOURCE_DIR}/src/alglib-cpp/gpl2.txt"
    "${PROJECT_SOURCE_DIR}/src/alglib-cpp/gpl3.txt" DESTINATION notices/alglib)
file(GLOB navier_eigen_notices "${PROJECT_SOURCE_DIR}/src/eigen-3.4.0/COPYING.*")
install(FILES ${navier_eigen_notices} DESTINATION notices/eigen)
get_target_property(navier_qt_qmake Qt6::qmake IMPORTED_LOCATION)
get_filename_component(navier_qt_bin "${navier_qt_qmake}" DIRECTORY)
install(DIRECTORY "${navier_qt_bin}/../sbom/" DESTINATION notices/qt-sbom)
install(FILES "${PROJECT_SOURCE_DIR}/docs/portable-package.md" DESTINATION . RENAME README.md)
