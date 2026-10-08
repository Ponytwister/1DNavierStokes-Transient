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

# Keep the portable ZIP free of this machine-local shortcut. The local desktop
# deployment gets its own relative shortcut and refreshes the ignored root link.
install(CODE "
set(navier_shortcut_tool \"$<TARGET_FILE:navier_shortcut>\")
set(navier_source_dir \"${PROJECT_SOURCE_DIR}\")
set(navier_install_root \"\${CMAKE_INSTALL_PREFIX}\")
file(REAL_PATH \"\${navier_install_root}\" navier_install_real)
file(REAL_PATH \"\${navier_source_dir}/out/desktop\" navier_desktop_real)
string(TOLOWER \"\${navier_install_real}\" navier_install_lower)
string(TOLOWER \"\${navier_desktop_real}\" navier_desktop_lower)
if(navier_install_lower STREQUAL navier_desktop_lower)
    execute_process(
        COMMAND \"\${navier_shortcut_tool}\" \"\${navier_install_root}/Navier.lnk\"
            \"\${navier_install_root}/bin/NavierGui.exe\" \"bin/NavierGui.exe\"
        WORKING_DIRECTORY \"\${navier_install_root}\"
        RESULT_VARIABLE navier_deployed_shortcut_result)
    if(NOT navier_deployed_shortcut_result EQUAL 0)
        message(FATAL_ERROR \"Failed to create the deployed Navier shortcut\")
    endif()
    execute_process(
        COMMAND \"\${navier_shortcut_tool}\" \"\${navier_source_dir}/Navier.lnk\"
            \"\${navier_install_root}/bin/NavierGui.exe\" \"out/desktop/bin/NavierGui.exe\"
        WORKING_DIRECTORY \"\${navier_source_dir}\"
        RESULT_VARIABLE navier_root_shortcut_result)
    if(NOT navier_root_shortcut_result EQUAL 0)
        message(FATAL_ERROR \"Failed to refresh the root Navier shortcut\")
    endif()
endif()
")
