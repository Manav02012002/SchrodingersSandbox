cmake_minimum_required(VERSION 3.20)

project(schrodingers_sandbox VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SBOX_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(SBOX_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

execute_process(
    COMMAND git describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT GIT_VERSION)
    set(GIT_VERSION "0.1.0-dev")
endif()

message(STATUS "Version: ${GIT_VERSION}")

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/generated)
configure_file(
    ${CMAKE_SOURCE_DIR}/src/core/version.h.in
    ${CMAKE_BINARY_DIR}/generated/version.h
    @ONLY
)

if(SBOX_ENABLE_ASAN)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
endif()

if(SBOX_ENABLE_UBSAN)
    add_compile_options(-fsanitize=undefined)
    add_link_options(-fsanitize=undefined)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    add_compile_options(-g -O2)
endif()

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(external/glfw)
add_subdirectory(external/nfd)

find_package(OpenGL REQUIRED)
find_package(Eigen3 3.4 QUIET)
if (NOT Eigen3_FOUND)
    find_package(Eigen3 REQUIRED)
endif()

add_library(glad STATIC external/glad/src/gl.c)
target_include_directories(glad PUBLIC external/glad/include)

add_library(imgui STATIC
    external/imgui/imgui.cpp
    external/imgui/imgui_demo.cpp
    external/imgui/imgui_draw.cpp
    external/imgui/imgui_tables.cpp
    external/imgui/imgui_widgets.cpp
    external/imgui/backends/imgui_impl_glfw.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC external/imgui external/imgui/backends)
target_link_libraries(imgui PUBLIC glfw glad OpenGL::GL)
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_CUSTOM)
target_compile_options(imgui PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-include glad/gl.h>)

add_library(implot STATIC
    external/implot/implot.cpp
    external/implot/implot_items.cpp
    external/implot/implot_demo.cpp
)
target_include_directories(implot PUBLIC external/implot)
target_link_libraries(implot PUBLIC imgui)

add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE external/nlohmann)

add_executable(schrodingers_sandbox
    src/cli.cpp
    src/main.cpp
    src/backend/backend_manager.cpp
    src/analysis/crystal_field.cpp
    src/analysis/nci.cpp
    src/analysis/orbital_composition.cpp
    src/backend/python_env.cpp
    src/chem/ligand_library.cpp
    src/core/crash_handler.cpp
    src/core/special_functions.cpp
    src/core/hydrogen.cpp
    src/core/logging.cpp
    src/core/slater.cpp
    src/core/spline.cpp
    src/core/paths.cpp
    src/core/settings.cpp
    src/core/update_checker.cpp
    src/chem/coordination.cpp
    src/core/elements.cpp
    src/core/gaussian_eval.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/symmetry.cpp
    src/core/valence.cpp
    src/core/zmatrix.cpp
    src/core/basis_set.cpp
    src/core/molden_parser.cpp
    src/io/cube_io.cpp
    src/io/fchk_io.cpp
    src/io/pdb_io.cpp
    src/io/project_io.cpp
    src/io/sdf_io.cpp
    src/io/trajectory_io.cpp
    src/io/xyz_io.cpp
    src/editor/command.cpp
    src/editor/picking.cpp
    src/editor/select_mode.cpp
    src/editor/draw_mode.cpp
    src/editor/erase_mode.cpp
    src/editor/measure_mode.cpp
    src/editor/fragment_mode.cpp
    src/editor/fragment_library.cpp
    src/renderer/window.cpp
    src/renderer/shader.cpp
    src/renderer/camera.cpp
    src/renderer/basis_texture.cpp
    src/renderer/esp_surface.cpp
    src/renderer/gbuffer.cpp
    src/renderer/lod_renderer.cpp
    src/renderer/mol_renderer.cpp
    src/renderer/post_process.cpp
    src/renderer/shadow_map.cpp
    src/renderer/ssao.cpp
    src/renderer/screenshot.cpp
    src/renderer/volume_texture.cpp
    src/renderer/video_export.cpp
    src/ui/app.cpp
    src/ui/about_dialog.cpp
    src/ui/annotations.cpp
    src/ui/annotation_editor.cpp
    src/ui/bond_order_panel.cpp
    src/ui/ui_utils.cpp
    src/ui/panels.cpp
    src/ui/status_bar.cpp
    src/ui/energy_diagram.cpp
    src/ui/mo_diagram.cpp
    src/ui/mo_diagram_plot.cpp
    src/ui/orbital_composition_panel.cpp
    src/ui/optimization_panel.cpp
    src/ui/pes_panel.cpp
    src/ui/dos_panel.cpp
    src/ui/export_dialog.cpp
    src/ui/computation_panel.cpp
    src/ui/complex_builder.cpp
    src/ui/constraint_editor.cpp
    src/ui/charge_overlay.cpp
    src/ui/crystal_field_panel.cpp
    src/ui/d_orbital_viewer.cpp
    src/ui/context_menu.cpp
    src/ui/file_dialog.cpp
    src/ui/molecule_info.cpp
    src/ui/nci_panel.cpp
    src/ui/population_panel.cpp
    src/ui/esp_controls.cpp
    src/ui/ir_spectrum_panel.cpp
    src/ui/property_dashboard.cpp
    src/ui/reaction_path_panel.cpp
    src/ui/results_panel.cpp
    src/ui/settings_panel.cpp
    src/ui/setup_wizard.cpp
    src/ui/solvent_panel.cpp
    src/ui/spectrochemical_panel.cpp
    src/ui/symmetry_overlay.cpp
    src/ui/plot_utils.cpp
    src/ui/periodic_table.cpp
    src/ui/orbital_browser.cpp
    src/ui/properties_panel.cpp
    src/ui/editor_toolbar.cpp
)
target_include_directories(schrodingers_sandbox PRIVATE src ${CMAKE_BINARY_DIR}/generated external/stb)
target_link_libraries(schrodingers_sandbox PRIVATE imgui implot glad glfw OpenGL::GL Eigen3::Eigen nlohmann_json nfd)
target_compile_options(schrodingers_sandbox PRIVATE -Wall -Wextra -Wpedantic)

if(UNIX AND NOT APPLE)
    target_link_options(schrodingers_sandbox PRIVATE -rdynamic)
endif()

if(APPLE)
    set(MACOSX_BUNDLE_ICON_FILE AppIcon)
    set(APP_ICON_MACOSX ${CMAKE_SOURCE_DIR}/resources/AppIcon.icns)
    if(EXISTS ${APP_ICON_MACOSX})
        set_source_files_properties(${APP_ICON_MACOSX} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
        target_sources(schrodingers_sandbox PRIVATE ${APP_ICON_MACOSX})
    endif()

    set_target_properties(schrodingers_sandbox PROPERTIES
        MACOSX_BUNDLE TRUE
        MACOSX_BUNDLE_BUNDLE_NAME "Schrödinger's Sandbox"
        MACOSX_BUNDLE_BUNDLE_VERSION "${GIT_VERSION}"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "${GIT_VERSION}"
        MACOSX_BUNDLE_GUI_IDENTIFIER "com.schrodingerssandbox.app"
        MACOSX_BUNDLE_ICON_FILE "AppIcon"
        MACOSX_BUNDLE_COPYRIGHT "Copyright 2025 Manav Rawal"
    )
endif()

add_custom_command(TARGET schrodingers_sandbox POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data/shaders
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data/shaders ${CMAKE_BINARY_DIR}/data/shaders
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data/scripts
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/scripts/neb_driver.py ${CMAKE_BINARY_DIR}/data/scripts/neb_driver.py
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/scripts/pes_scan_driver.py ${CMAKE_BINARY_DIR}/data/scripts/pes_scan_driver.py
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/scripts/pyscf_driver.py ${CMAKE_BINARY_DIR}/data/scripts/pyscf_driver.py
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/scripts/xtb_driver.py ${CMAKE_BINARY_DIR}/data/scripts/xtb_driver.py
)

if(APPLE)
    add_custom_command(TARGET schrodingers_sandbox POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_BUNDLE_DIR:schrodingers_sandbox>/Contents/Resources/data"
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data/shaders "$<TARGET_BUNDLE_DIR:schrodingers_sandbox>/Contents/Resources/data/shaders"
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data/scripts "$<TARGET_BUNDLE_DIR:schrodingers_sandbox>/Contents/Resources/data/scripts"
    )
endif()

install(TARGETS schrodingers_sandbox
    RUNTIME DESTINATION bin
    BUNDLE DESTINATION .
)

install(DIRECTORY data/shaders DESTINATION share/schrodingers_sandbox/data)
install(DIRECTORY data/scripts DESTINATION share/schrodingers_sandbox/data)

install(FILES external/glfw/LICENSE.md
    DESTINATION share/schrodingers_sandbox/licenses
    RENAME glfw-LICENSE.md
)
install(FILES external/imgui/LICENSE.txt
    DESTINATION share/schrodingers_sandbox/licenses
    RENAME imgui-LICENSE.txt
)

enable_testing()
add_subdirectory(external/googletest)

add_executable(test_constants tests/test_constants.cpp)
target_include_directories(test_constants PRIVATE src)
target_link_libraries(test_constants PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_constants PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_special_functions
    tests/test_special_functions.cpp
    src/core/special_functions.cpp
)
target_include_directories(test_special_functions PRIVATE src)
target_link_libraries(test_special_functions PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_special_functions PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_hydrogen
    tests/test_hydrogen.cpp
    src/core/hydrogen.cpp
    src/core/special_functions.cpp
)
target_include_directories(test_hydrogen PRIVATE src)
target_link_libraries(test_hydrogen PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_hydrogen PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_slater
    tests/test_slater.cpp
    src/core/slater.cpp
    src/core/elements.cpp
)
target_include_directories(test_slater PRIVATE src)
target_link_libraries(test_slater PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_slater PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_molden_parser
    tests/test_molden_parser.cpp
    src/core/basis_set.cpp
    src/core/molden_parser.cpp
)
target_include_directories(test_molden_parser PRIVATE src)
target_link_libraries(test_molden_parser PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_molden_parser PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_molecular_system
    tests/test_molecular_system.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
)
target_include_directories(test_molecular_system PRIVATE src)
target_link_libraries(test_molecular_system PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_molecular_system PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_bond_perception
    tests/test_bond_perception.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
)
target_include_directories(test_bond_perception PRIVATE src)
target_link_libraries(test_bond_perception PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_bond_perception PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_xyz_io
    tests/test_xyz_io.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
    src/io/xyz_io.cpp
)
target_include_directories(test_xyz_io PRIVATE src)
target_link_libraries(test_xyz_io PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_xyz_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_sdf_io
    tests/test_sdf_io.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
    src/io/sdf_io.cpp
)
target_include_directories(test_sdf_io PRIVATE src)
target_link_libraries(test_sdf_io PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_sdf_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_cube_io
    tests/test_cube_io.cpp
    src/io/cube_io.cpp
)
target_include_directories(test_cube_io PRIVATE src)
target_link_libraries(test_cube_io PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_cube_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_fchk_io
    tests/test_fchk_io.cpp
    src/io/fchk_io.cpp
)
target_include_directories(test_fchk_io PRIVATE src)
target_link_libraries(test_fchk_io PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_fchk_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_pdb_io
    tests/test_pdb_io.cpp
    src/io/pdb_io.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
)
target_include_directories(test_pdb_io PRIVATE src)
target_link_libraries(test_pdb_io PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_pdb_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_trajectory_io
    tests/test_trajectory_io.cpp
    src/io/trajectory_io.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
)
target_include_directories(test_trajectory_io PRIVATE src)
target_link_libraries(test_trajectory_io PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_trajectory_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_zmatrix
    tests/test_zmatrix.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
    src/core/zmatrix.cpp
)
target_include_directories(test_zmatrix PRIVATE src)
target_link_libraries(test_zmatrix PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_zmatrix PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_symmetry
    tests/test_symmetry.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
    src/core/symmetry.cpp
)
target_include_directories(test_symmetry PRIVATE src)
target_link_libraries(test_symmetry PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_symmetry PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_project_io
    tests/test_project_io.cpp
    src/core/covalent_radii.cpp
    src/core/molecular_system.cpp
    src/core/elements.cpp
    src/io/project_io.cpp
)
target_include_directories(test_project_io PRIVATE src)
target_link_libraries(test_project_io PRIVATE GTest::gtest_main Eigen3::Eigen nlohmann_json)
target_compile_options(test_project_io PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_gaussian_eval
    tests/test_gaussian_eval.cpp
    src/core/basis_set.cpp
    src/core/gaussian_eval.cpp
)
target_include_directories(test_gaussian_eval PRIVATE src)
target_link_libraries(test_gaussian_eval PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_gaussian_eval PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_gpu_crossval
    tests/test_gpu_crossval.cpp
    src/core/basis_set.cpp
    src/core/gaussian_eval.cpp
    src/renderer/basis_texture.cpp
)
target_include_directories(test_gpu_crossval PRIVATE src external/glad/include)
target_link_libraries(test_gpu_crossval PRIVATE GTest::gtest_main Eigen3::Eigen glad OpenGL::GL)
target_compile_options(test_gpu_crossval PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_python_env
    tests/test_python_env.cpp
    src/backend/python_env.cpp
)
target_include_directories(test_python_env PRIVATE src)
target_link_libraries(test_python_env PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_python_env PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_command_stack
    tests/test_command_stack.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
    src/core/valence.cpp
    src/editor/command.cpp
)
target_include_directories(test_command_stack PRIVATE src)
target_link_libraries(test_command_stack PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_command_stack PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_picking
    tests/test_picking.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
    src/editor/picking.cpp
)
target_include_directories(test_picking PRIVATE src external/imgui)
target_link_libraries(test_picking PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_picking PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_valence
    tests/test_valence.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
    src/core/valence.cpp
)
target_include_directories(test_valence PRIVATE src)
target_link_libraries(test_valence PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_valence PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_fragment_library
    tests/test_fragment_library.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
    src/editor/fragment_library.cpp
)
target_include_directories(test_fragment_library PRIVATE src)
target_link_libraries(test_fragment_library PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_fragment_library PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_orbital_composition
    tests/test_orbital_composition.cpp
    src/analysis/orbital_composition.cpp
    src/core/basis_set.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/slater.cpp
    src/core/molecular_system.cpp
)
target_include_directories(test_orbital_composition PRIVATE src)
target_link_libraries(test_orbital_composition PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_orbital_composition PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_coordination
    tests/test_coordination.cpp
    src/chem/coordination.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
)
target_include_directories(test_coordination PRIVATE src)
target_link_libraries(test_coordination PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_coordination PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_ligand_library
    tests/test_ligand_library.cpp
    src/chem/coordination.cpp
    src/chem/ligand_library.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
)
target_include_directories(test_ligand_library PRIVATE src)
target_link_libraries(test_ligand_library PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_ligand_library PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_crystal_field
    tests/test_crystal_field.cpp
    src/analysis/crystal_field.cpp
    src/core/basis_set.cpp
    src/core/covalent_radii.cpp
    src/core/elements.cpp
    src/core/molecular_system.cpp
)
target_include_directories(test_crystal_field PRIVATE src)
target_link_libraries(test_crystal_field PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_crystal_field PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_spline
    tests/test_spline.cpp
    src/core/spline.cpp
)
target_include_directories(test_spline PRIVATE src)
target_link_libraries(test_spline PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_spline PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_settings
    tests/test_settings.cpp
    src/core/logging.cpp
    src/core/settings.cpp
)
target_include_directories(test_settings PRIVATE src)
target_link_libraries(test_settings PRIVATE GTest::gtest_main Eigen3::Eigen nlohmann_json)
target_compile_options(test_settings PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_cli
    tests/test_cli.cpp
    src/cli.cpp
)
target_include_directories(test_cli PRIVATE src ${CMAKE_BINARY_DIR}/generated)
target_link_libraries(test_cli PRIVATE GTest::gtest_main Eigen3::Eigen)
target_compile_options(test_cli PRIVATE -Wall -Wextra -Wpedantic)

add_executable(test_update_checker
    tests/test_update_checker.cpp
    src/core/update_checker.cpp
    src/backend/python_env.cpp
)
target_include_directories(test_update_checker PRIVATE src ${CMAKE_BINARY_DIR}/generated)
target_link_libraries(test_update_checker PRIVATE GTest::gtest_main Eigen3::Eigen nlohmann_json)
target_compile_options(test_update_checker PRIVATE -Wall -Wextra -Wpedantic)

if(DEFINED ENV{SBOX_RUN_INTEGRATION_TESTS} AND "$ENV{SBOX_RUN_INTEGRATION_TESTS}" STREQUAL "1")
    add_executable(test_pyscf_integration
        tests/test_pyscf_integration.cpp
        src/backend/backend_manager.cpp
        src/backend/python_env.cpp
        src/core/basis_set.cpp
        src/core/covalent_radii.cpp
        src/core/elements.cpp
        src/core/molecular_system.cpp
        src/core/molden_parser.cpp
        src/io/cube_io.cpp
    )
    target_include_directories(test_pyscf_integration PRIVATE src)
    target_link_libraries(test_pyscf_integration PRIVATE GTest::gtest_main Eigen3::Eigen nlohmann_json)
    target_compile_options(test_pyscf_integration PRIVATE -Wall -Wextra -Wpedantic)

    add_executable(test_pes_integration
        tests/test_pes_integration.cpp
        src/backend/backend_manager.cpp
        src/backend/python_env.cpp
        src/core/basis_set.cpp
        src/core/covalent_radii.cpp
        src/core/elements.cpp
        src/core/molecular_system.cpp
        src/core/molden_parser.cpp
        src/io/cube_io.cpp
    )
    target_include_directories(test_pes_integration PRIVATE src)
    target_link_libraries(test_pes_integration PRIVATE GTest::gtest_main Eigen3::Eigen nlohmann_json)
    target_compile_options(test_pes_integration PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(GoogleTest)
gtest_discover_tests(test_constants)

add_test(NAME test_special_functions COMMAND test_special_functions)
add_test(NAME test_hydrogen COMMAND test_hydrogen)
add_test(NAME test_slater COMMAND test_slater)
add_test(NAME test_molden_parser COMMAND test_molden_parser)
add_test(NAME test_molecular_system COMMAND test_molecular_system)
add_test(NAME test_bond_perception COMMAND test_bond_perception)
add_test(NAME test_xyz_io COMMAND test_xyz_io)
add_test(NAME test_sdf_io COMMAND test_sdf_io)
add_test(NAME test_cube_io COMMAND test_cube_io)
add_test(NAME test_fchk_io COMMAND test_fchk_io)
add_test(NAME test_pdb_io COMMAND test_pdb_io)
add_test(NAME test_trajectory_io COMMAND test_trajectory_io)
add_test(NAME test_zmatrix COMMAND test_zmatrix)
add_test(NAME test_symmetry COMMAND test_symmetry)
add_test(NAME test_project_io COMMAND test_project_io)
add_test(NAME test_gaussian_eval COMMAND test_gaussian_eval)
add_test(NAME test_gpu_crossval COMMAND test_gpu_crossval)
add_test(NAME test_python_env COMMAND test_python_env)
add_test(NAME test_command_stack COMMAND test_command_stack)
add_test(NAME test_picking COMMAND test_picking)
add_test(NAME test_valence COMMAND test_valence)
add_test(NAME test_fragment_library COMMAND test_fragment_library)
add_test(NAME test_orbital_composition COMMAND test_orbital_composition)
add_test(NAME test_coordination COMMAND test_coordination)
add_test(NAME test_ligand_library COMMAND test_ligand_library)
add_test(NAME test_crystal_field COMMAND test_crystal_field)
add_test(NAME test_spline COMMAND test_spline)
add_test(NAME test_settings COMMAND test_settings)
add_test(NAME test_cli COMMAND test_cli)
add_test(NAME test_update_checker COMMAND test_update_checker)

if(TARGET test_pyscf_integration)
    add_test(NAME test_pyscf_integration COMMAND test_pyscf_integration)
endif()

if(TARGET test_pes_integration)
    add_test(NAME test_pes_integration COMMAND test_pes_integration)
endif()

set(SBOX_LICENSE_FILE "${CMAKE_SOURCE_DIR}/README.md")
if(EXISTS "${CMAKE_SOURCE_DIR}/LICENSE")
    set(SBOX_LICENSE_FILE "${CMAKE_SOURCE_DIR}/LICENSE")
endif()

if(APPLE)
    file(WRITE "${CMAKE_BINARY_DIR}/Info.plist" "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n<plist version=\"1.0\">\n<dict>\n    <key>CFBundleDisplayName</key>\n    <string>Schrödinger's Sandbox</string>\n    <key>CFBundleExecutable</key>\n    <string>schrodingers_sandbox</string>\n    <key>CFBundleIdentifier</key>\n    <string>com.schrodingerssandbox.app</string>\n    <key>CFBundleName</key>\n    <string>Schrödinger's Sandbox</string>\n    <key>CFBundlePackageType</key>\n    <string>APPL</string>\n    <key>CFBundleShortVersionString</key>\n    <string>${GIT_VERSION}</string>\n    <key>CFBundleVersion</key>\n    <string>${GIT_VERSION}</string>\n</dict>\n</plist>\n")
endif()

set(CPACK_PACKAGE_NAME "SchrodingersSandbox")
set(CPACK_PACKAGE_VERSION "${GIT_VERSION}")
set(CPACK_PACKAGE_VENDOR "Manav Rawal")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Interactive quantum chemistry visualisation")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/your-username/schrodingerssandbox")
set(CPACK_RESOURCE_FILE_LICENSE "${SBOX_LICENSE_FILE}")

if(APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_DMG_VOLUME_NAME "Schrödinger's Sandbox")
    set(CPACK_BUNDLE_NAME "Schrödinger's Sandbox")
    set(CPACK_BUNDLE_PLIST "${CMAKE_BINARY_DIR}/Info.plist")
elseif(WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_DISPLAY_NAME "Schrödinger's Sandbox")
    set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
    set(CPACK_NSIS_MENU_LINKS "bin/schrodingers_sandbox.exe" "Schrödinger's Sandbox")
    if(EXISTS "${CMAKE_SOURCE_DIR}/resources/icon.ico")
        set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/resources/icon.ico")
    endif()
else()
    set(CPACK_GENERATOR "TGZ;DEB")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Manav Rawal")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libeigen3-dev, libgl-dev")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
endif()

include(CPack)
