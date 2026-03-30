cmake_minimum_required(VERSION 3.20)

project(schrodingers_sandbox VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
    src/main.cpp
    src/backend/backend_manager.cpp
    src/analysis/orbital_composition.cpp
    src/backend/python_env.cpp
    src/core/special_functions.cpp
    src/core/hydrogen.cpp
    src/core/slater.cpp
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
    src/io/project_io.cpp
    src/io/sdf_io.cpp
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
    src/renderer/mol_renderer.cpp
    src/renderer/volume_texture.cpp
    src/ui/app.cpp
    src/ui/bond_order_panel.cpp
    src/ui/ui_utils.cpp
    src/ui/panels.cpp
    src/ui/status_bar.cpp
    src/ui/energy_diagram.cpp
    src/ui/mo_diagram.cpp
    src/ui/mo_diagram_plot.cpp
    src/ui/orbital_composition_panel.cpp
    src/ui/dos_panel.cpp
    src/ui/computation_panel.cpp
    src/ui/charge_overlay.cpp
    src/ui/context_menu.cpp
    src/ui/file_dialog.cpp
    src/ui/molecule_info.cpp
    src/ui/population_panel.cpp
    src/ui/esp_controls.cpp
    src/ui/ir_spectrum_panel.cpp
    src/ui/property_dashboard.cpp
    src/ui/results_panel.cpp
    src/ui/setup_wizard.cpp
    src/ui/plot_utils.cpp
    src/ui/periodic_table.cpp
    src/ui/orbital_browser.cpp
    src/ui/properties_panel.cpp
    src/ui/editor_toolbar.cpp
)
target_include_directories(schrodingers_sandbox PRIVATE src)
target_link_libraries(schrodingers_sandbox PRIVATE imgui implot glad glfw OpenGL::GL Eigen3::Eigen nlohmann_json nfd)
target_compile_options(schrodingers_sandbox PRIVATE -Wall -Wextra -Wpedantic)

add_custom_command(TARGET schrodingers_sandbox POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data/shaders
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data/shaders ${CMAKE_BINARY_DIR}/data/shaders
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data/scripts
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/scripts/pyscf_driver.py ${CMAKE_BINARY_DIR}/data/scripts/pyscf_driver.py
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/data/scripts/xtb_driver.py ${CMAKE_BINARY_DIR}/data/scripts/xtb_driver.py
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

if(TARGET test_pyscf_integration)
    add_test(NAME test_pyscf_integration COMMAND test_pyscf_integration)
endif()
