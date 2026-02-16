cmake_minimum_required(VERSION 3.20)

project(schrodingers_sandbox VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(external/glfw)

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

add_executable(schrodingers_sandbox
    src/main.cpp
    src/core/special_functions.cpp
    src/core/hydrogen.cpp
    src/core/slater.cpp
    src/core/elements.cpp
    src/renderer/window.cpp
    src/renderer/shader.cpp
    src/renderer/camera.cpp
    src/ui/app.cpp
    src/ui/panels.cpp
    src/ui/periodic_table.cpp
    src/ui/orbital_browser.cpp
    src/ui/properties_panel.cpp
)
target_include_directories(schrodingers_sandbox PRIVATE src)
target_link_libraries(schrodingers_sandbox PRIVATE imgui glad glfw OpenGL::GL Eigen3::Eigen)
target_compile_options(schrodingers_sandbox PRIVATE -Wall -Wextra -Wpedantic)

add_custom_command(TARGET schrodingers_sandbox POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data/shaders
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data/shaders ${CMAKE_BINARY_DIR}/data/shaders
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

include(GoogleTest)
gtest_discover_tests(test_constants)

add_test(NAME test_special_functions COMMAND test_special_functions)
add_test(NAME test_hydrogen COMMAND test_hydrogen)
add_test(NAME test_slater COMMAND test_slater)
