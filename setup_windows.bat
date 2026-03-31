@echo off
setlocal EnableExtensions EnableDelayedExpansion

echo === Schrödinger's Sandbox - Dependency Setup ===
echo.

cd /d "%~dp0"

set "PS=powershell -NoProfile -ExecutionPolicy Bypass -Command"

call :check_tool cmake "cmake not found. Install CMake and ensure it is on PATH."
call :check_tool git "git not found. Install Git for Windows."
call :check_tool python "python not found. Install Python 3 and ensure it is on PATH."

echo.
echo Setting up external dependencies...

if not exist "external\glfw\CMakeLists.txt" (
    echo   Cloning GLFW 3.4...
    if exist "external\glfw" rmdir /s /q "external\glfw"
    git clone --depth 1 --branch 3.4 https://github.com/glfw/glfw.git external\glfw || goto :error
    echo [OK] GLFW 3.4
) else (
    echo [OK] GLFW already present
)

if not exist "external\imgui\imgui.h" (
    echo   Cloning Dear ImGui ^(docking^)...
    if exist "external\imgui" rmdir /s /q "external\imgui"
    git clone --depth 1 --branch docking https://github.com/ocornut/imgui.git external\imgui || goto :error
    echo [OK] Dear ImGui ^(docking^)
) else (
    echo [OK] Dear ImGui already present
)

if not exist "external\implot\implot.h" (
    echo   Cloning ImPlot...
    if exist "external\implot" rmdir /s /q "external\implot"
    git clone --depth 1 https://github.com/epezent/implot.git external\implot || goto :error
    echo [OK] ImPlot
) else (
    echo [OK] ImPlot already present
)

if not exist "external\googletest\CMakeLists.txt" (
    echo   Cloning GoogleTest v1.14.0...
    if exist "external\googletest" rmdir /s /q "external\googletest"
    git clone --depth 1 --branch v1.14.0 https://github.com/google/googletest.git external\googletest || goto :error
    echo [OK] GoogleTest v1.14.0
) else (
    echo [OK] GoogleTest already present
)

if not exist "external\nfd" (
    echo   Cloning nativefiledialog-extended...
    if exist "external\nfd" rmdir /s /q "external\nfd"
    git clone --depth 1 https://github.com/btzy/nativefiledialog-extended.git external\nfd || goto :error
    echo [OK] nativefiledialog-extended
) else (
    echo [OK] nativefiledialog-extended already present
)

if not exist "external\nlohmann\json.hpp" (
    echo   Downloading nlohmann/json...
    if not exist "external\nlohmann" mkdir "external\nlohmann"
    %PS% "Invoke-WebRequest -UseBasicParsing 'https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp' -OutFile 'external/nlohmann/json.hpp'" || goto :error
    echo [OK] nlohmann/json v3.11.3
) else (
    echo [OK] nlohmann/json already present
)

if not exist "external\glad\src\gl.c" (
    echo   Generating glad loader...
    python -m pip install --quiet glad2 || goto :error
    if exist "external\glad" rmdir /s /q "external\glad"
    python -m glad --api gl:core=4.6 --out-path external\glad c || goto :error
    echo [OK] glad ^(OpenGL 4.6 core^)
) else (
    echo [OK] glad already present
)

echo.
echo Checking optional dependencies...

python -c "import pyscf" >nul 2>nul
if errorlevel 1 (
    echo [WARN] PySCF not found. Install with: python -m pip install pyscf
) else (
    echo [OK] PySCF available
)

python -c "import tblite" >nul 2>nul
if errorlevel 1 (
    where xtb >nul 2>nul
    if errorlevel 1 (
        echo [WARN] tblite/xtb not found. Install tblite-python or xtb.
    ) else (
        echo [OK] xtb CLI available
    )
) else (
    echo [OK] tblite available
)

python -c "import geometric" >nul 2>nul
if errorlevel 1 (
    echo [WARN] geomeTRIC not found. Install with: python -m pip install geometric
) else (
    echo [OK] geomeTRIC available
)

echo.
echo === Setup complete ===
echo.
echo To build:
echo   mkdir build ^&^& cd build
echo   cmake .. -DCMAKE_BUILD_TYPE=Release
echo   cmake --build . --config Release
echo.
echo To run tests:
echo   cd build ^&^& ctest --output-on-failure -C Release
echo.
echo To launch:
echo   build\Release\schrodingers_sandbox.exe
exit /b 0

:check_tool
where %1 >nul 2>nul
if errorlevel 1 (
    echo [FAIL] %~2
    exit /b 1
)
echo [OK] %1 found
exit /b 0

:error
echo [FAIL] Setup failed.
exit /b 1
