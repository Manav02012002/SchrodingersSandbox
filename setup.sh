#!/usr/bin/env bash
set -euo pipefail

echo "=== Schrödinger's Sandbox — Dependency Setup ==="
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok()   { echo -e "${GREEN}✓${NC} $1"; }
warn() { echo -e "${YELLOW}⚠${NC} $1"; }
fail() { echo -e "${RED}✗${NC} $1"; }

echo "Checking prerequisites..."

command -v cmake >/dev/null 2>&1 || {
    fail "cmake not found. Install with: brew install cmake (macOS) or apt install cmake (Linux)"
    exit 1
}
ok "cmake $(cmake --version | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"

command -v git >/dev/null 2>&1 || {
    fail "git not found"
    exit 1
}
ok "git $(git --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"

command -v python3 >/dev/null 2>&1 || {
    fail "python3 not found"
    exit 1
}
ok "python3 $(python3 --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"

EIGEN_FOUND=false
if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists eigen3 2>/dev/null; then
    ok "Eigen3 $(pkg-config --modversion eigen3)"
    EIGEN_FOUND=true
elif [ -d "/opt/homebrew/include/eigen3" ] || [ -d "/usr/include/eigen3" ] || [ -d "/usr/local/include/eigen3" ]; then
    ok "Eigen3 found in system includes"
    EIGEN_FOUND=true
fi

if [ "$EIGEN_FOUND" = false ]; then
    warn "Eigen3 not found. Install with: brew install eigen (macOS) or apt install libeigen3-dev (Linux)"
fi

echo ""
echo "Setting up external dependencies..."

if [ ! -f "external/glfw/CMakeLists.txt" ]; then
    echo "  Cloning GLFW 3.4..."
    rm -rf external/glfw
    git clone --depth 1 --branch 3.4 https://github.com/glfw/glfw.git external/glfw
    ok "GLFW 3.4"
else
    ok "GLFW (already present)"
fi

if [ ! -f "external/imgui/imgui.h" ]; then
    echo "  Cloning Dear ImGui (docking)..."
    rm -rf external/imgui
    git clone --depth 1 --branch docking https://github.com/ocornut/imgui.git external/imgui
    ok "Dear ImGui (docking)"
else
    ok "Dear ImGui (already present)"
fi

if [ ! -f "external/implot/implot.h" ]; then
    echo "  Cloning ImPlot..."
    rm -rf external/implot
    git clone --depth 1 https://github.com/epezent/implot.git external/implot
    ok "ImPlot"
else
    ok "ImPlot (already present)"
fi

if [ ! -f "external/googletest/CMakeLists.txt" ]; then
    echo "  Cloning GoogleTest v1.14.0..."
    rm -rf external/googletest
    git clone --depth 1 --branch v1.14.0 https://github.com/google/googletest.git external/googletest
    ok "GoogleTest v1.14.0"
else
    ok "GoogleTest (already present)"
fi

if [ ! -d "external/nfd" ]; then
    echo "  Cloning nativefiledialog-extended..."
    rm -rf external/nfd
    git clone --depth 1 https://github.com/btzy/nativefiledialog-extended.git external/nfd
    ok "nativefiledialog-extended"
else
    ok "nativefiledialog-extended (already present)"
fi

if [ ! -f "external/nlohmann/json.hpp" ]; then
    echo "  Downloading nlohmann/json..."
    mkdir -p external/nlohmann
    curl -sL https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp -o external/nlohmann/json.hpp
    ok "nlohmann/json v3.11.3"
else
    ok "nlohmann/json (already present)"
fi

if [ ! -f "external/stb/stb_image_write.h" ]; then
    mkdir -p external/stb
    curl -sL https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h -o external/stb/stb_image_write.h
    ok "stb_image_write"
else
    ok "stb_image_write (already present)"
fi

if [ ! -f "external/glad/src/gl.c" ]; then
    echo "  Generating glad loader..."
    python3 -m pip install --quiet glad2 2>/dev/null || python3 -m pip install --quiet --user glad2
    rm -rf external/glad
    python3 -m glad --api gl:core=4.6 --out-path external/glad c
    ok "glad (OpenGL 4.6 core)"
else
    ok "glad (already present)"
fi

echo ""
echo "Checking optional dependencies..."

if python3 -c "import pyscf; print(pyscf.__version__)" 2>/dev/null; then
    ok "PySCF $(python3 -c 'import pyscf; print(pyscf.__version__)')"
else
    warn "PySCF not found. Molecular calculations will be unavailable."
    echo "       Install with: pip3 install pyscf"
fi

if python3 -c "import tblite" 2>/dev/null; then
    ok "tblite available"
elif command -v xtb >/dev/null 2>&1; then
    ok "xtb CLI available"
else
    warn "tblite/xtb not found. Semi-empirical calculations will be unavailable."
    echo "       Install with: pip3 install tblite-python  OR  conda install xtb"
fi

if python3 -c "import geometric" 2>/dev/null; then
    ok "geomeTRIC available"
else
    warn "geomeTRIC not found. Geometry optimisation will be unavailable."
    echo "       Install with: pip3 install geometric"
fi

echo ""
echo "=== Setup complete ==="
echo ""
echo "To build:"
echo "  mkdir -p build && cd build"
echo "  cmake .. -DCMAKE_BUILD_TYPE=Release"
echo "  make -j\$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
echo ""
echo "To run tests:"
echo "  cd build && ctest --output-on-failure"
echo ""
echo "To launch:"
echo "  ./build/schrodingers_sandbox"
