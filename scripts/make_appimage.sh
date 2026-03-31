#!/usr/bin/env bash
set -euo pipefail

echo "=== Building Schrödinger's Sandbox AppImage ==="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_DIR}/build"
APPDIR="${BUILD_DIR}/AppDir"
VERSION=$(cd "$PROJECT_DIR" && git describe --tags --always 2>/dev/null || echo "0.1.0")

if ! command -v appimagetool &>/dev/null; then
    echo "Downloading appimagetool..."
    curl -sL https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage -o /tmp/appimagetool
    chmod +x /tmp/appimagetool
    APPIMAGETOOL=/tmp/appimagetool
else
    APPIMAGETOOL=appimagetool
fi

if [ ! -f "${BUILD_DIR}/schrodingers_sandbox" ]; then
    echo "Building project..."
    cd "$PROJECT_DIR"
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr
    cmake --build build -j"$(nproc)"
fi

echo "Creating AppDir..."
rm -rf "$APPDIR"
mkdir -p "${APPDIR}/usr/bin"
mkdir -p "${APPDIR}/usr/share/schrodingers_sandbox/data/shaders"
mkdir -p "${APPDIR}/usr/share/schrodingers_sandbox/data/scripts"
mkdir -p "${APPDIR}/usr/share/applications"
mkdir -p "${APPDIR}/usr/share/icons/hicolor/256x256/apps"
mkdir -p "${APPDIR}/usr/share/mime/packages"

cp "${BUILD_DIR}/schrodingers_sandbox" "${APPDIR}/usr/bin/"
cp -r "${PROJECT_DIR}/data/shaders/"* "${APPDIR}/usr/share/schrodingers_sandbox/data/shaders/"
cp -r "${PROJECT_DIR}/data/scripts/"* "${APPDIR}/usr/share/schrodingers_sandbox/data/scripts/"

cat > "${APPDIR}/usr/share/applications/schrodingers_sandbox.desktop" <<'EOF'
[Desktop Entry]
Type=Application
Name=Schrödinger's Sandbox
Comment=Interactive quantum chemistry visualisation
Exec=schrodingers_sandbox %f
Icon=schrodingers_sandbox
Categories=Education;Science;Chemistry;
MimeType=chemical/x-xyz;chemical/x-pdb;chemical/x-mol;application/x-sbox;
Terminal=false
EOF

cp "${APPDIR}/usr/share/applications/schrodingers_sandbox.desktop" "${APPDIR}/"

if [ -f "${PROJECT_DIR}/resources/icon_256x256.png" ]; then
    cp "${PROJECT_DIR}/resources/icon_256x256.png" "${APPDIR}/usr/share/icons/hicolor/256x256/apps/schrodingers_sandbox.png"
    cp "${PROJECT_DIR}/resources/icon_256x256.png" "${APPDIR}/schrodingers_sandbox.png"
else
    echo "Warning: No icon found, using placeholder"
    python3 -c "
import struct, zlib
def create_png(w, h, color):
    def chunk(t, d):
        c = t + d
        return struct.pack('>I', len(d)) + c + struct.pack('>I', zlib.crc32(c) & 0xFFFFFFFF)
    raw = b''
    for _ in range(h):
        raw += b'\x00' + bytes(color) * w
    return b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 6, 0, 0, 0)) + chunk(b'IDAT', zlib.compress(raw)) + chunk(b'IEND', b'')
with open('${APPDIR}/schrodingers_sandbox.png', 'wb') as f:
    f.write(create_png(256, 256, [20, 60, 90, 255]))
" 2>/dev/null || touch "${APPDIR}/schrodingers_sandbox.png"
    cp "${APPDIR}/schrodingers_sandbox.png" "${APPDIR}/usr/share/icons/hicolor/256x256/apps/schrodingers_sandbox.png"
fi

if [ -f "${PROJECT_DIR}/resources/schrodingers-sandbox.xml" ]; then
    cp "${PROJECT_DIR}/resources/schrodingers-sandbox.xml" "${APPDIR}/usr/share/mime/packages/"
fi

cat > "${APPDIR}/AppRun" <<'EOF'
#!/bin/bash
SELF=$(readlink -f "$0")
HERE=${SELF%/*}
export PATH="${HERE}/usr/bin:${PATH}"
export LD_LIBRARY_PATH="${HERE}/usr/lib:${LD_LIBRARY_PATH:-}"
export SBOX_DATA_DIR="${HERE}/usr/share/schrodingers_sandbox/data"
exec "${HERE}/usr/bin/schrodingers_sandbox" "$@"
EOF
chmod +x "${APPDIR}/AppRun"

echo "Bundling shared libraries..."
for lib in libstdc++.so.6 libgcc_s.so.1; do
    LIB_PATH=$(ldconfig -p 2>/dev/null | grep "$lib" | head -1 | awk '{print $NF}')
    if [ -n "${LIB_PATH:-}" ] && [ -f "$LIB_PATH" ]; then
        mkdir -p "${APPDIR}/usr/lib"
        cp "$LIB_PATH" "${APPDIR}/usr/lib/" 2>/dev/null || true
    fi
done

echo "Building AppImage..."
ARCH=x86_64 "$APPIMAGETOOL" "$APPDIR" "${BUILD_DIR}/SchrodingersSandbox-${VERSION}-x86_64.AppImage"

echo ""
echo "AppImage created: ${BUILD_DIR}/SchrodingersSandbox-${VERSION}-x86_64.AppImage"
echo "Run with: ./${BUILD_DIR}/SchrodingersSandbox-${VERSION}-x86_64.AppImage"
