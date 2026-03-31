#!/usr/bin/env bash
set -euo pipefail

APP_PATH="${1:-build/schrodingers_sandbox.app}"
DMG_NAME="SchrodingersSandbox-$(git describe --tags --always 2>/dev/null || echo dev)"
DMG_PATH="build/${DMG_NAME}.dmg"

if [ ! -d "$APP_PATH" ]; then
    echo "Error: App bundle not found at $APP_PATH"
    echo "Build with: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build"
    exit 1
fi

STAGING="build/dmg_staging"
rm -rf "$STAGING"
mkdir -p "$STAGING"

cp -R "$APP_PATH" "$STAGING/"
ln -s /Applications "$STAGING/Applications"

hdiutil create -volname "Schrödinger's Sandbox" \
    -srcfolder "$STAGING" \
    -ov -format UDZO \
    "$DMG_PATH"

echo "DMG created: $DMG_PATH"

rm -rf "$STAGING"
