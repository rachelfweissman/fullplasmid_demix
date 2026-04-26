#!/usr/bin/env bash
# Build a standalone PlasmidDemix app for the current OS using PyInstaller.
#
# macOS  ->  dist/PlasmidDemix.app   (also a .dmg if create-dmg is installed)
# Linux  ->  dist/PlasmidDemix/      (folder bundle; zip it for distribution)
#
# Requirements: Python 3.10+ and `pip install -r requirements-gui.txt`.
# The bioinformatics tools are NOT bundled here — they install on first launch
# via micromamba into ~/.plasmid_demix/bio-env (see setup_tools.py).

set -euo pipefail

cd "$(dirname "$0")"

APP_NAME="PlasmidDemix"
ENTRY="gui.py"

OS="$(uname -s)"

echo "[build] Cleaning previous build…"
rm -rf build dist "${APP_NAME}.spec"

EXTRA_ARGS=()
if [[ "$OS" == "Darwin" ]]; then
    EXTRA_ARGS+=(--windowed --osx-bundle-identifier com.plasmid.demix)
    if [[ -f "assets/icon.icns" ]]; then
        EXTRA_ARGS+=(--icon "assets/icon.icns")
        echo "[build] Using icon: assets/icon.icns"
    else
        echo "[build] No assets/icon.icns found — see assets/README.md to add one."
    fi
elif [[ "$OS" == "Linux" ]]; then
    EXTRA_ARGS+=(--windowed)
    if [[ -f "assets/icon.png" ]]; then
        EXTRA_ARGS+=(--icon "assets/icon.png")
        echo "[build] Using icon: assets/icon.png"
    fi
fi

# Bundle the assets folder if it has an icon (for the runtime window icon).
if [[ -f "assets/icon.png" ]]; then
    EXTRA_ARGS+=(--add-data "assets/icon.png:assets")
fi

echo "[build] Running PyInstaller…"
pyinstaller \
    --name "$APP_NAME" \
    --noconfirm \
    --clean \
    --add-data "demixEnv2.yml:." \
    --add-data "fullPlasmidSeq_demix_RFW.py:." \
    --hidden-import customtkinter \
    --hidden-import Bio \
    --hidden-import plotly \
    --hidden-import pandas \
    --hidden-import openpyxl \
    "${EXTRA_ARGS[@]}" \
    "$ENTRY"

echo
echo "[build] Done."
if [[ "$OS" == "Darwin" ]]; then
    echo "  → dist/${APP_NAME}.app"
    if command -v create-dmg >/dev/null 2>&1; then
        echo "[build] Creating .dmg…"
        rm -f "dist/${APP_NAME}.dmg"
        create-dmg \
            --volname "$APP_NAME" \
            --window-size 540 380 \
            --icon-size 96 \
            --app-drop-link 380 180 \
            "dist/${APP_NAME}.dmg" \
            "dist/${APP_NAME}.app" || true
        echo "  → dist/${APP_NAME}.dmg"
    else
        echo "  (Install \`brew install create-dmg\` to also produce a .dmg installer.)"
    fi
else
    echo "  → dist/${APP_NAME}/"
    echo "  Zip the folder for distribution: (cd dist && tar czf ${APP_NAME}-linux.tar.gz ${APP_NAME})"
fi
