# Plasmid Demix — Desktop GUI

A standalone desktop app that wraps `fullPlasmidSeq_demix_RFW.py` in a
point-and-click interface so non-Python users can run the pipeline.

## Installing on macOS (for end users)

1. Download `PlasmidDemix.dmg` from the [latest release](https://github.com/rachelfweissman/fullplasmid_demix/releases/latest).
2. Open the `.dmg` and drag **PlasmidDemix** into your `Applications` folder.
3. **Important — one-time fix for the "damaged" warning:** because the app is not codesigned with an Apple Developer ID, macOS marks any downloaded copy as quarantined and (on Apple Silicon) refuses to open it with a misleading *"PlasmidDemix is damaged and can't be opened"* message. This is a false alarm. Open **Terminal** once and run:
   ```bash
   xattr -cr /Applications/PlasmidDemix.app
   ```
   Then double-click the app normally. You only have to do this once per machine.
4. First launch takes a few minutes while the bioinformatics tools install into `~/.plasmid_demix/`. After that, the app opens instantly.

If you'd rather skip the Terminal step, you can also right-click the **`.dmg`** (not the app) and choose **Open**, copy the app to Applications from there, and on some macOS versions that bypass works. The `xattr -cr` command always works.

## Installing on Linux

1. Download the `PlasmidDemix-linux.tar.gz` from the latest release.
2. Extract: `tar xzf PlasmidDemix-linux.tar.gz`
3. Run: `./PlasmidDemix/PlasmidDemix`

## Windows

Not supported as a native `.exe` because `flye` and `medaka` don't ship
for native Windows. Windows users should run the app inside WSL2 from
source (`python gui.py`).

## How it's structured

| File | Purpose |
| --- | --- |
| `gui.py` | CustomTkinter GUI — file pickers, options, live log panel |
| `setup_tools.py` | First-launch helper. Downloads micromamba and creates a self-contained env with `samtools`, `minimap2`, `flye`, `racon`, `medaka`, etc. into `~/.plasmid_demix/bio-env`. |
| `demixEnv2.yml` | Conda env spec — used by `setup_tools.py` and bundled into the app |
| `build_app.sh` | PyInstaller build script for macOS / Linux |
| `requirements-gui.txt` | Python deps for building the GUI bundle |

The GUI bundle itself only contains Python + the GUI code. Heavy
bioinformatics tools (~2 GB once unpacked) install on first launch via
micromamba — that keeps the distributable small and avoids re-bundling
binaries for every OS/arch combination.

## Running from source (development)

```bash
pip install -r requirements-gui.txt
python gui.py
```

On first run the app will detect that the bio tools aren't installed and
download them automatically (this takes a few minutes).

## Building a standalone app

### macOS

```bash
pip install -r requirements-gui.txt
./build_app.sh
open dist/PlasmidDemix.app
```

Optional: `brew install create-dmg` first to also produce a `.dmg` installer.

### Linux

```bash
pip install -r requirements-gui.txt
./build_app.sh
./dist/PlasmidDemix/PlasmidDemix
```

Distribute by tarring `dist/PlasmidDemix/`.

### Windows

The GUI itself works on Windows, but `flye` and `medaka` do **not** run
natively on Windows. Windows users should run the app inside WSL2.
Building a true Windows `.exe` would require swapping or removing those
dependencies and is out of scope here.

## Cross-compiling

PyInstaller does not cross-compile. To produce a macOS app, run
`build_app.sh` on a Mac. To produce a Linux bundle, run it on Linux
(GitHub Actions with a matrix of `macos-latest` and `ubuntu-latest` is
the easiest way).

## What users see

1. Double-click `PlasmidDemix.app` (macOS) or `PlasmidDemix` (Linux).
2. On first launch only: a few minutes of "[setup] …" log lines while
   the bio env is created in `~/.plasmid_demix/`.
3. Pick the Excel file, FASTQ folder, optional reference folder, output
   folder. Toggle quick / keep-temp. Click **Run**.
4. Logs stream into the bottom panel. Outputs land in the chosen output
   folder, same layout as the CLI script produces.
