# Plasmid Demix — Desktop GUI

A standalone desktop app that wraps `fullPlasmidSeq_demix_RFW.py` in a
point-and-click interface so non-Python users can run the pipeline.

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
