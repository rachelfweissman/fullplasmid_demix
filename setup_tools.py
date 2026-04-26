"""
First-launch setup: check for required external bioinformatics tools and offer
to install them via micromamba/conda into a dedicated environment that the GUI
will then prepend to PATH.

The actual demix pipeline (fullPlasmidSeq_demix_RFW.py) shells out to:
    minimap2, samtools, flye, racon, medaka_consensus, circlator
This module makes sure those are reachable.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, Iterable

REQUIRED_TOOLS = ["minimap2", "samtools", "flye", "racon", "medaka_consensus"]
OPTIONAL_TOOLS = ["circlator"]

# Where we install the bundled bioinformatics environment on first launch.
APP_DIR = Path(os.environ.get("DEMIX_HOME", Path.home() / ".plasmid_demix"))
ENV_DIR = APP_DIR / "bio-env"
MICROMAMBA_BIN = APP_DIR / "bin" / ("micromamba.exe" if os.name == "nt" else "micromamba")


class ToolsMissingError(RuntimeError):
    pass


def _log(log: Callable[[str], None] | None, msg: str) -> None:
    if log:
        log(msg)
    else:
        print(msg)


def _which_in(paths: Iterable[Path], tool: str) -> str | None:
    for p in paths:
        candidate = p / tool
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)
    return shutil.which(tool)


def _missing_tools(extra_path: Path | None = None) -> list[str]:
    paths: list[Path] = []
    if extra_path and extra_path.exists():
        paths.append(extra_path)
    return [t for t in REQUIRED_TOOLS if not _which_in(paths, t)]


def _prepend_path(env_bin: Path) -> None:
    if env_bin.exists():
        os.environ["PATH"] = f"{env_bin}{os.pathsep}{os.environ.get('PATH', '')}"


def _download_micromamba(log: Callable[[str], None] | None) -> Path:
    """Download a static micromamba binary for this OS/arch."""
    import platform
    import tarfile
    import urllib.request

    APP_DIR.mkdir(parents=True, exist_ok=True)
    (APP_DIR / "bin").mkdir(parents=True, exist_ok=True)
    if MICROMAMBA_BIN.exists():
        return MICROMAMBA_BIN

    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "darwin":
        plat = "osx-arm64" if machine in ("arm64", "aarch64") else "osx-64"
    elif system == "linux":
        plat = "linux-aarch64" if machine in ("arm64", "aarch64") else "linux-64"
    elif system == "windows":
        plat = "win-64"
    else:
        raise ToolsMissingError(f"Unsupported platform: {system}/{machine}")

    url = f"https://micro.mamba.pm/api/micromamba/{plat}/latest"
    _log(log, f"[setup] Downloading micromamba from {url}")
    tar_path = APP_DIR / "micromamba.tar.bz2"
    urllib.request.urlretrieve(url, tar_path)

    with tarfile.open(tar_path, "r:bz2") as tar:
        member_name = "Library/bin/micromamba.exe" if os.name == "nt" else "bin/micromamba"
        member = tar.getmember(member_name)
        member.name = MICROMAMBA_BIN.name
        tar.extract(member, path=APP_DIR / "bin")
    tar_path.unlink(missing_ok=True)
    os.chmod(MICROMAMBA_BIN, 0o755)
    _log(log, f"[setup] micromamba installed at {MICROMAMBA_BIN}")
    return MICROMAMBA_BIN


def _create_env(micromamba: Path, env_yml: Path, log: Callable[[str], None] | None) -> None:
    if (ENV_DIR / "bin").exists():
        return
    _log(log, f"[setup] Creating bioinformatics env at {ENV_DIR} (this may take several minutes)…")
    cmd = [
        str(micromamba), "create", "-y",
        "-p", str(ENV_DIR),
        "-f", str(env_yml),
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        _log(log, "[setup] " + line.rstrip())
    proc.wait()
    if proc.returncode != 0:
        raise ToolsMissingError(f"Environment creation failed (exit {proc.returncode})")


def ensure_tools_available(log: Callable[[str], None] | None = None) -> None:
    """
    Make sure all required tools are reachable on PATH. If not, install them
    into a dedicated micromamba env on first launch.
    """
    env_bin = ENV_DIR / ("Scripts" if os.name == "nt" else "bin")
    _prepend_path(env_bin)

    missing = _missing_tools(env_bin)
    if not missing:
        _log(log, "[setup] All bioinformatics tools detected.")
        return

    _log(log, f"[setup] Missing tools: {', '.join(missing)}")

    env_yml = _resolve_env_yaml()
    if env_yml is None:
        raise ToolsMissingError(
            "Could not find demixEnv2.yml next to the app. "
            "Please install the conda env manually: `mamba env create -f demixEnv2.yml`"
        )

    micromamba = _download_micromamba(log)
    _create_env(micromamba, env_yml, log)
    _prepend_path(env_bin)

    still_missing = _missing_tools(env_bin)
    if still_missing:
        raise ToolsMissingError(
            f"Tools still missing after install: {', '.join(still_missing)}"
        )
    _log(log, "[setup] Setup complete.")


def _resolve_env_yaml() -> Path | None:
    """Find demixEnv2.yml — handles both source layout and PyInstaller bundle."""
    candidates = []
    if getattr(sys, "_MEIPASS", None):  # PyInstaller bundle
        candidates.append(Path(sys._MEIPASS) / "demixEnv2.yml")
    here = Path(__file__).resolve().parent
    candidates += [here / "demixEnv2.yml", here.parent / "demixEnv2.yml"]
    for c in candidates:
        if c.exists():
            return c
    return None
