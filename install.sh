#!/usr/bin/env bash
# SemiSim 3D - installer for Linux (Ubuntu / Debian / WSL2, Fedora, Arch) and macOS.
#
#   bash install.sh
#
#   1. installs what is missing: Python 3 with Tk, NumPy, Matplotlib, a C compiler
#      (Linux: the distribution's packages; macOS: a virtual environment in .venv)
#   2. builds the solver core for this computer (semibuild.py)
#   3. runs a short self-test
#   4. adds a launcher: the command 'semisim3d' (and on macOS SemiSim3D.command,
#      which opens with a double click in Finder)
#
# Options (environment variables):
#   PYTHON=/path/to/python3   use this Python
#   SEMISIM_VENV=1            Linux: NumPy/Matplotlib from pip in .venv instead of the distribution
#   NO_OPENMP=1               macOS: do not install libomp (the core is then single-threaded)
set -u
cd "$(dirname "$0")" || exit 1
DIR="$(pwd -P)"
say() { printf '== %s\n' "$*"; }
note() { printf '   %s\n' "$*"; }
die() { printf '!! %s\n' "$*" >&2; exit 1; }
OS="$(uname -s)"
VENV="$DIR/.venv"
say "SemiSim 3D: installing in $DIR ($OS $(uname -m))"

# a Python that is new enough and has a usable Tk (8.6 or newer)
py_ok() {
    "$1" - <<'EOF' >/dev/null 2>&1
import sys, tkinter
assert sys.version_info >= (3, 8)
assert tkinter.TkVersion >= 8.6
EOF
}
has_mods() { "$1" -c "import numpy, matplotlib" >/dev/null 2>&1; }

make_venv() {   # make_venv BASE_PYTHON -> sets PYRUN
    say "creating the virtual environment .venv (NumPy, Matplotlib from pip)"
    "$1" -m venv "$VENV" || die "could not create .venv with $1 (Ubuntu: sudo apt install python3-venv)"
    "$VENV/bin/python" -m pip install --upgrade pip >/dev/null 2>&1 || true
    "$VENV/bin/python" -m pip install -r "$DIR/requirements.txt" || die "pip could not install NumPy/Matplotlib"
    PYRUN="$VENV/bin/python"
}

# ──────────────────────────────────────────────────────────── macOS
mac_setup() {
    if ! xcode-select -p >/dev/null 2>&1; then
        say "Apple's command line tools (the C compiler) are missing - starting their installer"
        xcode-select --install >/dev/null 2>&1 || true
        die "finish the 'command line developer tools' installation in the window that opened, then run:  bash install.sh"
    fi
    BREW=""
    for b in /opt/homebrew/bin/brew /usr/local/bin/brew "$(command -v brew 2>/dev/null)"; do
        [ -n "$b" ] && [ -x "$b" ] && { BREW="$b"; break; }
    done
    # Python: $PYTHON, python.org, Homebrew, PATH (not Apple's /usr/bin/python3: its Tk 8.5 is broken)
    PY=""
    CANDS="${PYTHON:-}"
    for v in 3.14 3.13 3.12 3.11 3.10 3.9; do
        CANDS="$CANDS /Library/Frameworks/Python.framework/Versions/$v/bin/python$v"
    done
    CANDS="$CANDS /opt/homebrew/bin/python3 /usr/local/bin/python3 $(command -v python3 2>/dev/null)"
    for p in $CANDS; do
        [ -x "$p" ] || continue
        if py_ok "$p"; then PY="$p"; break; fi
    done
    if [ -z "$PY" ]; then
        if [ -n "$BREW" ]; then
            say "no Python with Tk found - installing Homebrew's python and python-tk"
            "$BREW" install python python-tk || die "brew install python python-tk failed"
            P="$("$BREW" --prefix)/bin/python3"
            py_ok "$P" && PY="$P"
        fi
        [ -n "$PY" ] || die "no Python 3 with Tk found. Install Python from https://www.python.org/downloads/macos/ (it includes Tk), then run:  bash install.sh"
    fi
    note "Python: $PY ($("$PY" -c 'import sys,tkinter; print(sys.version.split()[0], "with Tk", tkinter.TkVersion)'))"
    if [ -x "$VENV/bin/python" ] && py_ok "$VENV/bin/python" && has_mods "$VENV/bin/python"; then
        PYRUN="$VENV/bin/python"; note "using the existing .venv"
        "$PYRUN" -m pip install -q -r "$DIR/requirements.txt" >/dev/null 2>&1 || true
    else
        rm -rf "$VENV"; make_venv "$PY"
    fi
    # OpenMP for the solver's parallel loops: Homebrew's libomp
    if [ -z "${NO_OPENMP:-}" ] && [ -n "$BREW" ]; then
        if ! "$BREW" list libomp >/dev/null 2>&1; then
            say "installing libomp (multi-core solving) with Homebrew"
            "$BREW" install libomp || note "libomp could not be installed - the core will be single-threaded"
        fi
    elif [ -z "$BREW" ]; then
        note "Homebrew not found: the core is built single-threaded (for all cores: install Homebrew, then  brew install libomp  and  bash install.sh)"
    fi
}

# ──────────────────────────────────────────────────────────── Linux
linux_setup() {
    PY="${PYTHON:-python3}"
    command -v "$PY" >/dev/null 2>&1 || PY=""
    need=""
    if [ -z "$PY" ]; then need="python"
    else
        "$PY" -c "import tkinter" >/dev/null 2>&1 || need="$need tk"
        if [ -z "${SEMISIM_VENV:-}" ]; then
            "$PY" -c "import numpy" >/dev/null 2>&1 || need="$need numpy"
            "$PY" -c "import matplotlib" >/dev/null 2>&1 || need="$need matplotlib"
        fi
    fi
    command -v gcc >/dev/null 2>&1 || command -v cc >/dev/null 2>&1 || need="$need gcc"
    if [ -n "$need" ]; then
        say "installing:$need (sudo asks for your password)"
        if command -v apt-get >/dev/null 2>&1; then
            pk="python3"
            for n in $need; do case $n in
                tk) pk="$pk python3-tk";; numpy) pk="$pk python3-numpy";; matplotlib) pk="$pk python3-matplotlib";;
                gcc) pk="$pk gcc";; esac; done
            [ -n "${SEMISIM_VENV:-}" ] && pk="$pk python3-venv"
            sudo apt-get update || note "(apt-get update reported errors - continuing)"
            sudo apt-get install -y $pk
        elif command -v dnf >/dev/null 2>&1; then
            pk="python3"
            for n in $need; do case $n in
                tk) pk="$pk python3-tkinter";; numpy) pk="$pk python3-numpy";; matplotlib) pk="$pk python3-matplotlib";;
                gcc) pk="$pk gcc";; esac; done
            sudo dnf install -y $pk
        elif command -v pacman >/dev/null 2>&1; then
            pk="python"
            for n in $need; do case $n in
                tk) pk="$pk tk";; numpy) pk="$pk python-numpy";; matplotlib) pk="$pk python-matplotlib";;
                gcc) pk="$pk gcc";; esac; done
            sudo pacman -S --needed --noconfirm $pk
        elif command -v zypper >/dev/null 2>&1; then
            pk="python3"
            for n in $need; do case $n in
                tk) pk="$pk python3-tk";; numpy) pk="$pk python3-numpy";; matplotlib) pk="$pk python3-matplotlib";;
                gcc) pk="$pk gcc";; esac; done
            sudo zypper install -y $pk
        else
            die "install these with your package manager, then run install.sh again:$need"
        fi
        PY="${PYTHON:-python3}"
    fi
    command -v "$PY" >/dev/null 2>&1 || die "python3 not found"
    PY="$(command -v "$PY")"
    py_ok "$PY" || die "$PY has no usable Tk (Ubuntu: sudo apt install python3-tk) - or run: PYTHON=/usr/bin/python3 bash install.sh"
    if [ -n "${SEMISIM_VENV:-}" ] || ! has_mods "$PY"; then
        if [ -x "$VENV/bin/python" ] && has_mods "$VENV/bin/python"; then PYRUN="$VENV/bin/python"
        else rm -rf "$VENV"; make_venv "$PY"; fi
    else
        PYRUN="$PY"
    fi
}

case "$OS" in
    Darwin) mac_setup ;;
    Linux)  linux_setup ;;
    *) die "unsupported system '$OS' (Linux, WSL2 and macOS are supported; on Windows use WSL2)" ;;
esac
note "running with: $PYRUN"

# ──────────────────────────────────────────────────────────── core
say "building the solver core for this computer"
"$PYRUN" "$DIR/semibuild.py" || die "the core could not be built - see build.log"

say "self-test (a few seconds)"
if ! "$PYRUN" "$DIR/validate.py" --only 1,3,5; then
    note "retrying with a portable build"
    "$PYRUN" "$DIR/semibuild.py" --portable && "$PYRUN" "$DIR/validate.py" --only 1,3,5 || die "self-test failed"
fi

# ──────────────────────────────────────────────────────────── launchers
mkdir -p "$HOME/.local/bin"
L="$HOME/.local/bin/semisim3d"
printf '#!/bin/sh\nexec "%s" "%s/main.py" "$@"\n' "$PYRUN" "$DIR" > "$L"
chmod +x "$L"
if [ "$OS" = "Darwin" ]; then
    C="$DIR/SemiSim3D.command"
    printf '#!/bin/sh\ncd "%s" && exec "%s" "%s/main.py"\n' "$DIR" "$PYRUN" "$DIR" > "$C"
    chmod +x "$C"
else
    # a menu entry (on WSL2 it appears in the Windows Start menu)
    A="$HOME/.local/share/applications"; mkdir -p "$A"
    cat > "$A/semisim3d.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=SemiSim 3D
Comment=3-D drift-diffusion semiconductor device simulator
Exec=$L
Icon=$DIR/docs/icon.png
Terminal=false
Categories=Science;Education;
EOF
fi

echo
say "done. Start SemiSim 3D with:"
note "$PYRUN $DIR/main.py"
case ":$PATH:" in
    *":$HOME/.local/bin:"*) note "or just:  semisim3d" ;;
    *) if [ "$OS" = "Linux" ] && grep -qs '.local/bin' "$HOME/.profile"; then
           note "or 'semisim3d' in a new terminal"
       else
       note "or 'semisim3d' after adding ~/.local/bin to your PATH:"
       if [ "$OS" = "Darwin" ]; then note "  echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> ~/.zprofile   (then open a new terminal)"
       else note "  echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> ~/.bashrc   (then open a new terminal)"; fi
       fi ;;
esac
[ "$OS" = "Darwin" ] && note "or double-click SemiSim3D.command in this folder (Finder)"
exit 0
