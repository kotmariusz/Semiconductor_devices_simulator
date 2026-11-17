#!/usr/bin/env python3
"""Install SemiSim 3D on Windows (double-click install.bat, which runs this).

  1. checks this Python: 64-bit, 3.9 or newer, with Tk 8.6+;
  2. makes a private Python environment in .venv and installs the packages
     of requirements.txt into it (numpy, matplotlib);
  3. the solver core: keeps a ready-made semiconductor_core.dll that matches
     the sources (the release zip has one), else builds it with a MinGW-w64
     gcc (semibuild.py), else says how to get one;
  4. a short self-test (validate.py sections 1, 3, 5);
  5. SemiSim3D.bat (no console) and "SemiSim3D (console).bat" (shows error
     messages) in this folder, and shortcuts in the Start menu and on the
     Desktop.

Needs no administrator rights.  Run it again at any time (after an update,
for example); it keeps what is already there.

    py install_win.py [--no-shortcuts] [--portable]
"""
import os, sys, shutil, struct, platform, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
VENV = os.path.join(HERE, ".venv")
VPY = os.path.join(VENV, "Scripts", "python.exe")
VPYW = os.path.join(VENV, "Scripts", "pythonw.exe")
NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def say(msg):
    print(f"\n==> {msg}", flush=True)


def note(msg):
    print(f"    {msg}", flush=True)


def die(msg):
    print(f"\n!!  {msg}\n", flush=True)
    sys.exit(1)


# Not passed on to the programs started here: variables an IDE (Spyder, VS
# Code, PyCharm) sets for its own Python - the .venv Python would load the
# IDE's packages or its inline matplotlib backend - and SEMISIM_SO, which
# would make the self-test load another core than the one in this folder.
DROP_VARS = ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP", "MPLBACKEND", "SEMISIM_SO")

# Windows exit codes of a crashed program
CRASH = {0xC0000005: "access violation", 0xC00000FD: "stack overflow",
         0xC000001D: "illegal instruction", 0xC0000409: "stack buffer overrun",
         0xC0000135: "a DLL was not found", 0xC0000142: "a DLL failed to initialise"}


def child_env():
    env = {k: v for k, v in os.environ.items() if k.upper() not in DROP_VARS}
    env["PYTHONIOENCODING"] = "utf-8"
    return env


def run(cmd):
    """Run a command; True when it succeeds.  Its output is passed through this
    script's own output, so it also shows in an IDE's console (Spyder's
    console, for example, does not show what a child process prints)."""
    print("    $ " + " ".join(f'"{c}"' if " " in c else c for c in cmd), flush=True)
    try:
        p = subprocess.Popen(cmd, cwd=HERE, env=child_env(), stdin=subprocess.DEVNULL,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, encoding="utf-8", errors="replace", bufsize=1)
    except OSError as e:
        note(f"cannot run it: {e}")
        return False
    with p.stdout:
        for line in p.stdout:
            print("    | " + line.rstrip("\r\n"), flush=True)
    rc = p.wait() & 0xFFFFFFFF
    if rc and rc not in (1, 2):
        note(f"exit code {rc} (0x{rc:08X}{', ' + CRASH[rc] if rc in CRASH else ''})")
    return rc == 0


def check_python():
    say("checking Python")
    note(f"{sys.executable}  ({sys.version.split()[0]}, {struct.calcsize('P') * 8}-bit, {platform.machine()})")
    if struct.calcsize("P") != 8:
        die("this Python is 32-bit; install the 64-bit Python from https://www.python.org/downloads/windows/")
    if platform.machine().lower() in ("arm64", "aarch64"):
        die("this Python is built for ARM64; install the x64 (AMD64) Python from python.org - it runs on ARM\n"
            "    Windows too - and run install.bat again")
    if sys.version_info < (3, 9):
        die(f"Python {sys.version.split()[0]} is too old; install Python 3.9 or newer from python.org")
    try:
        import tkinter
        tkv = tkinter.TkVersion
    except Exception:
        die("this Python has no Tk (tkinter).  Install Python from python.org with 'tcl/tk and IDLE'\n"
            "    ticked, or: Settings -> Apps -> Python -> Modify, and tick it")
    if tkv < 8.6:
        die(f"this Python's Tk {tkv} is too old (8.6 or newer is needed)")
    note(f"Tk {tkv}: ok")
    exe = sys.executable.lower()
    if "spyder" in exe or os.path.isdir(os.path.join(sys.base_prefix, "conda-meta")):
        note("this Python belongs to " + ("Spyder" if "spyder" in exe else "a conda environment (Anaconda)")
             + "; .venv will depend on it, and an update of it can break SemiSim 3D.")
        note("Recommended: install Python from python.org, then double-click install.bat.")
    if os.environ.get("SEMISIM_SO"):
        note(f"SEMISIM_SO is set ({os.environ['SEMISIM_SO']}): the self-test uses the core in this folder;")
        note("SemiSim 3D itself would load that file - remove the variable when you are done testing it.")


def venv_home():
    """The folder of the Python that .venv was made from (pyvenv.cfg 'home')."""
    try:
        with open(os.path.join(VENV, "pyvenv.cfg"), encoding="utf-8") as f:
            for line in f:
                k, _, v = line.partition("=")
                if k.strip().lower() == "home":
                    return os.path.normcase(os.path.abspath(v.strip()))
    except OSError:
        pass
    return None


def make_venv():
    say("Python packages (in .venv, private to SemiSim 3D)")
    base = os.path.normcase(os.path.dirname(os.path.abspath(getattr(sys, "_base_executable", "") or sys.executable)))
    home = venv_home()
    if os.path.isfile(VPY) and home and home != base:
        # e.g. made from Spyder's Python, now run with python.org's
        note(f".venv was made from another Python ({home}); making it again")
        shutil.rmtree(VENV, ignore_errors=True)
        if os.path.exists(VPY):
            die("could not remove the old .venv - close SemiSim 3D and run install.bat again")
    if not os.path.isfile(VPY):
        if not run([sys.executable, "-m", "venv", VENV]):
            die("could not make the environment .venv (is this folder writable?)")
    ok = run([VPY, "-m", "pip", "install", "--disable-pip-version-check", "-q", "-r",
              os.path.join(HERE, "requirements.txt")])
    if not ok:
        die("installing numpy and matplotlib failed - is the computer online?  Then run install.bat again.")
    note("numpy, matplotlib: ok")


def core(portable):
    say("the solver core")
    sys.path.insert(0, HERE)
    import semibuild as B
    dll = os.path.join(HERE, B.lib_name(win=True))
    st = B.read_stamp(HERE)
    if (os.path.isfile(dll) and st and st.get("src") == B.src_hash()
            and (st.get("system"), st.get("machine")) == ("Windows", "amd64") and not portable):
        rep = B.dll_report(dll)
        if not rep["bad"] and rep["sentinel"]:
            note(f"ready-made core found ({'from a release' if st.get('prebuilt') else 'built here'}, "
                 f"C runtime {rep['crt'] or '?'}): kept")
            return
    args = [VPY, os.path.join(HERE, "semibuild.py")] + (["--portable"] if portable else [])
    if not run(args):
        die("there is no usable solver core.\n\n" + B._HINT_WIN.replace("\n", "\n    "))


def self_test(portable):
    say("self-test (a few seconds)")
    v = [VPY, os.path.join(HERE, "validate.py"), "--only", "1,3,5"]
    if run(v):
        return
    sys.path.insert(0, HERE)
    import semibuild as B
    if not portable and B.compilers():
        note("retrying with a portable build of the core")
        if run([VPY, os.path.join(HERE, "semibuild.py"), "--portable"]) and run(v):
            return
    die("the self-test failed - please report it (with the text above) at\n"
        "    https://github.com/kotmariusz/Semiconductor_devices_simulator/issues")


def launchers():
    say("launchers in this folder")
    a = os.path.join(HERE, "SemiSim3D.bat")
    b = os.path.join(HERE, "SemiSim3D (console).bat")
    with open(a, "w", encoding="ascii", newline="\r\n") as f:
        f.write('@echo off\nrem start SemiSim 3D without a console window (made by install.bat)\n'
                'start "" "%~dp0.venv\\Scripts\\pythonw.exe" "%~dp0main.py" %*\n')
    with open(b, "w", encoding="ascii", newline="\r\n") as f:
        f.write('@echo off\nrem start SemiSim 3D with a console that shows its messages (made by install.bat)\n'
                'chcp 65001 >nul\n"%~dp0.venv\\Scripts\\python.exe" "%~dp0main.py" %*\npause\n')
    note(os.path.basename(a)); note(os.path.basename(b))


PS_SHORTCUTS = r"""
$ErrorActionPreference = 'Stop'
$ws = New-Object -ComObject WScript.Shell
$start = Join-Path ([Environment]::GetFolderPath('Programs')) 'SemiSim 3D.lnk'
$desk  = Join-Path ([Environment]::GetFolderPath('Desktop')) 'SemiSim 3D.lnk'
foreach ($p in @($start, $desk)) {
    $s = $ws.CreateShortcut($p)
    $s.TargetPath = $env:SS_TARGET
    $s.Arguments = '"' + $env:SS_MAIN + '"'
    $s.WorkingDirectory = $env:SS_DIR
    $s.IconLocation = $env:SS_ICON + ',0'
    $s.Description = 'SemiSim 3D - semiconductor device simulator'
    $s.Save()
    Write-Output $p
}
"""


def shortcuts():
    say("shortcuts: Start menu and Desktop")
    env = dict(os.environ, SS_TARGET=VPYW, SS_MAIN=os.path.join(HERE, "main.py"), SS_DIR=HERE,
               SS_ICON=os.path.join(HERE, "docs", "icon.ico"))
    try:
        r = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass",
                            "-Command", PS_SHORTCUTS], env=env, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=120, creationflags=NO_WINDOW)
    except Exception as e:
        note(f"no shortcuts ({e}); start SemiSim3D.bat in this folder instead"); return
    if r.returncode == 0:
        for line in r.stdout.strip().splitlines():
            note(line)
    else:
        note("no shortcuts (PowerShell said: " + (r.stderr.strip().splitlines() or ["?"])[-1] + ")")
        note("start SemiSim3D.bat in this folder instead")


def main(argv):
    if sys.platform != "win32":
        die("install_win.py is for Windows; on Linux, WSL2 and macOS run:  bash install.sh")
    no_sc = "--no-shortcuts" in argv
    portable = "--portable" in argv
    print(f"SemiSim 3D - installing in {HERE}", flush=True)
    check_python()
    make_venv()
    core(portable)
    self_test(portable)
    launchers()
    if not no_sc:
        shortcuts()
    say("done.  Start SemiSim 3D from the Start menu, the Desktop shortcut or SemiSim3D.bat")
    note("if it does not start, use 'SemiSim3D (console).bat' to see why")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
