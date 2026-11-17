#!/usr/bin/env python3
"""Build and load the SemiSim 3D solver core (semiconductor_3d.c) on Linux,
WSL2, macOS and Windows.

    python3 semibuild.py              build for this computer (fastest code)
    python3 semibuild.py --portable   build for any CPU of this architecture
    python3 semibuild.py --check      show the compiler, OpenMP and the core in use
    python3 semibuild.py --cross CC -o DIR/semiconductor_core.dll
                                      compile a Windows core with a MinGW-w64
                                      cross-compiler (no load test, any host)

main.py calls load_core(), which builds the core by itself when it is
missing, was built from another version of semiconductor_3d.c, or cannot be
loaded on this computer (a Linux build copied to a Mac, for example).

Needs a C compiler: gcc or clang on Linux (Ubuntu: sudo apt install gcc),
the Xcode command line tools on macOS (xcode-select --install), MinGW-w64
gcc on Windows (winget install BrechtSanders.WinLibs.POSIX.UCRT, or MSYS2's
mingw-w64-ucrt-x86_64-gcc) - or, on Windows, the ready-made core of a
release zip.  OpenMP (multi-core solving) comes with gcc; on macOS it needs
Homebrew's libomp (brew install libomp) - without it the core is built
single-threaded.
"""
import os, sys, json, time, ctypes, struct, hashlib, platform, shutil, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "semiconductor_3d.c")
SENTINEL = "api3_set_caxis"             # newest entry point; an older core lacks it
IS_MAC = sys.platform == "darwin"
IS_WIN = sys.platform == "win32"
# DLLs a Windows core must not need: Python (3.8+) looks for a DLL's own
# dependencies only in system folders and the DLL's folder, never on PATH
WIN_BAD_DEPS = ("libgomp-1.dll", "libwinpthread-1.dll", "libgcc_s_seh-1.dll", "libgcc_s_dw2-1.dll",
                "msys-2.0.dll", "cygwin1.dll")


def lib_name(win=None):
    """File name of the core: .dll on Windows, .so elsewhere (also on macOS:
    ctypes loads any file name)."""
    return "semiconductor_core.dll" if (IS_WIN if win is None else win) else "semiconductor_core.so"


NAME = lib_name()


class CoreError(RuntimeError):
    pass


# ───────────────────────────────────────────────────────── paths, stamp
def _writable(d):
    try:
        fd, p = tempfile.mkstemp(dir=d, prefix=".wtest"); os.close(fd); os.remove(p); return True
    except Exception:
        return False


def cache_dir(win=None):
    """The user folder used when the program folder is read-only."""
    if IS_WIN if win is None else win:
        base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), "AppData", "Local")
        return os.path.join(base, "SemiSim3D")
    return os.path.join(os.path.expanduser("~"), ".cache", "semisim3d")


def out_dir():
    """Where the built core lives: next to the sources, or a user cache
    folder when the program folder is read-only."""
    if _writable(HERE):
        return HERE
    d = cache_dir()
    os.makedirs(d, exist_ok=True)
    return d


def lib_path():
    return os.path.join(out_dir(), NAME)


def stamp_path(folder=None):
    return os.path.join(folder or out_dir(), ".core_build.json")


def src_hash(path=None):
    """Hash of the core's source; line endings do not count, so a core built
    from a Windows (CRLF) checkout matches the same source anywhere."""
    with open(path or SRC, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()[:16]


def read_stamp(folder=None):
    try:
        with open(stamp_path(folder), encoding="utf-8") as f:
            d = json.load(f)
        return d if isinstance(d, dict) else None
    except Exception:
        return None


def _machine():
    return platform.system(), platform.machine().lower()


def dll_report(path):
    """What a Windows core needs: its C runtime ('ucrt' / 'msvcrt' / None),
    DLLs it must not need, and whether it has the newest entry point."""
    try:
        with open(path, "rb") as f:
            data = f.read()
    except OSError:
        return dict(crt=None, bad=[], sentinel=False)
    low = data.lower()
    crt = ("ucrt" if (b"ucrtbase.dll" in low or b"api-ms-win-crt" in low)
           else "msvcrt" if b"msvcrt.dll" in low else None)
    return dict(crt=crt, bad=[n for n in WIN_BAD_DEPS if n.encode() in low],
                sentinel=SENTINEL.encode() in data)


# ───────────────────────────────────────────────────────── running tools
def _run_kw(cc=None):
    """subprocess options: no console window flashing up on Windows (a
    start through pythonw), and the compiler's own folder first on PATH
    (MSYS2's gcc needs its DLLs)."""
    kw = dict(stdin=subprocess.DEVNULL)
    if IS_WIN:
        kw["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)
    if cc and os.path.dirname(cc):
        env = dict(os.environ)
        env["PATH"] = os.path.dirname(cc) + os.pathsep + env.get("PATH", "")
        kw["env"] = env
    return kw


def _python_console():
    """python.exe next to pythonw.exe: a windowless Python has no stdout for
    the load tests."""
    exe = sys.executable or "python"
    if IS_WIN and os.path.basename(exe).lower() == "pythonw.exe":
        p = os.path.join(os.path.dirname(exe), "python.exe")
        if os.path.isfile(p):
            return p
    return exe


# ───────────────────────────────────────────────────────── compiler probing
def _win_dirs():
    """Folders where Windows C compilers usually live but are not on PATH."""
    la = os.environ.get("LOCALAPPDATA", "")
    up = os.environ.get("USERPROFILE", "")
    d = [r"C:\msys64\ucrt64\bin", r"C:\msys64\mingw64\bin", r"C:\msys64\clang64\bin", r"C:\mingw64\bin",
         os.path.join(HERE, "tools", "mingw64", "bin")]
    if la:
        d.append(os.path.join(la, "Microsoft", "WinGet", "Links"))
    if up:
        d.append(os.path.join(up, "scoop", "shims"))
    return d


def target_ok(machine, win=None):
    """None when a compiler's target (gcc -dumpmachine) can build a core for
    this Python, else the reason it cannot."""
    m = (machine or "").strip().lower()
    if not (IS_WIN if win is None else win):
        return None
    if not m:
        return "it does not report its target (-dumpmachine)"
    if not m.startswith("x86_64"):
        return f"it builds for {m}, not 64-bit x86 (x86_64)"
    if "msys" in m or "cygwin" in m:
        return f"it is an MSYS/Cygwin compiler ({m}): its DLLs need msys-2.0.dll/cygwin1.dll - use the MinGW-w64 (UCRT64) gcc"
    if "windows-msvc" in m:
        return f"it targets the MSVC libraries ({m}) - use a MinGW-w64 gcc or clang"
    if "mingw" in m or "windows-gnu" in m:
        return None
    return f"unknown target {m}"


def _dumpmachine(cc):
    try:
        return subprocess.run([cc, "-dumpmachine"], capture_output=True, text=True, timeout=30,
                              **_run_kw(cc)).stdout.strip()
    except Exception:
        return ""


REJECTED = {}            # compiler -> why it was not used (for --check and build.log)


def compilers():
    """C compilers to try, best first."""
    c = []
    if os.environ.get("CC"):
        c.append(os.environ["CC"])
    if IS_MAC:   # Homebrew gcc has OpenMP built in; Apple clang needs libomp
        c += ["clang", "cc"] + [f"gcc-{v}" for v in range(16, 10, -1)]
    elif IS_WIN:
        c += ["gcc", "x86_64-w64-mingw32-gcc", "clang", "cc"]
    else:
        c += ["gcc", "cc", "clang"]
    found = []
    for x in c:
        p = shutil.which(x)
        if p and p not in found:
            found.append(p)
    if IS_WIN:
        for d in _win_dirs():
            for x in ("gcc.exe", "x86_64-w64-mingw32-gcc.exe", "clang.exe"):
                p = os.path.join(d, x)
                if os.path.isfile(p) and p not in found:
                    found.append(p)
        out = []
        for p in found:
            why = target_ok(_dumpmachine(p))
            if why is None:
                out.append(p)
            else:
                REJECTED[p] = why
        return out
    return found


def _is_gcc(cc):
    try:
        v = subprocess.run([cc, "--version"], capture_output=True, text=True, timeout=30, **_run_kw(cc)).stdout
    except Exception:
        return False
    return ("Free Software Foundation" in v or "gcc" in os.path.basename(cc)) and "clang" not in v.lower()


def _libomp_flags():
    """OpenMP flag sets for Apple clang + libomp (Homebrew / MacPorts)."""
    pre = []
    if shutil.which("brew"):
        try:
            p = subprocess.run(["brew", "--prefix", "libomp"], capture_output=True, text=True,
                               timeout=60).stdout.strip()
            if p:
                pre.append(p)
        except Exception:
            pass
    for p in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"):
        if p not in pre:
            pre.append(p)
    sets = []
    for p in pre:
        if os.path.isfile(os.path.join(p, "include", "omp.h")):
            sets.append((["-Xpreprocessor", "-fopenmp", f"-I{p}/include"],
                         [f"-L{p}/lib", "-lomp", f"-Wl,-rpath,{p}/lib"]))
    if os.path.isfile("/opt/local/include/libomp/omp.h"):          # MacPorts
        sets.append((["-Xpreprocessor", "-fopenmp", "-I/opt/local/include/libomp"],
                     ["-L/opt/local/lib/libomp", "-lomp", "-Wl,-rpath,/opt/local/lib/libomp"]))
    return sets


def _omp_sets(cc):
    """(compile flags, link flags) that switch OpenMP on, best first."""
    if IS_MAC and not _is_gcc(cc):
        return _libomp_flags() + [(["-fopenmp"], ["-fopenmp"])]   # LLVM clang from Homebrew
    return [(["-fopenmp"], ["-fopenmp"])]


def python_arch_problem():
    """None, or why no core can be built/loaded for this Python."""
    if struct.calcsize("P") != 8:
        return "this is a 32-bit Python; SemiSim 3D needs a 64-bit Python (python.org: 'Windows installer (64-bit)')"
    if IS_WIN and platform.machine().lower() in ("arm64", "aarch64"):
        return ("this Python is built for ARM64 Windows; install the x64 Python from python.org "
                "(it runs on ARM Windows too)")
    return None


def _arch_sets(native, cc):
    m = platform.machine().lower()
    arm = m in ("arm64", "aarch64")
    # build for the architecture of THIS Python (an x86-64 Python under Rosetta
    # needs an x86-64 core); GCC on macOS has no -arch option
    base = (["-arch", "arm64" if arm else "x86_64"] if IS_MAC and not _is_gcc(cc) else [])
    if native:
        tune = ["-mcpu=native"] if arm else ["-march=native"]
        return [base + tune, base]
    if m in ("x86_64", "amd64") and not IS_MAC:
        return [base + ["-march=x86-64-v2"], base]
    return [base]


PROBE = r"""
#ifdef _OPENMP
#include <omp.h>
#endif
int semisim_probe(void){
#ifdef _OPENMP
    int n=0;
    #pragma omp parallel
    {
        #pragma omp atomic
        n++;
    }
    return n>0 ? 1000+n : -1;
#else
    return 1;
#endif
}
"""


def opt_flags(win=None):
    """Compile flags.  Windows: no -fPIC (all DLL code is relocatable)."""
    return (["-O3", "-ffast-math", "-fno-finite-math-only", "-shared"] if (IS_WIN if win is None else win)
            else ["-O3", "-ffast-math", "-fno-finite-math-only", "-fPIC", "-shared"])


def link_extra(win=None):
    """Extra link flags.  Windows: static libgomp, winpthreads and libgcc
    inside the DLL, so it needs nothing but Windows' own DLLs."""
    return ["-static"] if (IS_WIN if win is None else win) else []


OPT = opt_flags()


def _try(cmd, log, cc=None):
    log.write("$ " + " ".join(cmd) + "\n"); log.flush()
    try:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=900, **_run_kw(cc or cmd[0]))
        return r.returncode == 0
    except Exception as e:
        log.write(f"   {e}\n")
        return False


def _probe(cc, arch, omp, tmp, log):
    """Compile and load a tiny OpenMP test library; returns 1000+threads
    with OpenMP, 1 without, None when the flags do not work."""
    src = os.path.join(tmp, "probe.c")
    with open(src, "w", encoding="utf-8") as f:
        f.write(PROBE)
    so = os.path.join(tmp, f"probe{abs(hash((cc, tuple(arch), tuple(omp[0])))) % 10**8}"
                           f"{os.path.splitext(NAME)[1]}")
    if not _try([cc] + OPT + arch + omp[0] + ["-o", so, src] + omp[1] + link_extra(), log, cc):
        return None
    try:
        r = subprocess.run([_python_console(), "-c",
                            "import ctypes,sys; print(ctypes.CDLL(sys.argv[1]).semisim_probe())", so],
                           capture_output=True, text=True, timeout=120, **_run_kw())
    except Exception as e:
        log.write(f"   {e}\n")
        return None
    log.write(r.stdout + r.stderr)
    try:
        return int(r.stdout.strip().splitlines()[-1])
    except Exception:
        return None


def _loads(path, log):
    try:
        r = subprocess.run([_python_console(), "-c",
                            "import ctypes,sys; L=ctypes.CDLL(sys.argv[1]); L.api3_init(); "
                            f"assert hasattr(L,'{SENTINEL}'); print('ok')", path],
                           capture_output=True, text=True, timeout=120, **_run_kw())
    except Exception as e:
        log.write(f"   {e}\n")
        return False
    log.write(r.stdout + r.stderr)
    return r.returncode == 0 and "ok" in r.stdout


def _replace(src, dst):
    """Copy the new core into place.  On Windows a loaded DLL cannot be
    replaced, and antivirus programs hold new files for a moment: retry."""
    shutil.copyfile(src, dst + ".new")
    for k in range(6 if IS_WIN else 1):
        try:
            os.replace(dst + ".new", dst)
            return
        except OSError:
            if k == (5 if IS_WIN else 0):
                try:
                    os.remove(dst + ".new")
                except OSError:
                    pass
                raise CoreError(f"{dst} is in use - close the other SemiSim 3D windows and start again")
            time.sleep(0.3)


# ───────────────────────────────────────────────────────── build
def build(native=True, out=None, verbose=True):
    """Compile semiconductor_3d.c.  Returns the stamp dict, or None."""
    out = out or lib_path()
    if not os.path.isfile(SRC):
        raise CoreError(f"{SRC} is missing")
    why = python_arch_problem()
    if why:
        raise CoreError(f"No solver core can be built: {why}.")
    logp = os.path.join(out_dir(), "build.log")
    say = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
    tmp = tempfile.mkdtemp(prefix="semisim_build_")
    try:
        with open(logp, "w", encoding="utf-8", errors="replace") as log:
            REJECTED.clear()
            ccs = compilers()
            for p, w in REJECTED.items():
                log.write(f"not used: {p} - {w}\n")
            if not ccs:
                log.write("no C compiler found\n")
                return None
            # best working (compiler, arch, OpenMP) combination, tested on a tiny
            # library: the first compiler that does OpenMP, else the first that works
            best = None
            for cc in ccs:
                found = None
                for arch in _arch_sets(native, cc):              # best tuning first
                    for omp in _omp_sets(cc) + [([], [])]:       # OpenMP first
                        r = _probe(cc, arch, omp, tmp, log)
                        if r is not None:
                            found = (cc, arch, omp, r > 1000)
                            break
                    if found:
                        break
                if found and (best is None or (found[3] and not best[3])):
                    best = found
                if best is not None and best[3]:
                    break
            if best is None:
                return None
            cc, arch, omp, has_omp = best
            tries = [(omp, has_omp)] + ([(([], []), False)] if has_omp else [])
            tmpout = os.path.join(tmp, NAME)
            for (cf, lf), om in tries:
                flags = OPT + arch + cf
                say(f"   {os.path.basename(cc)} {' '.join(arch + cf) or '(default target)'}"
                    f"{'' if om else '   [no OpenMP: single-threaded]'}")
                lf = lf + link_extra()
                if _try([cc] + flags + ["-o", tmpout, SRC] + lf + ["-lm"], log, cc) and _loads(tmpout, log):
                    rep = dll_report(tmpout) if IS_WIN else {}
                    if rep.get("bad"):
                        log.write(f"the core still needs {', '.join(rep['bad'])}: not used\n")
                        continue
                    _replace(tmpout, out)
                    st = dict(src=src_hash(), cc=cc, flags=flags + lf, openmp=om, native=native,
                              system=platform.system(), machine=platform.machine().lower(),
                              python=sys.version.split()[0], prebuilt=False)
                    if IS_WIN:
                        st["crt"] = rep.get("crt")
                    with open(stamp_path(os.path.dirname(os.path.abspath(out))), "w", encoding="utf-8") as f:
                        json.dump(st, f, indent=1)
                    return st
    finally:
        shutil.rmtree(tmp, ignore_errors=True)      # a file still being scanned must not stop us
    return None


def cross_build(cc, out, verbose=True):
    """Compile a portable Windows core with a MinGW-w64 cross-compiler (no
    load test: the host cannot run it) and write its stamp next to it."""
    if not os.path.isfile(SRC):
        raise CoreError(f"{SRC} is missing")
    why = target_ok(_dumpmachine(cc), win=True)
    if why:
        raise CoreError(f"{cc} cannot build a Windows core: {why}")
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    flags = opt_flags(win=True) + ["-march=x86-64-v2", "-fopenmp"]
    lf = ["-fopenmp"] + link_extra(win=True)
    logp = os.path.join(os.path.dirname(out), "build.log")
    with open(logp, "w", encoding="utf-8", errors="replace") as log:
        if not _try([cc] + flags + ["-o", out, SRC] + lf + ["-lm"], log, cc):
            raise CoreError(f"the cross build failed - see {logp}")
    rep = dll_report(out)
    if rep["bad"] or not rep["sentinel"]:
        raise CoreError(f"the cross-built core is not usable: needs {rep['bad']}, entry point "
                        f"{'found' if rep['sentinel'] else 'missing'}")
    st = dict(src=src_hash(), cc=cc, flags=flags + lf, openmp=True, native=False, system="Windows",
              machine="amd64", python="any 64-bit", prebuilt=True, crt=rep["crt"])
    with open(stamp_path(os.path.dirname(out)), "w", encoding="utf-8") as f:
        json.dump(st, f, indent=1)
    if verbose:
        print(f"built {out}  (C runtime: {rep['crt']}, OpenMP: yes, portable)")
    return st


# ───────────────────────────────────────────────────────── load
_HINT_LINUX = ("Install a C compiler and try again:\n"
               "    sudo apt install gcc        (Ubuntu / Debian / WSL)\n"
               "    sudo dnf install gcc        (Fedora)\n"
               "then run:  python3 semibuild.py")
_HINT_MAC = ("Install Apple's command line tools (C compiler) and try again:\n"
             "    xcode-select --install\n"
             "optional, for multi-core solving:  brew install libomp\n"
             "then run:  python3 semibuild.py   (or: bash install.sh)")
_HINT_WIN = ("Two ways to get the solver core:\n"
             "  1. download the Windows zip of this version (GitHub -> Releases); it contains\n"
             "     semiconductor_core.dll - unzip it over this folder and start again, or\n"
             "  2. install a C compiler and build it here:\n"
             "        winget install BrechtSanders.WinLibs.POSIX.UCRT\n"
             "     then open a NEW terminal in this folder and run:  py semibuild.py")


def hint():
    return _HINT_WIN if IS_WIN else _HINT_MAC if IS_MAC else _HINT_LINUX


def _dlclose(lib):
    try:
        import _ctypes
        (_ctypes.dlclose if hasattr(_ctypes, "dlclose") else _ctypes.FreeLibrary)(lib._handle)
    except Exception:
        pass


def load_core(verbose=True, progress=None):
    """The loaded core (ctypes.CDLL) and its path; builds it when needed.
    progress: optional callable(message) told when a build starts (a window
    without a console shows it - the build takes 10-60 s)."""
    env = os.environ.get("SEMISIM_SO")
    if env:
        lib = ctypes.CDLL(env)
        if not hasattr(lib, SENTINEL):
            raise CoreError(f"SEMISIM_SO points to a solver core of another version of SemiSim 3D:\n"
                            f"    {env}\nRemove the variable (or point it to a core built from this version).")
        return lib, env
    why = python_arch_problem()
    if why:
        raise CoreError(f"SemiSim 3D cannot run on this Python: {why}.")
    path = lib_path()
    st = read_stamp()
    have_src = os.path.isfile(SRC)
    reason, loaded = None, False
    if not os.path.isfile(path):
        reason = "first start: the solver core is not built yet"
    elif have_src and st is not None and st.get("src") != src_hash():
        reason = "semiconductor_3d.c has changed since the last build"
    elif st is not None and (st.get("system"), st.get("machine")) != _machine():
        reason = "the core was built on another kind of computer"
    if reason is None:
        try:
            lib = ctypes.CDLL(path); loaded = True
            if hasattr(lib, SENTINEL):
                return lib, path
            _dlclose(lib)
            reason = "the core was built from an older version"
        except OSError as e:
            reason = f"the core cannot be loaded on this computer ({e})"
    if have_src and compilers():
        if verbose:
            print(f"SemiSim 3D: building the solver core - {reason}.\n"
                  f"   This takes 10-60 s and happens once.", flush=True)
        if progress is not None:
            try:
                progress(f"Building the solver core ({reason}).\nThis takes 10-60 s and happens once.")
            except Exception:
                pass
        st = build(native=True, verbose=verbose)
        if st:
            if verbose:
                print("   done" + ("" if st["openmp"] else
                                   " (single-threaded: " + ("brew install libomp, then python3 semibuild.py"
                                                            if IS_MAC else "no OpenMP in this compiler") + ")"),
                      flush=True)
            if loaded:             # the old file may still be mapped: load a private copy
                fd, cp = tempfile.mkstemp(suffix=os.path.splitext(NAME)[1], prefix="semisim_core_"); os.close(fd)
                shutil.copyfile(path, cp)
                return ctypes.CDLL(cp), cp
            return ctypes.CDLL(path), path
        msg = f"The solver core could not be built ({reason}); see {os.path.join(out_dir(), 'build.log')}."
    elif not have_src:
        msg = f"{path} is missing or unusable and semiconductor_3d.c is not here to build it ({reason})."
    else:
        msg = f"No C compiler found to build the solver core ({reason})."
        if REJECTED:
            msg += "\nCompilers found but not usable:\n" + "\n".join(f"  {p}: {w}" for p, w in REJECTED.items())
    raise CoreError(msg + "\n\n" + hint())


# ───────────────────────────────────────────────────────── command line
def _utf8_streams():
    """A Windows console (or a redirected one) may use cp1250/cp1252: print
    anything without failing."""
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


def main(argv):
    import argparse
    _utf8_streams()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--portable", action="store_true", help="no CPU-specific instructions")
    ap.add_argument("--check", action="store_true", help="only report the current state")
    ap.add_argument("--cross", metavar="CC", default=None,
                    help="compile a Windows core with this MinGW-w64 cross-compiler (needs -o)")
    ap.add_argument("-o", "--output", default=None, help=f"output file (default: {NAME})")
    a = ap.parse_args(argv)
    if a.check:
        st = read_stamp()
        REJECTED.clear()
        ccs = compilers()
        print(f"source     : {SRC} ({src_hash() if os.path.isfile(SRC) else 'missing'})")
        print(f"core       : {lib_path()} ({'present' if os.path.isfile(lib_path()) else 'missing'})")
        print(f"build info : {json.dumps(st) if st else 'none'}")
        print(f"compilers  : {', '.join(ccs) or 'none'}")
        for p, w in REJECTED.items():
            print(f"   not used: {p} - {w}")
        if IS_WIN and os.path.isfile(lib_path()):
            rep = dll_report(lib_path())
            print(f"C runtime  : {rep['crt'] or 'unknown'}"
                  + ("" if rep["crt"] != "msvcrt" else "   (works; libgomp then ignores GOMP_SPINCOUNT)"))
            if rep["bad"]:
                print(f"!! the core needs {', '.join(rep['bad'])} - rebuild it (py semibuild.py)")
        if python_arch_problem():
            print(f"!! {python_arch_problem()}")
        if IS_MAC:
            print(f"libomp     : {', '.join(f[0][2][2:] for f in _libomp_flags()) or 'not found (brew install libomp)'}")
        return 0
    if a.cross:
        if not a.output:
            print("--cross needs -o PATH/semiconductor_core.dll"); return 2
        try:
            cross_build(a.cross, a.output)
        except CoreError as e:
            print(f"!! {e}"); return 1
        return 0
    print(f"Building the SemiSim 3D core ({'portable' if a.portable else 'tuned for this CPU'}) ...", flush=True)
    try:
        st = build(native=not a.portable, out=a.output)
    except CoreError as e:
        print(f"!! {e}"); return 1
    if not st:
        print(f"!! build failed - see {os.path.join(out_dir(), 'build.log')}\n\n" + hint())
        return 1
    print(f"built {a.output or lib_path()}  (OpenMP: {'yes' if st['openmp'] else 'no - single-threaded'})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
