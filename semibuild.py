#!/usr/bin/env python3
"""Build and load the SemiSim 3D solver core (semiconductor_3d.c) on Linux,
WSL2 and macOS.

    python3 semibuild.py              build for this computer (fastest code)
    python3 semibuild.py --portable   build for any CPU of this architecture
    python3 semibuild.py --check      show the compiler, OpenMP and the core in use

main.py calls load_core(), which builds the core by itself when it is
missing, was built from another version of semiconductor_3d.c, or cannot be
loaded on this computer (a Linux build copied to a Mac, for example).

Needs a C compiler: gcc or clang on Linux (Ubuntu: sudo apt install gcc),
the Xcode command line tools on macOS (xcode-select --install).  OpenMP
(multi-core solving) comes with gcc; on macOS it needs Homebrew's libomp
(brew install libomp) - without it the core is built single-threaded.
"""
import os, sys, json, ctypes, hashlib, platform, shutil, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "semiconductor_3d.c")
NAME = "semiconductor_core.so"          # also on macOS: ctypes loads any file name
SENTINEL = "api3_set_state"             # newest entry point; an older core lacks it
IS_MAC = sys.platform == "darwin"


class CoreError(RuntimeError):
    pass


# ───────────────────────────────────────────────────────── paths, stamp
def _writable(d):
    try:
        fd, p = tempfile.mkstemp(dir=d, prefix=".wtest"); os.close(fd); os.remove(p); return True
    except Exception:
        return False


def out_dir():
    """Where the built core lives: next to the sources, or a user cache
    folder when the program folder is read-only."""
    if _writable(HERE):
        return HERE
    d = os.path.join(os.path.expanduser("~"), ".cache", "semisim3d")
    os.makedirs(d, exist_ok=True)
    return d


def lib_path():
    return os.path.join(out_dir(), NAME)


def stamp_path():
    return os.path.join(out_dir(), ".core_build.json")


def src_hash():
    with open(SRC, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def read_stamp():
    try:
        with open(stamp_path()) as f:
            d = json.load(f)
        return d if isinstance(d, dict) else None
    except Exception:
        return None


def _machine():
    return platform.system(), platform.machine().lower()


# ───────────────────────────────────────────────────────── compiler probing
def compilers():
    """C compilers to try, best first."""
    c = []
    if os.environ.get("CC"):
        c.append(os.environ["CC"])
    if IS_MAC:   # Homebrew gcc has OpenMP built in; Apple clang needs libomp
        c += ["clang", "cc"] + [f"gcc-{v}" for v in range(16, 10, -1)]
    else:
        c += ["gcc", "cc", "clang"]
    out = []
    for x in c:
        p = shutil.which(x)
        if p and p not in out:
            out.append(p)
    return out


def _is_gcc(cc):
    try:
        v = subprocess.run([cc, "--version"], capture_output=True, text=True, timeout=30).stdout
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

OPT = ["-O3", "-ffast-math", "-fno-finite-math-only", "-fPIC", "-shared"]


def _try(cmd, log):
    log.write("$ " + " ".join(cmd) + "\n"); log.flush()
    try:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=900)
        return r.returncode == 0
    except Exception as e:
        log.write(f"   {e}\n")
        return False


def _probe(cc, arch, omp, tmp, log):
    """Compile and load a tiny OpenMP test library; returns 1000+threads
    with OpenMP, 1 without, None when the flags do not work."""
    src = os.path.join(tmp, "probe.c")
    with open(src, "w") as f:
        f.write(PROBE)
    so = os.path.join(tmp, f"probe{abs(hash((cc, tuple(arch), tuple(omp[0])))) % 10**8}.so")
    if not _try([cc] + OPT + arch + omp[0] + ["-o", so, src] + omp[1], log):
        return None
    r = subprocess.run([sys.executable, "-c",
                        "import ctypes,sys; print(ctypes.CDLL(sys.argv[1]).semisim_probe())", so],
                       capture_output=True, text=True, timeout=120)
    log.write(r.stdout + r.stderr)
    try:
        return int(r.stdout.strip().splitlines()[-1])
    except Exception:
        return None


def _loads(path, log):
    r = subprocess.run([sys.executable, "-c",
                        "import ctypes,sys; L=ctypes.CDLL(sys.argv[1]); L.api3_init(); "
                        f"assert hasattr(L,'{SENTINEL}'); print('ok')", path],
                       capture_output=True, text=True, timeout=120)
    log.write(r.stdout + r.stderr)
    return r.returncode == 0 and "ok" in r.stdout


# ───────────────────────────────────────────────────────── build
def build(native=True, out=None, verbose=True):
    """Compile semiconductor_3d.c.  Returns the stamp dict, or None."""
    out = out or lib_path()
    if not os.path.isfile(SRC):
        raise CoreError(f"{SRC} is missing")
    logp = os.path.join(out_dir(), "build.log")
    say = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
    with open(logp, "w") as log, tempfile.TemporaryDirectory() as tmp:
        ccs = compilers()
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
            if _try([cc] + flags + ["-o", tmpout, SRC] + lf + ["-lm"], log) and _loads(tmpout, log):
                shutil.copyfile(tmpout, out + ".new")
                os.replace(out + ".new", out)
                st = dict(src=src_hash(), cc=cc, flags=flags + lf, openmp=om, native=native,
                          system=platform.system(), machine=platform.machine().lower(),
                          python=sys.version.split()[0])
                with open(stamp_path(), "w") as f:
                    json.dump(st, f, indent=1)
                return st
    return None


# ───────────────────────────────────────────────────────── load
_HINT_LINUX = ("Install a C compiler and try again:\n"
               "    sudo apt install gcc        (Ubuntu / Debian / WSL)\n"
               "    sudo dnf install gcc        (Fedora)\n"
               "then run:  python3 semibuild.py")
_HINT_MAC = ("Install Apple's command line tools (C compiler) and try again:\n"
             "    xcode-select --install\n"
             "optional, for multi-core solving:  brew install libomp\n"
             "then run:  python3 semibuild.py   (or: bash install.sh)")


def _dlclose(lib):
    try:
        import _ctypes
        (_ctypes.dlclose if hasattr(_ctypes, "dlclose") else _ctypes.FreeLibrary)(lib._handle)
    except Exception:
        pass


def load_core(verbose=True):
    """The loaded core (ctypes.CDLL) and its path; builds it when needed."""
    env = os.environ.get("SEMISIM_SO")
    if env:
        return ctypes.CDLL(env), env
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
        st = build(native=True, verbose=verbose)
        if st:
            if verbose:
                print("   done" + ("" if st["openmp"] else
                                   " (single-threaded: " + ("brew install libomp, then python3 semibuild.py"
                                                            if IS_MAC else "no OpenMP in this compiler") + ")"),
                      flush=True)
            if loaded:             # the old file may still be mapped: load a private copy
                fd, cp = tempfile.mkstemp(suffix=".so", prefix="semisim_core_"); os.close(fd)
                shutil.copyfile(path, cp)
                return ctypes.CDLL(cp), cp
            return ctypes.CDLL(path), path
        msg = f"The solver core could not be built ({reason}); see {os.path.join(out_dir(), 'build.log')}."
    elif not have_src:
        msg = f"{path} is missing or unusable and semiconductor_3d.c is not here to build it ({reason})."
    else:
        msg = f"No C compiler found to build the solver core ({reason})."
    raise CoreError(msg + "\n\n" + (_HINT_MAC if IS_MAC else _HINT_LINUX))


# ───────────────────────────────────────────────────────── command line
def main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--portable", action="store_true", help="no CPU-specific instructions")
    ap.add_argument("--check", action="store_true", help="only report the current state")
    ap.add_argument("-o", "--output", default=None, help="output file (default: semiconductor_core.so)")
    a = ap.parse_args(argv)
    if a.check:
        st = read_stamp()
        print(f"source     : {SRC} ({src_hash() if os.path.isfile(SRC) else 'missing'})")
        print(f"core       : {lib_path()} ({'present' if os.path.isfile(lib_path()) else 'missing'})")
        print(f"build info : {json.dumps(st) if st else 'none'}")
        print(f"compilers  : {', '.join(compilers()) or 'none'}")
        if IS_MAC:
            print(f"libomp     : {', '.join(f[0][2][2:] for f in _libomp_flags()) or 'not found (brew install libomp)'}")
        return 0
    print(f"Building the SemiSim 3D core ({'portable' if a.portable else 'tuned for this CPU'}) ...", flush=True)
    st = build(native=not a.portable, out=a.output)
    if not st:
        print(f"!! build failed - see {os.path.join(out_dir(), 'build.log')}\n\n"
              + (_HINT_MAC if IS_MAC else _HINT_LINUX))
        return 1
    print(f"built {a.output or lib_path()}  (OpenMP: {'yes' if st['openmp'] else 'no - single-threaded'})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
