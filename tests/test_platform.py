#!/usr/bin/env python3
"""Platform test: the build-and-load logic of semibuild.py for Linux, macOS
and Windows, without a display and in seconds.

  - the source hash ignores line endings (a Windows checkout has CRLF);
  - the Windows rules, checked on any system with the platform injected:
    file name, flags, which compilers are accepted, the hint, the cache folder;
  - with a MinGW-w64 cross-compiler installed: a Windows core is built and
    needs nothing but Windows' own DLLs;
  - the core loads, and an environment knob set from Python reaches it
    (on Windows that needs the core's own env_get, not getenv);
  - on Windows itself: `semibuild.py --check` names the C runtime.

    python3 tests/test_platform.py
"""
import os, sys, json, shutil, ctypes, tempfile, subprocess

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import semibuild as B  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"   ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


print("1. source hash and line endings")
tmp = tempfile.mkdtemp(prefix="semisim_plat_")
with open(B.SRC, "rb") as f:
    raw = f.read()
lf, crlf = os.path.join(tmp, "lf.c"), os.path.join(tmp, "crlf.c")
with open(lf, "wb") as f:
    f.write(raw.replace(b"\r\n", b"\n"))
with open(crlf, "wb") as f:
    f.write(raw.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
check("CRLF and LF copies hash the same", B.src_hash(lf) == B.src_hash(crlf), B.src_hash(lf))
check("the hash of the LF source is the plain SHA-256 prefix",
      B.src_hash(lf) == __import__("hashlib").sha256(raw.replace(b"\r\n", b"\n")).hexdigest()[:16])

print("2. Windows rules (platform injected)")
check("core file name", B.lib_name(win=True) == "semiconductor_core.dll" and B.lib_name(win=False) == "semiconductor_core.so")
check("no -fPIC on Windows, kept elsewhere", "-fPIC" not in B.opt_flags(win=True) and "-fPIC" in B.opt_flags(win=False))
check("static runtime on Windows only", B.link_extra(win=True) == ["-static"] and B.link_extra(win=False) == [])
cases = {"x86_64-w64-mingw32": None, "x86_64-w64-windows-gnu": None,
         "x86_64-pc-msys": "msys", "x86_64-pc-cygwin": "cygwin", "i686-w64-mingw32": "x86_64",
         "x86_64-pc-windows-msvc": "MSVC", "aarch64-w64-mingw32": "x86_64", "": "target"}
for m, want in cases.items():
    got = B.target_ok(m, win=True)
    ok = (got is None) if want is None else (got is not None and want.lower() in got.lower())
    check(f"compiler target '{m or '(none)'}' {'accepted' if want is None else 'rejected'}", ok, got or "")
check("any target is fine off Windows", B.target_ok("x86_64-linux-gnu", win=False) is None)
old = os.environ.get("LOCALAPPDATA")
os.environ["LOCALAPPDATA"] = os.path.join(tmp, "AppData", "Local")
check("cache folder on Windows is %LOCALAPPDATA%\\SemiSim3D",
      B.cache_dir(win=True) == os.path.join(tmp, "AppData", "Local", "SemiSim3D"))
if old is None:
    del os.environ["LOCALAPPDATA"]
else:
    os.environ["LOCALAPPDATA"] = old
_w = B.IS_WIN
B.IS_WIN = True
check("the Windows hint names the release zip and winget", "Releases" in B.hint() and "winget" in B.hint())
B.IS_WIN = _w
rep = B.dll_report(B.SRC)            # not a DLL: nothing found, no crash
check("dll_report on a non-DLL", rep["crt"] is None and rep["bad"] == [])

print("3. Windows cross build")
cross = next((c for c in ("x86_64-w64-mingw32-gcc-posix", "x86_64-w64-mingw32-gcc") if shutil.which(c)), None)
if cross is None:
    print("   (no MinGW-w64 cross-compiler here - skipped; Ubuntu: sudo apt install gcc-mingw-w64-x86-64-posix)")
else:
    out = os.path.join(tmp, "dist", "semiconductor_core.dll")
    try:
        st = B.cross_build(cross, out, verbose=False)
        r = B.dll_report(out)
        check("cross-built core has the newest entry point", r["sentinel"])
        check("... needs no libgomp / winpthread / libgcc DLL", r["bad"] == [], ", ".join(r["bad"]))
        check("... its C runtime is known", r["crt"] in ("ucrt", "msvcrt"), r["crt"] or "")
        with open(B.stamp_path(os.path.dirname(out)), encoding="utf-8") as f:
            s2 = json.load(f)
        check("... its stamp: prebuilt, Windows, this source",
              s2.get("prebuilt") is True and s2.get("system") == "Windows" and s2.get("src") == B.src_hash(), "")
        check("... and the stamp carries the C runtime", s2.get("crt") == r["crt"])
    except B.CoreError as e:
        check("cross build", False, str(e))

print("4. the core on this computer")
os.environ["SEMISIM_AAM"] = "3"          # a knob read by api3_init: Anderson depth
try:
    lib, path = B.load_core(verbose=False)
    check("the core loads", hasattr(lib, B.SENTINEL), os.path.basename(path))
    lib.api3_init()
    lib.api3_set_grid.argtypes = [ctypes.c_int] * 3 + [ctypes.c_double] * 3
    lib.api3_set_grid(6, 6, 2, 1e-6, 1e-6, 1e-7)
    lib.api3_get_aa_depth.restype = ctypes.c_int
    d = lib.api3_get_aa_depth()
    check("an environment knob set from Python reaches the core (SEMISIM_AAM=3)", d == 3, f"depth {d}")
    lib.api3_get_stat.restype = ctypes.c_longlong
    lib.api3_get_stat.argtypes = [ctypes.c_int]
    check("the iteration counter is a 64-bit integer", lib.api3_get_stat(3) >= 0)
finally:
    del os.environ["SEMISIM_AAM"]
if B.IS_WIN:
    r = subprocess.run([sys.executable, os.path.join(HERE, "semibuild.py"), "--check"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
    check("semibuild --check names the C runtime", "C runtime" in r.stdout, r.stdout.strip().splitlines()[-1] if r.stdout else "")
    check("... and the core needs no extra DLL", "!! the core needs" not in r.stdout)

shutil.rmtree(tmp, ignore_errors=True)
print(f"\n{'all passed' if not FAILS else str(len(FAILS)) + ' FAILED: ' + '; '.join(FAILS)}")
sys.exit(1 if FAILS else 0)
