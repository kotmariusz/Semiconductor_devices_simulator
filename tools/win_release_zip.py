#!/usr/bin/env python3
"""Make the Windows release zip: the program (git archive of HEAD) plus a
ready-made portable solver core, so Windows users need no C compiler.

    python tools/win_release_zip.py CORE_DIR [OUT_DIR]

CORE_DIR holds semiconductor_core.dll and its .core_build.json, made by
    python semibuild.py --portable -o CORE_DIR/semiconductor_core.dll      (on Windows)
    python3 semibuild.py --cross x86_64-w64-mingw32-gcc-posix -o CORE_DIR/semiconductor_core.dll   (Linux/WSL)
The zip is OUT_DIR/SemiSim3D-<version>-windows-x64.zip (OUT_DIR: dist).  In
GitHub Actions its path is also written to $GITHUB_ENV as ZIP.
"""
import os, re, sys, json, zipfile, subprocess

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import semibuild as B  # noqa: E402


def main(argv):
    if not argv:
        print(__doc__); return 2
    core_dir = os.path.abspath(argv[0])
    out_dir = os.path.abspath(argv[1] if len(argv) > 1 else os.path.join(HERE, "dist"))
    dll = os.path.join(core_dir, "semiconductor_core.dll")
    st = B.read_stamp(core_dir)
    if not os.path.isfile(dll) or not st:
        print(f"!! {dll} or its .core_build.json is missing"); return 1
    if st.get("src") != B.src_hash():
        print(f"!! the core was built from another source ({st.get('src')} vs {B.src_hash()})"); return 1
    rep = B.dll_report(dll)
    if rep["bad"] or not rep["sentinel"]:
        print(f"!! the core is not usable: needs {rep['bad']}, entry point {rep['sentinel']}"); return 1
    st["prebuilt"] = True
    with open(os.path.join(HERE, "main.py"), encoding="utf-8") as f:
        ver = re.search(r'^VERSION="([^"]+)"', f.read(), re.M).group(1)
    name = f"SemiSim3D-{ver}-windows-x64"
    os.makedirs(out_dir, exist_ok=True)
    zp = os.path.join(out_dir, name + ".zip")
    subprocess.run(["git", "archive", "--format=zip", f"--prefix={name}/", "-o", zp, "HEAD"], cwd=HERE, check=True)
    with zipfile.ZipFile(zp, "a", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(dll, f"{name}/semiconductor_core.dll")
        z.writestr(f"{name}/.core_build.json", json.dumps(st, indent=1))
    print(f"{zp}  ({os.path.getsize(zp) / 1e6:.1f} MB, core: C runtime {rep['crt']}, source {st['src']})")
    if os.environ.get("GITHUB_ENV"):
        with open(os.environ["GITHUB_ENV"], "a", encoding="utf-8") as f:
            f.write(f"ZIP={zp}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
