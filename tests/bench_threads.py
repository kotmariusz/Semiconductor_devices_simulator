#!/usr/bin/env python3
"""How many solver threads are fastest on this computer?

Times one reference job - the drain sweep 0-20 V of the 55 V VDMOS cell
(`RW: Si MOSFET 55V (IRLZ44N)`, gate on, 6 points) - with 1, 2, 4, ... threads
up to all logical CPUs, and prints a table.  SemiSim uses all logical CPUs
by default (Apple-silicon Macs: the performance cores); if "all" is clearly
slower than the best count here, please send the table to the author.

    python3 tests/bench_threads.py          (Windows: py tests\\bench_threads.py)
Takes a few minutes; keep other programs closed.
"""
import os, sys, time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
os.environ.setdefault("MPLBACKEND", "Agg")
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import main as M  # noqa: E402

NAME = "RW: Si MOSFET 55V (IRLZ44N)"


def job(nthreads):
    M.api_threads(nthreads)
    t = M.TEMPLATES[NAME](); cons = t["cons"]
    t0 = time.perf_counter()
    M.setup_device(t, models=t.get("models", M.M_DEFAULT))
    gate = [i for i, c in enumerate(cons) if c.bc == 2]
    drain = [i for i, c in enumerate(cons) if c.label.lower().startswith("drain")]
    va = [c.V for c in cons]
    for i in gate:
        va[i] = 10.
    R = M.run_sweep(cons, drain, 0., 20., 6, vapp=va, refine=False)
    ok = all(R.get("conv", [True]))
    return time.perf_counter() - t0, ok


def main():
    nl = os.cpu_count() or 1
    counts = sorted({1, 2, 4, 6, 8, 12, 16, nl // 2, nl} - {0})
    counts = [c for c in counts if c <= nl]
    print(f"{NAME}: drain sweep 0-20 V at V_GS = 10 V, 6 points; {nl} logical CPUs\n")
    print(" threads   time (s)   vs best")
    res = {}
    job(nl)                                    # warm-up (first-touch, caches)
    for n in counts:
        dt, ok = job(n)
        res[n] = dt
        print(f"  {n:5d}   {dt:8.1f}" + ("" if ok else "   (not all points converged)"), flush=True)
    best = min(res, key=res.get)
    print(f"\nfastest: {best} threads ({res[best]:.1f} s); all {nl}: {res[nl]:.1f} s "
          f"({100 * (res[nl] / res[best] - 1):+.0f} %)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
