#!/usr/bin/env python3
"""GUI smoke test: opens the main window, solves two templates, visits every
tab, changes the arrow density of the Current tab, toggles the editor panel,
fill-screen, the text size, Help and About, then quits.  Fails on any
exception or error dialog.  Needs a display (Linux CI: xvfb-run).

    python3 tests/gui_smoke.py            # as the platform draws it
    python3 tests/gui_smoke.py --flat     # with the macOS-style label buttons
"""
import os, sys, time, tempfile, traceback

if "--flat" in sys.argv:
    os.environ["SEMISIM_FLAT_BUTTONS"] = "1"
os.environ["HOME"] = tempfile.mkdtemp(prefix="semisim_home_")   # fresh settings file
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import main as M  # noqa: E402

ERRORS = []
M.messagebox.showerror = lambda *a, **k: ERRORS.append(("error dialog",) + a)
M.messagebox.showwarning = lambda *a, **k: print("warning dialog:", a, flush=True)
M.messagebox.showinfo = lambda *a, **k: print("info dialog:", a, flush=True)
M.messagebox.askyesno = lambda *a, **k: True

app = M.App()
_orig_report = app.report_callback_exception
def report(exc, val, tb):
    ERRORS.append(("callback exception", "".join(traceback.format_exception(exc, val, tb))))
    _orig_report(exc, val, tb)
app.report_callback_exception = report


def pump(t=0.3):
    end = time.time() + t
    while time.time() < end:
        app.update(); time.sleep(0.02)


def wait(tmax=300):
    t0 = time.time()
    while app.busy and time.time() - t0 < tmax:
        pump(0.3)
    assert not app.busy, "solver did not finish"
    pump(0.8)


def step(name, fn):
    t0 = time.time()
    fn()
    print(f"  ok  {name} ({time.time() - t0:.1f} s)", flush=True)


try:
    pump(1.0)
    print(f"windowing system: {app.tk.call('tk', 'windowingsystem')}, Tk {app.tk.call('info', 'patchlevel')}, "
          f"flat buttons: {M.FLAT_BUTTONS}, core: {M._SO}", flush=True)

    def diode():
        app._load("PN Diode"); pump(0.5); app._run(); wait()
        assert "CONVERGED" in app.sv.get() and "NOT" not in app.sv.get(), app.sv.get()
    step("PN diode solves", diode)

    def tabs():
        for i in range(len(app._TABS)):
            app.nb.select(i); pump(0.6)
    step("every plot tab draws", tabs)

    def arrows():
        app.nb.select(app.tc); pump(0.5)
        for n2, n3 in ((40, 16), (90, 30), (6, 2), (22, 8)):
            app.cur_n2d.set(n2); app.cur_n3d.set(n3); app._arrows_changed(); pump(0.6)
            app._upd_curr_plot(); pump(0.6)
            q3 = [a for a in app.fc.axes if getattr(a, "name", "") == "3d"]
            assert q3 and "arrows" in q3[0].get_title(), "3-D arrow title"
        app.cur_n2d.set(999); app._arrows_changed(); pump(0.4)
        assert app._arrow_n()[0] == M.ARROWS_2D[1], app._arrow_n()
        app.cur_n2d.set(22); app._arrows_changed(); pump(0.4)
    step("Current tab arrow density", arrows)

    def nmos_sweep():
        app._load("NMOS"); pump(0.5)
        app.ivN.set(6); app._sweep(); wait()
        app.nb.select(app.ti); pump(0.8)
    step("NMOS short sweep", nmos_sweep)

    def window():
        app._toggle_drawer(); pump(0.4); app._toggle_drawer(); pump(0.4)
        app._text_step(1); pump(0.4); app._text_step(None); pump(0.4)
        app._toggle_fill(); pump(1.0); app._toggle_fill(); pump(1.0)
        app._reset_window(); pump(0.6)
    step("drawer, text size, fill screen", window)

    def buttons():
        b = app.b_run   # a toolbar button: disabled while busy, invokable when idle
        assert str(b.cget("state")) in ("normal", "active"), b.cget("state")
        app._set_busy(True); pump(0.2)
        assert str(b.cget("state")) == "disabled"
        app._set_busy(False); pump(0.2)
    step("toolbar button states", buttons)

    def helpwin():
        app._help("Keyboard shortcuts"); pump(0.8)
        app._help("Sources"); pump(0.5)
        txt = app._help_text.get("1.0", "end")
        assert "Scharfetter" in txt and "Status and terms" in [t for t, _ in app._help_marks], "help topics"
        app._about(); pump(0.8)
    step("help (incl. status and sources) and about", helpwin)
except Exception:
    ERRORS.append(("exception", traceback.format_exc()))

try:
    app._quit()
except Exception:
    ERRORS.append(("quit", traceback.format_exc()))

if ERRORS:
    print("\nFAILED:")
    for e in ERRORS:
        print(" ", *e)
    sys.exit(1)
print("GUI smoke test passed")
