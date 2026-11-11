#!/usr/bin/env python3
"""GUI smoke test: opens the main window, solves two templates, visits every
tab, changes the arrow density of the Current tab, hides and shows every
panel, uses the Layers window (per-part opacity) and the structure under the
current arrows, pins a sweep and runs an A/B compare of two template
variants, fills the parameter table (copy, CSV), toggles the editor panel,
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
import numpy as np  # noqa: E402

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

    def panels():
        redraw = {"td": app._ddev, "ts": app._upd_slice, "tb": app._upd_band_plot, "tc": app._upd_curr_plot,
                  "tf": app._upd_fi_plot, "ti": app._replot_iv}
        for tab, keys in M.App._PANELS.items():
            for k, _ in keys:                      # each panel off on its own, then all off, then all on
                app._panels_all(tab, True); app.pv[tab][k].set(False); redraw[tab](); pump(0.15)
            app._panels_all(tab, False); redraw[tab](); pump(0.15)
            app._panels_all(tab, True); redraw[tab](); pump(0.15)
        app._layers_win(); pump(0.3)
        for what in ("semi", "nocon", "glass", "all"):
            app._lay_set(what); pump(0.2)
        app.lay_ins3d.set(True); app.lay_sec.set(True); app._ddev(); pump(0.2)
        tv = app._lay_tv; parts = [i for i in tv.get_children() if i != "d"]
        tv.selection_set(parts[:2]); pump(0.2)                 # own opacity for two parts
        app.lay_sa.set(0.8); app._lay_sel_alpha(); pump(0.5)
        assert len(app.lay_alpha) == 2 and "0.80" in tv.set(parts[0], "opac"), tv.item(parts[0])
        app.nb.select(app.tc); pump(0.3)                       # the structure under the current arrows
        for bg, st3 in (("|J| + structure", True), ("structure", True), ("|J|", False)):
            app.cur_bg.set(bg); app.cur_st3.set(st3); app.cur_sa.set(0.4); app._upd_curr_plot(); pump(0.5)
        assert app._cur_sa_lbl.cget("text") == "0.40", app._cur_sa_lbl.cget("text")
        app._lay_sel_auto(); pump(0.3)
        assert not app.lay_alpha, app.lay_alpha
        app.lay_ins3d.set(False); app.lay_sec.set(False); app._lay_set("all"); app._lay_top.destroy()
    step("panel toggles of every tab, Layers window, structure under the current", panels)

    def compare():
        app._load("BJT effects: Kirk effect (f_T vs I_C)"); pump(0.5)
        app.ivS.set(0.75); app.ivE.set(0.95); app.ivN.set(3); app._sweep(); wait()
        assert app.ivd is not None and np.isfinite(app.ivd["R"]["Qp"]).all(), "stored charge per point"
        assert "f_T,max" in app.res_par.get("1.0", "end"), "BJT panel parameters"
        app._pin_sweep(); pump(0.3); app._pins_win(); pump(0.3)
        app._prof_choices(); em = [v for v in app.cb_prof.cget("values") if v.startswith("Emitter")][0]
        app.cur_prof.set(em); app._prof_pick(); pump(0.4)
        app._ab_win(); pump(0.3)
        v = list(app.ab_vcb.cget("values")); assert len(v) == 2, v
        app.ab_var.set(v[0]); app._ab_variant("A"); app.ab_var.set(v[1]); app._ab_variant("B")
        app.ivS.set(0.75); app.ivE.set(0.95); app.ivN.set(3); app._ab_run(); wait()
        assert len(app.pins) == 3 and app.pins[-1]["s"] is not None, [p["name"] for p in app.pins]
        app._pin_tv.selection_set("1"); app._pins_fields(); pump(0.4)
        app._ab_load("A"); pump(0.4)
        app._ab_top.destroy(); app._pin_top.destroy()
        app._ptab_win(); pump(0.4)                             # the parameter table of A and B
        t = app._ptab_table(); names = t[0][2:4]; rows = {r[0]: r for r in t[1:]}
        assert names == [p["name"] for p in app.pins[1:]], (names, [p["name"] for p in app.pins])
        assert "Current gain β,max" in rows and "Transit frequency f_T,max" in rows, list(rows)
        assert "%" in rows["Transit frequency f_T,max"][-1], rows["Transit frequency f_T,max"]
        app.ptab_all.set(True); app.ptab_ref.set(names[1]); app.ptab_mode.set("Differences only"); app._ptab_fill()
        pump(0.3); assert len(app._ptab_table()[0]) == 2 + 3 + 2, app._ptab_table()[0]
        app._ptab_copy(); assert app.clipboard_get().count("\t") > 20
        csv_path = os.path.join(os.environ["HOME"], "ptab.csv")
        M.filedialog.asksaveasfilename = lambda *a, **k: csv_path
        app._ptab_csv(); assert open(csv_path, encoding="utf-8").read().startswith("Parameter,Unit,")
        app._ptab_top.destroy()
    step("pinned sweeps, A/B compare of template variants, current under a contact, parameter table", compare)

    def bad_file():                      # a malformed device file: one error dialog, editor unchanged
        import json
        app._load("NMOS"); pump(0.2)
        d = app._device_json(); d["regs"][0]["x1"] = d["regs"][0]["x0"] - 1.
        p = os.path.join(os.environ["HOME"], "bad_device.json")
        with open(p, "w") as f: json.dump(d, f)
        M.filedialog.askopenfilename = lambda *a, **k: p
        before = (app.tname, len(app.regs), len(app.cons)); n0 = len(ERRORS)
        app._open_device(); pump(0.2)
        assert len(ERRORS) == n0 + 1 and "region 1" in str(ERRORS[-1]), ERRORS[n0:]
        assert (app.tname, len(app.regs), len(app.cons)) == before
        del ERRORS[n0:]
    step("a malformed device file is refused", bad_file)

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
