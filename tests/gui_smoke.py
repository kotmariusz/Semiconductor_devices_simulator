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

for _s in (sys.stdout, sys.stderr):        # a Windows console may not be UTF-8
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
if "--flat" in sys.argv:
    os.environ["SEMISIM_FLAT_BUTTONS"] = "1"
TMP = tempfile.mkdtemp(prefix="semisim_smoke_")                  # every file the test writes
os.environ["HOME"] = TMP                                         # fresh settings file (Linux, macOS)
os.environ["SEMISIM_SETTINGS"] = os.path.join(TMP, "settings.json")   # ... and on Windows
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
        for q in ("net doping", "ionised doping", "trapped charge"):     # 3-D slicer maps
            app.sl_qty.set(q); app._upd_slice(); pump(0.2)
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

    def flow():                                  # what flows, threshold, flowlines, log map, animation, GIF
        app.nb.select(app.tc); pump(0.4)
        for fld in M.CUR_FIELDS:
            app.cur_field.set(fld); app.cur_thr.set("1 %"); app.cur_lines.set(True); app.cur_log.set(True)
            app._upd_curr_plot(); pump(0.5)
            maps = [a for a in app.fc.axes if "plane" in a.get_title()]
            assert maps and any("flowlines" in t.get_text() or "streamlines" in t.get_text() for t in maps[0].texts), fld
        app.cur_field.set("current J"); app._anim_toggle(); pump(1.5)
        assert app._anim is not None and app._anim_art is not None and app._anim_bg is not None, "animation running"
        P0 = app._anim.P.copy(); pump(0.6)
        assert np.abs(app._anim.P - P0).max() > 0, "particles move"
        app._anim_toggle(); pump(0.3); assert not app._anim_on and app._anim is None
        gif = os.path.join(TMP, "flow.gif")
        app._anim_gif(path=gif, frames=8); pump(0.3)
        from PIL import Image
        with Image.open(gif) as im:
            assert im.n_frames == 8, im.n_frames
            cols = set()                         # one palette for all frames: the colour scale cannot drift
            for k in range(im.n_frames):
                im.seek(k); cols |= set(map(tuple, np.asarray(im.convert("RGB")).reshape(-1, 3)[::7].tolist()))
            assert len(cols) <= M.GIF_COLORS, len(cols)
        app.cur_thr.set("off"); app.cur_lines.set(False); app.cur_log.set(False); app._upd_curr_plot(); pump(0.3)
        # the 3-D view animated alone, turned while it runs, and saved as a GIF
        vis0 = {k: v.get() for k, v in app.pv["tc"].items()}
        for k in app.pv["tc"]: app.pv["tc"][k].set(k == "3d")
        app.cur_anim3.set(True); app._anim_toggle(); pump(1.5)
        ok2, ok3 = app._anim_ready(); assert ok3 and not ok2, (ok2, ok3)
        Q0 = app._anim3.P.copy(); pump(0.6)
        assert np.abs(app._anim3.P - Q0).max() > 0, "3-D particles move"
        app._anim3_ax.view_init(elev=30, azim=-40); app.cvc.draw(); pump(0.6)
        assert app._anim_ready()[1], "3-D background kept after a turn"
        app._anim_toggle(); pump(0.3); assert app._anim3 is None and app._anim3_art is None
        app._anim_gif(path=gif, frames=6); pump(0.3)
        with Image.open(gif) as im: assert im.n_frames == 6, im.n_frames
        for k, v in vis0.items(): app.pv["tc"][k].set(v)
        app._upd_curr_plot(); pump(0.3)
    step("Current tab: electron / hole flow, threshold, flowlines, animation (map and 3-D), GIF", flow)

    def defects3d():                             # the defects of the solved device in the 3-D current view
        app._load("Defects: leaky p-n junction (particle, void, dislocation)"); pump(0.4); app._run(); wait()
        app.nb.select(app.tc); pump(0.3)
        vis0 = {k: v.get() for k, v in app.pv["tc"].items()}
        for k in app.pv["tc"]: app.pv["tc"][k].set(k == "3d")
        app._upd_curr_plot(); pump(0.6)
        a3 = [a for a in app.fc.axes if hasattr(a, "get_zlim")][0]
        assert "3 defects" in a3.get_title(), a3.get_title()
        app.lay_def.set(False); app._lay_redraw(); pump(0.6)
        a3 = [a for a in app.fc.axes if hasattr(a, "get_zlim")][0]
        assert "defect" not in a3.get_title(), a3.get_title()
        app.lay_def.set(True)
        for k, v in vis0.items(): app.pv["tc"][k].set(v)
        app._load("PN Diode"); pump(0.4); app._run(); wait()
    step("defects shown in the 3-D current view (Layers → Defects)", defects3d)

    def mosfet():                                # C-V sweep, capacitances and on-resistance tools, R_DS(on)
        app._load("RW: Si MOSFET 55V (IRLZ44N)"); pump(0.5)
        disp = {v: k for k, v in app._gmap.items()}
        for c in app.cons:
            if c.bc == 2: c.V = 0.
        app.iv_tgt.set(disp["Drain"]); app.ivS.set(0.); app.ivE.set(20.); app.ivN.set(3); app.iv_cap.set(True)
        app._sweep(); wait()
        R = app.ivd["R"]; assert "K" in R and np.isfinite(R["K"]).all(), "capacitance matrices"
        assert any(a.get_title().startswith("Capacitances") for a in app.fiv.axes), "C-V panel"
        par = app.res_par.get("1.0", "end"); assert "C_iss" in par and "C_rss" in par and "C_o(tr)" in par, par
        app._cap_here(); wait(); assert "C_iss" in app._last_report[1], app._last_report
        app.iv_cap.set(False)
        for c in app.cons:
            if c.bc == 2: c.V = 5.
        app.ivS.set(0.); app.ivE.set(0.4); app.ivN.set(3); app._sweep(); wait()
        assert "R_DS(on)" in app.res_par.get("1.0", "end"), app.res_par.get("1.0", "end")
        app._ron_here(); pump(0.3); assert "R_i" in app._last_report[1] and "n epi" in app._last_report[1], app._last_report
    step("MOSFET: C-V sweep (C_iss, C_oss, C_rss), Analysis menu, R_DS(on)", mosfet)

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
        csv_path = os.path.join(TMP, "ptab.csv")
        M.filedialog.asksaveasfilename = lambda *a, **k: csv_path
        app._ptab_csv(); assert open(csv_path, encoding="utf-8").read().startswith("Parameter,Unit,")
        app._ptab_top.destroy()
    step("pinned sweeps, A/B compare of template variants, current under a contact, parameter table", compare)

    def bad_file():                      # a malformed device file: one error dialog, editor unchanged
        import json
        app._load("NMOS"); pump(0.2)
        d = app._device_json(); d["regs"][0]["x1"] = d["regs"][0]["x0"] - 1.
        p = os.path.join(TMP, "bad_device.json")
        with open(p, "w") as f: json.dump(d, f)
        M.filedialog.askopenfilename = lambda *a, **k: p
        before = (app.tname, len(app.regs), len(app.cons)); n0 = len(ERRORS)
        app._open_device(); pump(0.2)
        assert len(ERRORS) == n0 + 1 and "region 1" in str(ERRORS[-1]), ERRORS[n0:]
        assert (app.tname, len(app.regs), len(app.cons)) == before
        del ERRORS[n0:]
    step("a malformed device file is refused", bad_file)

    def caxis_file():                    # the 4H-SiC c-axis travels with a device file
        import json
        app._load("NMOS"); pump(0.2)
        app.caxis.set("x"); d = app._device_json()
        assert d["solver"]["caxis"] == "x" and app._device_dict()["caxis"] == "x", d["solver"]
        p = os.path.join(TMP, "caxis_device.json")
        with open(p, "w", encoding="utf-8") as f: json.dump(d, f)
        app.caxis.set("y"); M.filedialog.askopenfilename = lambda *a, **k: p
        app._open_device(); pump(0.2)
        assert app.caxis.get() == "x", app.caxis.get()
        app._load("NMOS"); pump(0.2)
        assert app.caxis.get() == "y", app.caxis.get()      # a template: the wafer normal
    step("the c-axis option travels with a device file", caxis_file)

    def profiles():                      # doping profiles: editor, structure view, device file (version 2)
        import json
        app._load("RW: Si MOSFET 55V (IRLZ44N), diffused"); pump(0.5)
        assert len(app.implants) == 4 and len(app.itree.get_children()) == 4
        app.nb.select(0); pump(0.8)                                  # Structure: junction contours, windows
        app.ib_dt.set("n"); app.ib_kind.set("erfc"); app.ib_face.set("y0"); app._fullimp()
        app.ib_peak.set("1e18"); app.ib_Lv.set("50"); app._addimp(); pump(0.3)
        assert len(app.implants) == 5 and app.implants[-1].kind == "erfc"
        app._select(app.itree, 4); pump(0.2); app._delimp(); pump(0.3)
        assert len(app.implants) == 4
        d = app._device_json()
        assert d["version"] == 2 and len(d["implants"]) == 4, (d["version"], len(d.get("implants", [])))
        p = os.path.join(TMP, "profiles_device.json")
        with open(p, "w", encoding="utf-8") as f: json.dump(d, f)
        app._load("NMOS"); pump(0.3)
        assert not app.implants and app._device_json()["version"] == 1 and "implants" not in app._device_json()
        M.filedialog.askopenfilename = lambda *a, **k: p
        app._open_device(); pump(0.5)
        assert len(app.implants) == 4 and app.implants[0].label == "p-body L (diffused)"
    step("doping profiles: editor, structure view, device file", profiles)

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
