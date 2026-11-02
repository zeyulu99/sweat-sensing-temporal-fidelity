"""Laplace-domain (analytical) solutions for the transport-sensor cascade.

Requested in the second review round: express the envelopes analytically so
that the confounding factors are connected and the numerical grid is bounded.

Transfer functions (s in min^-1, time constants in min)
  plug flow        Hf = exp(-s tau_f)
  well mixed       Hf = 1 / (1 + s tau_f)
  gamma, shape k   Hf = (1 + s tau_f / k)^-k
  first-order sensor Hs = 1 / (1 + s tau_s)
Mean delay and variance of the cascade add stage by stage:
  mean = tau_f + tau_s ; var = {0, tau_f^2, tau_f^2/k} + tau_s^2

Outputs (outputs/revision_2026-10-27/analytic/)
  A_periodic_exact_vs_numeric.csv  exact periodic envelope vs dense numeric grid
  B_convergence.csv                numeric metrics vs time step, against exact
  C_ramp_rule_vs_numeric.csv       tau_f,max = T/4 - tau_s vs numeric
  D_pulse_moment_vs_numeric.csv    Gaussian-moment approximation vs numeric
  E_design_rules.csv               small-argument rules vs exact boundaries
  figS_analytic_vs_numeric.{pdf,png,svg}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compare_figure2_strategies import strict_pass, transport_fast
from run_evidence_closure_analysis import canonical_metrics, canonical_waveform

PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT / "outputs" / "revision_2026-10-27" / "analytic"
OUT.mkdir(parents=True, exist_ok=True)
DENSE = PROJECT / "outputs" / "revision_2026-09-22" / "round2_analyses" / "A_envelope_summary.csv"
TOPOS = ["plug", "cstr", "gamma_rtd"]
NAMES = {"plug": "Plug flow", "cstr": "Well mixed", "gamma_rtd": "Gamma RTD (k = 3)"}
TS = [2.0, 5.0, 10.0, 30.0, 60.0]
K = 3.0


# ------------------------------------------------------------ transfer functions
def gain_phase(topo, tf, ts, w):
    """Exact gain and (unwrapped) phase lag of the cascade at angular frequency w."""
    if topo == "plug":
        gf, pf = 1.0, w * tf
    elif topo == "cstr":
        gf, pf = 1.0 / np.sqrt(1 + (w * tf) ** 2), np.arctan(w * tf)
    else:
        a = w * tf / K
        gf, pf = (1 + a ** 2) ** (-K / 2), K * np.arctan(a)
    gs, ps = 1.0 / np.sqrt(1 + (w * ts) ** 2), np.arctan(w * ts)
    return gf * gs, pf + ps


def periodic_metrics(topo, tf, ts, T):
    """Exact steady-state scores for a unit-range sinusoid of period T."""
    w = 2 * np.pi / T
    A, phi = gain_phase(topo, tf, ts, w)
    nrmse = np.sqrt(1 - 2 * A * np.cos(phi) + A ** 2) / (2 * np.sqrt(2))
    return {"amplitude_retention": A, "nrmse": nrmse, "pearson_r": np.cos(phi),
            "feature_delay_min": phi / w}


def boundary(passes, lo=1e-4, hi=15.0):
    """Largest tau_f of the contiguous passing range starting at lo (bisection)."""
    if not passes(lo):
        return np.nan
    if passes(hi):
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if passes(mid) else (lo, mid)
    return lo


# ------------------------------------------------------------ A. periodic, exact
def part_a():
    dense = pd.read_csv(DENSE)
    rows = []
    for ts in (0.0, 1.0):
        for T in TS:
            for topo in TOPOS:
                b = boundary(lambda tf: strict_pass("sinusoid", T, periodic_metrics(topo, tf, ts, T)))
                d = dense[(dense.tau_s_min == ts) & (dense.waveform == "sinusoid")
                          & (dense.T_min == T) & (dense.topology == topo)].iloc[0]
                rows.append({"tau_s_min": ts, "T_min": T, "topology": topo,
                             "exact_tau_f_max": b,
                             "numeric_tau_f_max": d.strict_max_tau_f,
                             "numeric_next_fail": d.strict_next_fail})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "A_periodic_exact_vs_numeric.csv", index=False)
    return df


# ------------------------------------------------------------ B. convergence
def part_b():
    rows = []
    for T, topo, tf, ts in [(30.0, "cstr", 0.335, 1.0), (10.0, "cstr", 0.30, 0.0),
                            (2.0, "cstr", 0.08, 0.0), (30.0, "gamma_rtd", 0.335, 1.0)]:
        ex = periodic_metrics(topo, tf, ts, T)
        for dt in (0.05, 0.01, 0.002):
            dur = max(180.0, 8 * T)
            tm = np.arange(0.0, dur + dt, dt)
            u = 0.5 + 0.5 * np.sin(2 * np.pi * tm / T)
            _, s = transport_fast(tm * 60.0, u, topo, tf, sensor_min=ts)
            m = canonical_metrics("sinusoid", T, tm * 60.0, u, s)
            rows.append({"T_min": T, "topology": topo, "tau_f_min": tf, "tau_s_min": ts,
                         "dt_min": dt, "numeric_nrmse": m["nrmse"],
                         "exact_nrmse": ex["nrmse"],
                         "numeric_delay_min": m["feature_delay_min"],
                         "exact_delay_min": ex["feature_delay_min"],
                         "numeric_amplitude": m["amplitude_retention"],
                         "exact_amplitude": ex["amplitude_retention"]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "B_convergence.csv", index=False)
    return df


# ------------------------------------------------------------ C. ramp
def part_c():
    """For an LTI cascade, a ramp is reproduced after the transient with a lag
    equal to the first moment, tau_f + tau_s, independent of topology. With
    full amplitude retention, the delay limit max(1, T/4) then binds."""
    dense = pd.read_csv(DENSE)
    rows = []
    for ts in (0.0, 1.0):
        for T in TS:
            rule = max(1.0, T / 4) - ts
            for topo in TOPOS:
                d = dense[(dense.tau_s_min == ts) & (dense.waveform == "ramp")
                          & (dense.T_min == T) & (dense.topology == topo)].iloc[0]
                rows.append({"tau_s_min": ts, "T_min": T, "topology": topo,
                             "rule_tau_f_max": min(rule, 15.0),
                             "numeric_tau_f_max": d.strict_max_tau_f,
                             "binding_numeric": d.binding_criterion})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "C_ramp_rule_vs_numeric.csv", index=False)
    return df


# ------------------------------------------------------------ D. pulse
def kernel_var(topo, tf, ts):
    vf = {"plug": 0.0, "cstr": tf ** 2, "gamma_rtd": tf ** 2 / K}[topo]
    return vf + ts ** 2


def pulse_moment_metrics(topo, tf, ts, T):
    """Gaussian-moment approximation: the output is a Gaussian with the input
    area, its centroid shifted by tau_f + tau_s and its variance increased by
    the kernel variance. Scored with the same metric code on the same window."""
    t, u = canonical_waveform("pulse", T)
    tm = t / 60.0
    s0 = T / 2.355
    centre = tm[np.argmax(u)]
    s1 = np.sqrt(s0 ** 2 + kernel_var(topo, tf, ts))
    y = (s0 / s1) * np.exp(-0.5 * ((tm - centre - (tf + ts)) / s1) ** 2)
    return canonical_metrics("pulse", T, t, u, y)


def part_d():
    dense = pd.read_csv(DENSE)
    rows = []
    for ts in (0.0, 1.0):
        for T in TS:
            for topo in TOPOS:
                b = boundary(lambda tf: strict_pass("pulse", T, pulse_moment_metrics(topo, tf, ts, T)))
                d = dense[(dense.tau_s_min == ts) & (dense.waveform == "pulse")
                          & (dense.T_min == T) & (dense.topology == topo)].iloc[0]
                rows.append({"tau_s_min": ts, "T_min": T, "topology": topo,
                             "moment_tau_f_max": b,
                             "numeric_tau_f_max": d.strict_max_tau_f})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "D_pulse_moment_vs_numeric.csv", index=False)
    return df


# ------------------------------------------------------------ E. design rules
def part_e(a):
    """Small-argument expansion of the exact periodic scores:
      NRMSE ~ phi / (2 sqrt 2), phi ~ w (tau_f + tau_s)  ->  tau_f + tau_s <= 0.045 T
      A     ~ 1 - w^2 sum(var) / 2                     ->  sqrt(sum var) <= 0.071 T
    The NRMSE rule binds first; compare it with the exact boundary."""
    c_nrmse = 2 * np.sqrt(2) * 0.10 / (2 * np.pi)         # 0.0450
    c_amp = np.sqrt(2 * (1 - 0.90)) / (2 * np.pi)          # 0.0712
    rows = []
    for _, r in a.iterrows():
        rows.append({"tau_s_min": r.tau_s_min, "T_min": r.T_min, "topology": r.topology,
                     "rule_tau_f_max": max(c_nrmse * r.T_min - r.tau_s_min, np.nan if c_nrmse * r.T_min <= r.tau_s_min else 0),
                     "exact_tau_f_max": r.exact_tau_f_max})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "E_design_rules.csv", index=False)
    return c_nrmse, c_amp, df


# ------------------------------------------------------------ figure S
def figure(a, c, d):
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator, FuncFormatter
    plt.rcParams.update({"font.family": "Arial", "font.size": 7.5,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    colors = {"plug": "#151515", "cstr": "#245DB8", "gamma_rtd": "#E87522"}
    marks = {"plug": "o", "cstr": "s", "gamma_rtd": "^"}
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.3), sharex=True, sharey=True)
    panels = [("Pulse", d, "moment_tau_f_max", "Gaussian-moment approximation"),
              ("Ramp", c, "rule_tau_f_max", r"$\tau_{\mathrm{f}}+\tau_{\mathrm{s}} = T/4$"),
              ("Periodic", a, "exact_tau_f_max", "exact (Laplace)")]
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    for row, ts in enumerate((0.0, 1.0)):
        for col, (title, df, colname, lab) in enumerate(panels):
            ax = axes[row, col]
            sub = df[df.tau_s_min == ts]
            for topo in TOPOS:
                q = sub[sub.topology == topo].sort_values("T_min")
                y = q[colname].where(q[colname] > 0)
                ax.plot(q.T_min, y, color=colors[topo], lw=1.1, zorder=2)
                ax.plot(q.T_min, q.numeric_tau_f_max, marks[topo], color=colors[topo],
                        ms=4.0, mfc="white", mew=1.0, zorder=3, ls="none")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.xaxis.set_major_locator(FixedLocator(TS)); ax.xaxis.set_major_formatter(fmt)
            ax.xaxis.set_minor_locator(NullLocator()); ax.xaxis.set_minor_formatter(NullFormatter())
            ax.yaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.5, 1, 5, 10]))
            ax.yaxis.set_major_formatter(fmt); ax.yaxis.set_minor_formatter(NullFormatter())
            ax.set_xlim(1.7, 70); ax.set_ylim(0.04, 20)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if row == 0:
                ax.set_title(title, fontsize=8.0, fontweight="bold")
                ax.text(0.04, 0.97, f"Lines: {lab}", transform=ax.transAxes, va="top",
                        fontsize=6.3, color="#444444")
            if row == 1:
                ax.set_xlabel(r"Feature timescale, $T$ (min)")
            if col == 0:
                ax.set_ylabel(r"Largest passing $\tau_{\mathrm{f}}$ (min)")
            ax.text(0.96, 0.05, rf"$\tau_{{\mathrm{{s}}}}$ = {ts:g} min", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=6.6)
    handles = [plt.Line2D([], [], color=colors[t], marker=marks[t], mfc="white", lw=1.1,
                          label=NAMES[t]) for t in TOPOS]
    handles.append(plt.Line2D([], [], color="#777777", marker="o", mfc="white", ls="none",
                              label="Numerical grid (markers)"))
    fig.legend(handles=handles, frameon=False, fontsize=6.6, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    for ext in ("pdf", "svg"):
        fig.savefig(OUT / f"figS_analytic_vs_numeric.{ext}")
    fig.savefig(OUT / "figS_analytic_vs_numeric.png", dpi=600)
    plt.close(fig)


def main():
    a = part_a(); b = part_b(); c = part_c(); d = part_d()
    cn, ca, e = part_e(a)
    figure(a, c, d)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 20)
    print("A periodic\n", a.round(3).to_string())
    print("B convergence\n", b.round(4).to_string())
    print("C ramp\n", c.round(3).to_string())
    print("D pulse\n", d.round(3).to_string())
    print(f"E rules: NRMSE-limited c = {cn:.4f}, amplitude-limited c = {ca:.4f}")
    print(e.round(3).to_string())


if __name__ == "__main__":
    main()
