"""Does the LMC magnetic field run along the rims of the supergiant shells?

Question (Prof. Alina, 11.11.2026: JWST "bubbly" galaxies; Inutsuka et al. 2015: the ISM is made of bubbles):
if expanding shells shape the field, the field should run tangential to their rims. Outcome use: tangential in both
H I VGT and Planck dust -> the shells are a structure the LMC-field work must treat explicitly; tangential in VGT
only -> VGT at the rims follows the swept-up gas geometry, which Planck does not confirm; neither -> shells do not
organise the field at our resolution.

Tracers, all on the H I grid binned 4x (80" pixels; the VGT angle is smooth on ~13.7', so nothing is lost):
- VGT: psi_crop_y1_x113_none.fits (LMC H I, 182-339 km/s), averaged as Stokes cos 2psi, sin 2psi per 4x4 block.
- Planck 353 GHz: B angle = IAU polarization angle + 90 deg, Q/U smoothed to 10' FWHM (from 4.8') because single
  Planck pixels are noise dominated in the diffuse LMC. No foreground subtraction: a uniform foreground adds the same
  Q/U everywhere, and a uniform orientation scores exactly 0 on a full ring (the tangent turns through 180 deg), so it
  dilutes but cannot fake tangency. Noise also only dilutes; the random-ring null below carries the same dilution.
- N_HI structure (control for geometry): orientation of the integrated H I map's iso-contours (gradient + 90 deg),
  smoothed to the VGT resolution (13.7' FWHM). A bright ring has tangential iso-contours by construction, so this
  shows how much tangency the gas geometry alone gives.

Caveat on what a tangential VGT field means here: this VGT sums channel-map intensity gradients, and an expanding
shell is a ring in every channel (of different radius per channel), so swept-up gas alone makes VGT tangential at a
rim, field or no field. The N_HI-contour control works on the integrated map and therefore does not capture that
(controls at the channel level: VGT on channels outside the expansion velocities, on high-pass-filtered channels, or on
a synthetic field-free expanding shell). Planck is the independent check; it is not foreground-subtracted here, and the
foreground varies on degree scales, so "a uniform foreground only dilutes" is approximate.

Statistic: tangency T = <cos 2(theta - theta_tangent)> over the rim annulus [A0 R, A1 R]; +1 along the rim, -1 radial,
0 unrelated. Computed for every possible centre at once by FFT convolution (cos 2theta_t = -cos 2theta_r, and the
annulus kernel weighted by cos/sin 2theta_r is symmetric under d -> -d, so convolution = correlation).
Null: the same ring (same radius) centred anywhere else in the covered LMC (annulus >= MIN_COVER valid). Per shell,
p = share of valid centres with T >= observed. Pooled: mean T over the shells vs 20000 draws placing every shell at
a random valid centre. Caveat: random rings can land on other real shells (Kim et al. 1999 list ~100 in the LMC),
which makes the null conservative.

Shells: Book et al. 2008 supergiant shells LMC 1-9 (centres rounded to 1 min RA / 10' Dec, diameters in pc at
50 kpc, 1 deg = 873 pc), the same list as demcels_quicklook.py.

Usage: python3 scripts/optical/shell_rim_orientation.py --out DIR [--halpha demcels_ha_on_hi_grid.fits]
"""
import argparse
import csv
import os
import sys

import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import compare_vgt_polarization as cvp  # noqa: E402
from demcels_quicklook import SHELLS, PC_PER_DEG  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PSI = f"{REPO}/data/sims/vgt_lmc_none_dn20/psi_crop_y1_x113_none.fits"
NHI = f"{REPO}/data/sims/lmc_nhi_full.npy"   # full 2174^2 H I grid; the VGT crop starts at x = 112
PLANCK = f"{REPO}/data/external/HFI_SkyMap_353-BPassCorrected_2048_R4.00_full.fits"
X0 = 112
BIN = 4
VGT_FWHM_ARCMIN = 13.7      # effective VGT resolution (fig_region_panel.py, README)
PLANCK_FWHM_ARCMIN = 10.0
PLANCK_NATIVE_ARCMIN = 4.8
MIN_COVER = 0.7
# Book et al. 2008 (sec. 4-5) conclude LMC 7 and LMC 9 are not physical shells ("false shells"): they are scored
# and shown, but kept out of the pooled value and the Holm family, and serve as negative controls.
FALSE_SHELLS = {"LMC 7", "LMC 9"}
NDRAW = 20000


def nanbin(a, n):
    ny, nx = (a.shape[0] // n) * n, (a.shape[1] // n) * n
    return np.nanmean(a[:ny, :nx].reshape(ny // n, n, nx // n, n), axis=(1, 3))


def nansmooth(a, s):
    m = np.isfinite(a)
    num = gaussian_filter(np.where(m, a, 0.0), s)
    den = gaussian_filter(m.astype(float), s)
    return np.where(m & (den > 0.3), num / np.maximum(den, 1e-9), np.nan)


def load_tracers():
    psi, wcs_full = cvp.load_psi(PSI)
    ny, nx = psi.shape
    with np.errstate(invalid="ignore"):
        c, s = nanbin(np.cos(2 * psi), BIN), nanbin(np.sin(2 * psi), BIN)
    wcs = wcs_full.slice((slice(None, None, BIN), slice(None, None, BIN)))
    pix_arcmin = abs(wcs.wcs.cdelt[0]) * 60
    tr = {"VGT (H I)": 0.5 * np.arctan2(s, c)}

    nhi = np.load(NHI)[:ny, X0:X0 + nx].astype(float)
    nhi[nhi <= 0] = np.nan
    nb = nansmooth(nanbin(nhi, BIN), VGT_FWHM_ARCMIN / 2.3548 / pix_arcmin)
    gy, gx = np.gradient(nb)
    tr["N_HI contours"] = np.arctan2(gy, gx) + np.pi / 2
    tr["N_HI contours"][~np.isfinite(tr["VGT (H I)"])] = np.nan

    planck = cvp.read_planck(PLANCK)
    p, l, b = cvp.planck_on_grid(planck, wcs, c.shape)
    sk = np.sqrt(PLANCK_FWHM_ARCMIN ** 2 - PLANCK_NATIVE_ARCMIN ** 2) / 2.3548 / pix_arcmin
    q, uu = gaussian_filter(p["Q"], sk), gaussian_filter(p["U"], sk)
    chi = 0.5 * np.arctan2(uu, q)                      # IAU polarization angle, Galactic east of north
    tr["Planck 353 GHz"] = cvp.galactic_angle_to_pixel(wcs, l.ravel(), b.ravel(), (chi + np.pi / 2).ravel()) \
        .reshape(c.shape)
    tr["Planck 353 GHz"][~np.isfinite(tr["VGT (H I)"])] = np.nan  # same footprint for all tracers
    return tr, wcs, nb


def ring_kernels(r_px, a0, a1):
    n = int(np.ceil(a1 * r_px)) + 1
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1].astype(float)
    rr = np.hypot(xx, yy)
    ann = ((rr >= a0 * r_px) & (rr <= a1 * r_px)).astype(float)
    th = np.arctan2(yy, xx)
    return ann, ann * np.cos(2 * th), ann * np.sin(2 * th)


def tangency_map(ang, r_px, a0, a1):
    valid = np.isfinite(ang)
    c = np.where(valid, np.cos(2 * ang), 0.0)
    s = np.where(valid, np.sin(2 * ang), 0.0)
    ann, kc, ks = ring_kernels(r_px, a0, a1)
    n = fftconvolve(valid.astype(float), ann, mode="same")
    t = -(fftconvolve(c, kc, mode="same") + fftconvolve(s, ks, mode="same")) / np.maximum(n, 1e-9)
    cover = n / ann.sum()
    return np.where(cover >= MIN_COVER, t, np.nan), cover


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--halpha", default=None)
    ap.add_argument("--annuli", default="0.8:1.2,0.7:1.3,0.9:1.1")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    tr, wcs, nb = load_tracers()
    pix_deg = abs(wcs.wcs.cdelt[0])
    rng = np.random.default_rng(20261111)
    rows, maps = [], {}
    for annulus in args.annuli.split(","):
        a0, a1 = map(float, annulus.split(":"))
        for name, ang in tr.items():
            obs, nulls, pvals = [], [], []
            for shell, (ra, dec, d) in SHELLS.items():
                r_px = d / 2 / PC_PER_DEG / pix_deg
                tmap, cover = tangency_map(ang, r_px, a0, a1)
                x, y = wcs.world_to_pixel(SkyCoord(ra, dec, frame="icrs"))
                ix, iy = int(round(float(x))), int(round(float(y)))
                inside = 0 <= iy < tmap.shape[0] and 0 <= ix < tmap.shape[1]
                t = tmap[iy, ix] if inside else np.nan
                pool = tmap[np.isfinite(tmap)]
                p = (1 + np.sum(pool >= t)) / (1 + pool.size) if np.isfinite(t) else np.nan
                rows.append(dict(annulus=annulus, tracer=name, shell=shell,
                                 kind="shell" if shell not in FALSE_SHELLS else "false shell (control)",
                                 r_deg=round(r_px * pix_deg, 3),
                                 T=round(float(t), 3) if np.isfinite(t) else "", p=round(float(p), 4) if np.isfinite(p) else "",
                                 p_holm="", null_mean=round(float(pool.mean()), 3),
                                 null_p95=round(float(np.percentile(pool, 95)), 3),
                                 n_centres=pool.size, cover=round(float(cover[iy, ix]), 2) if inside else 0))
                if np.isfinite(t) and shell not in FALSE_SHELLS:
                    obs.append(t)
                    nulls.append(pool)
                    pvals.append((p, len(rows) - 1))
                if annulus == "0.8:1.2":
                    maps[(name, shell)] = tmap
            # Holm correction over the physical shells (one family per tracer and annulus)
            order = sorted(pvals)
            k, running = len(order), 0.0
            for i, (pv, ri) in enumerate(order):
                running = max(running, min(1.0, (k - i) * pv))
                rows[ri]["p_holm"] = round(running, 4)
            draws = np.mean([rng.choice(pl, NDRAW) for pl in nulls], axis=0)
            m = float(np.mean(obs))
            rows.append(dict(annulus=annulus, tracer=name, shell=f"pooled ({len(obs)} shells)", kind="pooled",
                             r_deg="", T=round(m, 3), p=round(float((1 + np.sum(draws >= m)) / (1 + NDRAW)), 4),
                             p_holm="", null_mean=round(float(draws.mean()), 3),
                             null_p95=round(float(np.percentile(draws, 95)), 3), n_centres=NDRAW, cover=""))
        # across shells: does VGT tangency follow the geometry control or Planck?
        from scipy.stats import spearmanr
        T = {n: [float(r["T"]) for r in rows if r["annulus"] == annulus and r["tracer"] == n and r["kind"] == "shell"]
             for n in tr}
        for other in ("N_HI contours", "Planck 353 GHz"):
            rho, pr = spearmanr(T["VGT (H I)"], T[other])
            print(f"annulus {annulus}: Spearman over the {len(T[other])} shells, VGT T vs {other} T: {rho:+.2f} (p = {pr:.2f})")
    with open(os.path.join(args.out, "shell_rim_tangency.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    np.savez_compressed(os.path.join(args.out, "shell_rim_inputs.npz"), nhi=nb,
                        **{k: v for k, v in tr.items()}, header=wcs.to_header_string())
    for r in rows:
        if r["annulus"] == "0.8:1.2":
            print(r)
    plot(args.out, tr, wcs, nb, rows, args.halpha)


def plot(out, tr, wcs, nb, rows, halpha_path):
    """Fig 1: raw zooms, one per shell: N_HI (grey), Halpha contours where DeMCELS covers it, VGT (blue) and Planck
    (orange) field segments every 4 binned pixels (5.3'), rim annulus dashed. Fig 2: tangency per shell and tracer
    with the random-ring 95th percentile."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    BLUE, ORANGE, GREY = "#1f5fae", "#e3a21a", "#7a7a7a"
    ha = None
    if halpha_path:
        h = fits.getdata(halpha_path).astype(float)[:, X0:X0 + 2048][:2048]  # same full H I grid as N_HI
        ha = nansmooth(nanbin(h, BIN), 1.0)
    pix_deg = abs(wcs.wcs.cdelt[0])
    get = {(r["tracer"], r["shell"]): r for r in rows if r["annulus"] == "0.8:1.2"}
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    for ax, (shell, (ra, dec, d)) in zip(axes.ravel(), SHELLS.items()):
        r_px = d / 2 / PC_PER_DEG / pix_deg
        x, y = (float(v) for v in wcs.world_to_pixel(SkyCoord(ra, dec, frame="icrs")))
        h = 1.7 * r_px
        x0, x1, y0, y1 = int(x - h), int(x + h), int(y - h), int(y + h)
        sub = nb[max(y0, 0):y1, max(x0, 0):x1]
        ax.imshow(nb, origin="lower", cmap="Greys", vmin=np.nanpercentile(sub, 2), vmax=np.nanpercentile(sub, 99.5))
        if ha is not None and np.isfinite(ha[max(y0, 0):y1, max(x0, 0):x1]).mean() > 0.2:
            hs = ha[max(y0, 0):y1, max(x0, 0):x1]
            lv = np.nanpercentile(hs, [85, 95])
            ax.contour(np.where(np.isfinite(ha), ha, -1e9), levels=lv, colors="#c2185b", linewidths=0.7)
        th = np.linspace(0, 2 * np.pi, 200)
        for f in (0.8, 1.2):
            ax.plot(x + f * r_px * np.cos(th), y + f * r_px * np.sin(th), "--", color="k", lw=0.8)
        yy, xx = np.mgrid[y0:y1:4, x0:x1:4]
        ok = (yy >= 0) & (xx >= 0) & (yy < nb.shape[0]) & (xx < nb.shape[1])
        yy, xx = yy[ok], xx[ok]
        for name, col, off in (("VGT (H I)", BLUE, 0.0), ("Planck 353 GHz", ORANGE, 0.0)):
            a = tr[name][yy, xx]
            m = np.isfinite(a)
            L = 1.6
            seg = np.stack([np.c_[xx[m] - L * np.cos(a[m]), yy[m] - L * np.sin(a[m])],
                            np.c_[xx[m] + L * np.cos(a[m]), yy[m] + L * np.sin(a[m])]], axis=1)
            ax.add_collection(LineCollection(seg, colors=col, linewidths=1.4 if col == BLUE else 1.1,
                                             alpha=0.9 if col == BLUE else 0.8))
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_xticks([]); ax.set_yticks([])
        v, p = get[("VGT (H I)", shell)], get[("Planck 353 GHz", shell)]
        g = get[("N_HI contours", shell)]
        ax.set_title(f"{shell} (r = {d / 2 / PC_PER_DEG:.2f} deg)\nVGT T={v['T']} p={v['p']}   "
                     f"Planck T={p['T']} p={p['p']}   N_HI T={g['T']}", fontsize=10)
    fig.suptitle("Field at the rims of the LMC supergiant shells (Book et al. 2008)\nGrey: H I column; pink: DeMCELS "
                 "H-alpha contours (85th, 95th pct); blue: VGT field; orange: Planck 353 GHz field (10'); dashed: "
                 "rim annulus 0.8-1.2 R\nT = <cos 2(angle - rim tangent)>: +1 along the rim, -1 radial, 0 unrelated; "
                 "p = share of random ring centres with T at least as high", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(out, "shell_rim_zooms.png"), dpi=110)
    plt.close(fig)

    shells = list(SHELLS) + ["pooled (7 shells)"]
    fig, ax = plt.subplots(figsize=(11, 5))
    w = 0.27
    for k, (name, col) in enumerate((("VGT (H I)", BLUE), ("Planck 353 GHz", ORANGE), ("N_HI contours", GREY))):
        xs = np.arange(len(shells)) + (k - 1) * w
        vals = [float(get[(name, s)]["T"]) for s in shells]
        p95 = [float(get[(name, s)]["null_p95"]) for s in shells]
        ax.bar(xs, vals, w, color=col, label=name)
        ax.scatter(xs, p95, marker="_", s=120, color="k", zorder=3, label="random rings, 95th pct" if k == 0 else None)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(np.arange(len(shells)), [s.replace(" shells", "") + ("\n(false shell)" if s in FALSE_SHELLS else "")
                                           for s in shells])
    ax.set_ylabel("rim tangency T")
    ax.set_title("Is the field along the shell rims? Bars above the black tick beat 95% of random ring positions\n"
                 "(LMC 7 and LMC 9: not real shells per Book et al. 2008, kept as controls; pooled over the 7 real shells)",
                 fontsize=10)
    ax.legend(frameon=False, ncol=4, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    fig.savefig(os.path.join(out, "shell_rim_tangency.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
