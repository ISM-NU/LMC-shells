"""Quick-look PNGs of the DeMCELS mosaics on the H I grid: Halpha, [S II], [S II]/Halpha, with the
9 LMC supergiant shells of Book et al. 2008 (LMC 1-9) overlaid.

Inputs: OUT_DIR with demcels_{ha,s2}_on_hi_grid[_bgmatch].fits from demcels_mosaic.py.
Output: PNGs in PNG_DIR (log stretch for Halpha/[S II]; ratio panel diverging, centred on 0.4).
Ratio: ([S II] DN * 2.1e-18) / (Halpha DN * 2.4e-18) using the DR1 calibrations (Points et al. 2024), from maps
smoothed with a 1-pixel (20") Gaussian, shown only where smoothed Halpha > NSIG x noise (robust sigma from
adjacent-pixel differences). Ratio > ~0.4 marks shock-dominated gas (SNRs, superbubble walls).
Caveats: zero level of the continuum-subtracted images is uncertain by field (use *_bgmatch for display only);
the N662 filter may include [N II]; ratios near the Halpha threshold are noise- and background-limited.
Shells: centre J2000, diameter in pc at 50 kpc, 1 deg = 873 pc.

Usage: python demcels_quicklook.py OUT_DIR PNG_DIR [suffix: ''|_bgmatch|_10arcsec ...]
"""
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

SHELLS = {  # Book et al. 2008: name -> (RA hms, Dec dms, diameter pc)
    "LMC 1": ("05h00m", "-65d40m", 750), "LMC 2": ("05h44m", "-69d20m", 900), "LMC 3": ("05h30m", "-69d00m", 1000),
    "LMC 4": ("05h32m", "-66d40m", 1000), "LMC 5": ("05h22m", "-66d00m", 800), "LMC 6": ("04h59m", "-68d36m", 600),
    "LMC 7": ("04h54m", "-69d33m", 800), "LMC 8": ("05h04m", "-70d30m", 900), "LMC 9": ("05h26m", "-71d00m", 890)}
PC_PER_DEG = 873.0
CAL_HA, CAL_S2 = 2.4e-18, 2.1e-18


def nansmooth(a, s):
    m = np.isfinite(a)
    num = gaussian_filter(np.where(m, a, 0), s)
    den = gaussian_filter(m.astype(float), s)
    return np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)


def robust_sigma(a):
    d = (a[:, 1:] - a[:, :-1])[np.isfinite(a[:, 1:]) & np.isfinite(a[:, :-1])]
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)


def draw_shells(ax, wcs, pixscale_deg):
    th = np.linspace(0, 2 * np.pi, 181)
    for name, (ra, dec, d) in SHELLS.items():
        c = SkyCoord(ra, dec, frame="icrs")
        rad = d / 2 / PC_PER_DEG
        p = c.directional_offset_by(th * u.rad, rad * u.deg)
        x, y = wcs.world_to_pixel(p)
        ax.plot(x, y, color="w", lw=2.2, alpha=0.9); ax.plot(x, y, color="k", lw=1.0)
        cx, cy = wcs.world_to_pixel(c)
        ax.text(cx, cy, name.split()[1], clip_on=True, color="w", fontsize=9, ha="center", va="center", fontweight="bold",
                path_effects=[__import__("matplotlib.patheffects", fromlist=["x"]).withStroke(linewidth=2, foreground="k")])


def main(outdir, pngdir, suffixes):
    for suf in suffixes:
        base = "on_hi_grid" + suf
        ha, hh = fits.getdata(f"{outdir}/demcels_ha_{base}.fits"), fits.getheader(f"{outdir}/demcels_ha_{base}.fits")
        s2 = fits.getdata(f"{outdir}/demcels_s2_{base}.fits")
        w = WCS(hh)
        cov = np.isfinite(ha) | np.isfinite(s2)
        ys, xs = np.where(cov)
        pad = 20
        sl = (slice(max(ys.min() - pad, 0), ys.max() + pad), slice(max(xs.min() - pad, 0), xs.max() + pad))
        sig = robust_sigma(ha)
        hs, ss = nansmooth(ha, 1), nansmooth(s2, 1)
        sigs = sig / 2.0  # white-noise reduction of a sigma=1 px Gaussian is ~1/(2 sqrt(pi) sigma)^0.5 ~ 0.53
        ratio = (ss * CAL_S2) / (hs * CAL_HA)
        ratio[~(hs > 5 * sigs)] = np.nan
        panels = [("H$\\alpha$ (continuum-subtracted)", ha, "ha"), ("[S II] 6717,31", s2, "s2"), ("[S II]/H$\\alpha$ (flux ratio)", ratio, "r")]
        fig = plt.figure(figsize=(21, 7.6))
        for i, (title, a, k) in enumerate(panels):
            ax = fig.add_subplot(1, 3, i + 1, projection=w)
            if k == "r":
                im = ax.imshow(a, origin="lower", cmap="RdBu_r", norm=TwoSlopeNorm(0.4, 0.0, 0.9))
            else:
                v = a[np.isfinite(a)]
                lo = max(np.percentile(v, 5), 0) + 1.0
                im = ax.imshow(np.clip(a, 0, None) + 1.0, origin="lower", cmap="inferno",
                               norm=LogNorm(lo, np.percentile(v, 99.8) + 1.0))
            ax.set_facecolor("0.15")
            draw_shells(ax, w, abs(hh["CDELT1"]))
            ax.set_xlim(sl[1].start, sl[1].stop); ax.set_ylim(sl[0].start, sl[0].stop)
            ax.set_title(title); ax.coords[0].set_axislabel("RA (J2000)"); ax.coords[1].set_axislabel("Dec (J2000)")
            ax.coords.grid(color="w", alpha=0.15)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label({"ha": "DN / 0.27\" pix (log, +1)", "s2": "DN / 0.27\" pix (log, +1)", "r": "[S II]/H$\\alpha$ (0.4 = shock threshold)"}[k])
        fig.suptitle(f"DeMCELS DR1 on the H I grid ({suf or 'raw median mosaic'}); circles: Book et al. 2008 supergiant shells LMC 1-9", y=0.98)
        fig.tight_layout()
        out = f"{pngdir}/demcels_quicklook{suf or '_raw'}.png"
        fig.savefig(out, dpi=110); plt.close(fig)
        print("wrote", out, "sigma20", sig)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:] or [""])
