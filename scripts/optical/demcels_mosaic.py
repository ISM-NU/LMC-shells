"""Mosaic binned DeMCELS tiles onto the H I celestial grid (flux/intensity-preserving).

What: reads the per-tile .npz from demcels_fetch_bin.py, and reprojects + co-adds them
(reproject.reproject_exact, mean of overlaps) onto the celestial grid of the H I cube
(RA---TAN/DEC--TAN, 20" pixels), and onto a 2x finer 10" grid with the same centre.
Writes, per band (ha, s2), per grid (20", 10"), per statistic (median, mean of each 30x30 block):
  demcels_<band>_on_hi_grid[_10arcsec][_<stat>].fits       raw mosaic (no background matching)
  demcels_<band>_on_hi_grid[_10arcsec][_<stat>]_bgmatch.fits   with per-tile additive offset matching
  demcels_coverage_on_hi_grid[_10arcsec].fits              uint8, number of tiles with data (0 = hole), per band
The default product (no stat suffix) is the block *median* raw mosaic (robust to residual star cores);
the *_mean* versions are the plain block average.
Units: DN per 0.27" pixel, i.e. a surface brightness; DeMCELS DR1 calibration is 1 DN ~ 2.4e-18 erg/cm2/s
(Halpha) and ~2.1e-18 erg/cm2/s ([S II]) per 0.27" pixel (Points et al. 2024) -> multiply by
2.4e-18/0.27^2 = 3.29e-17 erg/s/cm2/arcsec2 (Halpha), 2.1e-18/0.27^2 = 2.88e-17 ([S II]).
Assumptions: tiles are continuum- and star-subtracted by the DeMCELS pipeline (ha_sub_r = N662 - scaled r;
s2_sub_r = N673 - scaled r). Field-to-field background offsets exist; bgmatch solves one additive offset
per tile from overlaps (reproject match_background) and is intended for display only.
NaN = not covered.

Usage: python demcels_mosaic.py BINNED_DIR HI_CUBE.fits OUT_DIR [bands=ha,s2] [grids=20,10] [stats=median,mean | cov (coverage only)] [erode=8]
(The 10" grid takes ~8 min per mosaic; run it with stats=median only.)
"""
import glob, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_exact
from reproject.mosaicking import reproject_and_coadd

HI_KEEP = ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "EQUINOX")


def target(hi_path, fine):
    h = fits.getheader(hi_path)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    w = WCS(h).celestial
    if fine:  # 2x finer, same centre/orientation: pixel edges stay aligned with the 20" grid
        hh = w.to_header()
        hh["CDELT1"] /= 2; hh["CDELT2"] /= 2
        hh["CRPIX1"] = (h["CRPIX1"] - 0.5) * 2 + 0.5; hh["CRPIX2"] = (h["CRPIX2"] - 0.5) * 2 + 0.5
        w, nx, ny = WCS(hh), nx * 2, ny * 2
    return w, (ny, nx)


def trim_edges(a, k):
    """NaN-out bins within k bins of a NaN (no-data) region of the tile: DR1 stacks have thin strips along the
    boundary of the covered area (1-2 exposures, imperfect continuum/sky subtraction) that show up as bright
    or dark stripes. The tile's own outer border is not treated as an edge."""
    if k <= 0:
        return a
    from scipy.ndimage import binary_dilation
    bad = binary_dilation(~np.isfinite(a), iterations=k)
    return np.where(bad, np.nan, a).astype(np.float32)


def tile_wcs(z):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval, w.wcs.crpix, w.wcs.cd = z["crval"], z["crpix"], z["cd"]
    return w


def main(binned, hi, outdir, bands="ha,s2", grids="20,10", stats="median,mean", erode="8"):
    os.makedirs(outdir, exist_ok=True)
    for band, key in (("ha", "ha_sub_r"), ("s2", "s2_sub_r")):
        if band not in bands.split(","):
            continue
        files = sorted(glob.glob(f"{binned}/*.{key}.dr01.npz"))
        tiles = [np.load(f) for f in files]
        print(band, len(tiles), "tiles", flush=True)
        for fine in [g == "10" for g in grids.split(",")]:
            wt, shp = target(hi, fine)
            tag = "_10arcsec" if fine else ""
            ncov = None
            for stat in stats.split(","):
                inp = [(trim_edges(t["median" if stat == "cov" else stat], int(erode)), tile_wcs(t)) for t in tiles]
                for bg in ((False, True) if stat != "cov" else ()):
                    arr, fp = reproject_and_coadd(
                        inp, wt, shape_out=shp, reproject_function=reproject_exact,
                        combine_function="mean", match_background=bg, block_size=(1024, 1024))
                    arr = np.where(fp >= 0.5, arr, np.nan).astype(np.float32)
                    hdr = wt.to_header()
                    hdr["BUNIT"] = "DN/pix0.27arcsec"
                    hdr["COMMENT"] = f"DeMCELS DR1 {key}, 30x30 block {stat}, reproject_exact mosaic"
                    hdr["COMMENT"] = "DN->erg/cm2/s per 0.27arcsec pix: 2.4e-18 (Ha), 2.1e-18 (SII)"
                    hdr["BGMATCH"] = (bg, "per-tile additive background matching applied")
                    hdr["NTILES"] = len(tiles)
                    hdr["ERODE"] = (int(erode), "bins (8.1 arcsec) trimmed next to no-data fill in each tile")
                    name = f"demcels_{band}_on_hi_grid{tag}" + ("" if stat == "median" else "_mean") + ("_bgmatch" if bg else "")
                    fits.writeto(f"{outdir}/{name}.fits", arr, hdr, overwrite=True)
                    print("wrote", name, np.isfinite(arr).mean(), flush=True)
                if stat in ("median", "cov"):  # coverage count (tiles contributing), once per band/grid
                    cnt = np.zeros(shp, np.uint8)
                    for a, w in inp:
                        f = reproject_exact((np.isfinite(a).astype(np.float32), w), wt, shape_out=shp, return_footprint=False)
                        cnt += (np.nan_to_num(f) >= 0.5).astype(np.uint8)
                    hdr = wt.to_header(); hdr["BUNIT"] = "n_tiles"
                    fits.writeto(f"{outdir}/demcels_coverage_{band}_on_hi_grid{tag}.fits", cnt, hdr, overwrite=True)


if __name__ == "__main__":
    main(*sys.argv[1:8])
