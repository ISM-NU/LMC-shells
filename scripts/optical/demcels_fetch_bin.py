"""Download DeMCELS DR1 continuum-subtracted tiles and block-reduce them immediately.

What: for each tile product (Halpha `ha_sub_r`, [S II] `s2_sub_r`) fetch the full tile (8700x8700 px of
0.27", ~300 MB) through the NOIRLab Data Lab cutout service, block-average (and block-median) it by
BIN x BIN (default 30 -> 8.1"), save a small .npz and delete the full-resolution file.
Inputs: tiles JSON from demcels_tiles.py; output dir; optional workers/bin factor.
Output per tile: OUTDIR/<siaRef>.npz with `mean`, `median` (float32, DN per 0.27" pixel, NaN-aware; constant no-data fill
masked by mask_fill(); bins with <50% valid pixels are NaN)
and the binned WCS (crval, crpix, cd, naxis) of the tile.
Assumptions: the products are already star-subtracted and continuum-subtracted by the DeMCELS pipeline
(PROCTYPE=StarSubtracted, SFILTER=r); flux units are DN per 0.27" pixel. Resumable (skips done tiles).
The cutout is requested at the SIA centre with SIZE 0.7 deg, larger than the tile, so the service returns the
whole tile; partial tiles are trimmed to a multiple of BIN pixels.

Usage: python demcels_fetch_bin.py tiles.json OUTDIR RAWDIR [workers=8] [bin=30]
"""
import json, os, subprocess, sys, time
import numpy as np
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor


def mask_fill(d, minfrac=2e-3, grow=4):
    """NaN-out the constant fill value that DR1 tiles carry where a tile has no data.

    The pipeline pads uncovered parts of a tile (chip gaps, field edges; up to ~30% of a tile) with ONE exactly
    repeated float (e.g. 53.37 DN in LMC_c01_T01 Halpha), which differs per tile and is not zero. Real float32 sky
    data never repeats a value in >0.2% of pixels, so the most frequent value is the fill when it exceeds that
    fraction. The mask is grown by `grow` px to drop the Lanczos ringing next to the boundary.
    Returns (data with NaN, fill value or None, masked fraction).
    """
    from scipy.ndimage import maximum_filter1d
    sub = d[::5, ::5]
    sub = sub[np.isfinite(sub)]
    v, c = np.unique(sub, return_counts=True)
    j = int(np.argmax(c))
    if c[j] < minfrac * sub.size:
        return d, None, 0.0
    m = (d == v[j]).astype(np.uint8)
    n = 2 * grow + 1
    m = maximum_filter1d(maximum_filter1d(m, n, axis=0), n, axis=1).astype(bool)
    d = d.copy()
    d[m] = np.nan
    return d, float(v[j]), float(m.mean())


def block(a, b, fn, minvalid=0.5):
    ny, nx = (a.shape[0] // b) * b, (a.shape[1] // b) * b
    v = a[:ny, :nx].reshape(ny // b, b, nx // b, b).swapaxes(1, 2).reshape(ny // b, nx // b, b * b)
    with np.errstate(all="ignore"), __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore")
        out = fn(v, axis=2).astype(np.float32)
        out[np.isfinite(v).mean(axis=2) < minvalid] = np.nan  # bins mostly outside the data
    return out


def work(job):
    name, ra, dec, outdir, rawdir, b = job
    outp = os.path.join(outdir, name.replace(".fits", ".npz"))
    if os.path.exists(outp):
        return name, "skip", 0
    raw = os.path.join(rawdir, name)
    url = (f"https://datalab.noirlab.edu/svc/cutout?col=demcels_dr1&siaRef={name}&extn=0"
           f"&POS={ra},{dec}&SIZE=0.7,0.7")
    t0 = time.time()
    for k in range(4):
        r = subprocess.run(["curl", "-s", "-f", "-m", "900", "-o", raw, url])
        if r.returncode == 0 and os.path.getsize(raw) > 1e6:
            break
        time.sleep(10 * (k + 1))
    else:
        if os.path.exists(raw): os.remove(raw)
        return name, "FAILED download", 0
    try:
        with fits.open(raw, memmap=False) as h:
            d, hd = h[0].data.astype(np.float32), h[0].header
        ny, nx = d.shape
        d, fillv, fillfrac = mask_fill(d)
        # bin pixel i' (1-based) covers old pixels b*(i'-1)+1 .. b*i'  -> centre at old b*(i'-1)+(b+1)/2
        crpix1 = (hd["CRPIX1"] - (b + 1) / 2) / b + 1
        crpix2 = (hd["CRPIX2"] - (b + 1) / 2) / b + 1
        np.savez_compressed(
            outp + ".tmp.npz", mean=block(d, b, np.nanmean), median=block(d, b, np.nanmedian),
            crval=[hd["CRVAL1"], hd["CRVAL2"]], crpix=[crpix1, crpix2],
            cd=np.array([[hd["CD1_1"], hd.get("CD1_2", 0)], [hd.get("CD2_1", 0), hd["CD2_2"]]]) * b,
            naxis=[nx, ny], fill=np.nan if fillv is None else fillv, fillfrac=fillfrac)
        os.replace(outp + ".tmp.npz", outp)
    finally:
        os.remove(raw)
    return name, f"ok fill={fillv} masked={fillfrac:.3f}", time.time() - t0


if __name__ == "__main__":
    tj, outdir, rawdir = sys.argv[1:4]
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    b = int(sys.argv[5]) if len(sys.argv) > 5 else 30
    os.makedirs(outdir, exist_ok=True); os.makedirs(rawdir, exist_ok=True)
    tiles = json.load(open(tj))["tiles"]
    jobs = [(n, v[1], v[2], outdir, rawdir, b) for n, v in sorted(tiles.items())]
    with ProcessPoolExecutor(nw) as ex:
        for i, (n, st, dt) in enumerate(ex.map(work, jobs)):
            print(f"[{i+1}/{len(jobs)}] {n} {st} {dt:.0f}s", flush=True)
