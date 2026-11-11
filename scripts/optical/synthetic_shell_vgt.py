"""How much rim tangency does an expanding shell give VGT with no magnetic field at all?

Why: shell_rim_orientation.py finds the H I VGT direction tangential to the rims of LMC supergiant shells (T = 0.33-0.41
at LMC 1, 2, 4, 5). A swept-up wall is a ring of bright gas in every channel map, so channel intensity gradients point
radially and VGT (gradient + 90 deg) reads tangential with no field. This script measures that geometry-only effect on
field-free synthetic cubes and calibrates it against an observable the real shells also have: the tangency of the
column-density (moment-0) contours, "contour T". If field-free shells with a given contour T give a VGT T well below
the observed shells at similar contour T, the observed excess is not shell geometry.

Model (no magnetic field anywhere; everything statistically isotropic except the shell):
- Grid like the real LMC H I cube: 20" pixels (CDELT 0.0055556 deg), 1 km/s channels, 1' beam (FWHM 3 px).
  Box 600 x 600 px, velocity -70..+70 km/s (141 channels).
- Background H I layer: one Gaussian line per pixel, T_B(v) = N / (sqrt(2 pi) sigma) exp(-(v - v_c)^2 / 2 sigma^2)
  with N in K km/s.
  * N lognormal, sigma_ln = 0.4, mean N0 = 526 K km/s, so the peak T_B is ~35 K (mean) / ~32 K (median), the typical
    LMC peak brightness of 30-40 K.
  * v_c Gaussian, rms 4 km/s; sigma = 6 km/s (most LMC flux has sigma 3-10 km/s).
  * ln N and v_c are independent isotropic Gaussian random fields with P(k) ~ k^(-11/3) (Kolmogorov, projected),
    periodic on the box, amplitude ~ k^(-11/6) on white noise.
- Shell, centred at pixel (301.5, 301.5) (0-based; this is the exact centre of binned pixel (75, 75) after the 4x
  block binning used for scoring, so no sub-pixel offset enters T), R = 90 px = 0.5 deg (Book et al. 2008
  supergiant shells: radii 0.34-0.57 deg). Wall thickness 0.15 R, i.e. the wall is R_in = 0.925 R .. R_out = 1.075 R.
  Wall gas line width sigma 5 km/s, systemic velocity 0 (the box centre).
  (a) "sphere": thin spherical wall expanding at v_exp = 25 km/s. Front and back caps each carry half the path length
      L(r) = 2 [sqrt(R_out^2 - r^2) - sqrt(max(R_in^2 - r^2, 0))] through the wall, at LOS velocities
      -/+ v_exp sqrt(1 - (r/R)^2) (0 for r > R), each modulated by its own lognormal field (sigma_ln 0.4,
      same spectrum as the background; the two caps are different gas).
  (b) "ring": blown-out shell = cylindrical wall in the disk plane at inclination i = 30 deg (the LMC disk). The wall
      is a Gaussian in projected radius centred on R (FWHM 0.15 R), LOS velocity v_exp sin(i) cos(phi - phi0) with
      phi0 random per seed, one lognormal modulation field. The ring is kept circular on the sky (not foreshortened
      to b/a = cos i = 0.87) so it sits on the scoring annulus exactly as the circular Book et al. shells are scored;
      the 13% ellipticity would stay inside the 0.8-1.2 R annulus anyway.
  Both: background column inside R_in reduced to 20% (the evacuated interior), hard edge, softened by the beam.
- Rim contrast C = mean wall column / mean background column over the rim annulus R_in <= r <= R_out, using the
  realised modulation fields (so the label C is exact per cube). C in {0.25, 0.5, 1, 2}, seeds {0, 1, 2}; the same
  background realisation is used for all C and both geometries at a given seed, so cubes differ only by the shell.
  Plus 3 background-only cubes (same 3 backgrounds) as the null.
- Every channel convolved with the 3 px FWHM beam; then noise: white noise convolved with the same beam, rescaled to
  0.24 K rms per channel (the real cube's noise).
- FITS: NAXIS1 = x, NAXIS2 = y, NAXIS3 = v as the real cube; RA---TAN / DEC--TAN, CDELT1 = -0.0055556,
  CDELT2 = 0.0055556, VELO-LSR in m/s (CDELT3 = 1000).

Real-background variant (same shells injected into real LMC H I, so the background turbulence is real): the
single-Gaussian synthetic background has far less small-scale velocity structure than real channel maps, so its
geometry-only T is likely an upper bound.
- Background: 600 x 600 patches of the unclipped LMC crop REAL_CROP (159 channels of 1 km/s from 181.69 km/s, 20"
  pixels, already beam-convolved and noisy, noise 0.2355 K = the production VGT noise of this crop).
- Patch rule (select_patches, written to real_patches.json): patch corners on a 16 px grid; the patch has emission
  (moment 0 > 5 sigma, sigma = noise sqrt(159) K km/s) over >= 90% of its area; no Book et al. 2008 shell centre
  (SHELLS in shell_rim_orientation / demcels_quicklook, converted to crop pixels) lies within 1.5 R_shell of the
  injected ring (| distance to the ring centre - R | > 1.5 R_shell); among those, the most typical patches first
  (patch mean moment 0 closest to the mean over all emission pixels of the crop, which keeps the brightest
  30 Dor-like fields out), greedily, with no two patches overlapping. 2338 patches pass; crop emission mean moment 0
  400.7 K km/s. Chosen (crop pixels, 0-based corner y0, x0; ring centre = corner + 301.5; nearest Book shell in units
  of its radius from the injected ring):
    p0: y0 = 1040, x0 = 1360, emission 99.3%, mean moment 0 401 K km/s, LMC 1 at 3.09 R_shell
    p1: y0 =    0, x0 =  768, emission 100%,  mean moment 0 401 K km/s, LMC 9 at 2.91 R_shell
    p2: y0 = 1008, x0 =  320, emission 92.4%, mean moment 0 466 K km/s, LMC 4 at 1.54 R_shell (just passes)
- Injected shell: the same sphere / ring models and seeds (patch index p = seed), contrast C relative to the patch's
  own mean moment 0 over the rim annulus, systemic velocity = the patch's moment 1 at the shell centre (intensity-
  weighted mean velocity of the spectrum averaged over the central 10 x 10 px, channels above 3 sigma), walls
  beam-convolved (3 px FWHM) before adding. The interior hole is applied as for the synthetic cubes (data inside
  R_in scaled to 20%, edge softened by the beam); scaling the data also scales its noise, so beam-convolved noise of
  rms sigma sqrt(1 - f^2) is added there to keep the noise level flat (f = the hole factor; this is the only noise
  added, zero outside the hole). The 3 uninjected patches are the null. Each cube keeps the crop WCS (CRPIX shifted
  to the cut-out). VGT noise = 0.2355 (production value for this crop).

Hole follow-up (background "realhole", real background only, same 3 patches): every shell above carries the same
interior hole (data inside R_in at 20%), whose hard edge sits inside the scoring annulus, so contrast C alone never
isolates a weak wall. Interior factor HOLES_EXTRA = 1.0 (no hole) and 0.6, both geometries, C in CONTRASTS_HOLE
(0.1-1), 3 patches: 48 more cubes, scored into the same CSV (column hole; the earlier runs are hole 0.2). The real
uninjected patches are their null.

VGT: the production LMC configuration (data/sims/vgt_lmc_none_dn20/run_metadata.toml, [config]: dn 20, ker 4,
min_valid 100, axial_vonmises, central_checked, scalar noise, ...) with noise = 0.24 and the paths changed. ny_chunks
is "auto" (execution only; chunk stitching is exact). Each cube is one Slurm job of hpc/irgetas/vgt.slurm (main.jl).

Scoring (T at the exact centre, and T averaged over centres offset like the rounding of the Book et al. 2008
coordinates: OFFSET_RA_ARCMIN / OFFSET_DEC_ARCMIN below): the same functions as the real-shell test (shell_rim_orientation.nanbin, nansmooth, tangency_map; BIN = 4,
VGT_FWHM_ARCMIN = 13.7, MIN_COVER). psi is block-binned 4x as Stokes (as load_tracers); contour angle = gradient of the
moment-0 map (binned 4x, smoothed to 13.7' FWHM) + 90 deg, blanked where VGT is blank; T at the shell centre for the
annulus 0.8-1.2 R. Null: T at every valid centre of the background-only cubes, and on the shell cubes at centres whose
annulus stays > 10 binned px (~beam + VGT window + smoothing) outside the wall.

Everything runs on the cluster (hpc/irgetas/README.md); --out holds cubes/ (27 x 0.2 GB), configs/ (absolute paths,
so generate them where VGT runs), vgt/<cube>/, the CSV and the figures. The real-shell CSV (shell_rim_orientation.py
output) is read from --real-csv for the overplot.

Usage (from the repo root on the cluster, each step inside a Slurm job or submitting them):
  python scripts/optical/synthetic_shell_vgt.py make   --background synthetic|real [--out DIR]  # cubes + configs
  python scripts/optical/synthetic_shell_vgt.py submit --background synthetic|real [--after ID] # one vgt.slurm job/cube
  python scripts/optical/synthetic_shell_vgt.py score [--out DIR] [--real-csv F]  # all finished backgrounds: CSV, figs
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import time

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter, map_coordinates

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
from shell_rim_orientation import nanbin, nansmooth, tangency_map, BIN, VGT_FWHM_ARCMIN, MIN_COVER  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEFAULT_OUT = f"{REPO}/validation_output/synshell"

NPIX = 600
PIX_DEG = 0.0055556
VEL = np.arange(-70.0, 70.0 + 0.5, 1.0)          # km/s, 141 channels
BEAM_FWHM_PX = 3.0
NOISE_K = 0.24
SIGMA_LN = 0.4
N0 = 526.0                                       # K km/s: mean peak T_B = N0 / (sqrt(2 pi) 6) = 35 K
VC_RMS, SIGMA_BG = 4.0, 6.0
R_PX, WALL = 90.0, 0.15
R_IN, R_OUT = R_PX * (1 - WALL / 2), R_PX * (1 + WALL / 2)
XC = YC = 301.5                                  # centre of binned pixel (75, 75)
V_EXP, SIGMA_WALL, INCL_DEG, HOLE = 25.0, 5.0, 30.0, 0.2
CONTRASTS = (0.25, 0.5, 1.0, 2.0)
HOLES_EXTRA, CONTRASTS_HOLE = (1.0, 0.6), (0.1, 0.25, 0.5, 1.0)
SEEDS = (0, 1, 2)
GEOMS = ("sphere", "ring")
ANNULUS = (0.8, 1.2)
EXCLUDE_MARGIN_BIN = 10.0
REAL_CROP = f"{REPO}/validation_output/rohsa/lmc2048_n16_smooth10_noclip/crop_y1_x113/crop_y1_x113.fits"
REAL_NOISE = 0.2355293228260705                  # production VGT noise of this crop (vgt_lmc_none_dn20 metadata)
PATCH_STEP, EMIT_SIGMA, EMIT_FRAC, BOOK_EXCL = 16, 5.0, 0.9, 1.5
# Book et al. 2008 centres are rounded to 1 min RA (+-2.7' at Dec -69) and 10' Dec (+-5'), i.e. 0.1-0.25 R off a
# true centre, so T at the exact synthetic centre is an upper bound for what the real-shell test can see. "offset" T =
# mean of the T map (bilinear) over N_OFFSET centres drawn uniformly in that box; the same draws for every cube.
OFFSET_RA_ARCMIN, OFFSET_DEC_ARCMIN, N_OFFSET = 2.7, 5.0, 200
VGT_KEYS = dict(dn=20, ker=4, min_valid=100, fitter_mode="axial_vonmises", finalize_mask_mode="strict_nonpositive",
                gradient_mode="central_checked", gradient_sample_mode="paired_finite",
                histogram_edge_mode="modern_range", zero_stokes_policy="nan_before_smoothing",
                output_support_mode="smoothed_stokes_support", rbar2_weight=False, rbar2_gate=0.0)
NOISE_KEYS = dict(noise=NOISE_K, noise_mode="scalar", noise_tile_size=128, noise_sigma_multiplier=1.0,
                  noise_spatial_exponent=1.0, noise_channel_exponent=1.0)


def powerlaw_field(rng, n=NPIX, slope=-11.0 / 3.0):
    """Unit-variance isotropic Gaussian random field with P(k) ~ k^slope (periodic box)."""
    ky = np.fft.fftfreq(n)[:, None]
    kx = np.fft.rfftfreq(n)[None, :]
    k = np.hypot(kx, ky)
    amp = np.zeros_like(k)
    amp[k > 0] = k[k > 0] ** (slope / 2)
    f = np.fft.irfft2(np.fft.rfft2(rng.standard_normal((n, n))) * amp, s=(n, n))
    return (f - f.mean()) / f.std()


def lognormal(rng, s=SIGMA_LN):
    """Mean-1 lognormal field."""
    return np.exp(s * powerlaw_field(rng) - s * s / 2)


def gauss_line(v, vc, s):
    return np.exp(-0.5 * ((v - vc) / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def components(geom, contrast, seed):
    """List of (column K km/s, centroid km/s, sigma km/s) 2D components making up the cube."""
    rng = np.random.default_rng(1000 + seed)
    nbg = N0 * lognormal(rng)
    vc = VC_RMS * powerlaw_field(rng)
    if geom == "none":
        return [(nbg, vc, SIGMA_BG)]
    yy, xx = np.mgrid[0:NPIX, 0:NPIX].astype(float)
    r = np.hypot(xx - XC, yy - YC)
    nbg = np.where(r < R_IN, HOLE * nbg, nbg)
    return [(nbg, vc, SIGMA_BG)] + shell_walls(geom, contrast * N0, seed, 0.0)


def shell_walls(geom, target, seed, v0):
    """Wall components; target = mean wall column over the rim annulus, v0 = systemic velocity."""
    yy, xx = np.mgrid[0:NPIX, 0:NPIX].astype(float)
    r = np.hypot(xx - XC, yy - YC)
    rim = (r >= R_IN) & (r <= R_OUT)
    srng = np.random.default_rng(2000 + seed)     # shell randomness independent of the background draw
    if geom == "sphere":
        path = 2 * (np.sqrt(np.clip(R_OUT ** 2 - r ** 2, 0, None)) - np.sqrt(np.clip(R_IN ** 2 - r ** 2, 0, None)))
        mf, mb = lognormal(srng), lognormal(srng)
        raw_f, raw_b = 0.5 * path * mf, 0.5 * path * mb
        a = target / np.mean((raw_f + raw_b)[rim])
        vlos = V_EXP * np.sqrt(np.clip(1 - (r / R_PX) ** 2, 0, None))
        return [(a * raw_f, v0 - vlos, SIGMA_WALL), (a * raw_b, v0 + vlos, SIGMA_WALL)]
    if geom == "ring":
        sw = WALL * R_PX / 2.3548
        raw = np.exp(-0.5 * ((r - R_PX) / sw) ** 2) * lognormal(srng)
        a = target / np.mean(raw[rim])
        phi0 = srng.uniform(0, 2 * np.pi)
        vlos = V_EXP * np.sin(np.radians(INCL_DEG)) * np.cos(np.arctan2(yy - YC, xx - XC) - phi0)
        return [(a * raw, v0 + vlos, SIGMA_WALL)]
    raise ValueError(geom)


def build_cube(geom, contrast, seed):
    comps = components(geom, contrast, seed)
    sb = BEAM_FWHM_PX / 2.3548
    nrng = np.random.default_rng(3000 + seed)
    cube = np.empty((VEL.size, NPIX, NPIX), np.float32)
    for j, v in enumerate(VEL):
        ch = sum(n * gauss_line(v, c, s) for n, c, s in comps)
        noise = gaussian_filter(nrng.standard_normal((NPIX, NPIX)), sb, mode="wrap")
        cube[j] = gaussian_filter(ch, sb, mode="wrap") + noise * (NOISE_K / noise.std())
    return cube


def header():
    h = fits.Header()
    for k, v in (("CTYPE1", "RA---TAN"), ("CRVAL1", 80.0), ("CRPIX1", NPIX / 2 + 0.5), ("CDELT1", -PIX_DEG),
                 ("CUNIT1", "deg"), ("CTYPE2", "DEC--TAN"), ("CRVAL2", -69.0), ("CRPIX2", NPIX / 2 + 0.5),
                 ("CDELT2", PIX_DEG), ("CUNIT2", "deg"), ("CTYPE3", "VELO-LSR"), ("CRVAL3", VEL[0] * 1000.0),
                 ("CRPIX3", 1.0), ("CDELT3", 1000.0), ("CUNIT3", "m/s"), ("BUNIT", "K"),
                 ("BMAJ", BEAM_FWHM_PX * PIX_DEG), ("BMIN", BEAM_FWHM_PX * PIX_DEG), ("BPA", 0.0),
                 ("EQUINOX", 2000.0)):
        h[k] = v
    h["COMMENT"] = "Field-free synthetic expanding-shell cube (scripts/optical/synthetic_shell_vgt.py)"
    return h


def real_velocity(h):
    return (h["CRVAL3"] + (np.arange(h["NAXIS3"]) + 1 - h["CRPIX3"]) * h["CDELT3"]) / 1000.0


def select_patches(out):
    """Pick 3 real-background patches by the rule in the module docstring; writes real_patches.json."""
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    from demcels_quicklook import SHELLS, PC_PER_DEG
    with fits.open(REAL_CROP, memmap=True) as hd:
        h, data = hd[0].header, hd[0].data
        m0 = np.zeros(data.shape[1:])
        for j in range(data.shape[0]):
            m0 += np.nan_to_num(data[j].astype(float))
    nv = h["NAXIS3"]
    emit = m0 > EMIT_SIGMA * REAL_NOISE * np.sqrt(nv) * abs(h["CDELT3"]) / 1000.0
    ii_e = np.pad(np.cumsum(np.cumsum(emit, 0), 1), ((1, 0), (1, 0)))
    ii_m = np.pad(np.cumsum(np.cumsum(m0, 0), 1), ((1, 0), (1, 0)))
    box = lambda ii, y, x: ii[y + NPIX, x + NPIX] - ii[y, x + NPIX] - ii[y + NPIX, x] + ii[y, x]  # noqa: E731
    wcs = WCS(h, naxis=2)
    book = []
    for name, (ra, dec, d) in SHELLS.items():
        bx, by = (float(v) for v in wcs.world_to_pixel(SkyCoord(ra, dec, frame="icrs")))
        book.append((name, by, bx, d / 2 / PC_PER_DEG / abs(h["CDELT1"])))
    ref = m0[emit].mean()
    cands = []
    for y0 in range(0, m0.shape[0] - NPIX + 1, PATCH_STEP):
        for x0 in range(0, m0.shape[1] - NPIX + 1, PATCH_STEP):
            frac = box(ii_e, y0, x0) / NPIX ** 2
            if frac < EMIT_FRAC:
                continue
            yc, xc = y0 + YC, x0 + XC
            if any(abs(np.hypot(yc - by, xc - bx) - R_PX) <= BOOK_EXCL * rs for _, by, bx, rs in book):
                continue
            mean = box(ii_m, y0, x0) / NPIX ** 2
            cands.append((abs(mean - ref), y0, x0, frac, mean))
    cands.sort()
    chosen = []
    for _, y0, x0, frac, mean in cands:
        if all(abs(y0 - c["y0"]) >= NPIX or abs(x0 - c["x0"]) >= NPIX for c in chosen):
            chosen.append(dict(y0=y0, x0=x0, emit_frac=round(float(frac), 3), mean_m0=round(float(mean), 1)))
        if len(chosen) == len(SEEDS):
            break
    if len(chosen) < len(SEEDS):
        raise RuntimeError(f"only {len(chosen)} patches satisfy the rule ({len(cands)} candidates)")
    for c in chosen:
        c["nearest_book"] = min(((round(float(abs(np.hypot(c["y0"] + YC - by, c["x0"] + XC - bx) - R_PX) / rs), 2), n)
                                 for n, by, bx, rs in book))
    print(f"{len(cands)} candidate patches; crop emission mean moment 0 {ref:.1f} K km/s; chosen: {chosen}")
    return chosen, ref


def build_real_cube(geom, contrast, p, info, hole=HOLE):
    """Patch p of the real crop, with the field-free shell injected (geom 'none': the patch as is)."""
    with fits.open(REAL_CROP, memmap=True) as hd:
        h = hd[0].header.copy()
        data = hd[0].data[:, info["y0"]:info["y0"] + NPIX, info["x0"]:info["x0"] + NPIX].astype(float)
    h["CRPIX1"] -= info["x0"]
    h["CRPIX2"] -= info["y0"]
    h["COMMENT"] = (f"LMC crop patch y0={info['y0']} x0={info['x0']}"
                    + ("" if geom == "none" else f", field-free {geom} shell injected, C={contrast:g}, hole={hole:g}"))
    if geom == "none":
        return data.astype(np.float32), h
    vel = real_velocity(h)
    sb = BEAM_FWHM_PX / 2.3548
    yy, xx = np.mgrid[0:NPIX, 0:NPIX].astype(float)
    f = gaussian_filter(np.where(np.hypot(xx - XC, yy - YC) < R_IN, hole, 1.0), sb, mode="nearest")
    comp = np.sqrt(np.clip(1 - f ** 2, 0, None))
    walls = shell_walls(geom, contrast * info["rim_m0"], p, info["v_sys"])
    nrng = np.random.default_rng(4000 + p)
    for j, v in enumerate(vel):
        ch = sum(n * gauss_line(v, c, s) for n, c, s in walls)
        noise = gaussian_filter(nrng.standard_normal((NPIX, NPIX)), sb, mode="wrap")
        data[j] = data[j] * f + gaussian_filter(ch, sb, mode="nearest") + comp * noise * (REAL_NOISE / noise.std())
    return data.astype(np.float32), h


def patch_info(out):
    with open(os.path.join(out, "real_patches.json")) as fh:
        return json.load(fh)


def cube_list(bg):
    """(background, geometry, C, seed, hole) per cube; 'realhole' = the hole follow-up on the real patches."""
    if bg == "realhole":
        return [("real", g, c, s, h) for h in HOLES_EXTRA for g in GEOMS for c in CONTRASTS_HOLE for s in SEEDS]
    out = [(bg, "none", 0.0, s, None) for s in SEEDS]
    out += [(bg, g, c, s, HOLE) for g in GEOMS for c in CONTRASTS for s in SEEDS]
    return out


def name_of(bg, geom, contrast, seed, hole=HOLE):
    if bg == "real":
        if geom == "none":
            return f"real_bg_p{seed}"
        return f"real_{geom}_C{contrast:g}_p{seed}" if hole == HOLE else f"real_h{hole:g}_{geom}_C{contrast:g}_p{seed}"
    return f"bg_s{seed}" if geom == "none" else f"{geom}_C{contrast:g}_s{seed}"


def toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return f'"{v}"'
    return repr(v)


def write_config(path, cube_path, out_dir, noise):
    lines = ["# Production LMC VGT settings (data/sims/vgt_lmc_none_dn20/run_metadata.toml) on a field-free shell test cube.",
             "[input]", f'input_file = "{cube_path}"', f'output_dir = "{out_dir}"', "keep_intermediate_chunks = false",
             "", "[noise]"]
    lines += [f"{k} = {toml_value(noise if k == 'noise' else v)}" for k, v in NOISE_KEYS.items()]
    lines += ["", "[vgt]"] + [f"{k} = {toml_value(v)}" for k, v in VGT_KEYS.items()]
    lines += ["", "[execution]", 'ny_chunks = "auto"', ""]
    with open(path, "w") as f:
        f.write("\n".join(lines))


def make(out, bg):
    for d in ("cubes", "configs", "vgt"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    is_real = bg != "synthetic"
    if is_real:
        jpath = os.path.join(out, "real_patches.json")
        if not os.path.exists(jpath):
            assert bg == "real", "pick the patches with --background real first"
            chosen, _ = select_patches(out)
            for p, info in enumerate(chosen):
                patch, h = build_real_cube("none", 0.0, p, info)
                yy, xx = np.mgrid[0:NPIX, 0:NPIX].astype(float)
                r = np.hypot(xx - XC, yy - YC)
                m0 = np.nansum(patch, axis=0) * abs(h["CDELT3"]) / 1000.0
                info["rim_m0"] = float(m0[(r >= R_IN) & (r <= R_OUT)].mean())
                spec = np.nanmean(patch[:, int(YC) - 4:int(YC) + 6, int(XC) - 4:int(XC) + 6], axis=(1, 2))
                vel, w = real_velocity(h), np.where(spec > 3 * REAL_NOISE, spec, 0.0)
                info["v_sys"] = float(np.sum(vel * w) / np.sum(w))
            with open(jpath, "w") as fh:
                json.dump(chosen, fh, indent=1)
        infos = patch_info(out)
        print("real patches:", json.dumps(infos))
    for b, geom, c, s, hole in cube_list(bg):
        name = name_of(b, geom, c, s, hole)
        cpath = os.path.join(out, "cubes", f"{name}.fits")
        if not os.path.exists(cpath):
            t0 = time.time()
            if is_real:
                data, h = build_real_cube(geom, c, s, infos[s], HOLE if hole is None else hole)
            else:
                data, h = build_cube(geom, c, s), header()
            fits.PrimaryHDU(data, h).writeto(cpath, overwrite=True)
            print(f"{name}: cube written in {time.time() - t0:.0f} s", flush=True)
        write_config(os.path.join(out, "configs", f"{name}.toml"), cpath, os.path.join(out, "vgt", name),
                     REAL_NOISE if is_real else NOISE_K)


def psi_path(out, name):
    return os.path.join(out, "vgt", name, f"psi_{name}.fits")


def submit(out, bg, after=None):
    """One hpc/irgetas/vgt.slurm job per cube without a finished run (run_metadata.toml is written last)."""
    ids = []
    for b, geom, c, s, hole in cube_list(bg):
        name = name_of(b, geom, c, s, hole)
        if os.path.exists(os.path.join(out, "vgt", name, "run_metadata.toml")):
            continue
        cmd = ["sbatch", "--parsable", "--job-name=vgt_synshell", "--cpus-per-task=8", "--mem=16G", "--time=01:00:00",
               f"--export=ALL,VGT_REPO={REPO},VGT_CONFIG={os.path.join(out, 'configs', f'{name}.toml')}"]
        if after:
            cmd.append(f"--dependency=afterok:{after}")
        ids.append(subprocess.check_output(cmd + ["hpc/irgetas/vgt.slurm"], cwd=REPO, text=True).strip().split(";")[0])
        print(f"{name}: job {ids[-1]}", flush=True)
    print("VGT jobs:", ":".join(ids))


# ---------------------------------------------------------------- scoring
def binned_angles(out, name):
    """VGT psi binned 4x as Stokes, and the moment-0 contour angle at 13.7' (as load_tracers)."""
    psi = fits.getdata(psi_path(out, name)).astype(float)
    with np.errstate(invalid="ignore"):
        c, s = nanbin(np.cos(2 * psi), BIN), nanbin(np.sin(2 * psi), BIN)
    vgt = 0.5 * np.arctan2(s, c)
    cube = fits.getdata(os.path.join(out, "cubes", f"{name}.fits")).astype(float)
    m0 = np.nansum(cube, axis=0) * 1.0            # K km/s, 1 km/s channels
    m0[m0 <= 0] = np.nan
    pix_arcmin = PIX_DEG * BIN * 60
    nb = nansmooth(nanbin(m0, BIN), VGT_FWHM_ARCMIN / 2.3548 / pix_arcmin)
    gy, gx = np.gradient(nb)
    cont = np.arctan2(gy, gx) + np.pi / 2
    cont[~np.isfinite(vgt)] = np.nan
    return vgt, cont, nb, cube


def centre_bin():
    return int(round((YC - 1.5) / BIN)), int(round((XC - 1.5) / BIN))


def far_mask(shape, r_bin):
    """Centres whose scoring annulus stays EXCLUDE_MARGIN_BIN outside the shell wall."""
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    iy, ix = centre_bin()
    return np.hypot(yy - iy, xx - ix) > ANNULUS[1] * r_bin + R_OUT / BIN + EXCLUDE_MARGIN_BIN


def sanity_scoring(r_bin):
    n = NPIX // BIN
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    iy, ix = centre_bin()
    th = np.arctan2(yy - iy, xx - ix)
    t_tan = tangency_map(th + np.pi / 2, r_bin, *ANNULUS)[0][iy, ix]
    t_rad = tangency_map(th, r_bin, *ANNULUS)[0][iy, ix]
    return float(t_tan), float(t_rad)


def backgrounds(out):
    """Backgrounds whose VGT runs are all finished."""
    done = lambda bg: all(os.path.exists(os.path.join(out, "vgt", name_of(*k), "run_metadata.toml"))  # noqa: E731
                          for k in cube_list(bg))
    return [bg for bg in ("synthetic", "real", "realhole") if done(bg)]


def offset_mean(tmap, iy, ix):
    """Mean of the T map over N_OFFSET centres offset uniformly within the Book et al. coordinate rounding."""
    pix_arcmin = PIX_DEG * BIN * 60
    rng = np.random.default_rng(20261112)
    dy = rng.uniform(-1, 1, N_OFFSET) * OFFSET_DEC_ARCMIN / pix_arcmin
    dx = rng.uniform(-1, 1, N_OFFSET) * OFFSET_RA_ARCMIN / pix_arcmin
    return float(np.nanmean(map_coordinates(tmap, [iy + dy, ix + dx], order=1, mode="constant", cval=np.nan)))


def score_cubes(out, bg, rows, null_v, null_c, r_bin, iy, ix):
    """T at the exact and offset-sampled centre for every cube of background bg; appends to rows and prints means."""
    far_v, new = [], []
    for b, geom, c, s, hole in cube_list(bg):
        vgt, cont, _, _ = binned_angles(out, name_of(b, geom, c, s, hole))
        tv, cov = tangency_map(vgt, r_bin, *ANNULUS)
        tc = tangency_map(cont, r_bin, *ANNULUS)[0]
        if geom != "none":
            far_v.append(tv[far_mask(tv.shape, r_bin) & np.isfinite(tv)])
        new.append(dict(background=b, hole="" if hole is None else hole, geometry=geom, C=c, seed=s,
                        T_vgt_exact=round(float(tv[iy, ix]), 3), T_vgt_offset=round(offset_mean(tv, iy, ix), 3),
                        T_contour_exact=round(float(tc[iy, ix]), 3),
                        T_contour_offset=round(offset_mean(tc, iy, ix), 3), cover=round(float(cov[iy, ix]), 2),
                        null_mean=round(float(null_v.mean()), 3),
                        null_p95=round(float(np.percentile(null_v, 95)), 3),
                        null_contour_mean=round(float(null_c.mean()), 3),
                        null_contour_p95=round(float(np.percentile(null_c, 95)), 3)))
    far_v = np.concatenate(far_v)
    print(f"[{bg}] null on shell cubes, centres far from the wall ({far_v.size} centres): mean {far_v.mean():+.3f}, "
          f"95th pct {np.percentile(far_v, 95):+.3f}")
    print(f"[{bg}] means over seeds, exact centre / offset-sampled centre:")
    keys = sorted({(r["hole"], r["geometry"], r["C"]) for r in new}, key=lambda k: (str(k[0]), k[1], k[2]))
    for hole, g, c in keys:
        sel = [r for r in new if (r["hole"], r["geometry"], r["C"]) == (hole, g, c)]
        mean = lambda k: np.mean([r[k] for r in sel])  # noqa: E731
        seeds = ", ".join(f"{r['T_vgt_offset']:+.2f}" for r in sel)
        print(f"  hole={hole!s:4s} {g:6s} C={c:<4g}  T_vgt {mean('T_vgt_exact'):+.3f} / {mean('T_vgt_offset'):+.3f} "
              f"(offset per seed {seeds})  T_contour {mean('T_contour_exact'):+.3f} / {mean('T_contour_offset'):+.3f}")
    rows.extend(new)


def score(out, real_csv):
    r_bin = R_PX / BIN
    iy, ix = centre_bin()
    t_tan, t_rad = sanity_scoring(r_bin)
    print(f"sanity (i): tangential map T = {t_tan:+.4f}, radial map T = {t_rad:+.4f}")
    rows, nulls, nulls_full = [], {}, {}
    for bg in backgrounds(out):
        if bg == "realhole":
            score_cubes(out, "realhole", rows, *nulls_full["real"], r_bin, iy, ix)
            continue
        null_v, null_c, align, rank = [], [], [], []
        for s in SEEDS:
            vgt, cont, _, _ = binned_angles(out, name_of(bg, "none", 0, s))
            tv, tc = tangency_map(vgt, r_bin, *ANNULUS)[0], tangency_map(cont, r_bin, *ANNULUS)[0]
            ok = np.isfinite(tc)
            rank.append(float(np.mean(tc[ok] < tc[iy, ix])))
            null_v.append(tv[np.isfinite(tv)])
            null_c.append(tc[np.isfinite(tc)])
            m = np.isfinite(vgt) & np.isfinite(cont)
            align.append(float(np.mean(np.cos(2 * (vgt[m] - cont[m])))))
        null_v, null_c = np.concatenate(null_v), np.concatenate(null_c)
        nulls[bg] = null_v
        nulls_full[bg] = (null_v, null_c)
        print(f"\n[{bg}] sanity (ii): background-only VGT T at {null_v.size} centres: mean {null_v.mean():+.3f}, "
              f"95th pct {np.percentile(null_v, 95):+.3f}; contour T mean {null_c.mean():+.3f}, "
              f"95th pct {np.percentile(null_c, 95):+.3f}")
        print(f"[{bg}] no-shell contour T at the box centre: percentile among all centres of the same cube = "
              f"{', '.join(f'{x:.2f}' for x in rank)} (uniform 0-1 expected: the centre is not special)")
        print(f"[{bg}] convention check: <cos 2(psi_VGT - contour angle)> on background-only cubes = "
              f"{', '.join(f'{a:+.2f}' for a in align)} (> 0: same pixel-frame convention)")
        score_cubes(out, bg, rows, null_v, null_c, r_bin, iy, ix)
    with open(os.path.join(out, "synthetic_shell_tangency.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    plot_relation(out, rows, nulls, real_csv)
    plot_look(out, list(nulls))


def real_shells(path):
    if not path or not os.path.exists(path):
        print(f"no real-shell CSV at {path}; figure without the observed shells")
        return {}
    got = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["annulus"] == "0.8:1.2" and r["kind"] == "shell":
                got.setdefault(r["shell"], {})[r["tracer"]] = float(r["T"])
    return {k: (v["N_HI contours"], v["VGT (H I)"]) for k, v in got.items()}


def _shells_and_axes(ax, real, title):
    ax.plot([-1, 1], [-1, 1], color="#999", lw=0.8, ls=":", zorder=1)
    ax.axhline(0, color="#bbb", lw=0.6)
    ax.axvline(0, color="#bbb", lw=0.6)
    if real:
        ax.scatter([v[0] for v in real.values()], [v[1] for v in real.values()], s=70, marker="D", color="black",
                   zorder=4)
        nudge = {"LMC 3": (-38, 6)}                     # LMC 3 and LMC 8 sit 0.03 apart
        for k, (x, y) in real.items():
            ax.annotate(k, (x, y), xytext=nudge.get(k, (5, 4)), textcoords="offset points", fontsize=8.5)
    ax.set_xlim(-0.35, 0.95)
    ax.set_ylim(-0.25, 1.0)
    ax.set_xlabel("contour T  (moment-0 iso-contours vs rim tangent)")
    ax.set_ylabel("VGT T  (VGT direction vs rim tangent)")
    ax.set_title(title, fontsize=9)


def _null_lines(ax, nulls, styles):
    for k, (bg, null_v) in enumerate(nulls.items()):
        p95 = np.percentile(null_v, 95)
        ax.axhline(p95, color="#777", lw=0.9, ls=styles[bg], zorder=1)
        ax.text(0.98, p95 + (0.012 if k == 0 else -0.012),
                f"no shell, {bg} background: 95% of single-ring T_VGT below {p95:.2f}",
                ha="right", va="bottom" if k == 0 else "top", fontsize=7.5, color="#555",
                transform=ax.get_yaxis_transform())


def _mean_curve(rows, cs):
    return [(np.mean([r["T_contour_offset"] for r in rows if r["C"] == c]),
             np.mean([r["T_vgt_offset"] for r in rows if r["C"] == c])) for c in cs]


def plot_relation(out, rows, nulls, real_csv):
    """Left: hole 0.2, synthetic vs real background. Right: real background, interior hole 0.2 / 0.6 / none."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    col = {"sphere": "#1f5fae", "ring": "#e3a21a"}
    size = {0.1: 14, 0.25: 25, 0.5: 50, 1.0: 100, 2.0: 200}
    mark = {"synthetic": ("o", "-", "x"), "real": ("s", "--", "+")}
    real = real_shells(real_csv)
    holes = sorted({r["hole"] for r in rows if r["background"] == "real" and r["hole"] != ""}, reverse=True)
    fig, axes = plt.subplots(1, 2 if len(holes) > 1 else 1, figsize=(15.5 if len(holes) > 1 else 7.8, 7.0),
                             squeeze=False)
    ax = axes[0, 0]
    _null_lines(ax, nulls, {bg: mark[bg][1] for bg in nulls})
    base = [r for r in rows if r["hole"] in (HOLE, "")]
    for bg in nulls:
        m, ls, nm = mark[bg]
        for g in GEOMS:
            for c in CONTRASTS:
                sel = [r for r in base if r["background"] == bg and r["geometry"] == g and r["C"] == c]
                ax.scatter([r["T_contour_exact"] for r in sel], [r["T_vgt_exact"] for r in sel], s=size[c], marker=m,
                           color=col[g], lw=0, alpha=0.18, zorder=2)
                ax.scatter([r["T_contour_offset"] for r in sel], [r["T_vgt_offset"] for r in sel], s=size[c],
                           marker=m, color=col[g], edgecolor="white", lw=1.0, alpha=0.9, zorder=3)
            ax.plot(*zip(*_mean_curve([r for r in base if r["background"] == bg and r["geometry"] == g], CONTRASTS)),
                    color=col[g], lw=1.2, ls=ls, zorder=2)
        sel = [r for r in base if r["background"] == bg and r["geometry"] == "none"]
        ax.scatter([r["T_contour_offset"] for r in sel], [r["T_vgt_offset"] for r in sel], s=45, marker=nm,
                   color="#555", zorder=3)
    handles = [Line2D([], [], color=col[g], lw=6, label=g) for g in GEOMS]
    handles += [Line2D([], [], marker=mark[bg][0], ls=mark[bg][1], color="#555", label=f"{bg} background")
                for bg in nulls]
    handles += [Line2D([], [], marker="o", ls="", color="#999", markersize=np.sqrt(size[c]), label=f"C = {c:g}")
                for c in CONTRASTS]
    handles += [Line2D([], [], marker=mark[bg][2], ls="", color="#555", label=f"no shell ({bg})") for bg in nulls]
    if real:
        handles.append(Line2D([], [], marker="D", ls="", color="black", label="LMC shells (observed)"))
    ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="upper left", ncol=2)
    _shells_and_axes(ax, real, "Field-free shells with the interior at 20% vs the observed LMC shells (annulus 0.8-1.2 R)"
                     "\nsolid: centre offset as the rounded Book et al. coordinates; faint: exact centre;"
                     "\n3 seeds / patches each, lines = means over them")
    if len(holes) > 1:
        ax = axes[0, 1]
        hcol = dict(zip(sorted(holes), ("#08463d", "#2a8f80", "#7fc4b8")))   # darker = more evacuated interior
        gm = {"sphere": ("o", "-"), "ring": ("^", "--")}
        _null_lines(ax, {"real": nulls["real"]}, {"real": "--"})
        rr = [r for r in rows if r["background"] == "real" and r["geometry"] != "none"]
        for h in holes:
            for g in GEOMS:
                sel = [r for r in rr if r["hole"] == h and r["geometry"] == g]
                cs = sorted({r["C"] for r in sel})
                ax.scatter([r["T_contour_offset"] for r in sel], [r["T_vgt_offset"] for r in sel],
                           s=[size[r["C"]] for r in sel], marker=gm[g][0], color=hcol[h], edgecolor="white", lw=0.8,
                           alpha=0.9, zorder=3)
                ax.plot(*zip(*_mean_curve(sel, cs)), color=hcol[h], lw=1.4, ls=gm[g][1], zorder=2)
        handles = [Line2D([], [], color=hcol[h], lw=6,
                          label=("no hole" if h == 1.0 else f"interior at {h:.0%}")) for h in holes]
        handles += [Line2D([], [], marker=gm[g][0], ls=gm[g][1], color="#555", label=g) for g in GEOMS]
        handles += [Line2D([], [], marker="o", ls="", color="#999", markersize=np.sqrt(size[c]), label=f"C = {c:g}")
                    for c in sorted({r["C"] for r in rr})]
        if real:
            handles.append(Line2D([], [], marker="D", ls="", color="black", label="LMC shells (observed)"))
        ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="lower right", ncol=2)
        _shells_and_axes(ax, real, "Real LMC background: how deep the interior hole is\n"
                         "centre offset as the rounded Book et al. coordinates; 3 patches each, lines = means")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "synthetic_shell_T_relation.png"), dpi=140)
    plt.close(fig)


def plot_look(out, bgs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    panels = [(bg, g, dv) for bg in bgs for g, dv in (("sphere", -15.0), ("ring", 8.0))]
    fig, axes = plt.subplots(len(panels), 3, figsize=(15, 4.8 * len(panels)), squeeze=False)
    th = np.linspace(0, 2 * np.pi, 200)
    for row, (bg, geom, dv) in zip(axes, panels):
        name = name_of(bg, geom, 1.0, 0)
        vgt, cont, nb, cube = binned_angles(out, name)
        if bg == "real":
            vel = real_velocity(fits.getheader(os.path.join(out, "cubes", f"{name}.fits")))
            v0 = patch_info(out)[0]["v_sys"]
        else:
            vel, v0 = VEL, 0.0
        j = int(np.argmin(np.abs(vel - (v0 + dv))))
        im = row[0].imshow(cube[j], origin="lower", cmap="Greys_r")
        fig.colorbar(im, ax=row[0], shrink=0.8, label="T_B (K)")
        row[0].set_title(f"{name}\nchannel {vel[j]:.0f} km/s (systemic {dv:+.0f} km/s)", fontsize=10)
        m0 = np.nansum(cube, axis=0)
        im = row[1].imshow(m0, origin="lower", cmap="Greys_r")
        fig.colorbar(im, ax=row[1], shrink=0.8, label="moment 0 (K km/s)")
        row[1].set_title("moment 0", fontsize=10)
        iy, ix = centre_bin()
        rb = R_PX / BIN
        h = int(1.7 * rb)                              # zoom like shell_rim_orientation.plot
        im = row[2].imshow(nb, origin="lower", cmap="Greys_r")
        fig.colorbar(im, ax=row[2], shrink=0.8, label="moment 0 at 13.7' (K km/s)")
        yy, xx = np.mgrid[iy - h:iy + h + 1:2, ix - h:ix + h + 1:2]
        a = vgt[yy, xx]
        ok = np.isfinite(a)
        L = 0.85
        seg = np.stack([np.c_[xx[ok] - L * np.cos(a[ok]), yy[ok] - L * np.sin(a[ok])],
                        np.c_[xx[ok] + L * np.cos(a[ok]), yy[ok] + L * np.sin(a[ok])]], axis=1)
        row[2].add_collection(LineCollection(seg, colors="#1f5fae", linewidths=1.2))
        for f in ANNULUS:
            row[2].plot(ix + f * rb * np.cos(th), iy + f * rb * np.sin(th), "--", color="#e3a21a", lw=1.2)
        row[2].set_xlim(ix - h, ix + h)
        row[2].set_ylim(iy - h, iy + h)
        row[2].set_title("VGT direction (binned 4x), zoom +-1.7 R\ndashed: scoring annulus 0.8-1.2 R", fontsize=10)
        for ax in row:
            ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Field-free shells, contrast C = 1, seed / patch 0", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "synthetic_shell_look.png"), dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("make", "submit", "score"))
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--background", choices=("synthetic", "real", "realhole"), default="synthetic")
    ap.add_argument("--after", default=None, help="Slurm job id the VGT jobs wait for (afterok), e.g. the make job")
    ap.add_argument("--real-csv", default=None, help="shell_rim_tangency.csv of the observed shells (overplot)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.step == "make":
        make(args.out, args.background)
    elif args.step == "submit":
        submit(args.out, args.background, args.after)
    else:
        score(args.out, args.real_csv)


if __name__ == "__main__":
    main()
