"""List DeMCELS DR1 stacked continuum-subtracted tiles (Halpha, [S II]) via the NOIRLab SIA service.

What: probes https://datalab.noirlab.edu/sia/demcels_dr1 on a regular sky grid and collects every
distinct tile product `<field>_<tile>.ha_sub_r.dr01.fits` (Halpha, N662 minus scaled r) and
`<field>_<tile>.s2_sub_r.dr01.fits` ([S II], N673 minus scaled r) with its SIA centre.
Inputs: none (network). Output: JSON {"tiles": {siaRef: [obs_id, ra, dec]}, "hit": {"ra,dec": n_rows, -1 = failed query}}.
Assumptions: grid spacing (0.25 deg) is below the tile size (0.65 deg), so every tile contains a probe
point; failed queries are recorded as -1 and must be re-run (rerun until none remain).

Usage: python demcels_tiles.py OUT.json [ra0 ra1 dec0 dec1 step]
"""
import json, re, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np


def query(p):
    ra, dec = p
    u = f"https://datalab.noirlab.edu/sia/demcels_dr1?POS={ra:.3f},{dec:.3f}&SIZE=0.1"
    s = None
    for _ in range(4):
        try:
            s = urllib.request.urlopen(u, timeout=90).read().decode()
            break
        except Exception:
            pass
    if s is None:
        return p, None
    out = []
    for r in re.findall(r"<TR>(.*?)</TR>", s, re.S):
        td = re.findall(r"<TD>(.*?)</TD>", r, re.S)
        if len(td) > 20 and ("ha_sub_r" in td[1] or "s2_sub_r" in td[1]):
            out.append((td[18], td[1].split("siaRef=")[1].split("&")[0], float(td[13]), float(td[14])))
    return p, out


def main():
    outp = sys.argv[1]
    ra0, ra1, dec0, dec1, step = (float(x) for x in sys.argv[2:7]) if len(sys.argv) > 6 else (64, 98, -74, -63.4, 0.25)
    pts = [(ra, dec) for dec in np.arange(dec0, dec1, step)
           for ra in np.arange(ra0, ra1, step / np.cos(np.radians(dec)))]
    tiles, hit = {}, {}
    with ThreadPoolExecutor(16) as ex:
        for p, o in ex.map(query, pts):
            hit[f"{p[0]:.2f},{p[1]:.2f}"] = -1 if o is None else len(o)
            for t, f, ra, dec in (o or []):
                tiles[f] = (t, ra, dec)
    json.dump({"tiles": tiles, "hit": hit}, open(outp, "w"))
    print(len(pts), "probes;", sum(v < 0 for v in hit.values()), "failed;", len(tiles), "tile products")


if __name__ == "__main__":
    main()
