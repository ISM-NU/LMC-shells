# LMC-shells

Experiments repo, shared with the team. It holds the research on the optical
supergiant shells of the LMC and on how the magnetic-field proxies from VGT
and Planck relate to the shell rims. It fetches and mosaics DeMCELS H-alpha and
[S II] tiles onto the H I grid, tests whether the VGT and Planck field
directions run tangential to the shell rims, and runs field-free synthetic
shells through VGT to measure the tangency that gas geometry alone produces.
Python only.

## Setup

Group members usually work from the LMC workspace (the private ISM-NU/LMC-workspace repository), whose `setup.sh` clones the related repositories side by side and links the shared `lmc-data` folder.

The scripts use two repositories cloned next to this one: ISM-NU/VGT.jl at `../VGT.jl` and ISM-NU/LMC-foreground at `../LMC-foreground`. Outside the group, LMC-foreground cannot be cloned yet (it is private until publication), so clone only VGT.jl:

```bash
git clone https://github.com/ISM-NU/LMC-shells.git
git clone https://github.com/ISM-NU/VGT.jl.git
cd LMC-shells
ln -s ../lmc-data data          # shared data folder next to the repo, ignored by git
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

Group members find the data in
`../LMC-foreground/docs/drive-index.md`.

## Run

From the repository root. Each script describes its purpose, inputs, outputs
and usage in its docstring.

```bash
python3 scripts/optical/demcels_tiles.py OUT.json [ra0 ra1 dec0 dec1 step]      # list DeMCELS tiles (network)
python3 scripts/optical/demcels_fetch_bin.py tiles.json OUTDIR RAWDIR [workers=8] [bin=30]
python3 scripts/optical/demcels_mosaic.py BINNED_DIR HI_CUBE.fits OUT_DIR [bands=ha,s2] [grids=20,10]
python3 scripts/optical/demcels_quicklook.py OUT_DIR PNG_DIR [suffix ...]
python3 scripts/optical/shell_rim_orientation.py --out DIR [--halpha demcels_ha_on_hi_grid.fits]
python3 scripts/optical/synthetic_shell_vgt.py {make|submit|score} --out DIR
```

`synthetic_shell_vgt.py submit` submits one `../VGT.jl/hpc/irgetas/vgt.slurm`
job per synthetic cube on Irgetas. Run the `make` and `score` steps as the
script describes. Do heavy steps on the cluster.

## Layout

- `scripts/optical/` - DeMCELS tile fetch, mosaic and quick-look; `shell_rim_orientation.py` (rim tangency test against VGT and Planck); `synthetic_shell_vgt.py` (field-free synthetic shells).

## Method versions

- VGT comes from `../VGT.jl` (normally `main`), LMC-foreground from `../LMC-foreground`, which provides `scripts/compare_vgt_polarization.py` and the data inventory. ISM-NU/LMC-foreground is private for now (unpublished results), so outside readers can read these scripts but not run the rim test until it opens.
- There is no pin: the scripts use whatever commit those checkouts are on. Record the VGT.jl commit with each result; versions used in a report are tagged in the method repos (`report-YYYY-MM-DD`).

## Related repos

- [VGT.jl](https://github.com/ISM-NU/VGT.jl) - the VGT method.
- [ROHSA-LMC](https://github.com/ISM-NU/ROHSA-LMC) - the ROHSA wrapper (method).
- LMC-foreground (private until publication) - foreground, Planck and starlight comparisons. Pinned here.
- LMC-cold-gas (private until publication) - ROHSA runs and cold H I work.

## License

MIT (`LICENSE`). DeMCELS images, Planck maps and VGT maps used here are not
redistributed; the scripts download or read them from `data/`.
