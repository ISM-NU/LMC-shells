# LMC-shells

Experiments repo, shared with the team. It holds the research on the optical
supergiant shells of the LMC and on how the magnetic-field proxies from VGT
and Planck relate to the shell rims. It fetches and mosaics DeMCELS H-alpha and
[S II] tiles onto the H I grid, tests whether the VGT and Planck field
directions run tangential to the shell rims, and runs field-free synthetic
shells through VGT to measure the tangency that gas geometry alone produces.
Python only.

## Setup

Usually this repository sits in the LMC workspace folder (`~/Projects/LMC`), whose `setup.sh` clones all related repositories next to each other and links the shared data; its AGENTS.md has the rules that span repositories.

```bash
git clone --recursive git@github.com:ISM-NU/LMC-shells.git
cd LMC-shells
ln -s ../lmc-data data          # shared data folder ~/Projects/LMC/lmc-data, ignored by git
git submodule foreach --recursive 'ln -sfn "$toplevel/data" data'   # same data inside submodules
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

Find the data in
`deps/LMC-foreground/docs/drive-index.md`.

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

`synthetic_shell_vgt.py submit` submits one `deps/VGT/hpc/irgetas/vgt.slurm`
job per synthetic cube on Irgetas. Run the `make` and `score` steps as the
script describes. Do heavy steps on the cluster.

## Layout

- `scripts/optical/` - DeMCELS tile fetch, mosaic and quick-look; `shell_rim_orientation.py` (rim tangency test against VGT and Planck); `synthetic_shell_vgt.py` (field-free synthetic shells).
- `deps/VGT` - submodule, the VGT method (ISM-NU/VGT-Denis, branch `split/vgt-method`).
- `deps/LMC-foreground` - submodule. Provides `scripts/compare_vgt_polarization.py` and the data inventory.

## Pinned method versions

- VGT: `deps/VGT`, branch `split/vgt-method`.
- LMC-foreground: `deps/LMC-foreground`.
- Update with `git -C deps/<name> pull`, then commit the new pointer.

## Related repos

- [VGT-Denis](https://github.com/ISM-NU/VGT-Denis) - the VGT method.
- [ROHSA-LMC](https://github.com/ISM-NU/ROHSA-LMC) - the ROHSA wrapper (method).
- [LMC-foreground](https://github.com/ISM-NU/LMC-foreground) - foreground, Planck and starlight comparisons. Pinned here.
- [LMC-cold-gas](https://github.com/ISM-NU/LMC-cold-gas) - ROHSA runs and cold H I work.
