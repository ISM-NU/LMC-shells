# AGENTS.md

## What this repo is

Research on the optical supergiant shells of the LMC and how the field proxies
(H I VGT, Planck) relate to the shell rims. Read `README.md` and the docstrings
of `scripts/optical/shell_rim_orientation.py` and `synthetic_shell_vgt.py`
first.

## Scope

- Belongs here: DeMCELS data handling, rim tests, synthetic shells, figures and notes for the shell work.
- Does not belong here: method code. VGT goes to ISM-NU/VGT-Denis, ROHSA to ISM-NU/ROHSA-LMC. Foreground and starlight comparison code goes to LMC-foreground. Cold-gas work goes to LMC-cold-gas.
- This repo imports `scripts/compare_vgt_polarization.py` from `deps/LMC-foreground`. Do not edit anything under `deps/`.

## Language

Use Julia where performance needs it. Simple scripts can be Python. This repo is Python only.

## Research rules

- Before a new analysis, state the research question and what each outcome would change. Look at maps before numbers. Do not conclude from 30 Doradus alone.
- Check the literature first and cite papers instead of re-deriving them.
- VGT is the group's method. Never test it as if it might be wrong. Explain disagreements physically (for example, swept-up shell gas makes channel gradients radial) and propose improvements.
- Planck is a comparison tracer with its own assumptions, not ground truth. It is not foreground-subtracted in the rim test.
- Every control (field-free shells, contour tangency, random-ring null) must be stated with its reason. Every bin or cut needs a reason, or comparable bins.
- Do not change a statistic or a mask to make a result look better.
- Show raw maps before derived figures.

## Data and Slurm

- Data live in `~/Projects/LMC/lmc-data`. `data` in the repo is a symlink made with `ln -s ../lmc-data data`. It is ignored by git. Never commit tiles, mosaics, cubes or run outputs.
- Heavy computation (VGT on synthetic cubes) runs on Irgetas via Slurm, not on laptops. Use `deps/VGT/hpc/irgetas/vgt.slurm`.
- Generated products go to ignored directories (`data/...`, `validation_output/`).

## Docs

- New scripts get a docstring with why, inputs, outputs and usage.
- There is no `docs/` yet. If durable evidence appears, add `docs/` with a status line per file and a `docs/README.md`.
