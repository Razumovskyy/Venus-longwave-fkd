# Venus Longwave FKDM K-Distribution

Venus longwave (10–6000 cm-1) absorption parameterization used here to speed up radiative transfer. CO2, H2O, and SO2 absorption is mapped to 32 effective cross-sections; A 32-point k-solution covers the band with errors targetted to be within a few percent versus the reference Monte Carlo model below 90 km.

## Layout
- `app/main.f90` — Reads `Atm_Description`, loops over atmospheres, loads profiles, runs KD/VAC/Planck pipelines, writes outputs.
- `src/Atmospheres.f90` — setup and utilities (initializes settings, reads profiles, prepares pressure/temperature/mixing ratio arrays).
- `src/KD.f90` — core KD/VAC/Planck logic (loads W_STAND and KD tables, optional CO2 T-correction, interpolates VAC, Planck lookup).
- `src/co2_corr.f90`, `src/vac.f90`, `src/planck.f90` — helpers factored from the legacy VAC_Planck path.
- Data tables: `data/F_LW/` (W_STAND, O.CO2, O.H2O, O.SO2, KD-Planck direct-access, WS(T)).
- Atmospheres: `data/atmospheres/HAUS.00`, `data/atmospheres/VIRA.*` (altitude, pressure, temperature, gas mixing ratios).

## Inputs
- `Atm_Description` in repo root: counts (atmospheres, gases, levels) and gas names.
- Atmosphere profiles in `data/atmospheres/`: read by Atmospheres routines when `atm_settings.ini` points to them.
- KD data in `data/F_LW/`: read-only tables for k-grid, gas cross-sections, Planck integrals, and CO2 temperature correction.

## Outputs
- `VAC` — VAC profiles per level (written where the executable is run).
- `PLANCK` — Planck-integrated values per level (same location).

## Run
```sh
fpm run
```
Ensure `data/F_LW/` and `ATM-REs/` are present relative to the run directory.
