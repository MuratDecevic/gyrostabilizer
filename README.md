# Gyrostabilizer Simulation

This MATLAB project simulates the roll motion of a small boat in waves, with
and without a gyrostabilizer.

## Requirements

- MATLAB
- No additional toolboxes are required

Open MATLAB in this project directory before running the scripts.

## Simulate the boat with a gyrostabilizer

Run:

```matlab
gyro
```

The script evaluates several combinations of gyroscope damping (`Bg`) and
restoring torque (`Cg`). It opens one figure for each combination. Every figure
shows:

- solid line: boat roll angle
- dashed line: gyroscope precession angle
- vertical axis: angle in degrees
- horizontal axis: time in seconds

Calculated results are saved in `solutions.mat`. Later runs reuse this file
when the model and simulation settings have not changed. If they have changed,
MATLAB prints `Cached solutions are stale. Recalculating.` and creates new
results. A complete first run may take some time and opens 49 figures.

## Simulate the boat without a gyrostabilizer

Run:

```matlab
noGyro
```

This opens one figure showing the boat's roll angle in degrees over time. Use
this result as the baseline for comparison with the gyrostabilizer results.

## Supporting files

- `dynamics.m` contains the boat and gyroscope equations of motion.
- `noGyrodynamics.m` contains the boat-only equation of motion.
- `waveForcing.m` creates the deterministic JONSWAP wave excitation used by
  both simulations.
- `solutions.mat` is the generated cache of gyrostabilizer simulation results.

Normally, run only `gyro.m` or `noGyro.m`; the other MATLAB files are called
automatically.
