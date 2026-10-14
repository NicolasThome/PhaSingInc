# Validation

Validated on Windows with Julia 1.12.1 on 13 October 2026.

## Numerical checks

- 60 test assertions cover model dynamics, deterministic symmetric frequency
  sampling, streaming diagnostics versus an independently saved trajectory,
  FFT amplitude/frequency conventions, restart continuity, CSV output and
  paper configuration, main/supplementary selection, rotating equilibria,
  the neutral phase mode, normalized eigenmode weights and spectrum filtering.
- The seeded N=5000 frequency vector and initial state match the archived
  K=1.0 restart exactly.
- A paper-sized **core** continuation step at K=1.0, with transient 3000 and
  recording time 600, matches the archived mean frequencies and final state
  exactly (maximum absolute differences both zero). This check does not include
  a full 151-point rerun or validation of every paper-sized cluster diagnostic.
- Figure 3 rotating-frame data were recomputed from the archived restarts.
- All 24 Figure 4 numerical FFT frequency/spectrum arrays match the original
  N=5000 FFT cache exactly (maximum absolute differences zero).

## Figure checks

All four main figures render as PDF and PNG, with the intended panel contents,
labels and dimensions. Figures 2, 3 and 5 are pixel-identical to the existing
paper PNGs on this Windows installation using Times New Roman.

Figure 4 preserves the current plotting source and its FFT arrays. Its PNG is
not pixel-identical to the older `Fig_3_FFT_grid.png` alternative output. The
current research tree has a modified Figure 4 script/PDF and no `Fig_4.png`;
the older alternative is not an exact current presentation reference. The new
Figure 4 has been visually inspected.

## Standalone pipeline

The code was copied to a separate temporary directory, with no sibling SLH or
FreqClusters folders. Environment instantiation, tests, small-data generation
and plotting were checked there. The manifest contains no sibling-package paths;
its own `AmpDefInc` entry uses the relative path `.`.

The small CLI run uses N=32 and produces a diagnostic PDF/PNG from saved CSV
data. Paper plotting loads saved derived data rather than integrating the model.
The two CLI entrypoints also validate arguments and report missing data files.
An Ubuntu GitHub Actions smoke workflow is supplied but was not run locally.

## Timing and limits

Measured derived-data generation times were approximately 23 seconds for Figure 3
and 32 seconds for Figure 4. The first paper-sized core continuation check took
approximately 60 seconds including Julia compilation and concurrent validation
activity. It reported about 564 MiB of cumulative allocations, **not peak RAM**.
These measurements are specific to this machine and are not full-run estimates.

A complete independent N=5000, 151-point continuation was not run during repository
preparation. Figures were reproduced from the archived continuation, with derived
Figure 3/4 calculations rerun. Linux rendering and numerical equivalence have not
been established. Different fonts, hardware or dependency versions may change
pixels or fine numerical details.

## Supplementary extension

All five supplementary outputs render as PDF and PNG and are pixel-identical
to the existing script outputs on this Windows installation. The `min(r)`
time-series keeps the correct I=0.89 label. Its previously misnamed worksheet
has been renamed to I=0.89 in the publication dataset, as confirmed by the author.
Only workbook-name metadata changed; measurement, fit and style XML entries are
byte-identical to the original workbook.

The archived supplementary spectrum contains N=5000 and 51 accepted K points.
Its frequency vector and first initial state match the seeded standalone code
exactly. A full generation/plotting smoke scan at N=12 completed all three
requested points, with equilibrium residuals below 3e−15. Its first-K FFT
stationarity check passed. Full N=5000 dense spectrum regeneration was not run.
