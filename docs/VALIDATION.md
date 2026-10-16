# Validation

Checks completed on Windows with Julia 1.12.1 on 13 October 2026:

- 60 test assertions covered model dynamics, seeded frequency sampling,
  trajectory diagnostics, FFT conventions, restart continuity, CSV output,
  figure selection, rotating equilibria and eigenmode analysis.
- The N=5000 seeded frequencies and initial state matched the archived first restart
  exactly. One paper-sized core continuation step at K=1.0 matched archived mean
  frequencies and the final state exactly.
- Figure 3 trajectories were recomputed from saved states. All 24 Figure 4
  numerical FFT arrays matched the archived spectra exactly.
- The four numerical main figures and supplementary plots rendered as PDF and
  PNG and were checked against the source figures.
- Installation, tests, small-data generation and plotting worked in a standalone
  copy. A supplementary stability check at N=12 completed three points with
  equilibrium residuals below 3e−15.

## PhaSingInc update, 16 October 2026

All 63 test assertions passed after the package rename. Main and supplementary
smoke pipelines completed. Figures 2, 3, 4 and 6 and all five supplementary
outputs rendered from the research data, with Figure 3/4 derived data recomputed.
The temporary validation workbook used the documented I=0.89 worksheet label.
Figure 6 now uses the same sampled phase-singularity marker as Figure 2 (K=0.947).

## Limits

The complete N=5000, 151-point continuation and the full N=5000 supplementary
stability scan were not independently rerun during preparation. The single
continuation-step check did not cover every paper-sized cluster diagnostic.
Small checks validate the pipeline, not the full paper results.

Linux rendering and numerical equivalence have not been established; the supplied
Ubuntu GitHub Actions workflow was not run locally. Fonts, hardware and dependency
versions may change pixels or fine numerical details.
