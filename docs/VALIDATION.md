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
- All main and supplementary figures rendered as PDF and PNG. Figures 2, 3, 5
  and the supplementary plots matched reference PNGs pixel-for-pixel on this
  installation; Figure 4 was visually checked.
- Installation, tests, small-data generation and plotting worked in a standalone
  copy. A supplementary stability check at N=12 completed three points with
  equilibrium residuals below 3e−15.

## Limits

The complete N=5000, 151-point continuation and the full N=5000 supplementary
stability scan were not independently rerun during preparation. The single
continuation-step check did not cover every paper-sized cluster diagnostic.
Small checks validate the pipeline, not the full paper results.

Linux rendering and numerical equivalence have not been established; the supplied
Ubuntu GitHub Actions workflow was not run locally. Fonts, hardware and dependency
versions may change pixels or fine numerical details.
