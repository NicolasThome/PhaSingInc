# Zenodo data layout

Extract the accompanying dataset into a directory of your choice and pass it
with `--data-dir`. The default is `data/` in the repository:

```text
data/
  experimental/
    minR_meanF_cont.xlsx
    9_points.xlsx
    FC_FFT_data.xlsx
  numerical/
    K_continuation_mean_frequency.csv
    K_continuation_min_radius.csv
    K_continuation_max_radius.csv
    K_continuation_order_parameters.csv
    K_continuation_restart/
      K_continuation_K_1p0.jld2
      ... (151 files down to K_continuation_K_0p85.jld2)
  derived/
    fig3_rotating_frame.jld2
    fig4_fft.jld2
  supplementary/
    experimental/
      KOP_Exp.xlsx
      min(r).xlsx
      Book2.xlsx
    numerical/
      locked_spectrum.jld2
      locked_spectrum_summary.csv
      locked_spectrum_recovery_diagnostics.csv
  SHA256SUMS
```

Include Figure 1's large experimental file in a separately documented directory
in the Zenodo record. No script in this repository reads that file.

## Figure inputs

| Figure | Data |
| --- | --- |
| 1 | Separately supplied raw experimental data; no plotting script |
| 2 | Continuation frequencies and radius extrema; `9_points.xlsx` and `minR_meanF_cont.xlsx` |
| 3 | Rotating-frame trajectories in `fig3_rotating_frame.jld2` |
| 4 | Numerical spectra in `fig4_fft.jld2`; experimental spectra in `FC_FFT_data.xlsx` |
| 5 | Continuation frequencies, minimum radii and cluster order parameters |

Preserve experimental workbook sheet order and measured values. Supplementary
inputs are described in [SUPPLEMENTARY.md](SUPPLEMENTARY.md).

## Main simulation parameters

| Setting | Value |
| --- | --- |
| Oscillators | N = 5000 |
| Coupling sweep | K = 1.0 to 0.85, 151 points, decreasing |
| C1, C2 | −0.5, −7.0 |
| Transient | 3000 at every K, including the first |
| Recording time / sampling interval | 600 / 0.1 |
| Relative / absolute tolerance | 1e−8 / 1e−8 |
| Frequency distribution | Gaussian, mean 0, standard deviation 0.3, symmetric pairs |
| Frequency seed / initial-state seed | 1234 / 5678 |
| Continuation noise | 0 |
| Amplitude-cluster threshold | 0.7 |
| Cluster diagnostic relative / absolute tolerance | 1e−4 / 1e−4 |

Each coupling point starts from the preceding final state. Preserve sweep order,
oscillator ordering and seeds to reproduce the selected branch. Figure 5 uses
the mean-frequency order parameter and transition markers K_HB=0.963 and K_def=0.948.

## Formats

CSV rows follow the decreasing K sweep. Radius files have K plus 5000 oscillator
columns. Mean-frequency CSV has K, system_mean_frequency, then 5000 oscillator
columns. Column names encode the oscillator index and rounded intrinsic frequency.
The order-parameter file includes overall, high-cluster and low-cluster diagnostics
and the high/low population counts.

Restarts store K, C1, C2, N, frequency vector, initial and final states, recording
settings, tolerances and seeds. States alternate real and imaginary components.
Frequency and state vectors are ordered consistently across K.

Figure 3 derived data contain the rotating-frame trajectory, unwrapped relative
phase, time vector, rotation rate, modulation period, winding diagnostic and
fixed-point status for K=0.970, 0.955, 0.948 and 0.880. Starting from saved final
states, the simulation uses transient 100, recording time 20, sampling 0.001 and
the restart tolerances. The low oscillator index is 199. The reference oscillator
is selected from high-amplitude oscillators with intrinsic frequency closest to zero.

Figure 4 uses K=1.000, 0.955 and 0.880, transient 150, recording time 700,
sampling 0.05 and tolerances 1e−7. Oscillator selection follows the original
Hopf/radius criteria. Spectra are one-sided amplitude FFTs, with the mean removed,
no window and no extra normalization. Numerical frequency axes are multiplied by
2π when plotted. Experimental frequencies retain their measured units.

Generated checksums describe the files in that particular archive. A newly
computed trajectory may differ slightly across hardware or solver versions;
scientific reproduction should compare numerical observables and figure structure,
rather than demand byte-identical numerical files.
