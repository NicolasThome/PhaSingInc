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

## Provenance

The source continuation dataset is
`C1_-0.5_C2_m7_Hopf_all_N_5000_threshold_0p7_Rfreq`.
All 151 restarts report N=5000, 5000 frequencies and 10000 real state components.
The old constant containing `N10000` pointed to this N=5000 dataset and was a
naming error. The current development script with N=2000 is not the paper configuration.

The three experimental workbooks come from the original Figure 2 and Figure 4
directories. Preserve workbook sheet order and the measured values.
`9_points.xlsx` supplies the selected experimental oscillator curves;
`minR_meanF_cont.xlsx` remains the companion measurement input read by Figure 2.

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
