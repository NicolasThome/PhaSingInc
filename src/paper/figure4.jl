module Figure4
using ..AmpDefInc
using CairoMakie
using FFTW
using JLD2
using LaTeXStrings
using OrdinaryDiffEqTsit5
using Printf
using SciMLBase
using Statistics


using ..Figure2: read_first_sheet_xlsx, read_mean_frequency_csv, read_numeric_csv
const PAPER_FONT = AmpDefInc.PAPER_FONT

const FIG3_DIR = OUTPUT_ROOT
const EXPERIMENT_XLSX = joinpath(DATA_ROOT, "experimental", "FC_FFT_data.xlsx")
const NUMERICAL_DIR = joinpath(DATA_ROOT, "numerical")
const RESTART_DIR = joinpath(NUMERICAL_DIR, "K_continuation_restart")
const MEAN_FREQUENCY_CSV = joinpath(NUMERICAL_DIR, "K_continuation_mean_frequency.csv")
const MIN_RADIUS_CSV = joinpath(NUMERICAL_DIR, "K_continuation_min_radius.csv")
const FFT_BASE_TYPOGRAPHY_TEST = lowercase(get(ENV, "BASE_TYPOGRAPHY_TEST", "false")) == "true"
const FFT_TYPOGRAPHY_TEST_SUFFIX = FFT_BASE_TYPOGRAPHY_TEST ? "_base_typography_test" : ""
const FIG4_OUTPUT_STEM = joinpath(FIG3_DIR, "Fig_4" * FFT_TYPOGRAPHY_TEST_SUFFIX)
const NUMERIC_CACHE = joinpath(DATA_ROOT, "derived", "fig4_fft.jld2")
const NUMERIC_SELECTION_REPORT = joinpath(FIG3_DIR, "Fig_3_FFT_numeric_oscillator_selection.txt")

# The workbook contains a, c, e1, e2. The user requested e2 and to ignore e1.
# Add :e1 here if you want a literal 4-column workbook mirror.
const COLUMN_SPECS = [
    (key = :a, label = "(a)", numerical_K = 1.0, numerical_role = :one_frequency),
    (key = :c, label = "(b)", numerical_K = 0.955, numerical_role = :post_hopf),
    (key = :e2, label = "(c)", numerical_K = 0.88, numerical_role = :post_zero),
]

const FORCE_RECOMPUTE_NUMERICS = false
const NUMERIC_TTR = 150.0
const NUMERIC_T = 700.0
const NUMERIC_DT = 0.05
const FFT_WINDOW = nothing
const FFT_SPECTRUM = :amplitude
const FFT_NORMALIZE = false
const FFT_DETREND = true
const NUMERIC_ANGULAR_FREQUENCY_SCALE = 2pi
const FFT_XLIMS = NUMERIC_ANGULAR_FREQUENCY_SCALE .* (0.0, 1.25)
const NUMERIC_FFT_YLIMS = (0.0, 0.3)
const EXPERIMENT_FFT_XLIMS = (0.0, 0.1)
const EXPERIMENT_FFT_YLIMS = (0.0, 1.0)
const EXPERIMENT_FFT_XTICKS = [0.0, 0.05]
const NUMERIC_FFT_XTICKS = (
    [0.0, pi, 2pi],
    [L"0", L"\pi", L"2\pi"],
)
const EXPERIMENT_FFT_YTICKS = ([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])
const NUMERIC_FFT_YTICKS = ([0.0, 0.15, 0.3], ["0.00", "0.15", "0.30"])
const XI_INSET_XLIMS = (0.02, 0.03)
const XI_INSET_YLIMS = (0.0, 0.4)
const XI_INSET_XTICKS = ([0.02, 0.03], ["0.02", "0.03"])
const A_INSET_XLIMS = (0.0, 0.01)
const A_INSET_YLIMS = (0.0, 0.25)
const A_INSET_XTICKS = ([0.0, 0.005, 0.01], ["0", "0.005", "0.01"])
const LOW_CLUSTER_R_MIN_THRESHOLD = 0.05
const HOPF_BASELINE_MIN_K = 0.975
const HOPF_PROXY_THRESHOLD = 0.015
const R_ZERO_THRESHOLD = 0.012

const FFT_LATEX_BASE_PT = 12.0
const LATEX_TEXT_WIDTH_IN = 6.7
const TEX_PT_PER_IN = 72.27
const FFT_WIDTH_PT = LATEX_TEXT_WIDTH_IN * TEX_PT_PER_IN
const FFT_FIGURE_SIZE_PT = (FFT_WIDTH_PT, 0.74 * FFT_WIDTH_PT)
const FFT_TICK_PT = 1.00 * FFT_LATEX_BASE_PT
const FFT_LABEL_PT = 1.00 * FFT_LATEX_BASE_PT
const FFT_TITLE_PT = 1.00 * FFT_LATEX_BASE_PT
const FFT_LINE_PT = 1.25
const NUMERIC_FFT_LINE_PT = 1.50
const FFT_ZOOM_BOX_LINE_PT = 0.70
const PEAK_LABEL_PT = 0.90 * FFT_LABEL_PT
const PEAK_LEADER_LINE_PT = 0.55
const BLUE = RGBf(0.0, 0.32, 0.62) # high-amplitude cluster
const RED = RGBf(0.84, 0.22, 0.13)  # low-amplitude cluster

function restart_file_for_K(K::Real)
    files = filter(name -> startswith(name, "K_continuation_K_") && endswith(name, ".jld2"), readdir(RESTART_DIR))
    isempty(files) && error("No restart files found in $RESTART_DIR.")
    values = Float64[]
    for name in files
        token = replace(replace(name, "K_continuation_K_" => ""), ".jld2" => "")
        push!(values, parse(Float64, replace(token, "p" => ".")))
    end
    idx = argmin(abs.(values .- Float64(K)))
    abs(values[idx] - Float64(K)) <= 1.0e-8 || @warn "Using nearest restart K=$(values[idx]) for requested K=$K"
    return joinpath(RESTART_DIR, files[idx])
end

function fft_style_axis!(ax)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.09)
    ax.ygridcolor = (:black, 0.09)
    ax.xticklabelsize = FFT_TICK_PT
    ax.yticklabelsize = FFT_TICK_PT
    ax.xlabelsize = FFT_LABEL_PT
    ax.ylabelsize = FFT_LABEL_PT
    ax.titlesize = FFT_TITLE_PT
    ax.spinewidth = 0.55
    ax.xtickwidth = 0.55
    ax.ytickwidth = 0.55
    ax.xticksize = 2.4
    ax.yticksize = 2.4
    return ax
end

function nonmissing_xy(x, y)
    xx = Float64[]
    yy = Float64[]
    for (a, b) in zip(x, y)
        a === missing && continue
        b === missing && continue
        push!(xx, Float64(a))
        push!(yy, Float64(b))
    end
    return xx, yy
end

function experiment_column_map(key::Symbol)
    if key === :a
        return (f = 1, fft_b = 2, fft_r = nothing, afft_b = 3, afft_r = nothing)
    elseif key === :c
        return (f = 6, fft_b = 7, fft_r = 8, afft_b = 9, afft_r = 10)
    elseif key === :e1
        return (f = 13, fft_b = 14, fft_r = 15, afft_b = 16, afft_r = 17)
    elseif key === :e2
        return (f = 13, fft_b = 20, fft_r = 21, afft_b = 22, afft_r = 23)
    else
        error("Unknown experimental key $(key).")
    end
end

function read_experiment_fft(path::AbstractString, key::Symbol)
    sheet = read_first_sheet_xlsx(path)
    cols = experiment_column_map(key)
    f = sheet[3:end, cols.f]
    b = sheet[3:end, cols.fft_b]
    ab = sheet[3:end, cols.afft_b]
    fr, r = cols.fft_r === nothing ? (Float64[], Float64[]) : nonmissing_xy(f, sheet[3:end, cols.fft_r])
    far, ar = cols.afft_r === nothing ? (Float64[], Float64[]) : nonmissing_xy(f, sheet[3:end, cols.afft_r])
    fb, sb = nonmissing_xy(f, b)
    fab, sab = nonmissing_xy(f, ab)
    return (emsi = [(f = fb, spectrum = sb, color = BLUE, label = "high")],
            radii = [(f = fab, spectrum = sab, color = BLUE, label = "high")],
            emsi_red = (f = fr, spectrum = r, color = RED, label = "low"),
            radii_red = (f = far, spectrum = ar, color = RED, label = "low"))
end

function add_optional_red!(series, red_series)
    isempty(red_series.f) || push!(series, red_series)
    return series
end

function first_along_decreasing_k(k_values, values, predicate)
    order = sortperm(k_values; rev = true)
    for idx in order
        predicate(values[idx]) && return k_values[idx]
    end
    return NaN
end

function baseline_residual(k_values, min_radius)
    baseline_rows = findall(>=(HOPF_BASELINE_MIN_K), k_values)
    X = hcat(ones(length(baseline_rows)), k_values[baseline_rows])
    residual = similar(min_radius)
    for j in axes(min_radius, 2)
        beta = X \ min_radius[baseline_rows, j]
        baseline = beta[1] .+ beta[2] .* k_values
        residual[:, j] .= max.(0.0, baseline .- min_radius[:, j])
    end
    return residual
end

function choose_numerical_indices()
    k_values, mean_frequency, omega, Ω = read_mean_frequency_csv(MEAN_FREQUENCY_CSV)
    k_min, min_radius, omega_min = read_numeric_csv(MIN_RADIUS_CSV; first_data_col = 2)
    k_values == k_min || error("Numerical K grids do not match.")
    omega == omega_min || error("Numerical omega columns do not match.")

    low_cluster = findall(vec(minimum(min_radius; dims = 1)) .<= LOW_CLUSTER_R_MIN_THRESHOLD)
    residual = baseline_residual(k_values, min_radius)
    k_h = fill(NaN, length(omega))
    k_0 = fill(NaN, length(omega))
    for idx in low_cluster
        k_h[idx] = first_along_decreasing_k(k_values, residual[:, idx], value -> value > HOPF_PROXY_THRESHOLD)
        k_0[idx] = first_along_decreasing_k(k_values, min_radius[:, idx], value -> value <= R_ZERO_THRESHOLD)
    end

    kc = 0.955
    ic = argmin(abs.(k_values .- kc))
    ke = 0.88
    ie = argmin(abs.(k_values .- ke))

    low_candidates = [
        idx for idx in low_cluster
        if isfinite(k_h[idx]) &&
           isfinite(k_0[idx]) &&
           k_h[idx] >= kc &&
           k_0[idx] >= ke &&
           min_radius[ic, idx] > R_ZERO_THRESHOLD
    ]
    isempty(low_candidates) && error("No low-amplitude oscillator satisfies Hopf-before-c and r_min=0-before-e2 criteria.")
    low_scores = [(-k_h[idx], -k_0[idx], idx) for idx in low_candidates]
    low_idx = low_candidates[argmin(low_scores)]

    high_candidates = setdiff(collect(eachindex(omega)), low_cluster)
    high_scores = vec(mean(min_radius[[1, ic, ie], high_candidates]; dims = 1))
    high_idx = high_candidates[argmax(high_scores)]

    return (;
        omega,
        k_values,
        min_radius,
        k_h,
        k_0,
        low_cluster,
        low_idx,
        high_idx,
        check_k_values = (1.0, kc, ke),
        check_k_indices = (argmin(abs.(k_values .- 1.0)), ic, ie),
    )
end

function fft_result_for_signal(signal, t)
    F = AmpDefInc.fft_matrix(
        reshape(signal, :, 1),
        t;
        detrend = FFT_DETREND,
        window = FFT_WINDOW,
        spectrum = FFT_SPECTRUM,
        onesided = true,
        normalize = FFT_NORMALIZE,
    )
    return F.freqs, vec(F.spectra[1, :])
end

function load_restart(path::AbstractString)
    data = JLD2.load(path)
    return (
        K = Float64(data["K"]),
        C1 = Float64(data["C1"]),
        C2 = Float64(data["C2"]),
        omega = Float64.(data["omega"]),
        u0 = Float64.(data["final_state"]),
    )
end

function simulate_selected_restart(K::Real, selected_indices)
    restart = load_restart(restart_file_for_K(K))
    save_idxs = Int[]
    for idx in selected_indices
        push!(save_idxs, 2idx - 1)
        push!(save_idxs, 2idx)
    end
    save_idxs = sort(unique(save_idxs))
    p = (restart.K, restart.C1, restart.C2, restart.omega)
    prob = SciMLBase.ODEProblem(AmpDefInc.Model.sl_harmo_hetero!, restart.u0, (0.0, NUMERIC_TTR + NUMERIC_T), p)
    saveat = NUMERIC_TTR:NUMERIC_DT:(NUMERIC_TTR + NUMERIC_T)
    sol = SciMLBase.solve(prob, OrdinaryDiffEqTsit5.Tsit5(); saveat, save_idxs, reltol = 1.0e-7, abstol = 1.0e-7)
    mat = Array(sol)'
    t = collect(sol.t) .- first(sol.t)
    colmap = Dict(save_idxs[i] => i for i in eachindex(save_idxs))
    return t, mat, colmap
end

function spectra_for_K(K::Real, entries)
    selected = unique(vcat([entry.indices for entry in entries]...))
    t, mat, colmap = simulate_selected_restart(K, selected)
    result = Dict{String, Any}()
    result["K"] = Float64(K)
    for entry in entries
        idx = entry.idx
        x = mat[:, colmap[2idx - 1]]
        y = mat[:, colmap[2idx]]
        r = hypot.(x, y)
        fx, sx = fft_result_for_signal(x, t)
        fr, sr = fft_result_for_signal(r, t)
        result["$(entry.name)_real_freqs"] = fx
        result["$(entry.name)_real_spectrum"] = sx
        result["$(entry.name)_radii_freqs"] = fr
        result["$(entry.name)_radii_spectrum"] = sr
        result["$(entry.name)_idx"] = idx
    end
    return result
end

function numerical_spectra(; generate::Bool = false)
    selectors = choose_numerical_indices()
    specs = Dict{Symbol, Any}(
        :a => [
            (name = "high", idx = selectors.high_idx, indices = [selectors.high_idx], color = BLUE, label = "high"),
            (name = "low", idx = selectors.low_idx, indices = [selectors.low_idx], color = RED, label = "low"),
        ],
        :c => [
            (name = "high", idx = selectors.high_idx, indices = [selectors.high_idx], color = BLUE, label = "high"),
            (name = "low", idx = selectors.low_idx, indices = [selectors.low_idx], color = RED, label = "low"),
        ],
        :e1 => [
            (name = "high", idx = selectors.high_idx, indices = [selectors.high_idx], color = BLUE, label = "high"),
            (name = "low", idx = selectors.low_idx, indices = [selectors.low_idx], color = RED, label = "low"),
        ],
        :e2 => [
            (name = "high", idx = selectors.high_idx, indices = [selectors.high_idx], color = BLUE, label = "high"),
            (name = "low", idx = selectors.low_idx, indices = [selectors.low_idx], color = RED, label = "low"),
        ],
    )
    selection_signature = Dict(
        "low_idx" => selectors.low_idx,
        "high_idx" => selectors.high_idx,
        "low_color" => "red",
        "high_color" => "blue",
        "fft_normalize" => FFT_NORMALIZE,
        "numeric_ttr" => NUMERIC_TTR,
        "numeric_t" => NUMERIC_T,
        "numeric_dt" => NUMERIC_DT,
        "fft_window" => string(FFT_WINDOW),
        "fft_detrend" => FFT_DETREND,
    )

    if isfile(NUMERIC_CACHE) && !generate
        cache = JLD2.load(NUMERIC_CACHE)
        cache_matches_specs = all(COLUMN_SPECS) do spec
            entries_key = string(spec.key, "_entries")
            haskey(cache, entries_key) &&
                length(cache[entries_key]) == length(specs[spec.key]) &&
                all(getindex.(cache[entries_key], :name) .== getindex.(specs[spec.key], :name))
        end
        cache_matches_selection = haskey(cache, "selection_signature") && cache["selection_signature"] == selection_signature
        cache_matches_specs && cache_matches_selection && return cache
        @info "Refreshing numerical FFT cache because selected oscillator set changed."
    end

    generate || error("Figure 4 FFT data are missing or incompatible. Run generate_data.jl --figures 4 --skip-continuation.")
    mkpath(dirname(NUMERIC_CACHE))
    cache = Dict{String, Any}()
    for spec in COLUMN_SPECS
        cache[string(spec.key)] = spectra_for_K(spec.numerical_K, specs[spec.key])
        cache[string(spec.key, "_entries")] = specs[spec.key]
    end
    cache["selectors"] = selectors
    cache["selection_signature"] = selection_signature
    cache["N"] = length(selectors.omega)
    jldopen(NUMERIC_CACHE, "w") do io
        for (key, value) in cache
            io[key] = value
        end
    end
    return cache
end

function add_fft_lines!(ax, series; linewidth = FFT_LINE_PT)
    for s in series
        isempty(s.f) && continue
        lines!(ax, s.f, s.spectrum; color = s.color, linewidth, label = s.label)
    end
end

function add_zoom_box!(ax, xlims, ylims)
    x1, x2 = xlims
    y1, y2 = ylims
    lines!(
        ax,
        [x1, x2, x2, x1, x1],
        [y1, y1, y2, y2, y1];
        color = :black,
        linewidth = FFT_ZOOM_BOX_LINE_PT,
    )
end

function add_peak_annotation!(ax, label, peak, labelpos; align = (:center, :center))
    leader_end = (
        labelpos[1] + 0.30 * (peak[1] - labelpos[1]),
        labelpos[2] + 0.30 * (peak[2] - labelpos[2]),
    )
    lines!(
        ax,
        [leader_end[1], peak[1]],
        [leader_end[2], peak[2]];
        color = :black,
        linewidth = PEAK_LEADER_LINE_PT,
    )
    text!(
        ax,
        labelpos[1],
        labelpos[2];
        text = label,
        fontsize = PEAK_LABEL_PT,
        align,
        color = :black,
    )
end

function add_top_right_inset!(location, series; xlims, ylims, xticks, linewidth = FFT_LINE_PT)
    inset = Axis(
        location;
        width = Relative(0.38),
        height = Relative(0.40),
        halign = 0.90,
        valign = 0.90,
        tellwidth = false,
        tellheight = false,
        backgroundcolor = (:white, 0.94),
    )
    fft_style_axis!(inset)
    inset.xticklabelsize = 0.72 * FFT_TICK_PT
    inset.yticklabelsize = 0.72 * FFT_TICK_PT
    inset.xlabelvisible = false
    inset.ylabelvisible = false
    inset.titlevisible = false
    inset.xticks = xticks
    inset.yticks = (Float64[], String[])
    inset.xticklabelpad = 1.0
    inset.ygridvisible = false
    add_fft_lines!(inset, series; linewidth)
    xlims!(inset, xlims...)
    ylims!(inset, ylims...)
    return inset
end

function numerical_series(cache, key::Symbol, signal::Symbol)
    data = cache[string(key)]
    entries = cache[string(key, "_entries")]
    series = []
    for entry in entries
        prefix = entry.name
        f = NUMERIC_ANGULAR_FREQUENCY_SCALE .* data["$(prefix)_$(signal)_freqs"]
        s = data["$(prefix)_$(signal)_spectrum"]
        push!(series, (f = f, spectrum = s, color = entry.color, label = entry.label))
    end
    return series
end

function write_numeric_selection_report(cache)
    selectors = cache["selectors"]
    open(NUMERIC_SELECTION_REPORT, "w") do io
        println(io, "Fig. 3 FFT numerical oscillator selection")
        println(io, "Generated by plot_fft_grid.jl")
        println(io)
        println(io, "Color convention")
        println(io, "  red  = low-amplitude cluster")
        println(io, "  blue = high-amplitude cluster")
        println(io)
        println(io, "The same omega_k values are used for all numerical rows.")
        println(io)
        println(io, @sprintf("red low-amplitude oscillator:  index = %d, omega_k = %.12g", selectors.low_idx, selectors.omega[selectors.low_idx]))
        println(io, @sprintf("blue high-amplitude oscillator: index = %d, omega_k = %.12g", selectors.high_idx, selectors.omega[selectors.high_idx]))
        println(io)
        println(io, "Low-amplitude oscillator checks")
        println(io, @sprintf("  K_H proxy = %.12g", selectors.k_h[selectors.low_idx]))
        println(io, @sprintf("  K_0 proxy = %.12g", selectors.k_0[selectors.low_idx]))
        println(io, "  Criteria: K_H >= 0.955, K_0 >= 0.88, r_min(K=0.955) > zero threshold")
        println(io)
        println(io, "r_min at plotted numerical rows")
        for (requested_k, row_idx) in zip(selectors.check_k_values, selectors.check_k_indices)
            actual_k = selectors.k_values[row_idx]
            low_r = selectors.min_radius[row_idx, selectors.low_idx]
            high_r = selectors.min_radius[row_idx, selectors.high_idx]
            println(io, @sprintf("  requested K = %.3f, grid K = %.12g: red r_min = %.12g, blue r_min = %.12g",
                                 requested_k, actual_k, low_r, high_r))
        end
        println(io)
        println(io, @sprintf("Low-cluster threshold used for selection: %.12g", LOW_CLUSTER_R_MIN_THRESHOLD))
        println(io, @sprintf("r_min zero threshold used for selection: %.12g", R_ZERO_THRESHOLD))
        println(io, @sprintf("Hopf baseline residual threshold used for selection: %.12g", HOPF_PROXY_THRESHOLD))
    end
    return NUMERIC_SELECTION_REPORT
end

function main()
    nrows = length(COLUMN_SPECS)
    cache = numerical_spectra()
    report_path = write_numeric_selection_report(cache)

    fig = Figure(
        size = FFT_FIGURE_SIZE_PT,
        backgroundcolor = :white,
        figure_padding = (0.45 * FFT_LATEX_BASE_PT, 0.55 * FFT_LATEX_BASE_PT, 0.72 * FFT_LATEX_BASE_PT, 0.70 * FFT_LATEX_BASE_PT),
        fonts = (; regular = PAPER_FONT, bold = PAPER_FONT, italic = PAPER_FONT, bold_italic = PAPER_FONT),
        fontsize = FFT_LATEX_BASE_PT,
    )

    column_labels = [
        L"\xi(t)",
        L"A(t)",
        L"\mathrm{Re}(W(t))",
        L"|W(t)|",
    ]
    axes = Matrix{Axis}(undef, nrows, 4)
    Label(fig[0, 1:2], "Electrochemical Silicon Etching"; fontsize = FFT_TITLE_PT, font = :bold)
    Label(fig[0, 3:4], "Heterogeneous Stuart-Landau"; fontsize = FFT_TITLE_PT, font = :bold)

    for (i, spec) in enumerate(COLUMN_SPECS)
        exp = read_experiment_fft(EXPERIMENT_XLSX, spec.key)
        exp_emsi = add_optional_red!(copy(exp.emsi), exp.emsi_red)
        exp_radii = add_optional_red!(copy(exp.radii), exp.radii_red)
        num_real = numerical_series(cache, spec.key, :real)
        num_radii = numerical_series(cache, spec.key, :radii)

        datasets = [exp_emsi, exp_radii, num_real, num_radii]
        for col in 1:4
            ax = Axis(
                fig[i, col];
                title = i == 1 ? column_labels[col] : "",
                xlabel = i == nrows ? (col <= 2 ? "frequency [Hz]" : "frequency [a.u.]") : "",
                ylabel = col == 1 ? spec.label : "",
            )
            fft_style_axis!(ax)
            line_width = col <= 2 ? FFT_LINE_PT : NUMERIC_FFT_LINE_PT
            add_fft_lines!(ax, datasets[col]; linewidth = line_width)
            if col <= 2
                xlims!(ax, EXPERIMENT_FFT_XLIMS...)
                ylims!(ax, EXPERIMENT_FFT_YLIMS...)
                ax.xticks = EXPERIMENT_FFT_XTICKS
                ax.yticks = EXPERIMENT_FFT_YTICKS
            else
                xlims!(ax, FFT_XLIMS...)
                ylims!(ax, NUMERIC_FFT_YLIMS...)
                ax.xticks = NUMERIC_FFT_XTICKS
                ax.yticks = NUMERIC_FFT_YTICKS
            end
            col in (2, 4) && hideydecorations!(ax; label = false, ticklabels = true, ticks = false, grid = false)
            axes[i, col] = ax
            i < nrows && hidexdecorations!(ax; grid = false)

            if spec.key === :c && col == 1
                add_peak_annotation!(ax, L"\Omega^{\mathrm{lock}}", (0.0251, 0.94), (0.0118, 0.78); align = (:center, :center))
                add_peak_annotation!(ax, L"\omega_{\mathrm{sb}}", (0.0234, 0.141), (0.066, 0.245); align = (:center, :center))
                add_zoom_box!(ax, XI_INSET_XLIMS, XI_INSET_YLIMS)
                add_top_right_inset!(fig[i, col], datasets[col]; xlims = XI_INSET_XLIMS, ylims = XI_INSET_YLIMS, xticks = XI_INSET_XTICKS)
            elseif spec.key === :c && col == 2
                add_zoom_box!(ax, A_INSET_XLIMS, A_INSET_YLIMS)
                inset = add_top_right_inset!(fig[i, col], datasets[col]; xlims = A_INSET_XLIMS, ylims = A_INSET_YLIMS, xticks = A_INSET_XTICKS)
                add_peak_annotation!(inset, L"\omega_{\mathrm{H}}", (0.0023, 0.145), (0.0064, 0.215); align = (:center, :center))
            elseif spec.key === :c && col == 3
                add_peak_annotation!(ax, L"\omega_{\mathrm{sb}}", (NUMERIC_ANGULAR_FREQUENCY_SCALE * 0.040, 0.072), (NUMERIC_ANGULAR_FREQUENCY_SCALE * 0.31, 0.112); align = (:center, :center))
                add_peak_annotation!(ax, L"\Omega^{\mathrm{lock}}", (NUMERIC_ANGULAR_FREQUENCY_SCALE * 0.983, 0.31), (NUMERIC_ANGULAR_FREQUENCY_SCALE * 0.78, 0.235); align = (:center, :center))
            elseif spec.key === :c && col == 4
                add_peak_annotation!(ax, L"\omega_{\mathrm{H}}", (NUMERIC_ANGULAR_FREQUENCY_SCALE * 1.021, 0.084), (NUMERIC_ANGULAR_FREQUENCY_SCALE * 0.76, 0.145); align = (:center, :center))
            end
        end
    end

    for col in 1:4
        linkyaxes!(axes[:, col]...)
    end
    for row in 1:nrows
        linkxaxes!(axes[row, 1:2]...)
        linkxaxes!(axes[row, 3:4]...)
    end

    rowsize!(fig.layout, 1, Relative(1 / nrows))
    rowsize!(fig.layout, 2, Relative(1 / nrows))
    rowsize!(fig.layout, 3, Relative(1 / nrows))
    rowgap!(fig.layout, 0.55 * FFT_LATEX_BASE_PT)
    rowgap!(fig.layout, 2, 1.05 * FFT_LATEX_BASE_PT)
    colgap!(fig.layout, 0.45 * FFT_LATEX_BASE_PT)

    png_path = FIG4_OUTPUT_STEM * ".png"
    pdf_path = FIG4_OUTPUT_STEM * ".pdf"
    fig4_png_path = FIG4_OUTPUT_STEM * ".png"
    fig4_pdf_path = FIG4_OUTPUT_STEM * ".pdf"
    save(png_path, fig; px_per_unit = 3)
    save(pdf_path, fig)
    @info "Saved FFT grid" png_path pdf_path fig4_png_path fig4_pdf_path cache = NUMERIC_CACHE report_path columns = [spec.key for spec in COLUMN_SPECS]
end

generate_data() = numerical_spectra(; generate = true)

end
