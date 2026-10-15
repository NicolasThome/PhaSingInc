module Figure2
using ..AmpDefInc
using XLSX
using CairoMakie
using Printf



const NINE_POINTS_FILE = joinpath(DATA_ROOT, "experimental", "9_points.xlsx")
const NUMERICAL_DIR = joinpath(DATA_ROOT, "numerical")
const MEAN_FREQUENCY_FILE = joinpath(NUMERICAL_DIR, "K_continuation_mean_frequency.csv")
const MIN_RADIUS_FILE = joinpath(NUMERICAL_DIR, "K_continuation_min_radius.csv")
const MAX_RADIUS_FILE = joinpath(NUMERICAL_DIR, "K_continuation_max_radius.csv")
const BASE_TYPOGRAPHY_TEST = lowercase(get(ENV, "BASE_TYPOGRAPHY_TEST", "false")) == "true"

# LaTeX-aware sizing. CairoMakie writes vector PDF dimensions in points,
# so set LATEX_FIGURE_WIDTH_PT to the width used in \includegraphics.
const LATEX_BASE_PT = 12.0
const LATEX_FIGURE_WIDTH_PT = 6.5 * 72.0
const FIGURE_SIZE_PT = (LATEX_FIGURE_WIDTH_PT, 0.70 * LATEX_FIGURE_WIDTH_PT)

const TICK_PT = 1.00 * LATEX_BASE_PT
const LABEL_PT = 1.00 * LATEX_BASE_PT
const TITLE_PT = 1.00 * LATEX_BASE_PT
const PANEL_LABEL_PT = 1.00 * LATEX_BASE_PT
const LINE_PT = 0.95
const EXPERIMENT_LINE_PT = 1.25
const MARKER_PT = 2.45
const EXPERIMENT_MARKER_PT = 1.85 * MARKER_PT
const COLORBAR_WIDTH_PT = 1.10 * LATEX_BASE_PT
const K_LIMITS = (0.85, 1.0)
const K_TICKS = ([0.85, 0.90, 0.95, 1.0], ["0.85", "0.90", "0.95", "1.00"])
const LOW_CLUSTER_R_MIN_THRESHOLD = 0.05
const HOPF_BASELINE_MIN_K = 0.975
const HOPF_PROXY_THRESHOLD = 0.015
const HIGHLIGHT_LINE_PT = 0.88
const HIGH_AMPLITUDE_LINE_PT = 1.33
const HIGH_AMPLITUDE_MARKER_PT = 1.45 * MARKER_PT
const ILLUMINATION_LIMITS = (0.77, 1.0)
const ILLUMINATION_MARGIN = 0.01
const ILLUMINATION_TICKS = ([0.77, 0.81, 0.85, 0.89, 0.94, 1.0], ["0.77", "0.81", "0.85", "0.89", "0.94", "1.00"])
const FREQUENCY_TRANSITION_THRESHOLD = -0.3
const DOTTED_EXPERIMENT_SHEETS = Set(["Sheet1", "Sheet3", "Sheet8"])

const PAPER_FONT = AmpDefInc.PAPER_FONT

function read_first_sheet_xlsx(path::AbstractString)
    book = XLSX.readxlsx(path)
    return book[1][:]
end

function read_all_sheets_xlsx(path::AbstractString)
    book = XLSX.readxlsx(path)
    return [name => book[name][:] for name in XLSX.sheetnames(book)]
end

function as_float(x)
    x isa Number && return Float64(x)
    x isa AbstractString && return parse(Float64, replace(x, "," => "."))
    return missing
end

function as_optional_float(x)
    x === missing && return missing
    x === nothing && return missing
    x isa Number && return Float64(x)
    if x isa AbstractString
        stripped = strip(x)
        isempty(stripped) && return missing
        stripped == "--" && return missing
        return parse(Float64, replace(stripped, "," => "."))
    end
    return missing
end

function read_experimental_data(path::AbstractString)
    sheet = read_first_sheet_xlsx(path)
    header = [lowercase(String(sheet[1, col])) for col in axes(sheet, 2)]

    illumination_col = findfirst(==("illumination"), header)
    omega_col = findfirst(==("omega"), header)
    min_col = findfirst(==("min"), header)
    any(isnothing, (illumination_col, omega_col, min_col)) &&
        error("Could not find illumination, omega, and min columns in $(basename(path)).")

    illumination = Float64[]
    omega_mhz = Union{Missing, Float64}[]
    emsi_min = Union{Missing, Float64}[]

    for row in 3:size(sheet, 1)
        illumination_value = as_float(sheet[row, illumination_col])
        illumination_value === missing && continue
        isapprox(illumination_value, 3.13; atol = 1e-8) && continue

        push!(illumination, illumination_value)
        push!(omega_mhz, as_float(sheet[row, omega_col]))
        push!(emsi_min, as_float(sheet[row, min_col]))
    end

    return illumination, omega_mhz, emsi_min
end

function read_nine_point_oscillators(path::AbstractString)
    sheets = read_all_sheets_xlsx(path)
    oscillators = NamedTuple[]

    for (sheet_idx, (sheet_name, sheet)) in enumerate(sheets)
        header = [cell === missing ? "" : lowercase(String(cell)) for cell in sheet[1, :]]
        illumination_col = findfirst(==("illumination"), header)
        omega_col = findfirst(==("omega"), header)
        min_col = findfirst(==("min"), header)
        any(isnothing, (illumination_col, omega_col, min_col)) &&
            error("Could not find illumination, omega, and min columns in sheet '$sheet_name'.")

        illumination = Float64[]
        omega_mhz = Union{Missing, Float64}[]
        amin = Union{Missing, Float64}[]

        for row in 2:size(sheet, 1)
            illumination_value = as_optional_float(sheet[row, illumination_col])
            illumination_value === missing && continue
            first(ILLUMINATION_LIMITS) <= illumination_value <= last(ILLUMINATION_LIMITS) || continue

            push!(illumination, illumination_value)
            push!(omega_mhz, as_optional_float(sheet[row, omega_col]))
            push!(amin, as_optional_float(sheet[row, min_col]))
        end

        order = sortperm(illumination; rev = true)
        illumination = illumination[order]
        omega_mhz = omega_mhz[order]
        amin = amin[order]

        transition_idx = typemax(Int)
        transition_illumination = -Inf
        transition_jump = NaN
        for idx in 1:(length(omega_mhz) - 1)
            ismissing(omega_mhz[idx]) && continue
            ismissing(omega_mhz[idx + 1]) && continue
            jump = omega_mhz[idx + 1] - omega_mhz[idx]
            if jump <= FREQUENCY_TRANSITION_THRESHOLD
                transition_idx = idx + 1
                transition_illumination = illumination[idx + 1]
                transition_jump = jump
                break
            end
        end

        push!(
            oscillators,
            (;
                sheet_idx,
                sheet_name,
                illumination,
                omega_mhz,
                amin,
                transition_idx,
                transition_illumination,
                transition_jump,
            ),
        )
    end

    return oscillators
end

function transition_order_colors(oscillators; colormap = :jet)
    transition_levels = sort(
        unique([osc.transition_illumination for osc in oscillators if isfinite(osc.transition_illumination)]);
        rev = true,
    )
    isempty(transition_levels) && error("No transitions detected with threshold $(FREQUENCY_TRANSITION_THRESHOLD).")

    colors = Dict{String, Any}()
    group_count = length(transition_levels)
    grad = CairoMakie.cgrad(colormap)

    for (group_idx, transition_level) in enumerate(transition_levels)
        group = sort(
            [osc for osc in oscillators if osc.transition_illumination == transition_level];
            by = osc -> osc.sheet_idx,
        )
        center = group_count == 1 ? 0.5 : (group_idx - 1) / (group_count - 1)
        span = min(0.20, group_count == 1 ? 0.20 : 0.70 / max(1, group_count - 1))
        lo = group_idx == 1 ? 0.0 : group_idx == group_count ? max(0.0, 1.0 - span) : center - 0.5 * span
        hi = group_idx == 1 ? min(1.0, span) : group_idx == group_count ? 1.0 : center + 0.5 * span
        fractions = length(group) == 1 ? [center] : collect(range(lo, hi; length = length(group)))
        for (osc, fraction) in zip(group, fractions)
            colors[osc.sheet_name] = grad[clamp(fraction, 0.0, 1.0)]
        end
    end

    for osc in oscillators
        haskey(colors, osc.sheet_name) && continue
        colors[osc.sheet_name] = grad[1.0]
    end
    return colors
end

function read_numeric_csv(path::AbstractString; first_data_col::Int)
    lines = collect(eachline(path))
    header = split(lines[1], ",")
    n_rows = length(lines) - 1
    n_cols = length(header)

    values = Matrix{Float64}(undef, n_rows, n_cols)
    for (row_idx, line) in enumerate(@view lines[2:end])
        fields = split(line, ",")
        length(fields) == n_cols || error("Unexpected number of columns in $(basename(path)) at row $(row_idx + 1).")
        @inbounds for col_idx in 1:n_cols
            values[row_idx, col_idx] = parse(Float64, fields[col_idx])
        end
    end

    omega = Float64[]
    for name in header[first_data_col:end]
        omega_match = match(r"omega_(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$", name)
        omega_match === nothing && error("Could not parse intrinsic frequency from column '$name'.")
        push!(omega, parse(Float64, omega_match.captures[1]))
    end

    return values[:, 1], values[:, first_data_col:end], omega
end

function read_mean_frequency_csv(path::AbstractString)
    lines = collect(eachline(path))
    header = split(lines[1], ",")
    header[2] == "system_mean_frequency" ||
        error("Expected system_mean_frequency in the second column of $(basename(path)).")

    k_values, mean_frequency, omega = read_numeric_csv(path; first_data_col = 3)
    system_mean_frequency = [parse(Float64, split(line, ",")[2]) for line in @view lines[2:end]]
    return k_values, mean_frequency, omega, system_mean_frequency
end

function k_continuation_color(value::Real, colormap, colorrange)
    lo, hi = Float64.(colorrange)
    fraction = hi == lo ? 0.5 : clamp((Float64(value) - lo) / (hi - lo), 0.0, 1.0)
    return CairoMakie.cgrad(colormap)[fraction]
end

function finite_points(x::Vector{Float64}, y::Vector{Union{Missing, Float64}})
    keep = .!ismissing.(y)
    return x[keep], Float64.(y[keep])
end

function finite_mean(values)
    total = 0.0
    count = 0
    for value in values
        isfinite(value) || continue
        total += Float64(value)
        count += 1
    end
    return count == 0 ? NaN : total / count
end

function row_finite_mean(values::AbstractMatrix)
    means = Vector{Float64}(undef, size(values, 1))
    for i in axes(values, 1)
        means[i] = finite_mean(@view values[i, :])
    end
    return means
end

function add_colored_curves!(ax, k_values, data, omega; colorrange, colormap, linewidth = 0.16 * LINE_PT)
    n_k, n_osc = size(data)
    points = Vector{Point2f}(undef, 2 * (n_k - 1) * n_osc)
    colors = Vector{Float32}(undef, length(points))

    idx = 1
    @inbounds for osc_idx in 1:n_osc
        curve_color = Float32(omega[osc_idx])
        for k_idx in 1:(n_k - 1)
            points[idx] = Point2f(k_values[k_idx], data[k_idx, osc_idx])
            colors[idx] = curve_color
            idx += 1
            points[idx] = Point2f(k_values[k_idx + 1], data[k_idx + 1, osc_idx])
            colors[idx] = curve_color
            idx += 1
        end
    end

    return linesegments!(
        ax,
        points;
        color = colors,
        colormap = colormap,
        colorrange = colorrange,
        linewidth = linewidth,
        transparency = true,
        alpha = 0.34,
    )
end

function style_axis!(ax)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.10)
    ax.ygridcolor = (:black, 0.10)
    ax.xticklabelsize = TICK_PT
    ax.yticklabelsize = TICK_PT
    ax.xlabelsize = LABEL_PT
    ax.ylabelsize = LABEL_PT
    ax.titlesize = TITLE_PT
    ax.titlefont = :bold
    ax.spinewidth = 0.65 * LINE_PT
    ax.xtickwidth = 0.65 * LINE_PT
    ax.ytickwidth = 0.65 * LINE_PT
    ax.xticksize = 3.0
    ax.yticksize = 3.0
    ax.xlabelpadding = 2.0
    ax.ylabelpadding = 2.0
    ax.xticklabelpad = 2.0
    ax.yticklabelpad = 2.0
    return ax
end

function add_panel_label!(ax, label; x = 0.018, y = 0.965, align = (:left, :top))
    text!(
        ax,
        x,
        y;
        text = label,
        space = :relative,
        align,
        fontsize = PANEL_LABEL_PT,
        font = :bold,
        color = :black,
    )
end

function in_k_window(k_values)
    lo, hi = K_LIMITS
    return findall(k -> lo <= k <= hi, k_values)
end

function baseline_residual(k_values, min_radius)
    baseline_rows = findall(>=(HOPF_BASELINE_MIN_K), k_values)
    length(baseline_rows) >= 2 || error("Need at least two high-K points for the Hopf fallback proxy.")

    residual = similar(min_radius)
    X = hcat(ones(length(baseline_rows)), k_values[baseline_rows])
    for j in axes(min_radius, 2)
        beta = X \ min_radius[baseline_rows, j]
        baseline = beta[1] .+ beta[2] .* k_values
        residual[:, j] .= max.(0.0, baseline .- min_radius[:, j])
    end
    return residual
end

function optional_max_radius(path, shape)
    isfile(path) || return nothing
    k_max, max_radius, omega_max = read_numeric_csv(path; first_data_col = 2)
    size(max_radius) == shape || error("Max-radius matrix shape does not match min-radius data.")
    return k_max, max_radius, omega_max
end

function hopf_amplitude_proxy(k_values, min_radius, max_radius)
    if max_radius !== nothing
        return 0.5 .* (max_radius .- min_radius), "A_r"
    end
    return baseline_residual(k_values, min_radius), "baseline residual"
end

function first_along_decreasing_k(k_values, values, predicate)
    order = sortperm(k_values; rev = true)
    for idx in order
        predicate(values[idx]) && return k_values[idx]
    end
    return NaN
end

function closest_to_zero_k(k_values, values)
    # Resolve equal minima in favor of the first point along decreasing K.
    order = sortperm(k_values; rev = true)
    finite_order = filter(idx -> isfinite(values[idx]), order)
    isempty(finite_order) && return NaN
    closest = finite_order[argmin(abs.(values[finite_order]))]
    return k_values[closest]
end

function transition_ticks(k_hb, k_ad)
    return K_TICKS
end

function add_transition_label!(
    ax,
    x,
    label;
    side::Symbol,
    y = 0.36,
    dx = 0.0018,
)
    isfinite(x) || return nothing
    x_label = side === :left ? x - dx : x + dx
    align = side === :left ? (:right, :bottom) : (:left, :bottom)
    text!(
        ax,
        x_label,
        y;
        text = label,
        align,
        fontsize = 0.86 * LABEL_PT,
        color = :black,
    )
    return nothing
end

function selected_low_cluster_oscillators(k_values, min_radius, max_radius, omega)
    window = in_k_window(k_values)
    k_window = k_values[window]
    min_window = min_radius[window, :]
    max_window = max_radius === nothing ? nothing : max_radius[window, :]
    hopf_proxy, hopf_source = hopf_amplitude_proxy(k_window, min_window, max_window)
    min_by_osc = vec(minimum(min_window; dims = 1))

    low_cluster = findall(min_by_osc .<= LOW_CLUSTER_R_MIN_THRESHOLD)
    isempty(low_cluster) && error("No low-amplitude cluster oscillators found with LOW_CLUSTER_R_MIN_THRESHOLD = $LOW_CLUSTER_R_MIN_THRESHOLD.")

    k_h = fill(NaN, size(min_window, 2))
    k_0 = fill(NaN, size(min_window, 2))
    for idx in low_cluster
        k_h[idx] = first_along_decreasing_k(k_window, hopf_proxy[:, idx], value -> value > HOPF_PROXY_THRESHOLD)
        k_0[idx] = closest_to_zero_k(k_window, min_window[:, idx])
    end

    finite_hopf = [idx for idx in low_cluster if isfinite(k_h[idx])]
    finite_zero = [idx for idx in low_cluster if isfinite(k_0[idx])]
    isempty(finite_hopf) && error("No Hopf-onset proxy crossings found in the low-amplitude cluster.")
    isempty(finite_zero) && error("No finite r_min minima found in the low-amplitude cluster.")

    first_hopf = finite_hopf[argmax(k_h[finite_hopf])]
    last_zero = finite_zero[argmin(k_0[finite_zero])]
    zero_order = sort(finite_zero; by = idx -> k_0[idx], rev = true)
    middle_zero = finite_zero[argmin(abs.(omega[finite_zero]))]
    high_amplitude = argmax(min_by_osc)

    selected = unique([first_hopf, last_zero, middle_zero])
    if length(selected) < 3
        for idx in zero_order
            idx in selected && continue
            push!(selected, idx)
            length(selected) == 3 && break
        end
    end
    if length(selected) < 3
        for idx in sort(finite_hopf; by = i -> k_h[i], rev = true)
            idx in selected && continue
            push!(selected, idx)
            length(selected) == 3 && break
        end
    end

    return (;
        selected,
        first_hopf,
        last_zero,
        middle_zero,
        high_amplitude,
        k_h_first = k_h[first_hopf],
        k_0_first_hopf = k_0[first_hopf],
        k_0_last_zero = k_0[last_zero],
        k_0_middle_zero = k_0[middle_zero],
        low_cluster,
        hopf_source,
        k_h,
        k_0,
    )
end

function plot_highlighted_oscillators!(ax, k_values, data, selected, omega; colorrange, colormap)
    for idx in selected
        color = k_continuation_color(omega[idx], colormap, colorrange)
        lines!(ax, k_values, data[:, idx]; color, linewidth = HIGHLIGHT_LINE_PT)
    end
end

function plot_high_amplitude_mean_frequency!(ax, k_values, mean_frequency, idx, omega; colorrange, colormap)
    color = k_continuation_color(omega[idx], colormap, colorrange)
    lines!(ax, k_values, mean_frequency[:, idx]; color, linewidth = HIGH_AMPLITUDE_LINE_PT)
    scatter!(
        ax,
        k_values,
        mean_frequency[:, idx];
        color,
        markersize = HIGH_AMPLITUDE_MARKER_PT,
        strokewidth = 0.25 * LINE_PT,
        strokecolor = :white,
    )
end

function plot_observed_oscillator!(ax, x::Vector{Float64}, y::Vector{Union{Missing, Float64}}, color; linestyle = :solid)
    segment_x = Float64[]
    segment_y = Float64[]

    function flush_segment!()
        if length(segment_x) >= 2
            lines!(ax, segment_x, segment_y; color, linewidth = EXPERIMENT_LINE_PT, linestyle)
        end
        empty!(segment_x)
        empty!(segment_y)
        return nothing
    end

    for idx in eachindex(x, y)
        if ismissing(y[idx])
            flush_segment!()
        else
            push!(segment_x, x[idx])
            push!(segment_y, y[idx])
        end
    end
    flush_segment!()

    finite_x, finite_y = finite_points(x, y)
    scatter!(
        ax,
        finite_x,
        finite_y;
        color,
        markersize = EXPERIMENT_MARKER_PT,
        strokewidth = 0.25 * LINE_PT,
        strokecolor = :white,
    )
    return nothing
end

function main()
    observed_oscillators = read_nine_point_oscillators(NINE_POINTS_FILE)
    observed_colors = transition_order_colors(observed_oscillators; colormap = :jet)
    k_mean, mean_frequency, omega_mean, system_mean_frequency = read_mean_frequency_csv(MEAN_FREQUENCY_FILE)
    k_minr, min_radius, omega_minr = read_numeric_csv(MIN_RADIUS_FILE; first_data_col = 2)
    max_radius_data = optional_max_radius(MAX_RADIUS_FILE, size(min_radius))

    k_mean == k_minr || error("The K grids in the numerical CSV files do not match.")
    omega_mean == omega_minr || error("The oscillator frequency columns in the numerical CSV files do not match.")
    max_radius = if max_radius_data === nothing
        nothing
    else
        k_maxr, max_radius_values, omega_maxr = max_radius_data
        k_mean == k_maxr || error("The K grid in the max-radius CSV does not match.")
        omega_mean == omega_maxr || error("The oscillator frequency columns in the max-radius CSV do not match.")
        max_radius_values
    end

    colorrange = extrema(omega_mean)
    colormap = :jet
    blue = RGBf(0.0, 0.32, 0.62)
    red = RGBf(0.84, 0.22, 0.13)
    highlight_info = selected_low_cluster_oscillators(k_mean, min_radius, max_radius, omega_mean)

    fig = Figure(
        size = FIGURE_SIZE_PT,
        backgroundcolor = :white,
        figure_padding = (0.45 * LATEX_BASE_PT, 0.55 * LATEX_BASE_PT, 0.45 * LATEX_BASE_PT, 0.45 * LATEX_BASE_PT),
        fonts = (; regular = PAPER_FONT, bold = PAPER_FONT, italic = PAPER_FONT, bold_italic = PAPER_FONT),
        fontsize = LATEX_BASE_PT,
    )

    ax_mean = Axis(
        fig[1, 1];
        ylabel = L"\Omega_k",
    )
    style_axis!(ax_mean)
    ax_mean.yticks = ([0.0, 3.0, 6.0], ["0", "3", "6"])
    add_colored_curves!(ax_mean, k_mean, mean_frequency, omega_mean; colorrange, colormap)
    plot_highlighted_oscillators!(ax_mean, k_mean, mean_frequency, highlight_info.selected, omega_mean; colorrange, colormap)
    plot_high_amplitude_mean_frequency!(ax_mean, k_mean, mean_frequency, highlight_info.high_amplitude, omega_mean; colorrange, colormap)

    ax_minr = Axis(
        fig[2, 1];
        xlabel = L"K",
        ylabel = L"r_{\min,k}",
    )
    style_axis!(ax_minr)
    ax_minr.yticks = ([0.0, 0.5, 1.0], ["0", "0.5", "1"])
    add_colored_curves!(ax_minr, k_mean, min_radius, omega_mean; colorrange, colormap)
    plot_highlighted_oscillators!(ax_minr, k_mean, min_radius, highlight_info.selected, omega_mean; colorrange, colormap)

    linkxaxes!(ax_mean, ax_minr)
    hidexdecorations!(ax_mean; grid = false)
    xlims!(ax_mean, K_LIMITS...)
    xlims!(ax_minr, K_LIMITS...)
    k_transition_ticks = transition_ticks(highlight_info.k_h_first, highlight_info.k_0_first_hopf)
    ax_mean.xticks = k_transition_ticks
    ax_minr.xticks = k_transition_ticks
    for ax in (ax_mean, ax_minr)
        vlines!(ax, [highlight_info.k_h_first]; color = :black, linestyle = :dot, linewidth = 1.15 * LINE_PT)
        vlines!(ax, [highlight_info.k_0_first_hopf]; color = :black, linestyle = :dashdot, linewidth = 1.15 * LINE_PT)
        vlines!(ax, [highlight_info.k_0_last_zero]; color = :black, linestyle = :dashdot, linewidth = 1.15 * LINE_PT)
    end
    add_transition_label!(ax_minr, highlight_info.k_0_first_hopf, L"K_{\mathrm{ps},\,1}"; side = :left)
    add_transition_label!(ax_minr, highlight_info.k_0_last_zero, L"K_{\mathrm{ps},\,N^\prime}"; side = :left)
    add_transition_label!(ax_minr, highlight_info.k_h_first, L"K_{\mathrm{HB}}"; side = :right)

    Colorbar(
        fig[1:2, 2];
        colormap = colormap,
        limits = colorrange,
        label = L"\omega_k",
        ticks = [-1.0, 0.0, 1.0],
        width = COLORBAR_WIDTH_PT,
        valign = :center,
        labelsize = LABEL_PT,
        ticklabelsize = TICK_PT,
    )

    ax_obs_omega = Axis(
        fig[1, 3];
        ylabel = L"\Omega_k\ [\mathrm{mHz}]",
    )
    style_axis!(ax_obs_omega)
    ax_obs_omega.yticks = ([23.0, 24.0, 25.0, 26.0], ["23", "24", "25", "26"])
    ax_obs_omega.xticks = ILLUMINATION_TICKS
    hidexdecorations!(ax_obs_omega; grid = false)

    ax_obs_min = Axis(
        fig[2, 3];
        xlabel = L"I\ [\mathrm{mW}/\mathrm{cm}^2]",
        ylabel = L"A_{\min,k}\ [\mathrm{a.u.}]",
    )
    style_axis!(ax_obs_min)
    ax_obs_min.xticks = ILLUMINATION_TICKS
    ax_obs_min.yticks = ([0.0, 0.4, 0.8], ["0", "0.4", "0.8"])

    for osc in observed_oscillators
        color = observed_colors[osc.sheet_name]
        linestyle = osc.sheet_name in DOTTED_EXPERIMENT_SHEETS ? :dot : :solid
        plot_observed_oscillator!(ax_obs_omega, osc.illumination, osc.omega_mhz, color; linestyle)
        plot_observed_oscillator!(ax_obs_min, osc.illumination, osc.amin, color; linestyle)
    end
    linkxaxes!(ax_obs_omega, ax_obs_min)
    xlims!(ax_obs_omega, first(ILLUMINATION_LIMITS) - ILLUMINATION_MARGIN, last(ILLUMINATION_LIMITS) + ILLUMINATION_MARGIN)
    xlims!(ax_obs_min, first(ILLUMINATION_LIMITS) - ILLUMINATION_MARGIN, last(ILLUMINATION_LIMITS) + ILLUMINATION_MARGIN)
    ylims!(ax_obs_omega, 22.9, 26.1)
    ylims!(ax_obs_min, -0.03, 0.83)

    add_panel_label!(ax_mean, "(a)"; x = 0.955, y = 0.855, align = (:right, :top))
    add_panel_label!(ax_minr, "(b)"; x = 0.955, y = 0.780, align = (:right, :top))
    add_panel_label!(ax_obs_omega, "(c)"; x = 0.965, y = 0.965, align = (:right, :top))
    add_panel_label!(ax_obs_min, "(d)"; x = 0.965, y = 0.965, align = (:right, :top))

    Label(fig[0, 1:2], "Heterogeneous Stuart-Landau"; font = :bold)
    Label(fig[0, 3], "Electrochemical Silicon Etching"; font = :bold)

    rowsize!(fig.layout, 1, Relative(0.50))
    rowsize!(fig.layout, 2, Relative(0.50))
    colsize!(fig.layout, 1, Relative(0.47))
    colsize!(fig.layout, 3, Relative(0.47))
    colgap!(fig.layout, 0.70 * LATEX_BASE_PT)
    rowgap!(fig.layout, 0.75 * LATEX_BASE_PT)

    save(joinpath(OUTPUT_ROOT, "Fig_2.png"), fig; px_per_unit = 3)
    save(joinpath(OUTPUT_ROOT, "Fig_2.pdf"), fig)
    @info "Saved Figure 2" K_HB=highlight_info.k_h_first K_ps_1=highlight_info.k_0_first_hopf K_ps_last=highlight_info.k_0_last_zero
end

end
