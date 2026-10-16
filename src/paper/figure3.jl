module Figure3
using ..PhaSingInc
using CairoMakie
using JLD2
using LaTeXStrings
using OrdinaryDiffEqTsit5
using Printf
using SciMLBase



const NUMERICAL_DIR = joinpath(DATA_ROOT, "numerical")
const RESTART_DIR = joinpath(NUMERICAL_DIR, "K_continuation_restart")
const MIN_RADIUS_FILE = joinpath(NUMERICAL_DIR, "K_continuation_min_radius.csv")

const K_VALUES = [0.970, 0.955, 0.948, 0.88]
const FIXED_POINT_K_VALUES = [0.970]
const TRAJECTORY_T = 30.0
const SECOND_PANEL_TTR = 100.0
const ROTATING_FRAME_TTR = 100.0
const ROTATING_FRAME_ANALYSIS_T = 20.0
const ROTATING_FRAME_DT = 0.001
const ROTATING_FRAME_PERIODS = 3
const TRAJECTORY_DT = 0.01
const LOW_DARK_BLUE_INDEX = 199
const HIGH_CLUSTER_THRESHOLD = 0.7

const LATEX_BASE_PT = 12.0
const PHASE_FIGURE_SIZE_PT = (6.5 * 72.0, 4.10 * 72.0)
const TICK_PT = 1.00 * LATEX_BASE_PT
const LABEL_PT = 1.00 * LATEX_BASE_PT
const TITLE_PT = 1.00 * LATEX_BASE_PT
const MARKER_PT = 6.5

const DARK_BLUE = RGBf(0.0, 0.0, 0.46)

const PAPER_FONT = PhaSingInc.PAPER_FONT

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

function read_numeric_csv(path::AbstractString; first_data_col::Int)
    lines = readlines(path)
    header = split(lines[1], ",")
    omega = Float64[]
    for label in header[first_data_col:end]
        m = match(r"omega_([-+0-9.eE]+)$", label)
        m === nothing && error("Could not parse omega from CSV header label: $label")
        push!(omega, parse(Float64, m.captures[1]))
    end

    k_values = Float64[]
    values = Matrix{Float64}(undef, length(lines) - 1, length(omega))
    for (i, line) in enumerate(lines[2:end])
        parts = split(line, ",")
        push!(k_values, parse(Float64, parts[1]))
        for j in eachindex(omega)
            values[i, j] = parse(Float64, parts[first_data_col + j - 1])
        end
    end
    return k_values, values, omega
end

function high_green_index()
    _, min_radius, omega = read_numeric_csv(MIN_RADIUS_FILE; first_data_col = 2)
    high_candidates = findall(vec(minimum(min_radius; dims = 1)) .>= HIGH_CLUSTER_THRESHOLD)
    isempty(high_candidates) && error("No high-amplitude oscillators found with threshold $HIGH_CLUSTER_THRESHOLD.")
    return high_candidates[argmin(abs.(omega[high_candidates]))], omega
end

function load_restart(path::AbstractString)
    data = JLD2.load(path)
    return (;
        K = Float64(data["K"]),
        C1 = Float64(data["C1"]),
        C2 = Float64(data["C2"]),
        omega = Float64.(data["omega"]),
        u0 = Float64.(data["final_state"]),
        reltol = haskey(data, "reltol") ? Float64(data["reltol"]) : 1.0e-9,
        abstol = haskey(data, "abstol") ? Float64(data["abstol"]) : 1.0e-9,
    )
end

function sample_times(T::Real; dt::Real = TRAJECTORY_DT)
    times = collect(0.0:Float64(dt):Float64(T))
    if isempty(times) || times[end] < Float64(T) - 1.0e-12
        push!(times, Float64(T))
    end
    return times
end

function simulate_complex_trajectories(
    K::Real,
    indices;
    trajectory_T::Real = TRAJECTORY_T,
    trajectory_dt::Real = TRAJECTORY_DT,
    trajectory_Ttr::Real = 0.0,
)
    restart = load_restart(restart_file_for_K(K))
    save_idxs = Int[]
    for idx in indices
        push!(save_idxs, 2idx - 1)
        push!(save_idxs, 2idx)
    end
    save_idxs = sort(unique(save_idxs))

    p = (restart.K, restart.C1, restart.C2, restart.omega)
    final_time = Float64(trajectory_Ttr) + Float64(trajectory_T)
    record_times = Float64(trajectory_Ttr) .+ sample_times(trajectory_T; dt = trajectory_dt)
    prob = SciMLBase.ODEProblem(PhaSingInc.Model.sl_harmo_hetero!, restart.u0, (0.0, final_time), p)
    sol = SciMLBase.solve(
        prob,
        OrdinaryDiffEqTsit5.Tsit5();
        saveat = record_times,
        save_idxs,
        reltol = restart.reltol,
        abstol = restart.abstol,
    )

    mat = Array(sol)'
    colmap = Dict(save_idxs[i] => i for i in eachindex(save_idxs))
    trajectories = Dict{Int, Vector{ComplexF64}}()
    for idx in indices
        x = mat[:, colmap[2idx - 1]]
        y = mat[:, colmap[2idx]]
        trajectories[idx] = ComplexF64.(x, y)
    end
    return restart.K, trajectories
end

panel_transient(K::Real) = isapprox(Float64(K), 0.948; atol = 1.0e-10, rtol = 0.0) ? SECOND_PANEL_TTR : 0.0

function phase_slope(times, phase)
    t_mean = sum(times) / length(times)
    phase_mean = sum(phase) / length(phase)
    numerator = sum((times .- t_mean) .* (phase .- phase_mean))
    denominator = sum(abs2, times .- t_mean)
    return numerator / denominator
end

function unwrap_continuous!(phase)
    for i in 2:length(phase)
        delta = phase[i] - phase[i - 1]
        delta -= 2pi * round(delta / (2pi))
        phase[i] = phase[i - 1] + delta
    end
    return phase
end

function modulation_segment(times, z; n_periods::Int = ROTATING_FRAME_PERIODS)
    radius = abs.(z)
    minima = [i for i in 2:(length(radius) - 1) if radius[i] <= radius[i - 1] && radius[i] < radius[i + 1]]
    valid_starts = 1:(length(minima) - n_periods)
    isempty(valid_starts) && error("Not enough radial minima to select $n_periods modulation periods.")
    center_offset = n_periods ÷ 2
    start_position = valid_starts[argmin([radius[minima[start + center_offset]] for start in valid_starts])]
    first_idx = minima[start_position]
    last_idx = minima[start_position + n_periods]
    closest_radius = minimum(radius[minima[start_position:(start_position + n_periods)]])
    return first_idx:last_idx, closest_radius
end

function fixed_point_diagnostic(K_used, Omega, rotating)
    fixed_point = sum(rotating) / length(rotating)
    fixed_point_aligned = abs(fixed_point) + 0im
    segment_times = collect(range(0.0, ROTATING_FRAME_PERIODS; length = ROTATING_FRAME_PERIODS + 1))
    return (;
        K = K_used,
        Omega,
        times = segment_times,
        rotating = fill(fixed_point_aligned, length(segment_times)),
        psi = zeros(length(segment_times)),
        q = 0.0,
        closest_radius = abs(fixed_point_aligned),
        modulation_period = 1.0,
        hopf_frequency = 2pi,
        fixed_point = true,
    )
end

function rotating_frame_diagnostics()
    high_idx, _ = high_green_index()
    indices = [LOW_DARK_BLUE_INDEX, high_idx]
    times = sample_times(ROTATING_FRAME_ANALYSIS_T; dt = ROTATING_FRAME_DT)
    diagnostics = Dict{Float64, NamedTuple}()

    for K in K_VALUES
        K_used, trajectories = simulate_complex_trajectories(
            K,
            indices;
            trajectory_T = ROTATING_FRAME_ANALYSIS_T,
            trajectory_dt = ROTATING_FRAME_DT,
            trajectory_Ttr = ROTATING_FRAME_TTR,
        )
        locked_phase = unwrap_continuous!(angle.(trajectories[high_idx]))
        Omega = phase_slope(times, locked_phase)
        blue = trajectories[LOW_DARK_BLUE_INDEX]
        rotating = blue .* cis.(-Omega .* times)
        if any(isapprox(K_used, K_fixed; atol = 5.0e-4, rtol = 0.0) for K_fixed in FIXED_POINT_K_VALUES)
            diagnostics[K] = fixed_point_diagnostic(K_used, Omega, rotating)
            @info "Rotating-frame diagnostic" K=K_used Omega q=diagnostics[K].q closest_radius=diagnostics[K].closest_radius modulation_period=diagnostics[K].modulation_period hopf_frequency=diagnostics[K].hopf_frequency fixed_point=true
            continue
        end
        segment, closest_radius = try
            modulation_segment(times, blue)
        catch err
            @info "Treating rotating-frame trajectory as a fixed point" K=K_used reason=sprint(showerror, err)
            diagnostics[K] = fixed_point_diagnostic(K_used, Omega, rotating)
            @info "Rotating-frame diagnostic" K=K_used Omega q=diagnostics[K].q closest_radius=diagnostics[K].closest_radius modulation_period=diagnostics[K].modulation_period hopf_frequency=diagnostics[K].hopf_frequency
            continue
        end
        segment_times = times[segment] .- times[first(segment)]
        segment_rotating = rotating[segment]
        psi = unwrap_continuous!(angle.(segment_rotating))
        psi .-= psi[1]
        phase_aligned_rotating = abs.(segment_rotating) .* cis.(psi)
        q = (psi[end] - psi[1]) / (2pi * ROTATING_FRAME_PERIODS)
        diagnostics[K] = (;
            K = K_used,
            Omega,
            times = segment_times,
            rotating = phase_aligned_rotating,
            psi,
            q,
            closest_radius,
            modulation_period = segment_times[end] / ROTATING_FRAME_PERIODS,
            hopf_frequency = 2pi * ROTATING_FRAME_PERIODS / segment_times[end],
            fixed_point = false,
        )
        @info "Rotating-frame diagnostic" K=K_used Omega q closest_radius modulation_period=diagnostics[K].modulation_period hopf_frequency=diagnostics[K].hopf_frequency
    end
    return diagnostics
end

function add_direction_markers!(ax, trajectory)
    n = length(trajectory)
    for fraction in (0.22, 0.50, 0.78)
        idx = clamp(round(Int, fraction * n), 3, n - 2)
        tangent = trajectory[idx + 2] - trajectory[idx - 2]
        abs(tangent) > 0 || continue
        scatter!(
            ax,
            [real(trajectory[idx])],
            [imag(trajectory[idx])];
            marker = :utriangle,
            rotation = angle(tangent) - pi / 2,
            color = DARK_BLUE,
            markersize = 6.0,
            strokecolor = :white,
            strokewidth = 0.30,
        )
    end
end

function save_low_only_with_phase_figure(path_stem::AbstractString, diagnostics)
    fig = Figure(
        size = PHASE_FIGURE_SIZE_PT,
        backgroundcolor = :white,
        figure_padding = (0.44 * LATEX_BASE_PT, 0.32 * LATEX_BASE_PT, 0.90 * LATEX_BASE_PT, 0.28 * LATEX_BASE_PT),
        fonts = (; regular = PAPER_FONT, bold = PAPER_FONT, italic = PAPER_FONT, bold_italic = PAPER_FONT),
        fontsize = LATEX_BASE_PT,
    )
    complex_axes = Axis[]
    phase_axes = Axis[]
    top_labels = ["(a)", "(b)", "(c)", "(d)"]
    bottom_labels = ["(e)", "(f)", "(g)", "(h)"]
    extent = maximum(
        max(maximum(abs.(real.(diagnostics[K].rotating))), maximum(abs.(imag.(diagnostics[K].rotating))))
        for K in K_VALUES
    )
    complex_limit = 1.08 * extent
    complex_ticks = [-0.3, 0.0, 0.3]

    for (col, K) in pairs(K_VALUES)
        diagnostic = diagnostics[K]
        complex_ax = Axis(
            fig[1, col];
            xlabel = L"\mathrm{Re}\,\tilde{W}",
            ylabel = col == 1 ? L"\mathrm{Im}\,\tilde{W}" : "",
            limits = (-complex_limit, complex_limit, -complex_limit, complex_limit),
        )
        style_axis!(complex_ax)
        complex_ax.xticks = complex_ticks
        complex_ax.yticks = complex_ticks
        text!(complex_ax, 0.035, 0.965; text = top_labels[col], space = :relative, align = (:left, :top), fontsize = LABEL_PT, font = :bold)
        text!(complex_ax, 0.965, 0.965; text = @sprintf("K = %.3f", K), space = :relative, align = (:right, :top), fontsize = 0.70 * LABEL_PT)
        push!(complex_axes, complex_ax)

        phase_ax = Axis(
            fig[2, col];
            xlabel = L"t/T_{\mathrm{H}}",
            ylabel = col == 1 ? L"\psi(t)-\psi(0)" : "",
        )
        style_phase_axis!(phase_ax)
        phase_ax.yticks = ([-6pi, -4pi, -2pi, 0.0], [L"-6\pi", L"-4\pi", L"-2\pi", L"0"])
        text!(phase_ax, 0.035, 0.065; text = bottom_labels[col], space = :relative, align = (:left, :bottom), fontsize = LABEL_PT, font = :bold)
        q_label = diagnostic.q < -0.5 ? L"q=-1" : L"q=0"
        text!(phase_ax, 0.965, 0.80; text = q_label, space = :relative, align = (:right, :top), fontsize = 0.78 * LABEL_PT)
        push!(phase_axes, phase_ax)

        z = diagnostic.rotating
        if diagnostic.fixed_point
            scatter!(complex_ax, [real(z[1])], [imag(z[1])]; color = DARK_BLUE, markersize = 1.15 * MARKER_PT, strokecolor = :white, strokewidth = 0.35)
        else
            lines!(complex_ax, real.(z), imag.(z); color = DARK_BLUE, linewidth = 0.82)
            add_direction_markers!(complex_ax, z)
        end
        scatter!(complex_ax, [0.0], [0.0]; marker = :xcross, color = :black, markersize = 5.0, strokewidth = 0.8)
        normalized_time = diagnostic.times ./ diagnostic.modulation_period
        lines!(phase_ax, normalized_time, diagnostic.psi; color = DARK_BLUE, linewidth = 0.82)
        xlims!(phase_ax, 0.0, ROTATING_FRAME_PERIODS)
        phase_ax.xticks = col == 1 ? ([0.0, 1.0, 2.0, 3.0], ["0", "1", "2", "3"]) : ([1.0, 2.0, 3.0], ["1", "2", "3"])

        if col != 1
            hideydecorations!(complex_ax; grid = false)
            hideydecorations!(phase_ax; grid = false)
        end
    end

    linkyaxes!(complex_axes...)
    linkyaxes!(phase_axes...)
    rowsize!(fig.layout, 1, Aspect(1, 1.0))
    rowsize!(fig.layout, 2, Fixed(1.10 * 72.0))
    rowgap!(fig.layout, 0.42 * LATEX_BASE_PT)
    colgap!(fig.layout, 0.50 * LATEX_BASE_PT)

    png_path = path_stem * ".png"
    pdf_path = path_stem * ".pdf"
    save(png_path, fig; px_per_unit = 3)
    save(pdf_path, fig)
    return png_path, pdf_path
end

function style_axis!(ax)
    ax.aspect = DataAspect()
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.09)
    ax.ygridcolor = (:black, 0.09)
    ax.xticklabelsize = TICK_PT
    ax.yticklabelsize = TICK_PT
    ax.xlabelsize = LABEL_PT
    ax.ylabelsize = LABEL_PT
    ax.xlabelpadding = 1.0
    ax.ylabelpadding = 1.0
    ax.titlesize = TITLE_PT
    ax.spinewidth = 0.65
    ax.xtickwidth = 0.65
    ax.ytickwidth = 0.65
    ax.xticksize = 2.8
    ax.yticksize = 2.8
    return ax
end

function style_phase_axis!(ax)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.09)
    ax.ygridcolor = (:black, 0.09)
    ax.xticklabelsize = TICK_PT
    ax.yticklabelsize = TICK_PT
    ax.xlabelsize = LABEL_PT
    ax.ylabelsize = LABEL_PT
    ax.xlabelpadding = 1.0
    ax.ylabelpadding = 1.0
    ax.spinewidth = 0.65
    ax.xtickwidth = 0.65
    ax.ytickwidth = 0.65
    ax.xticksize = 2.8
    ax.yticksize = 2.8
    return ax
end

const DIAGNOSTICS_FILE = joinpath(DATA_ROOT, "derived", "fig3_rotating_frame.jld2")

function generate_data()
    diagnostics = rotating_frame_diagnostics()
    mkpath(dirname(DIAGNOSTICS_FILE))
    JLD2.jldsave(DIAGNOSTICS_FILE; diagnostics, N = length(load_restart(restart_file_for_K(1.0)).omega))
end

function main()
    require_file(DIAGNOSTICS_FILE)
    diagnostics = JLD2.load(DIAGNOSTICS_FILE, "diagnostics")
    save_low_only_with_phase_figure(joinpath(OUTPUT_ROOT, "Fig_3"), diagnostics)
    @info "Saved Figure 3"
end

end
