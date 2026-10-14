module SupplementStability
using ..AmpDefInc
using CairoMakie
using JLD2
using LaTeXStrings

const FIG_DIR = joinpath(OUTPUT_ROOT, "supplementary")
const DATA_FILE = joinpath(DATA_ROOT, "supplementary", "numerical", "locked_spectrum.jld2")
const OUTPUT_STEM = joinpath(FIG_DIR, "Fig_linear_stability_localization")
const ZOOM_OUTPUT_STEM = joinpath(FIG_DIR, "Fig_linear_stability_localization_spectrum_xlim_m0p2_0p2")

const LATEX_BASE_PT = 12.0
const LATEX_FIGURE_WIDTH_PT = 6.5 * 72.0
const FIGURE_SIZE_PT = (LATEX_FIGURE_WIDTH_PT, 0.71 * LATEX_FIGURE_WIDTH_PT)
const TICK_PT = 1.00 * LATEX_BASE_PT
const LABEL_PT = 1.00 * LATEX_BASE_PT
const PANEL_PT = 1.00 * LATEX_BASE_PT
const RED = RGBf(0.86, 0.18, 0.12)

const PAPER_FONT = AmpDefInc.PAPER_FONT

function style_axis!(ax)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.08)
    ax.ygridcolor = (:black, 0.08)
    ax.xticklabelsize = TICK_PT
    ax.yticklabelsize = TICK_PT
    ax.xlabelsize = LABEL_PT
    ax.ylabelsize = LABEL_PT
    ax.spinewidth = 0.85
    ax.xtickwidth = 0.85
    ax.ytickwidth = 0.85
    ax.xticksize = 3.0
    ax.yticksize = 3.0
    return ax
end

function panel_label!(ax, label; position = (0.025, 0.97), align = (:left, :top))
    text!(
        ax,
        position...;
        text = label,
        space = :relative,
        align,
        fontsize = PANEL_PT,
        font = :bold,
        color = :black,
    )
    return ax
end

function finite_spectrum_data(all_eigenvalues, K_path, accepted)
    xs = Float64[]
    ys = Float64[]
    colors = Float64[]
    sizehint!(xs, count(accepted) * size(all_eigenvalues, 2))
    sizehint!(ys, count(accepted) * size(all_eigenvalues, 2))
    sizehint!(colors, count(accepted) * size(all_eigenvalues, 2))

    for row in findall(accepted)
        for value in @view all_eigenvalues[row, :]
            if isfinite(real(value)) && isfinite(imag(value))
                push!(xs, real(value))
                push!(ys, imag(value))
                push!(colors, K_path[row])
            end
        end
    end
    return xs, ys, colors
end

function main(; zoom::Bool = false, smoke::Bool = false)
    data = JLD2.load(DATA_FILE)
    K_path = Float64.(data["K_path"])
    omega = Float64.(data["omega"])
    accepted = Bool.(data["accepted"])
    all_eigenvalues = ComplexF64.(data["all_eigenvalues"])
    leading_eigenvalues = ComplexF64.(data["leading_eigenvalues"])
    leading_weights = Float64.(data["leading_weights"])
    locked_amplitudes = Float64.(data["locked_amplitudes"])
    critical_index = Int(data["critical_index"])
    dominant_indices = Int.(data["dominant_index"])
    K_hopf_values = Float64.(data["candidate_K_hopf"])

    critical_index > 0 || error("No critical Hopf index is stored in $(basename(DATA_FILE)).")
    dominant_index = dominant_indices[critical_index]
    dominant_index > 0 || error("No dominant oscillator is stored at the critical Hopf index.")

    spectrum_x, spectrum_y, spectrum_K = finite_spectrum_data(all_eigenvalues, K_path, accepted)
    if zoom
        visible = (-0.2 .<= spectrum_x) .& (spectrum_x .<= 0.2)
        spectrum_x = spectrum_x[visible]
        spectrum_y = spectrum_y[visible]
        spectrum_K = spectrum_K[visible]
    end
    critical_weights = vec(leading_weights[critical_index, :])
    critical_amplitudes = vec(locked_amplitudes[critical_index, :])
    omega_order = sortperm(omega)
    accepted_indices = findall(accepted)

    fig = Figure(
        size = FIGURE_SIZE_PT,
        backgroundcolor = :white,
        figure_padding = (0.45 * LATEX_BASE_PT, 0.30 * LATEX_BASE_PT, 0.35 * LATEX_BASE_PT, 0.30 * LATEX_BASE_PT),
        fonts = (; regular = PAPER_FONT, bold = PAPER_FONT, italic = PAPER_FONT, bold_italic = PAPER_FONT),
        fontsize = LATEX_BASE_PT,
    )

    top = GridLayout(fig[1, 1])
    bottom = GridLayout(fig[2, 1])

    ax_a = Axis(
        top[1, 1];
        xlabel = L"\mathrm{Re}(\lambda)",
        ylabel = L"\mathrm{Im}(\lambda)",
    )
    style_axis!(ax_a)
    vlines!(ax_a, [0.0]; color = (:black, 0.50), linestyle = :dash, linewidth = 0.85)
    spectrum_plot = scatter!(
        ax_a,
        spectrum_x,
        spectrum_y;
        color = spectrum_K,
        colormap = :viridis,
        colorrange = extrema(K_path[accepted_indices]),
        markersize = 1.20,
        rasterize = zoom ? false : 5,
    )
    Colorbar(
        top[1, 2],
        spectrum_plot;
        label = L"K",
        labelsize = LABEL_PT,
        ticklabelsize = TICK_PT,
        width = 0.60 * LATEX_BASE_PT,
        ticks = 0.95:0.01:1.00,
    )
    colgap!(top, 0.35 * LATEX_BASE_PT)

    ax_b = Axis(
        bottom[1, 1];
        xlabel = L"\omega_k",
        ylabel = L"p_k",
        limits = (nothing, nothing, -0.035, 1.035),
        xticks = ([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1", "", L"0", "", L"1"]),
        yticks = [0.0, 0.5, 1.0],
    )
    style_axis!(ax_b)
    scatter!(ax_b, omega[omega_order], critical_weights[omega_order]; color = :black, markersize = 1.8)
    scatter!(
        ax_b,
        [omega[dominant_index]],
        [critical_weights[dominant_index]];
        color = RED,
        strokecolor = :black,
        strokewidth = 0.55,
        markersize = 5.5,
    )

    amplitude_limits = extrema(critical_amplitudes[isfinite.(critical_amplitudes)])
    amplitude_padding = 0.035 * (amplitude_limits[2] - amplitude_limits[1])
    ax_c = Axis(
        bottom[1, 2];
        xlabel = L"\omega_k",
        ylabel = L"r_k",
        xticks = ([-1.0, -0.5, 0.0, 0.5, 1.0], [L"-1", "", L"0", "", L"1"]),
        limits = (
            nothing,
            nothing,
            amplitude_limits[1] - amplitude_padding,
            amplitude_limits[2] + amplitude_padding,
        ),
    )
    style_axis!(ax_c)
    scatter!(ax_c, omega[omega_order], critical_amplitudes[omega_order]; color = :black, markersize = 1.8)
    scatter!(
        ax_c,
        [omega[dominant_index]],
        [critical_amplitudes[dominant_index]];
        color = RED,
        strokecolor = :black,
        strokewidth = 0.55,
        markersize = 5.5,
    )

    ax_d = Axis(
        bottom[1, 3];
        xlabel = L"K",
        ylabel = L"\mathrm{Re}(\lambda_k)",
    )
    style_axis!(ax_d)
    hlines!(ax_d, [0.0]; color = (:black, 0.50), linestyle = :dash, linewidth = 0.80)
    lines!(
        ax_d,
        K_path[accepted_indices],
        real.(leading_eigenvalues[accepted_indices]);
        color = :black,
        linewidth = 0.95,
    )
    scatter!(
        ax_d,
        K_path[accepted_indices],
        real.(leading_eigenvalues[accepted_indices]);
        color = :black,
        markersize = 2.2,
    )
    if !isempty(K_hopf_values)
        scatter!(ax_d, [first(K_hopf_values)], [0.0]; color = RED, markersize = 5.5)
    end

    # Draw panel tags last so dense scatter layers cannot cover them.
    panel_label!(ax_a, "(a)")
    panel_label!(ax_b, "(b)"; position = (0.975, 0.995), align = (:right, :top))
    panel_label!(ax_c, "(c)"; position = (0.975, 0.995), align = (:right, :top))
    panel_label!(ax_d, "(d)"; position = (0.975, 0.995), align = (:right, :top))

    colgap!(bottom, 0.55 * LATEX_BASE_PT)
    rowgap!(fig.layout, 0.22 * LATEX_BASE_PT)
    rowsize!(fig.layout, 1, Relative(0.61))
    rowsize!(fig.layout, 2, Relative(0.39))

    zoom && xlims!(ax_a, -0.2, 0.2)
    output_stem = zoom ? ZOOM_OUTPUT_STEM : OUTPUT_STEM
    smoke && (output_stem *= "_smoke")
    png_path = output_stem * ".png"
    pdf_path = output_stem * ".pdf"
    save(png_path, fig; px_per_unit = 3)
    save(pdf_path, fig)
    @info "Saved linear-stability and mode-localization figure" png_path pdf_path critical_index dominant_index zoom_spectrum=zoom
end


end
