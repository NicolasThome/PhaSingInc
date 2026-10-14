module SupplementKOP
using ..AmpDefInc
using CairoMakie
using LaTeXStrings
using Printf

const FIG_DIR = joinpath(OUTPUT_ROOT, "supplementary")
const DATA_FILE = joinpath(DATA_ROOT, "supplementary", "experimental", "KOP_Exp.xlsx")
const OUTPUT_STEM = joinpath(FIG_DIR, "Fig_KOP_exp_broken_axis")

# Half-width figure: use doubled source typography, then include at 0.5\textwidth.
const LATEX_BASE_PT = 24.0
const LATEX_FIGURE_WIDTH_PT = 6.5 * 72.0
const FIGURE_SIZE_PT = (LATEX_FIGURE_WIDTH_PT, 0.62 * LATEX_FIGURE_WIDTH_PT)
const TICK_PT = 1.00 * LATEX_BASE_PT
const LABEL_PT = 1.00 * LATEX_BASE_PT
const PANEL_PT = 1.00 * LATEX_BASE_PT
const LINE_PT = 2.20
const MARKER_PT = 9.5
const HIGH_ILLUMINATION_DISPLAY = 1.075

const PAPER_FONT = AmpDefInc.PAPER_FONT

using ..Figure2: read_first_sheet_xlsx

function as_float(x)
    x isa Number && return Float64(x)
    x isa AbstractString && return parse(Float64, replace(x, "," => "."))
    return missing
end

function read_kop_data(path::AbstractString)
    sheet = read_first_sheet_xlsx(path)
    header = [lowercase(String(sheet[1, col])) for col in axes(sheet, 2)]

    illumination_col = findfirst(==("illumination"), header)
    label_col = findfirst(==("number"), header)
    order_col = findfirst(==("order parameter"), header)
    any(isnothing, (illumination_col, label_col, order_col)) &&
        error("Could not find illumination, number, and Order parameter columns in $(basename(path)).")

    rows = NamedTuple[]
    for row in 2:size(sheet, 1)
        illumination = as_float(sheet[row, illumination_col])
        order_parameter = as_float(sheet[row, order_col])
        label = sheet[row, label_col]
        (illumination === missing || order_parameter === missing || label === missing) && continue
        push!(rows, (; I = illumination, R = order_parameter, label = String(label)))
    end
    return rows
end

function style_axis!(ax)
    ax.xgridvisible = false
    ax.ygridvisible = true
    ax.xgridcolor = (:black, 0.10)
    ax.ygridcolor = (:black, 0.10)
    ax.xticklabelsize = TICK_PT
    ax.yticklabelsize = TICK_PT
    ax.xlabelsize = LABEL_PT
    ax.ylabelsize = LABEL_PT
    ax.spinewidth = 0.65 * LINE_PT
    ax.xtickwidth = 0.65 * LINE_PT
    ax.ytickwidth = 0.65 * LINE_PT
    ax.xticksize = 3.0
    ax.yticksize = 3.0
    ax.xlabelpadding = 3.0
    ax.ylabelpadding = 3.0
    ax.xticklabelpad = 2.0
    ax.yticklabelpad = 2.0
    return ax
end

function plot_segment!(ax, rows; color)
    sorted = sort(rows; by = row -> row.I)
    xs = [row.I for row in sorted]
    ys = [row.R for row in sorted]
    lines!(ax, xs, ys; color, linewidth = LINE_PT)
    scatter!(ax, xs, ys; color, marker = :circle, markersize = MARKER_PT, strokecolor = :white, strokewidth = 0.35 * LINE_PT)
end

display_illumination(row) = row.I >= 1.5 ? HIGH_ILLUMINATION_DISPLAY : row.I

function annotate_selected_points!(ax, rows, labels_by_illumination)
    for row in rows
        label = get(labels_by_illumination, round(row.I; digits = 3), nothing)
        label === nothing && continue
        x = display_illumination(row)
        offset = row.I > 2 ? (-0.40 * LABEL_PT, -0.18 * LABEL_PT) : (0, 0.42 * LABEL_PT)
        align = row.I > 2 ? (:right, :top) : (:center, :bottom)
        text!(
            ax,
            x,
            row.R;
            text = label,
            align,
            offset,
            fontsize = 0.78 * LABEL_PT,
            color = :black,
        )
    end
end

function add_bottom_axis_break_marks!(ax)
    ymin, ymax = 0.30, 1.03
    x_break = 1.035
    dx = 0.0038
    slash_gap = 0.010
    dy = 0.015 * (ymax - ymin)
    mark_linewidth = 0.55 * LINE_PT

    for x0 in (x_break, x_break + slash_gap)
        lines!(
            ax,
            [x0 - dx, x0 + dx],
            [ymin, ymin + 2dy];
            color = :black,
            linewidth = mark_linewidth,
            linecap = :butt,
        )
    end
    return nothing
end

function main()
    rows = read_kop_data(DATA_FILE)
    left_rows = [row for row in rows if row.I < 1.5]
    right_rows = [row for row in rows if row.I >= 1.5]
    high_row = only(right_rows)

    fig = Figure(
        size = FIGURE_SIZE_PT,
        backgroundcolor = :white,
        figure_padding = (0.30 * LATEX_BASE_PT, 0.25 * LATEX_BASE_PT, 0.28 * LATEX_BASE_PT, 0.18 * LATEX_BASE_PT),
        fonts = (; regular = PAPER_FONT, bold = PAPER_FONT, italic = PAPER_FONT, bold_italic = PAPER_FONT),
        fontsize = LATEX_BASE_PT,
    )

    ax = Axis(fig[1, 1]; xlabel = L"I\ [\mathrm{mW}/\mathrm{cm}^2]", ylabel = L"R")
    style_axis!(ax)

    color = RGBf(0.0, 0.32, 0.62)
    plot_segment!(ax, left_rows; color)
    scatter!(
        ax,
        [HIGH_ILLUMINATION_DISPLAY],
        [high_row.R];
        color,
        marker = :circle,
        markersize = MARKER_PT,
        strokecolor = :white,
        strokewidth = 0.35 * LINE_PT,
    )

    selected_labels = Dict(3.13 => "(a)", 0.94 => "(b)", 0.85 => "(c)")
    annotate_selected_points!(ax, rows, selected_labels)

    ylims!(ax, 0.30, 1.03)
    xlims!(ax, 0.75, 1.105)

    ax.xticks = ([0.80, 0.90, 1.00, HIGH_ILLUMINATION_DISPLAY], ["0.80", "0.90", "1.00", "3.13"])
    ax.yticks = ([0.4, 0.6, 0.8, 1.0], ["0.4", "0.6", "0.8", "1.0"])
    add_bottom_axis_break_marks!(ax)

    png_path = OUTPUT_STEM * ".png"
    pdf_path = OUTPUT_STEM * ".pdf"
    save(png_path, fig; px_per_unit = 3)
    save(pdf_path, fig)
    @info "Saved experimental KOP broken-axis figure" png_path pdf_path n_points = length(rows)
end

end
