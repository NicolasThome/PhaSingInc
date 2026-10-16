module SupplementTimeseries
using ..PhaSingInc
using CairoMakie
using LaTeXStrings

const FIG_DIR = joinpath(OUTPUT_ROOT, "supplementary")
const DATA_FILE = joinpath(DATA_ROOT, "supplementary", "experimental", "min(r).xlsx")
const OUTPUT_STEM = joinpath(FIG_DIR, "Fig_min_r_timeseries")

const LATEX_BASE_PT = 12.0
const LATEX_FIGURE_WIDTH_PT = 6.5 * 72.0
const FIGURE_SIZE_PT = (LATEX_FIGURE_WIDTH_PT, 0.62 * LATEX_FIGURE_WIDTH_PT)
const TICK_PT = 0.95 * LATEX_BASE_PT
const LABEL_PT = 1.00 * LATEX_BASE_PT
const PANEL_PT = 1.00 * LATEX_BASE_PT
const LINE_PT = 0.95
const FIT_LINE_PT = 1.45
const T_LIMITS = (0.0, 1250.0)
const T_AXIS_MARGIN = 25.0
const T_TICKS = ([0, 250, 500, 750, 1000, 1250], ["0", "250", "500", "750", "1000", "1250"])

const PAPER_FONT = PhaSingInc.PAPER_FONT

using ..Figure2
read_all_sheets_xlsx(path) = Dict(Figure2.read_all_sheets_xlsx(path))

function as_float(x)
    x isa Number && return Float64(x)
    x isa AbstractString && return parse(Float64, replace(strip(x), "," => "."))
    return missing
end

function fit_parameter_value(x)
    x isa Number && return Float64(x)
    if x isa AbstractString
        number_match = match(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", x)
        number_match === nothing && return missing
        return parse(Float64, number_match.match)
    end
    return missing
end

function read_timeseries_sheet(sheet)
    header = [cell === missing ? "" : lowercase(String(cell)) for cell in sheet[1, :]]
    t_col = findfirst(==("t"), header)
    emsi_col = findfirst(==("emsi"), header)
    any(isnothing, (t_col, emsi_col)) && error("Could not find t and EMSI columns.")

    t = Float64[]
    emsi = Float64[]
    for row in 2:size(sheet, 1)
        t_value = as_float(sheet[row, t_col])
        emsi_value = as_float(sheet[row, emsi_col])
        (t_value === missing || emsi_value === missing) && continue
        push!(t, t_value)
        push!(emsi, emsi_value)
    end

    params = Dict{String, Float64}()
    for row in axes(sheet, 1)
        key = sheet[row, 5]
        value = sheet[row, 6]
        key isa AbstractString || continue
        key_string = strip(key)
        key_string in ("y0", "xc", "w", "A") || continue
        parsed = fit_parameter_value(value)
        parsed === missing && error("Could not parse fit parameter '$key_string'.")
        params[key_string] = parsed
    end

    for key in ("y0", "xc", "w", "A")
        haskey(params, key) || error("Missing fit parameter '$key'.")
    end

    return (; t, emsi, y0 = params["y0"], xc = params["xc"], w = params["w"], A = params["A"])
end

fit_curve(data, t) = data.y0 .+ data.A .* sin.(pi .* (t .- data.xc) ./ data.w)

function style_axis!(ax)
    ax.xgridvisible = true
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

function panel_label!(ax, label, illumination)
    text!(
        ax,
        0.020,
        0.945;
        text = label,
        space = :relative,
        align = (:left, :top),
        fontsize = PANEL_PT,
        font = :bold,
        color = :black,
    )
    text!(
        ax,
        0.985,
        0.945;
        text = illumination,
        space = :relative,
        align = (:right, :top),
        fontsize = LABEL_PT,
        color = :black,
    )
    return ax
end

function plot_timeseries_panel!(ax, data)
    data_color = RGBf(0.0, 0.20, 0.45)
    fit_color = RGBf(0.86, 0.18, 0.12)
    lines!(ax, data.t, data.emsi; color = data_color, linewidth = LINE_PT)
    t_fit = collect(range(first(T_LIMITS), last(T_LIMITS), length = 1800))
    lines!(ax, t_fit, fit_curve(data, t_fit); color = fit_color, linewidth = FIT_LINE_PT)
    return ax
end

function main()
    sheets = read_all_sheets_xlsx(DATA_FILE)
    sheet_names = ("I=0.94", "I=0.89")
    for name in sheet_names
        haskey(sheets, name) || error("Could not find sheet '$name' in $(basename(DATA_FILE)).")
    end

    data_top = read_timeseries_sheet(sheets["I=0.94"])
    data_bottom = read_timeseries_sheet(sheets["I=0.89"])

    fig = Figure(
        size = FIGURE_SIZE_PT,
        backgroundcolor = :white,
        figure_padding = (0.45 * LATEX_BASE_PT, 0.80 * LATEX_BASE_PT, 0.35 * LATEX_BASE_PT, 0.25 * LATEX_BASE_PT),
        fonts = (; regular = PAPER_FONT, bold = PAPER_FONT, italic = PAPER_FONT, bold_italic = PAPER_FONT),
        fontsize = LATEX_BASE_PT,
    )

    ax_a = Axis(fig[1, 1]; ylabel = L"A(t)\ [\mathrm{a.u.}]")
    ax_b = Axis(fig[2, 1]; xlabel = L"t\ [\mathrm{s}]", ylabel = L"A(t)\ [\mathrm{a.u.}]")
    style_axis!(ax_a)
    style_axis!(ax_b)

    plot_timeseries_panel!(ax_a, data_top)
    plot_timeseries_panel!(ax_b, data_bottom)
    panel_label!(ax_a, "(a)", L"I = 0.94\ \mathrm{mW/cm^2}")
    panel_label!(ax_b, "(b)", L"I = 0.89\ \mathrm{mW/cm^2}")

    ax_a.xlabelvisible = false
    ylims!(ax_a, 0.06, 0.68)
    ylims!(ax_b, 0.06, 0.68)
    xlims!(ax_a, first(T_LIMITS) - T_AXIS_MARGIN, last(T_LIMITS) + T_AXIS_MARGIN)
    xlims!(ax_b, first(T_LIMITS) - T_AXIS_MARGIN, last(T_LIMITS) + T_AXIS_MARGIN)
    ax_a.xticks = T_TICKS
    ax_b.xticks = T_TICKS
    ax_a.yticks = ([0.1, 0.3, 0.5], ["0.1", "0.3", "0.5"])
    ax_b.yticks = ([0.1, 0.3, 0.5], ["0.1", "0.3", "0.5"])

    rowgap!(fig.layout, 0.70 * LATEX_BASE_PT)

    png_path = OUTPUT_STEM * ".png"
    pdf_path = OUTPUT_STEM * ".pdf"
    save(png_path, fig; px_per_unit = 3)
    save(pdf_path, fig)
    @info "Saved min-r time-series figure" png_path pdf_path top_sheet=sheet_names[1] bottom_sheet=sheet_names[2]
end


end
