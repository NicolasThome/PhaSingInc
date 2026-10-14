module SupplementTimeseriesBook2
using ..AmpDefInc
using CairoMakie
using LaTeXStrings

using ..SupplementTimeseries: FIG_DIR, FIGURE_SIZE_PT, LATEX_BASE_PT, PAPER_FONT,
    T_LIMITS, T_AXIS_MARGIN, T_TICKS, LINE_PT, FIT_LINE_PT,
    read_all_sheets_xlsx, as_float, fit_parameter_value, fit_curve,
    style_axis!, panel_label!

const BOOK2_DATA_FILE = joinpath(DATA_ROOT, "supplementary", "experimental", "Book2.xlsx")
const BOOK2_OUTPUT_STEM = joinpath(FIG_DIR, "Fig_min_r_timeseries_book2")

function read_timeseries_block(sheet, data_col::Int, fit_col::Int)
    t = Float64[]
    amp = Float64[]
    for row in 3:size(sheet, 1)
        t_value = as_float(sheet[row, data_col])
        amp_value = as_float(sheet[row, data_col + 1])
        (t_value === missing || amp_value === missing) && continue
        push!(t, t_value)
        push!(amp, amp_value)
    end

    params = Dict{String, Float64}()
    for row in axes(sheet, 1)
        key = sheet[row, fit_col]
        value = sheet[row, fit_col + 1]
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

    return (; t, emsi = amp, y0 = params["y0"], xc = params["xc"], w = params["w"], A = params["A"])
end

function main()
    sheets = read_all_sheets_xlsx(BOOK2_DATA_FILE)
    haskey(sheets, "Sheet1") || error("Could not find Sheet1 in $(basename(BOOK2_DATA_FILE)).")
    sheet = sheets["Sheet1"]

    data_top = read_timeseries_block(sheet, 1, 4)
    data_bottom = read_timeseries_block(sheet, 7, 10)

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

    plot_book2_timeseries_panel!(ax_a, data_top)
    plot_book2_timeseries_panel!(ax_b, data_bottom)
    panel_label!(ax_a, "(a)", L"I = 0.94\ \mathrm{mW/cm^2}")
    panel_label!(ax_b, "(b)", L"I = 0.85\ \mathrm{mW/cm^2}")

    ax_a.xlabelvisible = false
    ylims!(ax_a, 0.06, 0.95)
    ylims!(ax_b, 0.06, 0.95)
    xlims!(ax_a, first(T_LIMITS) - T_AXIS_MARGIN, last(T_LIMITS) + T_AXIS_MARGIN)
    xlims!(ax_b, first(T_LIMITS) - T_AXIS_MARGIN, last(T_LIMITS) + T_AXIS_MARGIN)
    ax_a.xticks = T_TICKS
    ax_b.xticks = T_TICKS
    ax_a.yticks = ([0.2, 0.5, 0.8], ["0.2", "0.5", "0.8"])
    ax_b.yticks = ([0.2, 0.5, 0.8], ["0.2", "0.5", "0.8"])

    rowgap!(fig.layout, 0.70 * LATEX_BASE_PT)

    png_path = BOOK2_OUTPUT_STEM * ".png"
    pdf_path = BOOK2_OUTPUT_STEM * ".pdf"
    save(png_path, fig; px_per_unit = 3)
    save(pdf_path, fig)
    @info "Saved Book2 min-r time-series figure" png_path pdf_path
end

function plot_book2_timeseries_panel!(ax, data)
    data_color = RGBf(0.86, 0.18, 0.12)
    fit_color = RGBf(0.0, 0.20, 0.65)
    lines!(ax, data.t, data.emsi; color = data_color, linewidth = LINE_PT)
    t_fit = collect(range(first(T_LIMITS), last(T_LIMITS), length = 1800))
    lines!(ax, t_fit, fit_curve(data, t_fit); color = fit_color, linewidth = FIT_LINE_PT)
    return ax
end


end
