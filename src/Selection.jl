const MAIN_FIGURES = [2, 3, 4, 6]
const SUPPLEMENT_IDS = ["supp-kop", "supp-timeseries", "supp-timeseries-book2",
                        "supp-stability", "supp-stability-zoom"]

function parse_figure_selection(value::AbstractString)
    tokens = strip.(split(value, ','))
    tokens == ["all"] && (tokens = vcat(string.(MAIN_FIGURES), SUPPLEMENT_IDS))
    tokens == ["supplementary"] && (tokens = copy(SUPPLEMENT_IDS))
    all(token -> token in vcat(string.(MAIN_FIGURES), SUPPLEMENT_IDS), tokens) ||
        error("Choose 2,3,4,6, all, supplementary, or supplementary IDs: $(join(SUPPLEMENT_IDS, ',')). Raw data for Figures 1 and 5 are provided on Zenodo; no plotting scripts are included.")
    main = [parse(Int, token) for token in unique(tokens) if token in string.(MAIN_FIGURES)]
    supplements = [token for token in unique(tokens) if token in SUPPLEMENT_IDS]
    return main, supplements
end
