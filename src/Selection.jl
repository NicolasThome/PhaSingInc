const SUPPLEMENT_IDS = ["supp-kop", "supp-timeseries", "supp-timeseries-book2",
                        "supp-stability", "supp-stability-zoom"]

function parse_figure_selection(value::AbstractString)
    tokens = strip.(split(value, ','))
    tokens == ["all"] && (tokens = vcat(string.(2:5), SUPPLEMENT_IDS))
    tokens == ["supplementary"] && (tokens = copy(SUPPLEMENT_IDS))
    all(token -> token in vcat(string.(2:5), SUPPLEMENT_IDS), tokens) ||
        error("Choose 2,3,4,5, all, supplementary, or supplementary IDs: $(join(SUPPLEMENT_IDS, ','))")
    main = [parse(Int, token) for token in unique(tokens) if token in string.(2:5)]
    supplements = [token for token in unique(tokens) if token in SUPPLEMENT_IDS]
    return main, supplements
end
