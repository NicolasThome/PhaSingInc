include(joinpath(@__DIR__, "..", "src", "Selection.jl"))

function parse_options(args; generation::Bool)
    options = Dict("data-dir" => joinpath(@__DIR__, "..", "data"),
                   "output-dir" => joinpath(@__DIR__, "..", "figures"),
                   "figures" => "2,3,4,5")
    flags = Set{String}()
    allowed = generation ? Set(["help", "smoke", "skip-continuation", "continuation-only", "skip-stability", "force"]) : Set(["help", "smoke"])
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unexpected argument: $arg")
        key = arg[3:end]
        if key in allowed
            push!(flags, key)
        elseif haskey(options, key)
            i < length(args) || error("Missing value for --$key")
            options[key] = args[i + 1]
            i += 1
        else
            error("Unknown option: $arg")
        end
        i += 1
    end
    figures, supplements = parse_figure_selection(options["figures"])
    ENV["AMPDEFINC_DATA"] = abspath(options["data-dir"])
    ENV["AMPDEFINC_OUTPUT"] = abspath(options["output-dir"])
    return (; options, flags, figures, supplements)
end

function print_help(generation)
    script = generation ? "generate_data.jl" : "plot_figures.jl"
    println("Usage: julia --project=. scripts/$script [options]")
    println("  --data-dir PATH       data root (default: repository/data)")
    println("  --output-dir PATH     figure root (default: repository/figures)")
    println("  --figures SELECTION   default: 2,3,4,5; also all or supplementary")
    println("                        supplementary IDs: $(join(SUPPLEMENT_IDS, ','))")
    println("  --smoke               small pipeline check, not paper reproduction")
    if generation
        println("  --skip-continuation   use existing numerical CSVs and restart states")
        println("  --continuation-only   generate continuation CSVs and restart states only")
        println("  --skip-stability      retain archived supplementary stability data")
        println("  --force               replace an existing continuation dataset")
    end
end
