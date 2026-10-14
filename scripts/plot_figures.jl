include("cli.jl")
opts = parse_options(ARGS; generation = false)
if "help" in opts.flags
    print_help(false)
    exit()
end

using AmpDefInc

check_supplement_inputs(opts.supplements; plotting = true, smoke = "smoke" in opts.flags)

if "smoke" in opts.flags && !isempty(opts.figures)
    require_file(joinpath(DATA_ROOT, "numerical", "K_continuation_order_parameters.csv"))
    using CairoMakie
    lines = readlines(joinpath(DATA_ROOT, "numerical", "K_continuation_order_parameters.csv"))
    values = [parse.(Float64, split(line, ',')) for line in lines[2:end]]
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1]; xlabel = "K", ylabel = "R", title = "Pipeline smoke check (N=32)")
    lines!(ax, getindex.(values, 1), getindex.(values, 2))
    mkpath(OUTPUT_ROOT)
    save(joinpath(OUTPUT_ROOT, "smoke.png"), fig)
    save(joinpath(OUTPUT_ROOT, "smoke.pdf"), fig)
else
    "smoke" in opts.flags || check_inputs(opts.figures; plotting = true)
    mkpath(OUTPUT_ROOT)
    for n in opts.figures
        figure = load_figure(n)
        Base.invokelatest(Base.invokelatest(getfield, figure, :main))
    end
end
if !isempty(opts.supplements)
    mkpath(joinpath(OUTPUT_ROOT, "supplementary"))
    for id in opts.supplements
        figure = load_supplement(id)
        fn = Base.invokelatest(getfield, figure, :main)
        if startswith(id, "supp-stability")
            Base.invokelatest(fn; zoom = endswith(id, "-zoom"), smoke = "smoke" in opts.flags)
        else
            Base.invokelatest(fn)
        end
    end
end
println("Figures ready in $OUTPUT_ROOT")
