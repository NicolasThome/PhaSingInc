include("cli.jl")
opts = parse_options(ARGS; generation = true)
if "help" in opts.flags
    print_help(true)
    exit()
end

using AmpDefInc

if "smoke" in opts.flags && "skip-continuation" in opts.flags
    error("--smoke cannot be combined with --skip-continuation.")
end
if "smoke" in opts.flags && any(id -> !(id in ("supp-stability", "supp-stability-zoom")), opts.supplements)
    error("Experimental supplements use measured workbooks. For a supplementary smoke check select --figures supp-stability.")
end
if !("continuation-only" in opts.flags)
    check_supplement_inputs(opts.supplements; plotting = false,
        smoke = "smoke" in opts.flags, skip_stability = "skip-stability" in opts.flags)
end
needs_stability = !("continuation-only" in opts.flags) &&
    any(id -> startswith(id, "supp-stability"), opts.supplements) && !("skip-stability" in opts.flags)
if needs_stability
    stability_path = joinpath(DATA_ROOT, "supplementary", "numerical", "locked_spectrum.jld2")
    if isfile(stability_path)
        "force" in opts.flags || error("Supplementary stability data already exist. Use --skip-stability to retain them, or --force to replace matching-N data.")
        expected_n = "smoke" in opts.flags ? 12 : 5000
        old_n = AmpDefInc.JLD2.load(stability_path, "N")
        old_n == expected_n || error("Existing supplementary data have N=$old_n; requested N=$expected_n. Use a different --data-dir.")
    end
end
needs_primary = !isempty(opts.figures) || "continuation-only" in opts.flags
if needs_primary && !("skip-continuation" in opts.flags)
    numerical = joinpath(DATA_ROOT, "numerical")
    cfg = "smoke" in opts.flags ? smoke_config() : paper_config()
    existing = isdir(numerical) && !isempty(readdir(numerical))
    existing && !("force" in opts.flags) && error("Numerical data already exist at $numerical. Use --skip-continuation to retain them or --force to replace them.")
    restart_dir = joinpath(numerical, "K_continuation_restart")
    if existing && isdir(restart_dir)
        restarts = filter(name -> endswith(name, ".jld2"), readdir(restart_dir; join=true))
        if !isempty(restarts)
            old_n = AmpDefInc.JLD2.load(first(restarts), "N")
            old_n == cfg.N || error("Existing data have N=$old_n; requested N=$(cfg.N). Use a different --data-dir to keep paper and smoke data separate.")
        end
    end
    # Check required experimental inputs before starting an expensive run.
    if !("smoke" in opts.flags || "continuation-only" in opts.flags)
        for name in (2 in opts.figures ? ["minR_meanF_cont.xlsx", "9_points.xlsx"] : String[])
            require_file(joinpath(DATA_ROOT, "experimental", name))
        end
        4 in opts.figures && require_file(joinpath(DATA_ROOT, "experimental", "FC_FFT_data.xlsx"))
    end
    timing = @timed run_k_continuation(cfg)
    println("Continuation elapsed: $(round(timing.time; digits=2)) s; Julia allocations: $(round(timing.bytes / 2.0^20; digits=1)) MiB (not peak RAM).")
end

if !("smoke" in opts.flags || "continuation-only" in opts.flags)
    check_inputs(opts.figures; plotting = false)
    for n in intersect(opts.figures, [3, 4])
        figure = load_figure(n)
        elapsed = @elapsed Base.invokelatest(Base.invokelatest(getfield, figure, :generate_data))
        println("Figure $n derived-data elapsed: $(round(elapsed; digits=2)) s.")
    end
end
if needs_stability
    load_stability_engine()
    config_fn = Base.invokelatest(getfield, AmpDefInc, :supplementary_stability_config)
    cfg = Base.invokelatest(config_fn; smoke = "smoke" in opts.flags)
    generation_fn = Base.invokelatest(getfield, AmpDefInc, :run_hopf_spectrum_detection)
    elapsed = @elapsed Base.invokelatest(generation_fn, cfg)
    println("Supplementary stability generation elapsed: $(round(elapsed; digits=2)) s.")
end
mkpath(DATA_ROOT)
write_checksums(DATA_ROOT, joinpath(DATA_ROOT, "SHA256SUMS"))
println("Data ready in $DATA_ROOT")
