__precompile__(false)
module PhaSingInc

using Distributions, Random, Statistics, JLD2, FFTW, SHA
import OrdinaryDiffEqTsit5, SciMLBase

const ROOT = dirname(@__DIR__)
const DATA_ROOT = abspath(get(ENV, "PHASINGINC_DATA", joinpath(ROOT, "data")))
const OUTPUT_ROOT = abspath(get(ENV, "PHASINGINC_OUTPUT", joinpath(ROOT, "figures")))
const PAPER_FONT = get(ENV, "PHASINGINC_FONT", let
    font = joinpath(get(ENV, "WINDIR", ""), "Fonts", "times.ttf")
    isfile(font) ? font : "TeX Gyre Heros Makie"
end)

export DATA_ROOT, OUTPUT_ROOT, require_file, paper_config, smoke_config,
       run_k_continuation, load_figure, check_inputs, write_checksums,
       load_supplement, load_stability_engine, check_supplement_inputs

require_file(path) = isfile(path) ? path : error("Missing input: $path. See docs/DATA.md for the Zenodo layout.")

include("Model.jl")
include("Continuation.jl")
include("FFT.jl")
include("Selection.jl")

normalize_diagnostic_mode(value) = value === :online ? value : error("Only online diagnostics are supported.")
normalize_omega_symmetry(value) = value in (:sym, :nosym) ? value : error("Unknown frequency symmetry: $value")

function paper_config(; output_dir = joinpath(DATA_ROOT, "numerical"))
    return KContinuationConfig(; output_dir,
        K_start = 1.0, K_end = 0.85, K_steps = 151,
        C1 = -0.5, C2 = -7.0, N = 5000,
        Ttr_first = 3000, Ttr = 3000, T = 600.0, dt = 0.1,
        reltol = 1.0e-8, abstol = 1.0e-8, diagnostic_mode = :online,
        omega_mean = 0.0, omega_sigma = 0.3, omega_seed = 1234,
        omega_symmetry = :sym, recenter_omega = false,
        initial_seed = 5678, continuation_noise = 0.0,
        amplitude_cluster_threshold = 0.7,
        amplitude_cluster_reltol = 1.0e-4, amplitude_cluster_abstol = 1.0e-4)
end

function smoke_config(; output_dir = joinpath(DATA_ROOT, "numerical"))
    cfg = paper_config(; output_dir)
    fields = (; (name => getfield(cfg, name) for name in fieldnames(typeof(cfg)))...)
    return KContinuationConfig(; merge(fields, (; N = 32, K_steps = 3,
        Ttr_first = 2, Ttr = 2, T = 2.0, dt = 0.05))...)
end

# Figure modules are loaded only after CLI input validation. They do not run
# simulations or create figures when included.
function load_figure(number::Integer)
    number in MAIN_FIGURES || error("Figure must be 2, 3, 4, or 6. Raw data for Figures 1 and 5 are provided on Zenodo; no plotting scripts are included.")
    name = Symbol("Figure", number)
    if !isdefined(@__MODULE__, name)
        number == 4 && load_figure(2)
        Base.include(@__MODULE__, joinpath(@__DIR__, "paper", lowercase(string(name)) * ".jl"))
    end
    return Base.invokelatest(getfield, @__MODULE__, name)
end

function check_inputs(figures; plotting::Bool)
    isempty(figures) && return nothing
    numerical = joinpath(DATA_ROOT, "numerical")
    for name in ("mean_frequency", "min_radius", "order_parameters")
        require_file(joinpath(numerical, "K_continuation_" * name * ".csv"))
    end
    header = first(eachline(joinpath(numerical, "K_continuation_min_radius.csv")))
    length(split(header, ',')) == 5001 || error("Paper figures require the N=5000 dataset. Use --smoke to plot a smoke run.")
    2 in figures && require_file(joinpath(numerical, "K_continuation_max_radius.csv"))
    for name in (2 in figures ? ["minR_meanF_cont.xlsx", "9_points.xlsx"] : String[])
        require_file(joinpath(DATA_ROOT, "experimental", name))
    end
    4 in figures && require_file(joinpath(DATA_ROOT, "experimental", "FC_FFT_data.xlsx"))
    if plotting
        3 in figures && require_file(joinpath(DATA_ROOT, "derived", "fig3_rotating_frame.jld2"))
        4 in figures && require_file(joinpath(DATA_ROOT, "derived", "fig4_fft.jld2"))
    else
        ks = union(3 in figures ? [0.970, 0.955, 0.948, 0.88] : Float64[],
                   4 in figures ? [1.0, 0.955, 0.88] : Float64[])
        for k in ks
            path = joinpath(numerical, "K_continuation_restart", "K_continuation_K_$(safe_float_token(k)).jld2")
            require_file(path)
            jldopen(path, "r") do f
                f["N"] == 5000 || error("Paper restart $path must have N=5000.")
                length(f["omega"]) == 5000 || error("Invalid frequency vector in $path.")
                length(f["final_state"]) == 10000 || error("Invalid state vector in $path.")
            end
        end
    end
    return nothing
end

function load_stability_engine()
    if !isdefined(@__MODULE__, :HopfSpectrumConfig)
        Base.include(@__MODULE__, joinpath(@__DIR__, "LockedSpectrum.jl"))
        Base.include(@__MODULE__, joinpath(@__DIR__, "HopfSpectrum.jl"))
    end
    return nothing
end

function load_supplement(id::AbstractString)
    modules = Dict("supp-kop" => :SupplementKOP,
        "supp-timeseries" => :SupplementTimeseries,
        "supp-timeseries-book2" => :SupplementTimeseriesBook2,
        "supp-stability" => :SupplementStability,
        "supp-stability-zoom" => :SupplementStability)
    haskey(modules, id) || error("Unknown supplementary figure: $id")
    name = modules[id]
    if !isdefined(@__MODULE__, name)
        id in ("supp-kop", "supp-timeseries", "supp-timeseries-book2") && load_figure(2)
        id == "supp-timeseries-book2" && load_supplement("supp-timeseries")
        Base.include(@__MODULE__, joinpath(@__DIR__, "supplementary", lowercase(string(name)) * ".jl"))
    end
    return Base.invokelatest(getfield, @__MODULE__, name)
end

function check_supplement_inputs(ids; plotting::Bool, smoke::Bool = false, skip_stability::Bool = false)
    files = Dict("supp-kop" => "KOP_Exp.xlsx", "supp-timeseries" => "min(r).xlsx",
                 "supp-timeseries-book2" => "Book2.xlsx")
    for id in ids
        if haskey(files, id)
            require_file(joinpath(DATA_ROOT, "supplementary", "experimental", files[id]))
        elseif plotting || skip_stability
            path = require_file(joinpath(DATA_ROOT, "supplementary", "numerical", "locked_spectrum.jld2"))
            jldopen(path, "r") do f
                n = f["N"]
                smoke || n == 5000 || error("The paper supplementary stability dataset must have N=5000.")
                size(f["all_eigenvalues"], 2) == 2n || error("Invalid supplementary eigenvalue matrix.")
                length(f["omega"]) == n || error("Invalid supplementary frequency vector.")
                all(f["accepted"]) || error("Supplementary stability data contain unaccepted continuation points.")
            end
        end
    end
    return nothing
end

function write_checksums(data_dir::AbstractString, path::AbstractString)
    files = sort([joinpath(dir, name) for (dir, _, names) in walkdir(data_dir)
                  for name in names if endswith(name, ".csv") || endswith(name, ".xlsx") || endswith(name, ".jld2")])
    open(path, "w") do io
        for file in files
            digest = open(SHA.sha256, file)
            println(io, bytes2hex(digest), "  ", replace(relpath(file, data_dir), '\\' => '/'))
        end
    end
    return path
end

end
