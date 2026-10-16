using Test, PhaSingInc, Distributions, Random, Statistics, JLD2
import SciMLBase, OrdinaryDiffEqTsit5

@testset "Model and seeded frequency sampling" begin
    u = [0.4, -0.3]
    du = zeros(2)
    PhaSingInc.Model.sl_harmo_hetero!(du, u, (2.0, -0.5, -7.0, [0.2]), 0.0)
    r2 = sum(abs2, u)
    expected = [(1-r2)*u[1] - 0.2*u[2] - 7*r2*u[2],
                (1-r2)*u[2] + 0.2*u[1] + 7*r2*u[1]]
    @test du ≈ expected
    for n in (31, 32)
        omega = PhaSingInc.Model.sample_eigenfrequencies(MersenneTwister(1234), Normal(0, 0.3), n; omega_symmetry=:sym)
        @test length(omega) == n
        @test sort(omega) ≈ sort(-omega)
        @test abs(mean(omega)) < 1e-15
        @test omega == PhaSingInc.Model.sample_eigenfrequencies(MersenneTwister(1234), Normal(0, 0.3), n; omega_symmetry=:sym)
    end
end

@testset "Streaming diagnostics agree with saved trajectories" begin
    cfg = smoke_config()
    omega = PhaSingInc.Model.sample_eigenfrequencies(MersenneTwister(1234), Normal(0, 0.3), 32; omega_symmetry=:sym)
    u0 = rand(MersenneTwister(17), 64) .- 0.5
    online = PhaSingInc.run_k_step_online_core(cfg, 0.95, omega, u0, 2)
    times = PhaSingInc.sample_times(2, cfg.T, cfg.dt)
    prob = SciMLBase.ODEProblem(PhaSingInc.Model.sl_harmo_hetero!, u0, (0.0, 4.0), (0.95, cfg.C1, cfg.C2, omega))
    sol = SciMLBase.solve(prob, OrdinaryDiffEqTsit5.Tsit5(); saveat=times,
        reltol=cfg.reltol, abstol=cfg.abstol)
    z = Array(sol)[1:2:end, :] .+ im .* Array(sol)[2:2:end, :]
    @test online.min_radius ≈ vec(minimum(abs.(z); dims=2)) atol=2e-6
    @test online.max_radius ≈ vec(maximum(abs.(z); dims=2)) atol=2e-6
    @test online.R_mean ≈ mean(abs.(vec(mean(z ./ abs.(z); dims=1)))) atol=2e-6
    phase_steps = angle.(z[:, 2:end] ./ z[:, 1:end-1])
    @test online.mean_frequency ≈ vec(sum(phase_steps; dims=2)) ./ cfg.T atol=2e-6
    @test online.final_state ≈ sol.u[end] atol=2e-6
end

@testset "FFT normalization and frequency units" begin
    t = collect(0:0.01:9.99)
    signal = 3 .+ 2 .* sin.(2pi .* 1.0 .* t)
    spectrum = PhaSingInc.fft_matrix(reshape(signal, :, 1), t;
        detrend=true, window=nothing, spectrum=:amplitude, onesided=true, normalize=false)
    peak = argmax(vec(spectrum.spectra))
    @test spectrum.freqs[peak] ≈ 1.0 atol=1e-12
    @test spectrum.spectra[1, peak] ≈ 2.0 atol=1e-12
    @test spectrum.spectra[1, 1] < 1e-12
end

@testset "Continuation output and restart provenance" begin
    mktempdir() do dir
        result = run_k_continuation(smoke_config(; output_dir=dir))
        @test result.K_list == [1.0, 0.925, 0.85]
        @test length(result.restart_files) == 3
        @test all(isfile, [result.mean_csv, result.min_r_csv, result.max_r_csv, result.order_parameter_csv])
        @test size(result.mean_frequency) == (3, 32)
        @test all(isfinite, result.mean_frequency)
        @test all(result.min_radius .<= result.max_radius)
        @test all(result.cluster_order_parameters.N_high .+ result.cluster_order_parameters.N_low .== 32)
        initial = JLD2.load(first(result.restart_files))
        second = JLD2.load(result.restart_files[2])
        @test initial["N"] == 32
        @test second["u0"] == initial["final_state"]
        @test initial["omega_seed"] == 1234
        @test initial["initial_seed"] == 5678
        @test length(initial["omega"]) == 32
        @test length(initial["final_state"]) == 64
        manifest = write_checksums(dir, joinpath(dir, "SHA256SUMS"))
        @test length(readlines(manifest)) == 7
        @test all(length(first(split(line))) == 64 for line in readlines(manifest))
    end
end

@testset "Paper settings" begin
    cfg = paper_config()
    @test cfg.N == 5000
    @test cfg.K_steps == 151
    @test cfg.Ttr_first == cfg.Ttr == 3000
    @test cfg.T == 600
    @test cfg.dt == 0.1
    @test cfg.continuation_noise == 0
end

@testset "Main and supplementary selection" begin
    @test PhaSingInc.parse_figure_selection("2,supp-kop,6") == ([2, 6], ["supp-kop"])
    main, supp = PhaSingInc.parse_figure_selection("all")
    @test main == [2, 3, 4, 6]
    @test supp == PhaSingInc.SUPPLEMENT_IDS
    @test PhaSingInc.parse_figure_selection("supplementary")[1] == Int[]
    @test PhaSingInc.parse_figure_selection("6") == ([6], String[])
    @test_throws ErrorException PhaSingInc.parse_figure_selection("1")
    @test_throws ErrorException PhaSingInc.parse_figure_selection("5")
    @test_throws ErrorException load_figure(5)
end

load_stability_engine()

@testset "Rotating equilibrium and neutral phase mode" begin
    n = 4
    state = repeat([1.0, 0.0], n)
    omega = zeros(n)
    cfg = PhaSingInc.LockedSpectrumConfig(; K=1.0, C1=-0.5, C2=-7.0, N=n)
    field = PhaSingInc.rotating_vector_field(state, cfg.K, cfg.C1, cfg.C2, omega, 7.0)
    @test field == zeros(2n)
    equilibrium = PhaSingInc.solve_locked_equilibrium(cfg, omega, vcat(state, 7.0))
    @test equilibrium.residual_norm < 1e-10
    @test equilibrium.Omega ≈ 7.0
    spectrum = PhaSingInc.spectrum_diagnostics(cfg, state, 7.0, omega)
    @test spectrum.J * spectrum.phase_tangent ≈ zeros(2n) atol=1e-12
    @test abs(spectrum.eigenvalues[spectrum.neutral_index]) < 1e-10
    @test spectrum.leading_index != spectrum.neutral_index
    @test vec(sum(spectrum.weights; dims=2)) ≈ ones(2n)
    @test all(spectrum.participation_ratio .>= 1 - 1e-12)
    @test all(spectrum.participation_ratio .<= n + 1e-12)
end

@testset "Supplementary stability settings and spectra filtering" begin
    cfg = PhaSingInc.supplementary_stability_config()
    @test cfg.N == 5000
    @test cfg.K_steps == 51
    @test cfg.K_start == 1.0 && cfg.K_end == 0.95
    @test cfg.Ttr_first == 2000 && cfg.T_first == 1000
    @test cfg.analysis_window == 100
    @test cfg.continuation_method == :fixed_K
    @test !cfg.stop_at_first_hopf
    module_ref = load_supplement("supp-stability")
    filter_fn = Base.invokelatest(getfield, module_ref, :finite_spectrum_data)
    x, y, k = Base.invokelatest(filter_fn, [1+2im NaN+0im; 3+4im 5+6im], [1.0, 0.99], [true, false])
    @test x == [1.0] && y == [2.0] && k == [1.0]
end
