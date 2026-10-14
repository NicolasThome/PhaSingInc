import ADTypes, ForwardDiff, NonlinearSolve, LinearAlgebra
using DynamicalSystemsBase: CoupledODEs, trajectory

function simulate_locked_trajectory(K, C1, C2, omega, u0; Ttr, T, Δt, reltol, abstol)
    ds = CoupledODEs(Model.sl_harmo_hetero!, u0, (K, C1, C2, omega); diffeq=(; reltol, abstol))
    return trajectory(ds, T; Ttr, Δt)
end

amplitude_matrix(X, indices) = hcat([hypot.(X[:, 2j-1], X[:, 2j]) for j in indices]...)

Base.@kwdef struct LockedSpectrumConfig
    output_dir::String = joinpath(DATA_ROOT, "supplementary", "numerical")
    filename_base::String = "locked_spectrum"

    K::Float64 = 1.2
    C1::Float64 = -2.1
    C2::Float64 = 2.0
    N::Int = 20
    Ttr::Int = 500
    T::Float64 = 300.0
    dt::Float64 = 0.1
    analysis_window::Float64 = 150.0

    omega_mean::Float64 = 0.0
    omega_sigma::Float64 = 0.3
    omega_seed::Int = 1234
    omega_symmetry::Any = :sym
    recenter_omega::Bool = false

    initial_seed::Int = 5678
    use_slh_initial_state::Bool = true

    simulation_reltol::Float64 = 1.0e-9
    simulation_abstol::Float64 = 1.0e-9
    nonlinear_reltol::Float64 = 1.0e-10
    nonlinear_abstol::Float64 = 1.0e-10
    nonlinear_maxiters::Int = 100
    residual_tolerance::Float64 = 1.0e-8

    phase_frequency_tolerance::Float64 = 1.0e-5
    fft_min_frequency::Float64 = 0.0
    fft_ignore_first_bins::Int = 1
    amplitude_fft_threshold::Float64 = 1.0e-4
    mean_field_fft_threshold::Float64 = 1.0e-4
    fft_window::Union{Nothing, Symbol} = :hann
    fft_spectrum::Symbol = :amplitude
    fft_detrend::Bool = true
    fft_normalize::Bool = true

    neutral_overlap_threshold::Float64 = 0.5
    neutral_abs_tolerance::Float64 = 1.0e-6

    figure_size::Tuple{Int, Int} = (1100, 760)
    px_per_unit::Float64 = 2.0
    save_pdf::Bool = false
    marker_size::Float64 = 8.0
    line_width::Float64 = 1.8
    colormap::Any = :balance
    omega_colorrange::Union{Nothing, Tuple{Float64, Float64}} = nothing
end

function run_locked_simulation(cfg::LockedSpectrumConfig)
    omega = Model.sample_eigenfrequencies(Random.MersenneTwister(cfg.omega_seed),
        Normal(cfg.omega_mean, cfg.omega_sigma), cfg.N;
        omega_symmetry = cfg.omega_symmetry, recenter_nosym = cfg.recenter_omega)
    Random.seed!(cfg.initial_seed)
    u0 = Model.generate_u0(cfg.N)
    simulation_kwargs = (;
        Ttr = cfg.Ttr,
        T = cfg.T,
        Symbol("\u0394t") => cfg.dt,
        reltol = cfg.simulation_reltol,
        abstol = cfg.simulation_abstol,
    )
    X, t = simulate_locked_trajectory(
        cfg.K,
        cfg.C1,
        cfg.C2,
        omega,
        u0;
        simulation_kwargs...,
    )
    return Matrix(X), collect(t), omega, u0
end

function locked_window_indices(t::AbstractVector{<:Real}, analysis_window::Real)
    isempty(t) && error("Time vector is empty.")
    start_time = t[end] - Float64(analysis_window)
    first_idx = searchsortedfirst(t, start_time)
    first_idx = clamp(first_idx, firstindex(t), lastindex(t))
    return collect(first_idx:lastindex(t))
end

function linear_slope(t::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(t) == length(y) || error("t and y must have the same length.")
    length(t) >= 2 || return NaN
    mt = Statistics.mean(t)
    my = Statistics.mean(y)
    denom = sum((Float64(ti) - mt)^2 for ti in t)
    denom == 0.0 && return NaN
    return sum((Float64(t[i]) - mt) * (Float64(y[i]) - my) for i in eachindex(t)) / denom
end

function complex_state_matrix(Xmat::AbstractMatrix)
    N = div(size(Xmat, 2), 2)
    Z = Matrix{ComplexF64}(undef, size(Xmat, 1), N)
    for j in 1:N
        Z[:, j] .= Xmat[:, 2j - 1] .+ im .* Xmat[:, 2j]
    end
    return Z
end

function mean_field_series(Xmat::AbstractMatrix)
    Z = complex_state_matrix(Xmat)
    return vec(sum(Z; dims = 2)) ./ size(Z, 2)
end

function phase_diagnostics(cfg::LockedSpectrumConfig, Xmat::AbstractMatrix, t, idx)
    Xw = @view Xmat[idx, :]
    tw = t[idx]
    N = div(size(Xw, 2), 2)
    mean_frequency = Vector{Float64}(undef, N)
    phases = Matrix{Float64}(undef, length(idx), N)

    for j in 1:N
        z = Xw[:, 2j - 1] .+ im .* Xw[:, 2j]
        phases[:, j] .= Model.unwrap(angle.(z))
        mean_frequency[j] = linear_slope(tw, @view phases[:, j])
    end

    mean_field = mean_field_series(Xw)
    mean_field_phase = Model.unwrap(angle.(mean_field))
    Omega = linear_slope(tw, mean_field_phase)
    frequency_spread = maximum(mean_frequency) - minimum(mean_frequency)
    frequency_std = Statistics.std(mean_frequency)
    frequency_locked = frequency_spread <= cfg.phase_frequency_tolerance

    relative_drift = 0.0
    mean_phase = Model.unwrap(angle.(mean_field))
    for j in 1:N
        rel = phases[:, j] .- mean_phase
        drel = diff(rel)
        Model.unwrap_diff!(drel)
        relative_drift = max(relative_drift, abs(sum(drel)))
    end

    return (;
        mean_frequency,
        frequency_spread,
        frequency_std,
        mean_field_phase,
        Omega,
        relative_phase_drift = relative_drift,
        frequency_locked,
    )
end

function spectral_peak_indices(freqs, min_frequency::Real, ignore_first_bins::Integer)
    candidates = collect(eachindex(freqs))
    nignore = min(max(Int(ignore_first_bins), 0) + 1, length(candidates))
    candidates = nignore >= length(candidates) ? Int[] : candidates[(nignore + 1):end]
    return [i for i in candidates if freqs[i] >= Float64(min_frequency)]
end

function max_spectral_peak(freqs, spectra; min_frequency::Real, ignore_first_bins::Integer)
    candidates = spectral_peak_indices(freqs, min_frequency, ignore_first_bins)
    isempty(candidates) && return NaN
    return maximum(spectra[:, candidates])
end

function stationary_fft_diagnostics(cfg::LockedSpectrumConfig, Xmat::AbstractMatrix, t, idx)
    Xw = Xmat[idx, :]
    tw = t[idx]
    indices = collect(1:div(size(Xw, 2), 2))
    amps = amplitude_matrix(Xw, indices)
    amp_fft = fft_matrix(
        amps,
        tw;
        detrend = cfg.fft_detrend,
        window = cfg.fft_window,
        spectrum = cfg.fft_spectrum,
        onesided = true,
        normalize = cfg.fft_normalize,
    )

    mean_field_amp = abs.(mean_field_series(Xw))
    mean_fft = fft_matrix(
        reshape(mean_field_amp, :, 1),
        tw;
        detrend = cfg.fft_detrend,
        window = cfg.fft_window,
        spectrum = cfg.fft_spectrum,
        onesided = true,
        normalize = cfg.fft_normalize,
    )

    amp_peak = max_spectral_peak(
        amp_fft.freqs,
        amp_fft.spectra;
        min_frequency = cfg.fft_min_frequency,
        ignore_first_bins = cfg.fft_ignore_first_bins,
    )
    mean_peak = max_spectral_peak(
        mean_fft.freqs,
        mean_fft.spectra;
        min_frequency = cfg.fft_min_frequency,
        ignore_first_bins = cfg.fft_ignore_first_bins,
    )

    return (;
        freqs = amp_fft.freqs,
        amplitude_spectra = amp_fft.spectra,
        amplitude_meta = amp_fft.meta,
        mean_field_freqs = mean_fft.freqs,
        mean_field_spectrum = vec(mean_fft.spectra[1, :]),
        mean_field_meta = mean_fft.meta,
        max_amplitude_fft_peak = amp_peak,
        max_mean_field_fft_peak = mean_peak,
        stationary = amp_peak <= cfg.amplitude_fft_threshold && mean_peak <= cfg.mean_field_fft_threshold,
    )
end

function rotating_initial_guess(Xmat::AbstractMatrix, t, idx, Omega::Real)
    Xw = @view Xmat[idx, :]
    tw = t[idx]
    N = div(size(Xw, 2), 2)
    z0 = Vector{ComplexF64}(undef, N)
    for j in 1:N
        z = Xw[:, 2j - 1] .+ im .* Xw[:, 2j]
        z_rot = z .* exp.(-im .* Float64(Omega) .* tw)
        z0[j] = Statistics.mean(z_rot)
    end
    mean_z = Statistics.mean(z0)
    if abs(mean_z) > 0.0
        z0 .*= exp(-im * angle(mean_z))
    end

    u0 = Vector{Float64}(undef, 2N + 1)
    for j in 1:N
        u0[2j - 1] = real(z0[j])
        u0[2j] = imag(z0[j])
    end
    u0[end] = Float64(Omega)
    return u0
end

function rotating_vector_field(u::AbstractVector, K::Real, C1::Real, C2::Real, omega::AbstractVector, Omega::Real)
    T = promote_type(eltype(u), typeof(K), typeof(C1), typeof(C2), eltype(omega), typeof(Omega))
    N = length(omega)
    du = Vector{T}(undef, 2N)

    sx = zero(T)
    sy = zero(T)
    for j in 1:N
        sx += u[2j - 1]
        sy += u[2j]
    end
    mx = sx / T(N)
    my = sy / T(N)

    for j in 1:N
        i = 2j - 1
        x = u[i]
        y = u[i + 1]
        dx = mx - x
        dy = my - y
        r2 = x * x + y * y
        delta = omega[j] - Omega
        du[i] = x - delta * y - r2 * x + C2 * r2 * y + K * (dx - C1 * dy)
        du[i + 1] = y + delta * x - r2 * y - C2 * r2 * x + K * (C1 * dx + dy)
    end
    return du
end

function locked_equilibrium_residual!(F, q, p)
    N = length(p.omega)
    u = @view q[1:2N]
    Omega = q[end]
    du = rotating_vector_field(u, p.K, p.C1, p.C2, p.omega, Omega)
    F[1:2N] .= du
    F[end] = sum(u[2:2:end]) / N
    return nothing
end

function solve_locked_equilibrium(cfg::LockedSpectrumConfig, omega, guess)
    p = (K = cfg.K, C1 = cfg.C1, C2 = cfg.C2, omega = Float64.(omega))
    prob = NonlinearSolve.NonlinearProblem(locked_equilibrium_residual!, guess, p)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.NewtonRaphson(; autodiff = ADTypes.AutoForwardDiff());
        abstol = cfg.nonlinear_abstol,
        reltol = cfg.nonlinear_reltol,
        maxiters = cfg.nonlinear_maxiters,
    )

    residual = similar(Float64.(sol.u))
    locked_equilibrium_residual!(residual, sol.u, p)
    return (;
        sol,
        q = Float64.(sol.u),
        u = Float64.(sol.u[1:(end - 1)]),
        Omega = Float64(sol.u[end]),
        residual = residual,
        residual_norm = LinearAlgebra.norm(residual, Inf),
        retcode = string(sol.retcode),
    )
end

function mode_weights(eigenvectors::AbstractMatrix, N::Integer)
    nmodes = size(eigenvectors, 2)
    weights = Matrix{Float64}(undef, nmodes, N)
    for j in 1:nmodes
        denom = 0.0
        for k in 1:N
            value = abs2(eigenvectors[2k - 1, j]) + abs2(eigenvectors[2k, j])
            weights[j, k] = value
            denom += value
        end
        denom == 0.0 || (weights[j, :] ./= denom)
    end
    return weights
end

function spectrum_diagnostics(cfg::LockedSpectrumConfig, u, Omega, omega)
    N = length(omega)
    vf = x -> rotating_vector_field(x, cfg.K, cfg.C1, cfg.C2, omega, Omega)
    J = ForwardDiff.jacobian(vf, u)
    eig = LinearAlgebra.eigen(J)
    eigenvalues = eig.values
    eigenvectors = eig.vectors

    tangent = Vector{Float64}(undef, 2N)
    for j in 1:N
        tangent[2j - 1] = -u[2j]
        tangent[2j] = u[2j - 1]
    end
    tangent_norm = LinearAlgebra.norm(tangent)
    overlaps = Vector{Float64}(undef, length(eigenvalues))
    for j in eachindex(eigenvalues)
        v = @view eigenvectors[:, j]
        denom = LinearAlgebra.norm(v) * tangent_norm
        overlaps[j] = denom == 0.0 ? 0.0 : abs(LinearAlgebra.dot(v, tangent)) / denom
    end

    neutral_index = argmax(overlaps)
    nonneutral = [j for j in eachindex(eigenvalues) if j != neutral_index]
    leading_index = nonneutral[argmax(real.(eigenvalues[nonneutral]))]

    weights = mode_weights(eigenvectors, N)
    omega_eff = vec(weights * Float64.(omega))
    participation_ratio = [1.0 / sum(weights[j, :] .^ 2) for j in axes(weights, 1)]
    dominant_index = [argmax(weights[j, :]) for j in axes(weights, 1)]
    dominant_weight = [weights[j, dominant_index[j]] for j in axes(weights, 1)]

    neutral_candidates = findall(j ->
        overlaps[j] >= cfg.neutral_overlap_threshold &&
        abs(eigenvalues[j]) <= cfg.neutral_abs_tolerance,
        eachindex(eigenvalues),
    )

    return (;
        J,
        eigenvalues,
        eigenvectors,
        neutral_index,
        leading_index,
        phase_tangent = tangent,
        phase_overlaps = overlaps,
        neutral_candidates,
        weights,
        omega_eff,
        participation_ratio,
        dominant_index,
        dominant_weight,
    )
end

function locked_amplitudes(u)
    N = div(length(u), 2)
    return [hypot(u[2j - 1], u[2j]) for j in 1:N]
end
