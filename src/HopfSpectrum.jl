Base.@kwdef struct HopfSpectrumConfig
    output_dir::String = joinpath(DATA_ROOT, "supplementary", "numerical")
    filename_base::String = "hopf_spectrum"

    K_values::Union{Nothing, Vector{Float64}} = nothing
    K_start::Float64 = 1.2
    K_end::Float64 = 0.9
    K_steps::Int = 5

    C1::Float64 = -2.1
    C2::Float64 = 2.0
    N::Int = 20
    Ttr_first::Int = 500
    T_first::Float64 = 300.0
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

    hopf_imag_tolerance::Float64 = 1.0e-6
    real_tolerance::Float64 = 1.0e-6
    stability_tolerance::Float64 = 1.0e-7
    stay_on_stable_branch::Bool = true
    stop_at_first_hopf::Bool = true
    steady_bifurcation_max_bisections::Int = 10
    steady_bifurcation_K_tolerance::Float64 = 1.0e-5
    localized_pr_threshold::Float64 = 3.0
    collective_pr_fraction::Float64 = 0.25

    neutral_overlap_threshold::Float64 = 0.5
    neutral_abs_tolerance::Float64 = 1.0e-6
    stop_on_failure::Bool = true
    enable_recovery::Bool = true
    continuation_method::Symbol = :hybrid
    recovery_max_halvings::Int = 6
    branch_switch_epsilons::Vector{Float64} = [1.0e-4, 1.0e-3, 1.0e-2]
    enable_pseudo_arclength::Bool = true
    pseudo_arclength_epsilons::Vector{Float64} = [1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3]
    pseudo_arclength_dK_factors::Vector{Float64} = [0.25, 0.5, 1.0, 2.0]
    pseudo_arclength_max_steps::Int = 12
    pseudo_arclength_maxiters_factor::Int = 8
    pseudo_arclength_tolerance::Float64 = 1.0e-8
    enable_simulation_recovery::Bool = true
    simulation_recovery_Ttr::Union{Nothing, Int} = nothing
    simulation_recovery_T::Union{Nothing, Float64} = nothing
    simulation_recovery_analysis_window::Union{Nothing, Float64} = nothing
    recovery_retry_maxiters_factor::Int = 5
    store_recovery_points::Bool = true

    figure_size::Tuple{Int, Int} = (1100, 760)
    px_per_unit::Float64 = 2.0
    save_pdf::Bool = false
    marker_size::Float64 = 8.0
    line_width::Float64 = 1.8
    colormap::Any = :balance
    omega_colorrange::Union{Nothing, Tuple{Float64, Float64}} = nothing
end

function hopf_k_values(cfg::HopfSpectrumConfig)
    cfg.K_values !== nothing && return Float64.(cfg.K_values)
    return collect(range(cfg.K_start, cfg.K_end; length = cfg.K_steps))
end

function hopf_continuation_method(cfg::HopfSpectrumConfig)
    cfg.continuation_method in (:fixed_K, :hybrid) && return cfg.continuation_method
    error("continuation_method must be :fixed_K or :hybrid. Got $(repr(cfg.continuation_method)).")
end

function pseudo_arclength_enabled(cfg::HopfSpectrumConfig)
    return hopf_continuation_method(cfg) == :hybrid && cfg.enable_pseudo_arclength
end

function locked_config_for_hopf(
        cfg::HopfSpectrumConfig,
        K::Real;
        Ttr::Integer = cfg.Ttr_first,
        T::Real = cfg.T_first,
        analysis_window::Real = cfg.analysis_window,
        nonlinear_maxiters::Integer = cfg.nonlinear_maxiters)
    return LockedSpectrumConfig(
        output_dir = cfg.output_dir,
        filename_base = cfg.filename_base,
        K = Float64(K),
        C1 = cfg.C1,
        C2 = cfg.C2,
        N = cfg.N,
        Ttr = Int(Ttr),
        T = Float64(T),
        dt = cfg.dt,
        analysis_window = Float64(analysis_window),
        omega_mean = cfg.omega_mean,
        omega_sigma = cfg.omega_sigma,
        omega_seed = cfg.omega_seed,
        omega_symmetry = cfg.omega_symmetry,
        recenter_omega = cfg.recenter_omega,
        initial_seed = cfg.initial_seed,
        use_slh_initial_state = cfg.use_slh_initial_state,
        simulation_reltol = cfg.simulation_reltol,
        simulation_abstol = cfg.simulation_abstol,
        nonlinear_reltol = cfg.nonlinear_reltol,
        nonlinear_abstol = cfg.nonlinear_abstol,
        nonlinear_maxiters = Int(nonlinear_maxiters),
        residual_tolerance = cfg.residual_tolerance,
        phase_frequency_tolerance = cfg.phase_frequency_tolerance,
        fft_min_frequency = cfg.fft_min_frequency,
        fft_ignore_first_bins = cfg.fft_ignore_first_bins,
        amplitude_fft_threshold = cfg.amplitude_fft_threshold,
        mean_field_fft_threshold = cfg.mean_field_fft_threshold,
        fft_window = cfg.fft_window,
        fft_spectrum = cfg.fft_spectrum,
        fft_detrend = cfg.fft_detrend,
        fft_normalize = cfg.fft_normalize,
        neutral_overlap_threshold = cfg.neutral_overlap_threshold,
        neutral_abs_tolerance = cfg.neutral_abs_tolerance,
        figure_size = cfg.figure_size,
        px_per_unit = cfg.px_per_unit,
        save_pdf = cfg.save_pdf,
        marker_size = cfg.marker_size,
        line_width = cfg.line_width,
        colormap = cfg.colormap,
        omega_colorrange = cfg.omega_colorrange,
    )
end

function classify_hopf_mode(cfg::HopfSpectrumConfig, pr::Real)
    !isfinite(pr) && return "failed"
    pr <= cfg.localized_pr_threshold && return "localized"
    pr >= cfg.collective_pr_fraction * cfg.N && return "collective"
    return "mixed"
end

function first_hopf_guess_and_diagnostics(cfg::HopfSpectrumConfig, locked_cfg::LockedSpectrumConfig)
    Xmat, t, omega, u0 = run_locked_simulation(locked_cfg)
    idx = locked_window_indices(t, cfg.analysis_window)
    phase = phase_diagnostics(locked_cfg, Xmat, t, idx)
    fftdiag = stationary_fft_diagnostics(locked_cfg, Xmat, t, idx)
    guess = rotating_initial_guess(Xmat, t, idx, phase.Omega)
    return (; Xmat, t, omega, u0, idx, phase, fftdiag, guess)
end

function initialize_hopf_storage(N::Integer)
    return (;
        K_path = Float64[],
        residuals = Float64[],
        Omegas = Float64[],
        accepted = Bool[],
        retcodes = String[],
        leading_eigenvalues = ComplexF64[],
        neutral_values = ComplexF64[],
        neutral_indices = Int[],
        leading_indices = Int[],
        neutral_candidate_counts = Int[],
        participation_ratio = Float64[],
        omega_eff = Float64[],
        dominant_index = Int[],
        dominant_weight = Float64[],
        classifications = String[],
        all_eigenvalues = Vector{ComplexF64}[],
        leading_weights = Vector{Float64}[],
        locked_amplitudes = Vector{Float64}[],
        final_equilibria = Vector{Float64}[],
        recovery_method = String[],
        parent_K = Float64[],
        target_K = Float64[],
        is_inserted_recovery_point = Bool[],
        is_requested_K = Bool[],
        recovery_attempts = String[],
    )
end

function record_hopf_step!(store, K::Real, eq, spectrum, amplitudes, cfg::HopfSpectrumConfig;
        recovery_method::AbstractString = "direct",
        parent_K::Real = NaN,
        target_K::Real = K,
        is_inserted_recovery_point::Bool = false,
        is_requested_K::Bool = true,
        recovery_attempts::AbstractString = "direct")
    li = spectrum.leading_index
    ni = spectrum.neutral_index
    push!(store.K_path, Float64(K))
    push!(store.residuals, eq.residual_norm)
    push!(store.Omegas, eq.Omega)
    push!(store.accepted, eq.residual_norm <= cfg.residual_tolerance)
    push!(store.retcodes, eq.retcode)
    push!(store.leading_eigenvalues, spectrum.eigenvalues[li])
    push!(store.neutral_values, spectrum.eigenvalues[ni])
    push!(store.neutral_indices, ni)
    push!(store.leading_indices, li)
    push!(store.neutral_candidate_counts, length(spectrum.neutral_candidates))
    push!(store.participation_ratio, spectrum.participation_ratio[li])
    push!(store.omega_eff, spectrum.omega_eff[li])
    push!(store.dominant_index, spectrum.dominant_index[li])
    push!(store.dominant_weight, spectrum.dominant_weight[li])
    push!(store.classifications, classify_hopf_mode(cfg, spectrum.participation_ratio[li]))
    push!(store.all_eigenvalues, ComplexF64.(spectrum.eigenvalues))
    push!(store.leading_weights, Float64.(spectrum.weights[li, :]))
    push!(store.locked_amplitudes, Float64.(amplitudes))
    push!(store.final_equilibria, Float64.(eq.q))
    push!(store.recovery_method, String(recovery_method))
    push!(store.parent_K, Float64(parent_K))
    push!(store.target_K, Float64(target_K))
    push!(store.is_inserted_recovery_point, is_inserted_recovery_point)
    push!(store.is_requested_K, is_requested_K)
    push!(store.recovery_attempts, String(recovery_attempts))
    return store
end

function attempt_hopf_solve(cfg::HopfSpectrumConfig, K::Real, omega, guess; maxiters::Integer = cfg.nonlinear_maxiters)
    locked_cfg = locked_config_for_hopf(cfg, K; nonlinear_maxiters = maxiters)
    try
        eq = solve_locked_equilibrium(locked_cfg, omega, guess)
        spectrum = spectrum_diagnostics(locked_cfg, eq.u, eq.Omega, omega)
        amplitudes = locked_amplitudes(eq.u)
        return (;
            success = eq.residual_norm <= cfg.residual_tolerance,
            eq,
            spectrum,
            amplitudes,
            residual_norm = eq.residual_norm,
            retcode = eq.retcode,
            error = "",
        )
    catch err
        return (;
            success = false,
            eq = nothing,
            spectrum = nothing,
            amplitudes = Float64[],
            residual_norm = Inf,
            retcode = "Exception",
            error = sprint(showerror, err),
        )
    end
end

function pseudo_arclength_residual!(F, y, p)
    nq = length(y) - 1
    q = @view y[1:nq]
    K = y[end]
    locked_equilibrium_residual!(@view(F[1:nq]), q, (K = K, C1 = p.C1, C2 = p.C2, omega = p.omega))
    arc = zero(eltype(y))
    for j in eachindex(y)
        arc += (y[j] - p.y_pred[j]) * p.tangent[j]
    end
    F[end] = arc
    return nothing
end

function attempt_pseudo_arclength_solve(cfg::HopfSpectrumConfig, omega, y_guess, y_pred, tangent; maxiters::Integer = cfg.nonlinear_maxiters)
    p = (C1 = cfg.C1, C2 = cfg.C2, omega = Float64.(omega), y_pred = Float64.(y_pred), tangent = Float64.(tangent))
    try
        prob = NonlinearSolve.NonlinearProblem(pseudo_arclength_residual!, Float64.(y_guess), p)
        sol = NonlinearSolve.solve(
            prob,
            NonlinearSolve.NewtonRaphson(; autodiff = ADTypes.AutoForwardDiff());
            abstol = cfg.nonlinear_abstol,
            reltol = cfg.nonlinear_reltol,
            maxiters = Int(maxiters),
        )
        y = Float64.(sol.u)
        q = y[1:(end - 1)]
        K = y[end]
        residual = similar(q)
        locked_equilibrium_residual!(residual, q, (K = K, C1 = cfg.C1, C2 = cfg.C2, omega = Float64.(omega)))
        arclength_residual = LinearAlgebra.dot(y .- Float64.(y_pred), Float64.(tangent))
        locked_cfg = locked_config_for_hopf(cfg, K; nonlinear_maxiters = maxiters)
        eq = (;
            sol,
            q,
            u = Float64.(q[1:(end - 1)]),
            Omega = Float64(q[end]),
            residual,
            residual_norm = LinearAlgebra.norm(residual, Inf),
            retcode = string(sol.retcode),
            K = Float64(K),
            arclength_residual,
        )
        spectrum = spectrum_diagnostics(locked_cfg, eq.u, eq.Omega, omega)
        amplitudes = locked_amplitudes(eq.u)
        success = eq.residual_norm <= cfg.residual_tolerance &&
                  abs(arclength_residual) <= cfg.pseudo_arclength_tolerance
        return (;
            success,
            K = Float64(K),
            eq,
            spectrum,
            amplitudes,
            residual_norm = eq.residual_norm,
            arclength_residual,
            retcode = eq.retcode,
            error = "",
        )
    catch err
        return (;
            success = false,
            K = NaN,
            eq = nothing,
            spectrum = nothing,
            amplitudes = Float64[],
            residual_norm = Inf,
            arclength_residual = Inf,
            retcode = "Exception",
            error = sprint(showerror, err),
        )
    end
end

function hopf_step_result(K, attempt; method, parent_K, target_K, is_inserted, is_requested, attempts)
    return (;
        K = Float64(K),
        eq = attempt.eq,
        spectrum = attempt.spectrum,
        amplitudes = attempt.amplitudes,
        recovery_method = String(method),
        parent_K = Float64(parent_K),
        target_K = Float64(target_K),
        is_inserted_recovery_point = Bool(is_inserted),
        is_requested_K = Bool(is_requested),
        recovery_attempts = join(attempts, " | "),
    )
end

function append_hopf_step!(store, step, cfg::HopfSpectrumConfig)
    return record_hopf_step!(
        store,
        step.K,
        step.eq,
        step.spectrum,
        step.amplitudes,
        cfg;
        recovery_method = step.recovery_method,
        parent_K = step.parent_K,
        target_K = step.target_K,
        is_inserted_recovery_point = step.is_inserted_recovery_point,
        is_requested_K = step.is_requested_K,
        recovery_attempts = step.recovery_attempts,
    )
end

function rephase_equilibrium_guess(q::AbstractVector{<:Real})
    out = Float64.(q)
    N = div(length(out) - 1, 2)
    mean_z = sum(out[2j - 1] + im * out[2j] for j in 1:N) / N
    abs(mean_z) == 0.0 && return out
    rot = exp(-im * angle(mean_z))
    for j in 1:N
        z = (out[2j - 1] + im * out[2j]) * rot
        out[2j - 1] = real(z)
        out[2j] = imag(z)
    end
    return out
end

function branch_switch_guesses(cfg::HopfSpectrumConfig, previous_step)
    value = previous_step.spectrum.eigenvalues[previous_step.spectrum.leading_index]
    abs(imag(value)) <= cfg.hopf_imag_tolerance || return Vector{Float64}[]

    v = real.(previous_step.spectrum.eigenvectors[:, previous_step.spectrum.leading_index])
    nv = LinearAlgebra.norm(v)
    nv == 0.0 && return Vector{Float64}[]

    base = Float64.(previous_step.eq.q)
    guesses = Vector{Float64}[]
    N2 = length(base) - 1
    for eps in cfg.branch_switch_epsilons
        for sign in (-1.0, 1.0)
            q = copy(base)
            q[1:N2] .+= sign * eps .* v ./ nv
            q[end] = previous_step.eq.Omega
            push!(guesses, rephase_equilibrium_guess(q))
        end
    end
    return guesses
end

function simulation_recovery_config_values(cfg::HopfSpectrumConfig)
    Ttr = cfg.simulation_recovery_Ttr === nothing ? cfg.Ttr_first : cfg.simulation_recovery_Ttr
    T = cfg.simulation_recovery_T === nothing ? cfg.T_first : cfg.simulation_recovery_T
    window = cfg.simulation_recovery_analysis_window === nothing ? cfg.analysis_window : cfg.simulation_recovery_analysis_window
    return (Ttr = Int(Ttr), T = Float64(T), analysis_window = Float64(window))
end

function simulation_recovery_attempt(cfg::HopfSpectrumConfig, K::Real, omega)
    values = simulation_recovery_config_values(cfg)
    locked_cfg = locked_config_for_hopf(
        cfg,
        K;
        Ttr = values.Ttr,
        T = values.T,
        analysis_window = values.analysis_window,
        nonlinear_maxiters = cfg.nonlinear_maxiters * cfg.recovery_retry_maxiters_factor,
    )
    try
        Xmat, t, omega_recovered, _ = run_locked_simulation(locked_cfg)
        maximum(abs.(omega_recovered .- omega)) <= 1.0e-12 || error("Simulation recovery generated a different omega sample.")
        idx = locked_window_indices(t, locked_cfg.analysis_window)
        fftdiag = stationary_fft_diagnostics(locked_cfg, Xmat, t, idx)
        phase = phase_diagnostics(locked_cfg, Xmat, t, idx)
        guess = rotating_initial_guess(Xmat, t, idx, phase.Omega)
        attempt = attempt_hopf_solve(
            cfg,
            K,
            omega,
            guess;
            maxiters = cfg.nonlinear_maxiters * cfg.recovery_retry_maxiters_factor,
        )
        success = attempt.success && fftdiag.stationary
        return (; success, attempt, fft_stationary = fftdiag.stationary, error = "")
    catch err
        return (;
            success = false,
            attempt = nothing,
            fft_stationary = false,
            error = sprint(showerror, err),
        )
    end
end

function leading_real_part(step)
    value = step.spectrum.eigenvalues[step.spectrum.leading_index]
    return real(value)
end

function leading_imag_abs(step)
    value = step.spectrum.eigenvalues[step.spectrum.leading_index]
    return abs(imag(value))
end

function leading_value(step)
    return step.spectrum.eigenvalues[step.spectrum.leading_index]
end

function is_real_leading_mode(cfg::HopfSpectrumConfig, step)
    return leading_imag_abs(step) <= cfg.hopf_imag_tolerance
end

function is_stable_step(cfg::HopfSpectrumConfig, step)
    return leading_real_part(step) <= cfg.stability_tolerance
end

function has_real_stability_crossing(cfg::HopfSpectrumConfig, left_step, right_step)
    is_real_leading_mode(cfg, left_step) && is_real_leading_mode(cfg, right_step) || return false
    return is_stable_step(cfg, left_step) && !is_stable_step(cfg, right_step)
end

function has_hopf_crossing(cfg::HopfSpectrumConfig, left_step, right_step)
    return has_hopf_crossing(leading_value(left_step), leading_value(right_step), cfg)
end

function real_branch_obstruction(cfg::HopfSpectrumConfig, previous_step, recovered_step)
    leading_imag_abs(previous_step) <= cfg.hopf_imag_tolerance || return false
    leading_imag_abs(recovered_step) <= cfg.hopf_imag_tolerance || return false
    return leading_real_part(previous_step) <= cfg.real_tolerance &&
           leading_real_part(recovered_step) > cfg.real_tolerance
end

function localize_real_stability_crossing(cfg::HopfSpectrumConfig, stable_step, unstable_step, omega, attempts)
    retry_maxiters = cfg.nonlinear_maxiters * cfg.recovery_retry_maxiters_factor
    left = stable_step
    right = unstable_step
    inserted_stable_steps = Any[]

    for level in 1:cfg.steady_bifurcation_max_bisections
        abs(Float64(right.K) - Float64(left.K)) <= cfg.steady_bifurcation_K_tolerance && break
        K_mid = 0.5 * (Float64(left.K) + Float64(right.K))
        a = attempt_hopf_solve(cfg, K_mid, omega, left.eq.q; maxiters = retry_maxiters)
        push!(attempts, "steady_midpoint(level=$(level), K=$(K_mid), retcode=$(a.retcode), residual=$(a.residual_norm), error=$(a.error))")
        if !a.success
            push!(attempts, "steady_midpoint_failed_using_current_bracket(K_left=$(left.K), K_right=$(right.K))")
            return (; success = true, stable_step = left, unstable_step = right, inserted_stable_steps, attempts)
        end

        probe_step = hopf_step_result(
            K_mid,
            a;
            method = "steady_midpoint_probe",
            parent_K = left.K,
            target_K = unstable_step.K,
            is_inserted = true,
            is_requested = false,
            attempts,
        )
        mid_step = hopf_step_result(
            K_mid,
            a;
            method = is_stable_step(cfg, probe_step) ? "steady_midpoint_stable" : "steady_midpoint_unstable",
            parent_K = left.K,
            target_K = unstable_step.K,
            is_inserted = true,
            is_requested = false,
            attempts,
        )

        if is_stable_step(cfg, mid_step)
            push!(inserted_stable_steps, mid_step)
            left = mid_step
        else
            right = mid_step
        end
    end

    return (; success = true, stable_step = left, unstable_step = right, inserted_stable_steps, attempts)
end

function normalized_critical_vector(step)
    v = real.(step.spectrum.eigenvectors[:, step.spectrum.leading_index])
    nv = LinearAlgebra.norm(v)
    nv == 0.0 && return nothing
    return v ./ nv
end

function pseudo_arclength_step_result(cfg::HopfSpectrumConfig, attempt, method, parent_K, target_K, attempts)
    is_target = abs(Float64(attempt.K) - Float64(target_K)) <= cfg.steady_bifurcation_K_tolerance
    return hopf_step_result(
        attempt.K,
        attempt;
        method,
        parent_K,
        target_K,
        is_inserted = !is_target,
        is_requested = is_target,
        attempts,
    )
end

function try_finish_pseudo_arclength_to_target(cfg::HopfSpectrumConfig, source_step, switch_step, K_target::Real, omega, attempts)
    retry_maxiters = cfg.nonlinear_maxiters * cfg.pseudo_arclength_maxiters_factor
    steps = Any[switch_step]

    final_attempt = attempt_hopf_solve(cfg, K_target, omega, switch_step.eq.q; maxiters = retry_maxiters)
    push!(attempts, "pseudo_arclength_finish_fixed_K(K=$(K_target), retcode=$(final_attempt.retcode), residual=$(final_attempt.residual_norm), error=$(final_attempt.error))")
    if final_attempt.success
        final_step = hopf_step_result(
            K_target,
            final_attempt;
            method = "pseudo_arclength_finish",
            parent_K = switch_step.K,
            target_K = K_target,
            is_inserted = false,
            is_requested = true,
            attempts,
        )
        if is_stable_step(cfg, final_step) || has_hopf_crossing(cfg, switch_step, final_step)
            push!(steps, final_step)
            return (; success = true, steps, attempts)
        end
        push!(attempts, "pseudo_arclength_finish_rejected(K=$(K_target), lead=$(leading_value(final_step)))")
    end

    scan_direction = sign(Float64(K_target) - Float64(source_step.K))
    scan_direction == 0.0 && return (; success = false, steps = Any[], attempts)
    previous = source_step
    current = switch_step

    for istep in 1:cfg.pseudo_arclength_max_steps
        y_prev = vcat(previous.eq.q, previous.K)
        y_current = vcat(current.eq.q, current.K)
        tangent = y_current .- y_prev
        nt = LinearAlgebra.norm(tangent)
        nt == 0.0 && return (; success = false, steps = Any[], attempts)
        tangent ./= nt
        tangent[end] * scan_direction < 0.0 && (tangent .*= -1.0)

        dK_current = abs(Float64(current.K) - Float64(previous.K))
        dK_target = abs(Float64(K_target) - Float64(current.K))
        dK_predict = dK_current == 0.0 ? dK_target : min(max(dK_current, cfg.steady_bifurcation_K_tolerance), dK_target)
        abs(tangent[end]) <= eps(Float64) && return (; success = false, steps = Any[], attempts)
        ds = dK_predict / abs(tangent[end])
        y_pred = y_current .+ ds .* tangent
        y_pred[end] = Float64(current.K) + scan_direction * dK_predict
        y_guess = copy(y_pred)

        a = attempt_pseudo_arclength_solve(cfg, omega, y_guess, y_pred, tangent; maxiters = retry_maxiters)
        push!(attempts, "pseudo_arclength_continue(step=$(istep), K=$(a.K), retcode=$(a.retcode), residual=$(a.residual_norm), arc=$(a.arclength_residual), error=$(a.error))")
        a.success || return (; success = false, steps = Any[], attempts)

        next_step = pseudo_arclength_step_result(
            cfg,
            a,
            "pseudo_arclength_continue",
            current.K,
            K_target,
            attempts,
        )
        if !(is_stable_step(cfg, next_step) || has_hopf_crossing(cfg, current, next_step))
            push!(attempts, "pseudo_arclength_continue_rejected(K=$(next_step.K), lead=$(leading_value(next_step)))")
            return (; success = false, steps = Any[], attempts)
        end
        push!(steps, next_step)
        if abs(Float64(next_step.K) - Float64(K_target)) <= cfg.steady_bifurcation_K_tolerance ||
           (Float64(next_step.K) - Float64(K_target)) * scan_direction >= 0.0
            final_attempt = attempt_hopf_solve(cfg, K_target, omega, next_step.eq.q; maxiters = retry_maxiters)
            push!(attempts, "pseudo_arclength_final_fixed_K(K=$(K_target), retcode=$(final_attempt.retcode), residual=$(final_attempt.residual_norm), error=$(final_attempt.error))")
            if final_attempt.success
                final_step = hopf_step_result(
                    K_target,
                    final_attempt;
                    method = "pseudo_arclength_finish",
                    parent_K = next_step.K,
                    target_K = K_target,
                    is_inserted = false,
                    is_requested = true,
                    attempts,
                )
                if is_stable_step(cfg, final_step) || has_hopf_crossing(cfg, next_step, final_step)
                    push!(steps, final_step)
                    return (; success = true, steps, attempts)
                end
                push!(attempts, "pseudo_arclength_final_rejected(K=$(K_target), lead=$(leading_value(final_step)))")
            end
            return (; success = false, steps = Any[], attempts)
        end

        previous = current
        current = next_step
    end

    return (; success = false, steps = Any[], attempts)
end

function pseudo_arclength_branch_switch_recovery(cfg::HopfSpectrumConfig, bracket, K_target::Real, omega)
    attempts = bracket.attempts
    pseudo_arclength_enabled(cfg) || return (; success = false, steps = Any[], attempts)

    source = bracket.stable_step
    vcrit = normalized_critical_vector(source)
    vcrit === nothing && return (; success = false, steps = Any[], attempts)

    retry_maxiters = cfg.nonlinear_maxiters * cfg.pseudo_arclength_maxiters_factor
    base_q = Float64.(source.eq.q)
    base_y = vcat(base_q, Float64(source.K))
    target_dK = Float64(K_target) - Float64(source.K)
    target_dK == 0.0 && (target_dK = Float64(bracket.unstable_step.K) - Float64(source.K))
    scan_direction = sign(target_dK)
    scan_direction == 0.0 && return (; success = false, steps = Any[], attempts)

    for factor in cfg.pseudo_arclength_dK_factors
        dK = factor * target_dK
        for eps in cfg.pseudo_arclength_epsilons
            for sign_eps in (-1.0, 1.0)
                q_pred = copy(base_q)
                q_pred[1:(end - 1)] .+= sign_eps * eps .* vcrit
                q_pred[end] = source.eq.Omega
                q_pred = rephase_equilibrium_guess(q_pred)
                y_pred = vcat(q_pred, Float64(source.K) + dK)
                tangent = y_pred .- base_y
                nt = LinearAlgebra.norm(tangent)
                nt == 0.0 && continue
                tangent ./= nt
                y_guess = copy(y_pred)

                a = attempt_pseudo_arclength_solve(cfg, omega, y_guess, y_pred, tangent; maxiters = retry_maxiters)
                push!(attempts, "pseudo_arclength_switch(factor=$(factor), eps=$(eps), sign=$(sign_eps), K=$(a.K), retcode=$(a.retcode), residual=$(a.residual_norm), arc=$(a.arclength_residual), error=$(a.error))")
                a.success || continue
                (Float64(a.K) - Float64(source.K)) * scan_direction >= -cfg.steady_bifurcation_K_tolerance || continue

                switch_step = pseudo_arclength_step_result(
                    cfg,
                    a,
                    "pseudo_arclength_switch",
                    source.K,
                    K_target,
                    attempts,
                )
                if !(is_stable_step(cfg, switch_step) || has_hopf_crossing(cfg, source, switch_step))
                    push!(attempts, "pseudo_arclength_switch_rejected(K=$(switch_step.K), lead=$(leading_value(switch_step)))")
                    continue
                end

                steps = Any[]
                append!(steps, bracket.inserted_stable_steps)
                if abs(Float64(switch_step.K) - Float64(K_target)) <= cfg.steady_bifurcation_K_tolerance
                    push!(steps, switch_step)
                    return (; success = true, steps, attempts)
                end

                finish = try_finish_pseudo_arclength_to_target(cfg, source, switch_step, K_target, omega, attempts)
                if finish.success
                    append!(steps, finish.steps)
                    return (; success = true, steps, attempts)
                end
            end
        end
    end

    return (; success = false, steps = Any[], attempts)
end

function stable_branch_switch_recovery(cfg::HopfSpectrumConfig, stable_step, unstable_step, K_target::Real, omega; initial_attempts = String[])
    attempts = copy(initial_attempts)
    push!(attempts, "real_stability_crossing(K_left=$(stable_step.K), K_right=$(unstable_step.K), lead_left=$(leading_value(stable_step)), lead_right=$(leading_value(unstable_step)))")
    bracket = localize_real_stability_crossing(cfg, stable_step, unstable_step, omega, attempts)
    bracket.success || return (; success = false, steps = Any[], attempts = bracket.attempts)

    retry_maxiters = cfg.nonlinear_maxiters * cfg.recovery_retry_maxiters_factor
    switch_source = bracket.stable_step
    guesses = branch_switch_guesses(cfg, switch_source)

    switch_targets = Float64[]
    push!(switch_targets, Float64(bracket.unstable_step.K))
    if abs(Float64(K_target) - Float64(bracket.unstable_step.K)) > cfg.steady_bifurcation_K_tolerance
        push!(switch_targets, Float64(K_target))
    end

    if !isempty(guesses)
        for K_switch in switch_targets
            for (iguess, guess) in pairs(guesses)
                a = attempt_hopf_solve(cfg, K_switch, omega, guess; maxiters = retry_maxiters)
                push!(bracket.attempts, "stable_branch_switch(guess=$(iguess), K=$(K_switch), retcode=$(a.retcode), residual=$(a.residual_norm), error=$(a.error))")
                a.success || continue

                switch_step = hopf_step_result(
                    K_switch,
                    a;
                    method = "stable_branch_switch",
                    parent_K = switch_source.K,
                    target_K = K_target,
                    is_inserted = abs(K_switch - Float64(K_target)) > cfg.steady_bifurcation_K_tolerance,
                    is_requested = abs(K_switch - Float64(K_target)) <= cfg.steady_bifurcation_K_tolerance,
                    attempts = bracket.attempts,
                )
                is_stable_step(cfg, switch_step) || continue

                steps = Any[]
                append!(steps, bracket.inserted_stable_steps)
                push!(steps, switch_step)

                if abs(K_switch - Float64(K_target)) <= cfg.steady_bifurcation_K_tolerance
                    return (; success = true, steps, attempts = bracket.attempts)
                end

                final_attempt = attempt_hopf_solve(cfg, K_target, omega, switch_step.eq.q; maxiters = retry_maxiters)
                push!(bracket.attempts, "stable_branch_continue(K=$(K_target), retcode=$(final_attempt.retcode), residual=$(final_attempt.residual_norm), error=$(final_attempt.error))")
                if final_attempt.success
                    final_step = hopf_step_result(
                        K_target,
                        final_attempt;
                        method = "stable_branch_continue",
                        parent_K = switch_step.K,
                        target_K = K_target,
                        is_inserted = false,
                        is_requested = true,
                        attempts = bracket.attempts,
                    )
                    if is_stable_step(cfg, final_step) || has_hopf_crossing(cfg, switch_step, final_step)
                        push!(steps, final_step)
                        return (; success = true, steps, attempts = bracket.attempts)
                    end
                    push!(bracket.attempts, "stable_branch_continue_rejected(K=$(K_target), lead=$(leading_value(final_step)))")
                else
                    continuation_recovery = recover_hopf_step(cfg, switch_step.K, switch_step, K_target, omega, final_attempt)
                    append!(bracket.attempts, continuation_recovery.attempts)
                    if continuation_recovery.success && !isempty(continuation_recovery.steps)
                        recovered_steps = continuation_recovery.steps
                        final_recovered = recovered_steps[end]
                        if is_stable_step(cfg, final_recovered) || has_hopf_crossing(cfg, switch_step, final_recovered)
                            append!(steps, recovered_steps)
                            return (; success = true, steps, attempts = bracket.attempts)
                        end
                        push!(bracket.attempts, "stable_branch_recovery_rejected(K=$(final_recovered.K), lead=$(leading_value(final_recovered)))")
                    end
                end
            end
        end
    end

    pseudo = pseudo_arclength_branch_switch_recovery(cfg, bracket, K_target, omega)
    pseudo.success && return pseudo

    return (; success = false, steps = Any[], attempts = bracket.attempts)
end

function recover_hopf_step(cfg::HopfSpectrumConfig, K_prev::Real, previous_step, K_target::Real, omega, direct_attempt)
    attempts = String[
        "direct(retcode=$(direct_attempt.retcode), residual=$(direct_attempt.residual_norm), error=$(direct_attempt.error))",
    ]
    retry_maxiters = cfg.nonlinear_maxiters * cfg.recovery_retry_maxiters_factor

    retry = attempt_hopf_solve(cfg, K_target, omega, previous_step.eq.q; maxiters = retry_maxiters)
    push!(attempts, "retry_newton(maxiters=$(retry_maxiters), retcode=$(retry.retcode), residual=$(retry.residual_norm), error=$(retry.error))")
    if retry.success
        step = hopf_step_result(
            K_target,
            retry;
            method = "retry_newton",
            parent_K = K_prev,
            target_K = K_target,
            is_inserted = false,
            is_requested = true,
            attempts,
        )
        return (; success = true, steps = [step], attempts)
    end

    for level in 1:cfg.recovery_max_halvings
        nsegments = 2^level
        K_segments = collect(range(K_prev, K_target; length = nsegments + 1))[2:end]
        state = previous_step
        steps = Any[]
        ok = true
        local_attempts = copy(attempts)
        for (jseg, Kseg) in pairs(K_segments)
            a = attempt_hopf_solve(cfg, Kseg, omega, state.eq.q; maxiters = retry_maxiters)
            push!(local_attempts, "step_halving(level=$(level), segment=$(jseg)/$(length(K_segments)), K=$(Kseg), retcode=$(a.retcode), residual=$(a.residual_norm), error=$(a.error))")
            if a.success
                is_target = jseg == lastindex(K_segments)
                step = hopf_step_result(
                    Kseg,
                    a;
                    method = is_target ? "step_halving" : "step_halving_inserted",
                    parent_K = state.K,
                    target_K = K_target,
                    is_inserted = !is_target,
                    is_requested = is_target,
                    attempts = local_attempts,
                )
                push!(steps, step)
                state = step
            else
                ok = false
                break
            end
        end
        if ok
            final_step = steps[end]
            if cfg.enable_simulation_recovery && real_branch_obstruction(cfg, previous_step, final_step)
                sim = simulation_recovery_attempt(cfg, K_target, omega)
                residual = sim.attempt === nothing ? Inf : sim.attempt.residual_norm
                retcode = sim.attempt === nothing ? "Exception" : sim.attempt.retcode
                push!(local_attempts, "simulation_recovery_after_step_halving(fft_stationary=$(sim.fft_stationary), retcode=$(retcode), residual=$(residual), error=$(sim.error))")
                if sim.success
                    sim_step = hopf_step_result(
                        K_target,
                        sim.attempt;
                        method = "simulation_recovery",
                        parent_K = K_prev,
                        target_K = K_target,
                        is_inserted = false,
                        is_requested = true,
                        attempts = local_attempts,
                    )
                    if leading_real_part(sim_step) < leading_real_part(final_step) - cfg.real_tolerance
                        return (; success = true, steps = [sim_step], attempts = local_attempts)
                    end
                end
            end
            return (; success = true, steps, attempts = local_attempts)
        end
        append!(attempts, local_attempts[(length(attempts) + 1):end])
    end

    guesses = branch_switch_guesses(cfg, previous_step)
    for (iguess, guess) in pairs(guesses)
        a = attempt_hopf_solve(cfg, K_target, omega, guess; maxiters = retry_maxiters)
        push!(attempts, "branch_switch(guess=$(iguess), K=$(K_target), retcode=$(a.retcode), residual=$(a.residual_norm), error=$(a.error))")
        if a.success
            step = hopf_step_result(
                K_target,
                a;
                method = "branch_switch",
                parent_K = K_prev,
                target_K = K_target,
                is_inserted = false,
                is_requested = true,
                attempts,
            )
            return (; success = true, steps = [step], attempts)
        end

        K_mid = 0.5 * (Float64(K_prev) + Float64(K_target))
        amid = attempt_hopf_solve(cfg, K_mid, omega, guess; maxiters = retry_maxiters)
        push!(attempts, "branch_switch_midpoint(guess=$(iguess), K=$(K_mid), retcode=$(amid.retcode), residual=$(amid.residual_norm), error=$(amid.error))")
        if amid.success
            mid_step = hopf_step_result(
                K_mid,
                amid;
                method = "branch_switch_inserted",
                parent_K = K_prev,
                target_K = K_target,
                is_inserted = true,
                is_requested = false,
                attempts,
            )
            afinal = attempt_hopf_solve(cfg, K_target, omega, mid_step.eq.q; maxiters = retry_maxiters)
            push!(attempts, "branch_switch_from_midpoint(guess=$(iguess), K=$(K_target), retcode=$(afinal.retcode), residual=$(afinal.residual_norm), error=$(afinal.error))")
            if afinal.success
                final_step = hopf_step_result(
                    K_target,
                    afinal;
                    method = "branch_switch",
                    parent_K = K_mid,
                    target_K = K_target,
                    is_inserted = false,
                    is_requested = true,
                    attempts,
                )
                return (; success = true, steps = [mid_step, final_step], attempts)
            end
        end
    end

    if cfg.enable_simulation_recovery
        sim = simulation_recovery_attempt(cfg, K_target, omega)
        residual = sim.attempt === nothing ? Inf : sim.attempt.residual_norm
        retcode = sim.attempt === nothing ? "Exception" : sim.attempt.retcode
        push!(attempts, "simulation_recovery(fft_stationary=$(sim.fft_stationary), retcode=$(retcode), residual=$(residual), error=$(sim.error))")
        if sim.success
            step = hopf_step_result(
                K_target,
                sim.attempt;
                method = "simulation_recovery",
                parent_K = K_prev,
                target_K = K_target,
                is_inserted = false,
                is_requested = true,
                attempts,
            )
            return (; success = true, steps = [step], attempts)
        end
    end

    return (; success = false, steps = Any[], attempts)
end

function is_hopf_like(value::Complex, cfg::HopfSpectrumConfig)
    return isfinite(real(value)) && isfinite(imag(value)) && abs(imag(value)) >= cfg.hopf_imag_tolerance
end

function has_hopf_crossing(left::Complex, right::Complex, cfg::HopfSpectrumConfig)
    is_hopf_like(left, cfg) && is_hopf_like(right, cfg) || return false
    lr = real(left)
    rr = real(right)
    return lr * rr <= 0.0 || min(abs(lr), abs(rr)) <= cfg.real_tolerance
end

function interpolate_hopf_K(K_left::Real, K_right::Real, left::Complex, right::Complex)
    lr = real(left)
    rr = real(right)
    denom = rr - lr
    denom == 0.0 && return 0.5 * (Float64(K_left) + Float64(K_right))
    alpha = clamp(-lr / denom, 0.0, 1.0)
    return Float64(K_left) + alpha * (Float64(K_right) - Float64(K_left))
end

function detect_hopf_candidates(cfg::HopfSpectrumConfig, K_values, store, omega)
    left_indices = Int[]
    right_indices = Int[]
    K_hopf = Float64[]
    imag_hopf = Float64[]
    pr = Float64[]
    omega_eff = Float64[]
    dominant_index = Int[]
    dominant_omega = Float64[]
    dominant_weight = Float64[]
    classification = String[]

    for j in 2:length(K_values)
        store.accepted[j - 1] && store.accepted[j] || continue
        left = store.leading_eigenvalues[j - 1]
        right = store.leading_eigenvalues[j]
        has_hopf_crossing(left, right, cfg) || continue

        use_idx = abs(real(left)) <= abs(real(right)) ? j - 1 : j
        dom = store.dominant_index[use_idx]
        push!(left_indices, j - 1)
        push!(right_indices, j)
        push!(K_hopf, interpolate_hopf_K(K_values[j - 1], K_values[j], left, right))
        push!(imag_hopf, abs(imag(store.leading_eigenvalues[use_idx])))
        push!(pr, store.participation_ratio[use_idx])
        push!(omega_eff, store.omega_eff[use_idx])
        push!(dominant_index, dom)
        push!(dominant_omega, dom == 0 ? NaN : omega[dom])
        push!(dominant_weight, store.dominant_weight[use_idx])
        push!(classification, store.classifications[use_idx])
    end

    return (;
        left_indices,
        right_indices,
        K_hopf,
        imag_hopf,
        participation_ratio = pr,
        omega_eff,
        dominant_index,
        dominant_omega,
        dominant_weight,
        classification,
    )
end

function selected_critical_index(candidates, store)
    if !isempty(candidates.K_hopf)
        pairs = candidates.right_indices
        return pairs[argmin(abs.(real.(store.leading_eigenvalues[pairs])))]
    end
    accepted_indices = findall(store.accepted)
    isempty(accepted_indices) && return 0
    return accepted_indices[argmax(real.(store.leading_eigenvalues[accepted_indices]))]
end

function append_hopf_path_steps!(store, cfg::HopfSpectrumConfig, steps, previous_step)
    stopped_for_hopf = false
    for step in steps
        hopf_crossed = previous_step !== nothing && has_hopf_crossing(cfg, previous_step, step)
        if cfg.store_recovery_points || !step.is_inserted_recovery_point || hopf_crossed
            append_hopf_step!(store, step, cfg)
            println("  K = ", step.K,
                " residual = ", step.eq.residual_norm,
                ", lead = ", leading_value(step),
                ", PR = ", step.spectrum.participation_ratio[step.spectrum.leading_index],
                ", method = ", step.recovery_method)
        else
            println("  used unrecorded K = ", step.K,
                " residual = ", step.eq.residual_norm,
                ", lead = ", leading_value(step),
                ", method = ", step.recovery_method)
        end
        if cfg.stop_at_first_hopf && hopf_crossed
            println("  candidate Hopf crossing bracketed between K = ", previous_step.K, " and K = ", step.K)
            previous_step = step
            stopped_for_hopf = true
            break
        end
        previous_step = step
    end
    previous_guess = previous_step === nothing ? nothing : previous_step.eq.q
    return (; previous_step, previous_guess, stopped_for_hopf)
end

function write_hopf_spectrum_summary_csv(path::AbstractString, K_values, store, omega)
    open(path, "w") do io
        println(io, "K,residual,accepted,Re_lambda_lead,Im_lambda_lead,PR,omega_eff,dominant_index,dominant_omega,dominant_weight,classification,neutral_index,leading_index,neutral_candidate_count,recovery_method,parent_K,target_K,is_inserted_recovery_point,is_requested_K,recovery_attempts")
        for j in eachindex(K_values)
            dom = store.dominant_index[j]
            domomega = dom == 0 ? NaN : omega[dom]
            println(io, join((
                K_values[j],
                store.residuals[j],
                store.accepted[j],
                real(store.leading_eigenvalues[j]),
                imag(store.leading_eigenvalues[j]),
                store.participation_ratio[j],
                store.omega_eff[j],
                dom,
                domomega,
                store.dominant_weight[j],
                store.classifications[j],
                store.neutral_indices[j],
                store.leading_indices[j],
                store.neutral_candidate_counts[j],
                store.recovery_method[j],
                store.parent_K[j],
                store.target_K[j],
                store.is_inserted_recovery_point[j],
                store.is_requested_K[j],
                replace(store.recovery_attempts[j], "," => ";"),
            ), ","))
        end
    end
    return path
end

function write_recovery_diagnostics_csv(path::AbstractString, recovery_rows)
    open(path, "w") do io
        println(io, "parent_K,target_K,success,accepted_method,final_residual,attempts")
        for row in recovery_rows
            println(io, join((
                row.parent_K,
                row.target_K,
                row.success,
                row.accepted_method,
                row.final_residual,
                replace(row.attempts, "," => ";"),
            ), ","))
        end
    end
    return path
end

function vector_matrix(values::AbstractVector, ncols::Integer, fill_value)
    out = fill(fill_value, length(values), ncols)
    for (i, row) in pairs(values)
        n = min(length(row), ncols)
        n == 0 || (out[i, 1:n] .= row[1:n])
    end
    return out
end

function save_hopf_spectrum_data(path::AbstractString, cfg::HopfSpectrumConfig; kwargs...)
    JLD2.jldsave(path;
        C1 = cfg.C1,
        C2 = cfg.C2,
        N = cfg.N,
        Ttr_first = cfg.Ttr_first,
        T_first = cfg.T_first,
        dt = cfg.dt,
        analysis_window = cfg.analysis_window,
        omega_mean = cfg.omega_mean,
        omega_sigma = cfg.omega_sigma,
        omega_seed = cfg.omega_seed,
        omega_symmetry = normalize_omega_symmetry(cfg.omega_symmetry),
        recenter_omega = cfg.recenter_omega,
        initial_seed = cfg.initial_seed,
        residual_tolerance = cfg.residual_tolerance,
        continuation_method = hopf_continuation_method(cfg),
        hopf_imag_tolerance = cfg.hopf_imag_tolerance,
        real_tolerance = cfg.real_tolerance,
        stability_tolerance = cfg.stability_tolerance,
        stay_on_stable_branch = cfg.stay_on_stable_branch,
        stop_at_first_hopf = cfg.stop_at_first_hopf,
        steady_bifurcation_max_bisections = cfg.steady_bifurcation_max_bisections,
        steady_bifurcation_K_tolerance = cfg.steady_bifurcation_K_tolerance,
        branch_switch_epsilons = cfg.branch_switch_epsilons,
        enable_pseudo_arclength = cfg.enable_pseudo_arclength,
        pseudo_arclength_epsilons = cfg.pseudo_arclength_epsilons,
        pseudo_arclength_dK_factors = cfg.pseudo_arclength_dK_factors,
        pseudo_arclength_max_steps = cfg.pseudo_arclength_max_steps,
        pseudo_arclength_maxiters_factor = cfg.pseudo_arclength_maxiters_factor,
        pseudo_arclength_tolerance = cfg.pseudo_arclength_tolerance,
        localized_pr_threshold = cfg.localized_pr_threshold,
        collective_pr_fraction = cfg.collective_pr_fraction,
        kwargs...,
    )
    return path
end

function run_hopf_spectrum_detection(cfg::HopfSpectrumConfig)
    mkpath(cfg.output_dir)

    continuation_method = hopf_continuation_method(cfg)
    requested_K_values = hopf_k_values(cfg)
    isempty(requested_K_values) && error("K_values must not be empty.")
    store = initialize_hopf_storage(cfg.N)
    recovery_rows = NamedTuple[]

    println("Starting eigenvalue-based Hopf scan with ", length(requested_K_values), " requested K values.")
    println("Continuation method: ", continuation_method)
    first_locked_cfg = locked_config_for_hopf(cfg, first(requested_K_values))
    first_data = first_hopf_guess_and_diagnostics(cfg, first_locked_cfg)
    omega = first_data.omega
    previous_guess = first_data.guess
    previous_step = nothing
    stop_scan = false

    for (jK, K) in pairs(requested_K_values)
        println("Processing requested K = ", K, " (", jK, " / ", length(requested_K_values), ")")
        attempt = attempt_hopf_solve(cfg, K, omega, previous_guess)
        if attempt.success
            method = previous_step === nothing ? "initial_direct" : "direct"
            step = hopf_step_result(
                K,
                attempt;
                method,
                parent_K = previous_step === nothing ? NaN : previous_step.K,
                target_K = K,
                is_inserted = false,
                is_requested = true,
                attempts = ["$(method)(retcode=$(attempt.retcode), residual=$(attempt.residual_norm))"],
            )
            if cfg.stay_on_stable_branch && previous_step !== nothing && has_real_stability_crossing(cfg, previous_step, step)
                stable_recovery = stable_branch_switch_recovery(
                    cfg,
                    previous_step,
                    step,
                    K,
                    omega;
                    initial_attempts = ["direct_success_real_crossing(retcode=$(attempt.retcode), residual=$(attempt.residual_norm))"],
                )
                final_residual = stable_recovery.success && !isempty(stable_recovery.steps) ? stable_recovery.steps[end].eq.residual_norm : attempt.residual_norm
                accepted_method = stable_recovery.success && !isempty(stable_recovery.steps) ? stable_recovery.steps[end].recovery_method : "none"
                push!(recovery_rows, (;
                    parent_K = previous_step.K,
                    target_K = Float64(K),
                    success = stable_recovery.success,
                    accepted_method,
                    final_residual,
                    attempts = join(stable_recovery.attempts, " | "),
                ))
                if stable_recovery.success
                    appended = append_hopf_path_steps!(store, cfg, stable_recovery.steps, previous_step)
                    previous_step = appended.previous_step
                    previous_guess = appended.previous_guess
                    stop_scan = appended.stopped_for_hopf
                    stop_scan && break
                    continue
                else
                    @warn "Stopping Hopf scan because stable-branch switching failed after a real eigenvalue crossing." K=K attempts=join(stable_recovery.attempts, " | ")
                    cfg.stop_on_failure && break
                end
            else
                appended = append_hopf_path_steps!(store, cfg, [step], previous_step)
                previous_step = appended.previous_step
                previous_guess = appended.previous_guess
                stop_scan = appended.stopped_for_hopf
                stop_scan && break
                continue
            end
            continue
        end

        if previous_step === nothing || !cfg.enable_recovery
            push!(recovery_rows, (;
                parent_K = previous_step === nothing ? NaN : previous_step.K,
                target_K = Float64(K),
                success = false,
                accepted_method = "none",
                final_residual = attempt.residual_norm,
                attempts = "direct(retcode=$(attempt.retcode), residual=$(attempt.residual_norm), error=$(attempt.error))",
            ))
            @warn "Stopping Hopf scan because direct continuation failed and recovery is unavailable." K=K residual=attempt.residual_norm error=attempt.error
            cfg.stop_on_failure && break
            continue
        end

        recovery = recover_hopf_step(cfg, previous_step.K, previous_step, K, omega, attempt)
        final_residual = recovery.success && !isempty(recovery.steps) ? recovery.steps[end].eq.residual_norm : attempt.residual_norm
        accepted_method = recovery.success && !isempty(recovery.steps) ? recovery.steps[end].recovery_method : "none"
        push!(recovery_rows, (;
            parent_K = previous_step.K,
            target_K = Float64(K),
            success = recovery.success,
            accepted_method,
            final_residual,
            attempts = join(recovery.attempts, " | "),
        ))

        if recovery.success
            for step in recovery.steps
                if cfg.stay_on_stable_branch && previous_step !== nothing && has_real_stability_crossing(cfg, previous_step, step)
                    stable_recovery = stable_branch_switch_recovery(
                        cfg,
                        previous_step,
                        step,
                        K,
                        omega;
                        initial_attempts = ["recovery_real_crossing(method=$(step.recovery_method), K=$(step.K), lead=$(leading_value(step)))"],
                    )
                    stable_final_residual = stable_recovery.success && !isempty(stable_recovery.steps) ? stable_recovery.steps[end].eq.residual_norm : step.eq.residual_norm
                    stable_accepted_method = stable_recovery.success && !isempty(stable_recovery.steps) ? stable_recovery.steps[end].recovery_method : "none"
                    push!(recovery_rows, (;
                        parent_K = previous_step.K,
                        target_K = Float64(K),
                        success = stable_recovery.success,
                        accepted_method = stable_accepted_method,
                        final_residual = stable_final_residual,
                        attempts = join(stable_recovery.attempts, " | "),
                    ))
                    if stable_recovery.success
                        appended = append_hopf_path_steps!(store, cfg, stable_recovery.steps, previous_step)
                        previous_step = appended.previous_step
                        previous_guess = appended.previous_guess
                        stop_scan = appended.stopped_for_hopf
                    else
                        @warn "Stopping Hopf scan because stable-branch switching failed after recovery crossed a real eigenvalue." K=K attempts=join(stable_recovery.attempts, " | ")
                        stop_scan = cfg.stop_on_failure
                    end
                    break
                end

                appended = append_hopf_path_steps!(store, cfg, [step], previous_step)
                previous_step = appended.previous_step
                previous_guess = appended.previous_guess
                stop_scan = appended.stopped_for_hopf
                stop_scan && break
            end
            stop_scan && break
        else
            @warn "Stopping Hopf scan because recovery failed." K=K attempts=join(recovery.attempts, " | ")
            if cfg.stop_on_failure
                break
            end
        end
    end

    first_fft_stationary = first_data.fftdiag.stationary
    first_fft_stationary || @warn "First-K FFT stationarity safety check failed." max_amplitude_fft_peak=first_data.fftdiag.max_amplitude_fft_peak max_mean_field_fft_peak=first_data.fftdiag.max_mean_field_fft_peak

    K_path = store.K_path
    candidates = detect_hopf_candidates(cfg, K_path, store, omega)
    critical_index = selected_critical_index(candidates, store)
    requested_count = count(store.is_requested_K)
    requested_count == length(requested_K_values) && all(store.accepted) ||
        error("Stability continuation incomplete: $requested_count/$(length(requested_K_values)) requested points. No complete dataset was saved.")
    base_path = joinpath(cfg.output_dir, cfg.filename_base)

    csv_path = write_hopf_spectrum_summary_csv(base_path * "_summary.csv", K_path, store, omega)
    recovery_csv_path = write_recovery_diagnostics_csv(base_path * "_recovery_diagnostics.csv", recovery_rows)

    data_path = save_hopf_spectrum_data(base_path * ".jld2", cfg;
        K_values = requested_K_values,
        K_path,
        omega = Float64.(omega),
        first_u0 = Float64.(first_data.u0),
        first_t = Float64.(first_data.t),
        first_analysis_indices = Int.(first_data.idx),
        first_mean_frequency = first_data.phase.mean_frequency,
        first_frequency_spread = first_data.phase.frequency_spread,
        first_Omega_phase_fit = first_data.phase.Omega,
        first_fft_freqs = first_data.fftdiag.freqs,
        first_amplitude_spectra = first_data.fftdiag.amplitude_spectra,
        first_mean_field_fft_freqs = first_data.fftdiag.mean_field_freqs,
        first_mean_field_spectrum = first_data.fftdiag.mean_field_spectrum,
        first_max_amplitude_fft_peak = first_data.fftdiag.max_amplitude_fft_peak,
        first_max_mean_field_fft_peak = first_data.fftdiag.max_mean_field_fft_peak,
        first_fft_stationary,
        residuals = store.residuals,
        Omegas = store.Omegas,
        accepted = store.accepted,
        retcodes = store.retcodes,
        leading_eigenvalues = store.leading_eigenvalues,
        neutral_values = store.neutral_values,
        neutral_indices = store.neutral_indices,
        leading_indices = store.leading_indices,
        neutral_candidate_counts = store.neutral_candidate_counts,
        participation_ratio = store.participation_ratio,
        omega_eff = store.omega_eff,
        dominant_index = store.dominant_index,
        dominant_weight = store.dominant_weight,
        classifications = store.classifications,
        all_eigenvalues = vector_matrix(store.all_eigenvalues, 2cfg.N, ComplexF64(NaN, NaN)),
        leading_weights = vector_matrix(store.leading_weights, cfg.N, NaN),
        locked_amplitudes = vector_matrix(store.locked_amplitudes, cfg.N, NaN),
        final_equilibria = store.final_equilibria,
        recovery_method = store.recovery_method,
        parent_K = store.parent_K,
        target_K = store.target_K,
        is_inserted_recovery_point = store.is_inserted_recovery_point,
        is_requested_K = store.is_requested_K,
        recovery_attempts = store.recovery_attempts,
        recovery_rows,
        candidate_left_indices = candidates.left_indices,
        candidate_right_indices = candidates.right_indices,
        candidate_K_hopf = candidates.K_hopf,
        candidate_imag_hopf = candidates.imag_hopf,
        candidate_participation_ratio = candidates.participation_ratio,
        candidate_omega_eff = candidates.omega_eff,
        candidate_dominant_index = candidates.dominant_index,
        candidate_dominant_omega = candidates.dominant_omega,
        candidate_dominant_weight = candidates.dominant_weight,
        candidate_classification = candidates.classification,
        critical_index,
        summary_csv = csv_path,
        recovery_diagnostics_csv = recovery_csv_path,
    )

    println("Saved Hopf spectrum JLD2: ", data_path)
    println("Saved Hopf summary CSV: ", csv_path)
    println("Saved recovery diagnostics CSV: ", recovery_csv_path)
    println("First-K FFT stationary: ", first_fft_stationary)
    println("Accepted path points: ", count(store.accepted), " / ", length(store.accepted))
    println("Candidate Hopf crossings: ", length(candidates.K_hopf))

    return (;
        data_path,
        csv_path,
        recovery_csv_path,
        K_values = requested_K_values,
        K_path,
        omega,
        first_fft_stationary,
        store,
        candidates,
        critical_index,
    )
end

function supplementary_stability_config(; smoke::Bool = false)
    return HopfSpectrumConfig(
        output_dir = joinpath(DATA_ROOT, "supplementary", "numerical"),
        filename_base = "locked_spectrum",
        K_start = 1.0, K_end = smoke ? 0.98 : 0.95, K_steps = smoke ? 3 : 51,
        C1 = -0.5, C2 = -7.0, N = smoke ? 12 : 5000,
        Ttr_first = smoke ? 100 : 2000, T_first = smoke ? 30.0 : 1000.0,
        dt = 0.1, analysis_window = smoke ? 10.0 : 100.0,
        omega_mean = 0.0, omega_sigma = 0.3, omega_seed = 1234,
        omega_symmetry = :sym, recenter_omega = false, initial_seed = 5678,
        residual_tolerance = 1.0e-8, phase_frequency_tolerance = 1.0e-5,
        amplitude_fft_threshold = 1.0e-4, mean_field_fft_threshold = 1.0e-4,
        hopf_imag_tolerance = 1.0e-6, real_tolerance = 1.0e-6,
        stability_tolerance = 1.0e-7, stay_on_stable_branch = true,
        stop_at_first_hopf = false, continuation_method = :fixed_K,
        steady_bifurcation_max_bisections = 10, steady_bifurcation_K_tolerance = 1.0e-5,
        localized_pr_threshold = 3.0, collective_pr_fraction = 0.25)
end
