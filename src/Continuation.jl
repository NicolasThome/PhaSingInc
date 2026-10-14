Base.@kwdef struct KContinuationConfig
    output_dir::String = joinpath(DATA_ROOT, "numerical")
    filename_base::String = "K_continuation"

    K_start::Float64 = 1.0
    K_end::Float64 = 0.85
    K_steps::Int = 151
    C1::Float64 = -0.5
    C2::Float64 = -7.0
    N::Int = 5000
    Ttr_first::Int = 3000
    Ttr::Int = 3000
    T::Float64 = 600.0
    dt::Float64 = 0.1
    reltol::Float64 = 1.0e-8
    abstol::Float64 = 1.0e-8
    diagnostic_mode::Symbol = :online

    omega_mean::Float64 = 0.0
    omega_sigma::Float64 = 0.3
    omega_seed::Int = 1234
    omega_symmetry::Any = :sym
    recenter_omega::Bool = false

    initial_seed::Int = 5678
    continuation_noise::Float64 = 0.0
    savefiles::Bool = true

    amplitude_cluster_threshold::Union{Nothing, Float64} = 0.7
    amplitude_cluster_reltol::Float64 = 1.0e-4
    amplitude_cluster_abstol::Float64 = 1.0e-4

end

function mean_without_nan(values)
    finite_values = [value for value in values if isfinite(value)]
    isempty(finite_values) && return NaN
    return sum(finite_values) / length(finite_values)
end

@inline function wrapped_phase_difference(theta::Real, theta_prev::Real)
    dtheta = Float64(theta) - Float64(theta_prev)
    if dtheta > pi
        dtheta -= 2pi
    elseif dtheta < -pi
        dtheta += 2pi
    end
    return dtheta
end

k_continuation_alg(cfg) = OrdinaryDiffEqTsit5.Tsit5()

sample_times(Ttr::Integer, T::Real, dt::Real) = begin
    relative = collect(0.0:Float64(dt):Float64(T))
    if isempty(relative) || relative[end] < Float64(T) - 1.0e-10
        push!(relative, Float64(T))
    end
    Float64(Ttr) .+ relative
end

mutable struct OnlineStepDiagnostics
    omega::Vector{Float64}
    omega_var::Float64
    min_radius::Vector{Float64}
    max_radius::Vector{Float64}
    phase_delta_sum::Vector{Float64}
    prev_phase::Vector{Float64}
    prev_t::Float64
    elapsed::Float64
    sample_count::Int
    has_prev::Bool
    R_sum::Float64
    R_end::Float64
    Z_sum::Float64
    Z_end::Float64
    Rfreq_sum::Float64
    Rfreq_end::Float64
    Rfreq_count::Int
end

function OnlineStepDiagnostics(omega::AbstractVector{<:Real})
    n = length(omega)
    omega_float = Float64.(omega)
    return OnlineStepDiagnostics(
        omega_float,
        var(omega_float; corrected = false),
        fill(Inf, n),
        fill(-Inf, n),
        zeros(Float64, n),
        zeros(Float64, n),
        NaN,
        0.0,
        0,
        false,
        0.0,
        NaN,
        0.0,
        NaN,
        0.0,
        NaN,
        0,
    )
end

function update_online_step_diagnostics!(diag::OnlineStepDiagnostics, u::AbstractVector{<:Real}, t::Real)
    n = length(diag.omega)
    phase_sum = 0.0 + 0.0im
    amplitude_sum = 0.0 + 0.0im
    freq_mean = 0.0
    freq_m2 = 0.0
    nseen = 0

    @inbounds for j in 1:n
        x = Float64(u[2j - 1])
        y = Float64(u[2j])
        r = hypot(x, y)
        theta = atan(y, x)
        diag.min_radius[j] = min(diag.min_radius[j], r)
        diag.max_radius[j] = max(diag.max_radius[j], r)
        diag.prev_phase[j] = diag.has_prev ? diag.prev_phase[j] : theta
        phase_sum += cis(theta)
        amplitude_sum += x + im * y

        if diag.has_prev
            dtheta = wrapped_phase_difference(theta, diag.prev_phase[j])
            diag.phase_delta_sum[j] += dtheta
            omega_inst = dtheta / (Float64(t) - diag.prev_t)
            nseen += 1
            delta = omega_inst - freq_mean
            freq_mean += delta / nseen
            freq_m2 += delta * (omega_inst - freq_mean)
        end
        diag.prev_phase[j] = theta
    end

    diag.sample_count += 1
    diag.R_end = abs(phase_sum / n)
    diag.Z_end = abs(amplitude_sum / n)
    diag.R_sum += diag.R_end
    diag.Z_sum += diag.Z_end

    if diag.has_prev
        diag.elapsed += Float64(t) - diag.prev_t
        freq_var = freq_m2 / n
        diag.Rfreq_end = iszero(diag.omega_var) ? (freq_var <= 1.0e-14 ? 1.0 : 0.0) :
                         max(1.0 - freq_var / diag.omega_var, 0.0)
        diag.Rfreq_sum += diag.Rfreq_end
        diag.Rfreq_count += 1
    end

    diag.prev_t = Float64(t)
    diag.has_prev = true
    return diag
end

function finish_online_step_diagnostics(diag::OnlineStepDiagnostics)
    elapsed = diag.elapsed
    mean_frequency = elapsed <= 0.0 ? fill(NaN, length(diag.omega)) : diag.phase_delta_sum ./ elapsed
    return (;
        mean_frequency,
        min_radius = replace(diag.min_radius, Inf => NaN),
        max_radius = replace(diag.max_radius, -Inf => NaN),
        R_mean = diag.sample_count == 0 ? NaN : diag.R_sum / diag.sample_count,
        R_end = diag.R_end,
        Z_mean = diag.sample_count == 0 ? NaN : diag.Z_sum / diag.sample_count,
        Z_end = diag.Z_end,
        Rfreq_mean = diag.Rfreq_count == 0 ? NaN : diag.Rfreq_sum / diag.Rfreq_count,
        Rfreq_end = diag.Rfreq_end,
    )
end

mutable struct OnlineMinRadiusDiagnostics
    min_radius::Vector{Float64}
end

OnlineMinRadiusDiagnostics(n::Integer) = OnlineMinRadiusDiagnostics(fill(Inf, Int(n)))

function update_online_min_radius_diagnostics!(diag::OnlineMinRadiusDiagnostics, u::AbstractVector{<:Real})
    n = length(diag.min_radius)
    @inbounds for j in 1:n
        r = hypot(Float64(u[2j - 1]), Float64(u[2j]))
        diag.min_radius[j] = min(diag.min_radius[j], r)
    end
    return diag
end

function finish_online_min_radius_diagnostics(diag::OnlineMinRadiusDiagnostics)
    return replace(diag.min_radius, Inf => NaN)
end

mutable struct OnlineGroupDiagnostics
    indices::Vector{Int}
    omega_var::Float64
    prev_phase::Vector{Float64}
    prev_t::Float64
    has_prev::Bool
    sample_count::Int
    R_sum::Float64
    R_end::Float64
    Z_sum::Float64
    Z_end::Float64
    Rfreq_sum::Float64
    Rfreq_end::Float64
    Rfreq_count::Int
end

function OnlineGroupDiagnostics(indices::AbstractVector{<:Integer}, omega::AbstractVector{<:Real})
    idx = Int.(collect(indices))
    return OnlineGroupDiagnostics(
        idx,
        isempty(idx) ? NaN : var(Float64.(omega[idx]); corrected = false),
        zeros(Float64, length(idx)),
        NaN,
        false,
        0,
        0.0,
        NaN,
        0.0,
        NaN,
        0.0,
        NaN,
        0,
    )
end

function update_online_group_diagnostics!(diag::OnlineGroupDiagnostics, u::AbstractVector{<:Real}, t::Real)
    n = length(diag.indices)
    if n == 0
        diag.sample_count += 1
        diag.prev_t = Float64(t)
        diag.has_prev = true
        return diag
    end

    phase_sum = 0.0 + 0.0im
    amplitude_sum = 0.0 + 0.0im
    freq_mean = 0.0
    freq_m2 = 0.0
    nseen = 0

    @inbounds for (local_j, idx) in pairs(diag.indices)
        x = Float64(u[2idx - 1])
        y = Float64(u[2idx])
        theta = atan(y, x)
        diag.prev_phase[local_j] = diag.has_prev ? diag.prev_phase[local_j] : theta
        phase_sum += cis(theta)
        amplitude_sum += x + im * y

        if diag.has_prev
            dtheta = wrapped_phase_difference(theta, diag.prev_phase[local_j])
            omega_inst = dtheta / (Float64(t) - diag.prev_t)
            nseen += 1
            delta = omega_inst - freq_mean
            freq_mean += delta / nseen
            freq_m2 += delta * (omega_inst - freq_mean)
        end
        diag.prev_phase[local_j] = theta
    end

    diag.sample_count += 1
    diag.R_end = abs(phase_sum / n)
    diag.Z_end = abs(amplitude_sum / n)
    diag.R_sum += diag.R_end
    diag.Z_sum += diag.Z_end

    if diag.has_prev
        freq_var = freq_m2 / n
        diag.Rfreq_end = iszero(diag.omega_var) ? (freq_var <= 1.0e-14 ? 1.0 : 0.0) :
                         max(1.0 - freq_var / diag.omega_var, 0.0)
        diag.Rfreq_sum += diag.Rfreq_end
        diag.Rfreq_count += 1
    end

    diag.prev_t = Float64(t)
    diag.has_prev = true
    return diag
end

function finish_online_group_diagnostics(diag::OnlineGroupDiagnostics)
    isempty(diag.indices) && return (R_mean = NaN, R_end = NaN, Z_mean = NaN, Z_end = NaN, Rfreq_mean = NaN, Rfreq_end = NaN)
    return (;
        R_mean = diag.sample_count == 0 ? NaN : diag.R_sum / diag.sample_count,
        R_end = diag.R_end,
        Z_mean = diag.sample_count == 0 ? NaN : diag.Z_sum / diag.sample_count,
        Z_end = diag.Z_end,
        Rfreq_mean = diag.Rfreq_count == 0 ? NaN : diag.Rfreq_sum / diag.Rfreq_count,
        Rfreq_end = diag.Rfreq_end,
    )
end

function initialize_k_integrator(cfg::KContinuationConfig, K::Real, omega::AbstractVector{<:Real}, u0::AbstractVector{<:Real}, Ttr::Integer; reltol = cfg.reltol, abstol = cfg.abstol)
    total_time = Float64(Ttr) + cfg.T
    p = (Float64(K), cfg.C1, cfg.C2, Float64.(omega))
    prob = SciMLBase.ODEProblem(Model.sl_harmo_hetero!, Float64.(u0), (0.0, total_time), p)
    return SciMLBase.init(
        prob,
        k_continuation_alg(cfg);
        reltol = reltol,
        abstol = abstol,
        save_everystep = false,
        save_start = false,
    )
end

function step_to_sample!(integrator, t::Real)
    target = Float64(t)
    current = Float64(integrator.t)
    if target > current + 1.0e-12
        SciMLBase.step!(integrator, target - current, true)
    end
    return integrator
end

function run_k_step_online_core(cfg::KContinuationConfig, K::Real, omega::AbstractVector{<:Real}, u0::AbstractVector{<:Real}, Ttr::Integer)
    integrator = initialize_k_integrator(cfg, K, omega, u0, Ttr)
    diag = OnlineStepDiagnostics(omega)
    times = sample_times(Ttr, cfg.T, cfg.dt)

    for t in times
        step_to_sample!(integrator, t)
        update_online_step_diagnostics!(diag, integrator.u, t)
    end

    final_time = Float64(Ttr) + cfg.T
    step_to_sample!(integrator, final_time)
    return merge(finish_online_step_diagnostics(diag), (; final_state = Float64.(integrator.u)))
end

function run_k_step_online_cluster(cfg::KContinuationConfig, K::Real, omega::AbstractVector{<:Real}, u0::AbstractVector{<:Real}, Ttr::Integer, high_indices, low_indices)
    integrator = initialize_k_integrator(cfg, K, omega, u0, Ttr; reltol = cfg.amplitude_cluster_reltol, abstol = cfg.amplitude_cluster_abstol)
    high_diag = OnlineGroupDiagnostics(high_indices, omega)
    low_diag = OnlineGroupDiagnostics(low_indices, omega)
    times = sample_times(Ttr, cfg.T, cfg.dt)

    for t in times
        step_to_sample!(integrator, t)
        update_online_group_diagnostics!(high_diag, integrator.u, t)
        update_online_group_diagnostics!(low_diag, integrator.u, t)
    end

    return (high = finish_online_group_diagnostics(high_diag), low = finish_online_group_diagnostics(low_diag))
end

function run_k_step_cluster_classification(cfg::KContinuationConfig, K::Real, omega::AbstractVector{<:Real}, u0::AbstractVector{<:Real}, Ttr::Integer)
    integrator = initialize_k_integrator(cfg, K, omega, u0, Ttr; reltol = cfg.amplitude_cluster_reltol, abstol = cfg.amplitude_cluster_abstol)
    diag = OnlineMinRadiusDiagnostics(length(omega))
    times = sample_times(Ttr, cfg.T, cfg.dt)

    for t in times
        step_to_sample!(integrator, t)
        update_online_min_radius_diagnostics!(diag, integrator.u)
    end

    rough_min_radius = finish_online_min_radius_diagnostics(diag)
    all_indices = collect(eachindex(omega))
    high_indices = all_indices[rough_min_radius .>= Float64(cfg.amplitude_cluster_threshold)]
    low_indices = all_indices[rough_min_radius .< Float64(cfg.amplitude_cluster_threshold)]
    return (; high_indices, low_indices, min_radii = rough_min_radius)
end

function run_k_step_online(cfg::KContinuationConfig, K::Real, omega::AbstractVector{<:Real}, u0::AbstractVector{<:Real}, Ttr::Integer; cluster_order_parameters = nothing)
    core = run_k_step_online_core(cfg, K, omega, u0, Ttr)
    if cluster_order_parameters === nothing
        return merge(core, (; clusters = nothing, cluster_stats = nothing))
    end

    clusters = run_k_step_cluster_classification(cfg, K, omega, u0, Ttr)
    cluster_stats = run_k_step_online_cluster(cfg, K, omega, u0, Ttr, clusters.high_indices, clusters.low_indices)
    return merge(core, (; clusters, cluster_stats))
end

function write_continuation_csv(path::AbstractString, K_list, values::AbstractMatrix, indices, omega; system_mean_frequency = nothing)
    open(path, "w") do io
        labels = ["osc_$(idx)_omega_$(round(Float64(omega[idx]); digits = 6))" for idx in indices]
        if system_mean_frequency !== nothing
            labels = vcat("system_mean_frequency", labels)
        end
        println(io, join(vcat("K", labels), ","))
        for i in eachindex(K_list)
            row = [string(K_list[i])]
            if system_mean_frequency !== nothing
                push!(row, string(system_mean_frequency[i]))
            end
            append!(row, string.(values[i, :]))
            println(io, join(row, ","))
        end
    end
    return path
end

function write_order_parameter_csv(path::AbstractString, K_list, R_mean, R_end, Z_mean, Z_end, Rfreq_mean, Rfreq_end; cluster_order_parameters = nothing)
    open(path, "w") do io
        if cluster_order_parameters === nothing
            println(io, "K,R_mean,R_end,Z_mean,Z_end,Rfreq_mean,Rfreq_end")
        else
            println(io, "K,R_mean,R_end,Z_mean,Z_end,Rfreq_mean,Rfreq_end,N_high,N_low,R_high_mean,R_high_end,Z_high_mean,Z_high_end,Rfreq_high_mean,Rfreq_high_end,R_low_mean,R_low_end,Z_low_mean,Z_low_end,Rfreq_low_mean,Rfreq_low_end")
        end
        for i in eachindex(K_list)
            row = Any[
                K_list[i],
                R_mean[i],
                R_end[i],
                Z_mean[i],
                Z_end[i],
                Rfreq_mean[i],
                Rfreq_end[i],
            ]
            if cluster_order_parameters !== nothing
                append!(row, Any[
                    cluster_order_parameters.N_high[i],
                    cluster_order_parameters.N_low[i],
                    cluster_order_parameters.R_high_mean[i],
                    cluster_order_parameters.R_high_end[i],
                    cluster_order_parameters.Z_high_mean[i],
                    cluster_order_parameters.Z_high_end[i],
                    cluster_order_parameters.Rfreq_high_mean[i],
                    cluster_order_parameters.Rfreq_high_end[i],
                    cluster_order_parameters.R_low_mean[i],
                    cluster_order_parameters.R_low_end[i],
                    cluster_order_parameters.Z_low_mean[i],
                    cluster_order_parameters.Z_low_end[i],
                    cluster_order_parameters.Rfreq_low_mean[i],
                    cluster_order_parameters.Rfreq_low_end[i],
                ])
            end
            println(io, join(row, ","))
        end
    end
    return path
end

function safe_float_token(value::Real; digits::Integer = 8)
    rounded = round(Float64(value); digits = Int(digits))
    return replace(string(rounded), "-" => "m", "." => "p")
end

function k_step_restart_path(output_dir::AbstractString, filename_base::AbstractString, K::Real)
    restart_dir = joinpath(output_dir, filename_base * "_restart")
    mkpath(restart_dir)
    return joinpath(restart_dir, filename_base * "_K_$(safe_float_token(K)).jld2")
end

function save_k_step_restart(path::AbstractString, cfg::KContinuationConfig; jK, K, local_Ttr, omega, u0, final_state)
    JLD2.jldsave(path;
        K_index = Int(jK),
        K = Float64(K),
        C1 = cfg.C1,
        C2 = cfg.C2,
        N = cfg.N,
        Ttr = Int(local_Ttr),
        T = cfg.T,
        dt = cfg.dt,
        reltol = cfg.reltol,
        abstol = cfg.abstol,
        diagnostic_mode = normalize_diagnostic_mode(cfg.diagnostic_mode),
        amplitude_cluster_threshold = cfg.amplitude_cluster_threshold,
        amplitude_cluster_reltol = cfg.amplitude_cluster_reltol,
        amplitude_cluster_abstol = cfg.amplitude_cluster_abstol,
        omega = Float64.(omega),
        u0 = Float64.(u0),
        final_state = Float64.(final_state),
        omega_mean = cfg.omega_mean,
        omega_sigma = cfg.omega_sigma,
        omega_seed = cfg.omega_seed,
        omega_symmetry = normalize_omega_symmetry(cfg.omega_symmetry),
        recenter_omega = cfg.recenter_omega,
        initial_seed = cfg.initial_seed,
        continuation_noise = cfg.continuation_noise,
    )
    return path
end

function run_k_continuation(cfg::KContinuationConfig)
    mkpath(cfg.output_dir)
    cfg.diagnostic_mode === :online || error("Only online diagnostics are supported.")

    K_list = collect(range(cfg.K_start, cfg.K_end; length = cfg.K_steps))
    omega = Model.sample_eigenfrequencies(Random.MersenneTwister(cfg.omega_seed),
        Distributions.Normal(cfg.omega_mean, cfg.omega_sigma), cfg.N;
        omega_symmetry = cfg.omega_symmetry, recenter_nosym = cfg.recenter_omega)
    selected_indices = collect(eachindex(omega))
    Random.seed!(cfg.initial_seed)
    u0 = Model.generate_u0(cfg.N)
    noise_rng = Random.MersenneTwister(cfg.initial_seed + 1)

    system_mean_frequency = fill(NaN, length(K_list))
    mean_frequency = fill(NaN, length(K_list), length(selected_indices))
    min_radius = fill(NaN, length(K_list), length(selected_indices))
    max_radius = fill(NaN, length(K_list), length(selected_indices))
    R_mean = fill(NaN, length(K_list))
    R_end = fill(NaN, length(K_list))
    Z_mean = fill(NaN, length(K_list))
    Z_end = fill(NaN, length(K_list))
    Rfreq_mean = fill(NaN, length(K_list))
    Rfreq_end = fill(NaN, length(K_list))
    cluster_order_parameters = if cfg.amplitude_cluster_threshold === nothing
        nothing
    else
        (;
            N_high = fill(0, length(K_list)),
            N_low = fill(0, length(K_list)),
            R_high_mean = fill(NaN, length(K_list)),
            R_high_end = fill(NaN, length(K_list)),
            Z_high_mean = fill(NaN, length(K_list)),
            Z_high_end = fill(NaN, length(K_list)),
            Rfreq_high_mean = fill(NaN, length(K_list)),
            Rfreq_high_end = fill(NaN, length(K_list)),
            R_low_mean = fill(NaN, length(K_list)),
            R_low_end = fill(NaN, length(K_list)),
            Z_low_mean = fill(NaN, length(K_list)),
            Z_low_end = fill(NaN, length(K_list)),
            Rfreq_low_mean = fill(NaN, length(K_list)),
            Rfreq_low_end = fill(NaN, length(K_list)),
        )
    end
    if cluster_order_parameters !== nothing
        println(
            "Amplitude-cluster classification uses a light r_min prerun per K ",
            "(reltol = ", cfg.amplitude_cluster_reltol,
            ", abstol = ", cfg.amplitude_cluster_abstol, ").",
        )
        if cfg.diagnostic_mode == :online
            println("Online diagnostics: amplitude-cluster order parameters require one additional light streaming diagnostic pass per K.")
        end
    end
    final_states = Vector{Vector{Float64}}(undef, length(K_list))
    restart_files = String[]

    for (jK, K) in pairs(K_list)
        local_Ttr = jK == firstindex(K_list) ? cfg.Ttr_first : cfg.Ttr
        println("Processing K = ", K, " (", jK, " / ", length(K_list), ")")
        u0_for_step = copy(u0)

        step = run_k_step_online(cfg, K, omega, u0, local_Ttr; cluster_order_parameters = cluster_order_parameters)
        final_state = step.final_state
        system_mean_frequency[jK] = mean_without_nan(step.mean_frequency)
        mean_frequency[jK, :] .= step.mean_frequency[selected_indices]
        min_radius[jK, :] .= step.min_radius[selected_indices]
        max_radius[jK, :] .= step.max_radius[selected_indices]
        R_mean[jK] = step.R_mean
        R_end[jK] = step.R_end
        Z_mean[jK] = step.Z_mean
        Z_end[jK] = step.Z_end
        Rfreq_mean[jK] = step.Rfreq_mean
        Rfreq_end[jK] = step.Rfreq_end
        if cluster_order_parameters !== nothing
            clusters = step.clusters
            cluster_order_parameters.N_high[jK] = length(clusters.high_indices)
            cluster_order_parameters.N_low[jK] = length(clusters.low_indices)
            cluster_order_parameters.R_high_mean[jK] = step.cluster_stats.high.R_mean
            cluster_order_parameters.R_high_end[jK] = step.cluster_stats.high.R_end
            cluster_order_parameters.Z_high_mean[jK] = step.cluster_stats.high.Z_mean
            cluster_order_parameters.Z_high_end[jK] = step.cluster_stats.high.Z_end
            cluster_order_parameters.Rfreq_high_mean[jK] = step.cluster_stats.high.Rfreq_mean
            cluster_order_parameters.Rfreq_high_end[jK] = step.cluster_stats.high.Rfreq_end
            cluster_order_parameters.R_low_mean[jK] = step.cluster_stats.low.R_mean
            cluster_order_parameters.R_low_end[jK] = step.cluster_stats.low.R_end
            cluster_order_parameters.Z_low_mean[jK] = step.cluster_stats.low.Z_mean
            cluster_order_parameters.Z_low_end[jK] = step.cluster_stats.low.Z_end
            cluster_order_parameters.Rfreq_low_mean[jK] = step.cluster_stats.low.Rfreq_mean
            cluster_order_parameters.Rfreq_low_end[jK] = step.cluster_stats.low.Rfreq_end
        end

        final_states[jK] = final_state
        if cfg.savefiles
            push!(
                restart_files,
                save_k_step_restart(
                    k_step_restart_path(cfg.output_dir, cfg.filename_base, K),
                    cfg;
                    jK = jK,
                    K = K,
                    local_Ttr = local_Ttr,
                    omega = omega,
                    u0 = u0_for_step,
                    final_state = final_state,
                ),
            )
        end
        u0 = final_state .+ cfg.continuation_noise .* Random.randn(noise_rng, length(final_state))
    end

    base_path = joinpath(cfg.output_dir, cfg.filename_base)
    mean_csv = write_continuation_csv(base_path * "_mean_frequency.csv", K_list, mean_frequency, selected_indices, omega; system_mean_frequency = system_mean_frequency)
    min_r_csv = write_continuation_csv(base_path * "_min_radius.csv", K_list, min_radius, selected_indices, omega)
    max_r_csv = write_continuation_csv(base_path * "_max_radius.csv", K_list, max_radius, selected_indices, omega)
    order_parameter_csv = write_order_parameter_csv(
        base_path * "_order_parameters.csv",
        K_list,
        R_mean,
        R_end,
        Z_mean,
        Z_end,
        Rfreq_mean,
        Rfreq_end;
        cluster_order_parameters = cluster_order_parameters,
    )
    return (; K_list, omega, mean_frequency, min_radius, max_radius,
        R_mean, R_end, Z_mean, Z_end, Rfreq_mean, Rfreq_end,
        cluster_order_parameters, final_states, mean_csv, min_r_csv,
        max_r_csv, order_parameter_csv, restart_files)
end
