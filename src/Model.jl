module Model
using Distributions, Random, Statistics, LoopVectorization
using DSP: unwrap

@inline function unwrap_diff!(values::Vector{Float64})
    for k in eachindex(values)
        if values[k] > pi
            values[k] -= 2pi
        elseif values[k] < -pi
            values[k] += 2pi
        end
    end
    return values
end

function sl_harmo_hetero!(du, u, p, t)
    K, C1, C2, ωi = p
    len2 = length(u)
    N    = len2 >> 1

    sx = zero(eltype(u))
    sy = zero(eltype(u))
    @avx for i ∈ 1:2:len2
        sx += u[i]
        sy += u[i+1]
    end
    invN = one(eltype(u)) / N
    mx   = sx * invN
    my   = sy * invN

    @avx for j ∈ 1:N
        i1 = (j << 1) - 1
        xi = u[i1]
        yi = u[i1+1]

        Dx = mx - xi
        Dy = my - yi
        r2 = xi*xi + yi*yi
        c  = C2*r2
        ω = ωi[j]

        du[i1]     = xi - ω*yi - r2*xi + c*yi + K*(Dx - C1*Dy)
        du[i1 + 1] = yi + ω*xi- r2*yi - c*xi + K*(C1*Dx + Dy)
    end
    return nothing
end

generate_u0(N::Integer) = 2 .* rand(2*N) .- 1.1

function shuffle_eigenfrequencies!(rng, ωi::AbstractVector)
    for i in length(ωi):-1:2
        j = rand(rng, 1:i)
        ωi[i], ωi[j] = ωi[j], ωi[i]
    end
    return ωi
end

function sample_eigenfrequencies(rng, ω_dist::Distribution{Univariate,Continuous}, N::Integer;
        omega_symmetry::Symbol=:nosym,
        recenter_nosym::Bool=false)
    n = Int(N)
    n >= 0 || error("N must be non-negative")
    ω_mean = mean(ω_dist)

    if omega_symmetry === :nosym
        ωi = rand(rng, ω_dist, n)
        if recenter_nosym && n > 0
            ωi .-= mean(ωi)
            ωi .+= ω_mean
        end
        return ωi
    elseif omega_symmetry === :sym
        npairs = n ÷ 2
        offsets = abs.(rand(rng, ω_dist, npairs) .- ω_mean)
        ωi = Vector{Float64}(undef, n)
        ωi[1:npairs] .= ω_mean .+ offsets
        ωi[(npairs + 1):(2 * npairs)] .= ω_mean .- offsets
        if isodd(n)
            ωi[end] = ω_mean
        end
        return shuffle_eigenfrequencies!(rng, ωi)
    else
        error("Unknown omega_symmetry=$(omega_symmetry). Use :sym or :nosym.")
    end
end


end
