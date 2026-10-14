function fft_window(window, nt::Integer)
    window === nothing && return nothing
    idx = 0:(nt - 1)
    if window === :hann
        return 0.5 .- 0.5 .* cos.(2pi .* idx ./ (nt - 1))
    elseif window === :hamming
        return 0.54 .- 0.46 .* cos.(2pi .* idx ./ (nt - 1))
    elseif window === :blackman
        return 0.42 .- 0.5 .* cos.(2pi .* idx ./ (nt - 1)) .+ 0.08 .* cos.(4pi .* idx ./ (nt - 1))
    else
        error("Unknown window $(repr(window)). Use nothing, :hann, :hamming, or :blackman.")
    end
end

function fft_matrix(signals::AbstractMatrix, t; detrend::Bool, window, spectrum::Symbol, onesided::Bool, normalize::Bool)
    nt, nsignals = size(signals)
    length(t) == nt || error("Time vector length must match signal rows.")
    dtv = diff(t)
    dt_mean = sum(dtv) / length(dtv)
    w = fft_window(window, nt)
    wnorm = (w === nothing || !normalize) ? 1.0 : sum(w) / nt

    nfreqs = onesided ? fld(nt, 2) + 1 : nt
    freqs = onesided ? collect(0:(nfreqs - 1)) ./ (dt_mean * nt) :
            vcat(0:floor(Int, (nt - 1) / 2), ceil(Int, -nt / 2):-1) ./ (dt_mean * nt)
    spectra = zeros(Float64, nsignals, nfreqs)

    for j in 1:nsignals
        sig = collect(signals[:, j])
        detrend && (sig .-= sum(sig) / length(sig))
        w === nothing || (sig .*= w)
        F = FFTW.fft(sig)
        Fsel = onesided ? @view(F[1:nfreqs]) : F

        values = if spectrum === :amplitude
            abs.(Fsel) ./ (normalize ? nt * wnorm : nt)
        elseif spectrum === :power
            abs2.(Fsel) ./ (normalize ? nt^2 * wnorm^2 : nt^2)
        else
            error("spectrum must be :amplitude or :power. Got $(repr(spectrum)).")
        end

        if onesided
            if nt % 2 == 0
                length(values) > 2 && (values[2:end-1] .*= 2)
            else
                length(values) > 1 && (values[2:end] .*= 2)
            end
        end
        spectra[j, :] .= values
    end

    return (freqs = freqs, spectra = spectra, meta = (; detrend, window, spectrum, onesided, normalize, dt = dt_mean))
end
