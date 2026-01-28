import Base: +, -, *, /


"""
    Spectrum(x::Vector{N}, I::Vector{T}, unit::Symbol) where {N,T}

create a `Spectrum` from spectral positions `x` with intensitiy `I`.

# Arguments
- `x::Vector{N}`: spectral positions at which the intensities are given. The units of `x` have to be specified by `unit=:wavelength` (meter) or `unit=:wavenumber` (cm^-1) to allow for spectra to be created
  from both.
- `I::Vector{T}`: intensities as spectral positions `x`.
- `unit::Symbol`: Can be `wavelength` or `wavenumber`. Specifies the units of `x`

# Notes
- Internally, the spectral positions are stored in wavenumber sorted from lowest to highest values.

# Example
```Julia
Spectrum([404e-9, 405e-9, 406e-9], [0.0, 1.0, 0.0], :wavelength)
Spectrum([20000, 21000], [0.5, 1.0], :wavenumber)
```
"""
mutable struct Spectrum{T}
    ν::Vector{Float64}
    I::Vector{T}

    function Spectrum(x::Vector{N}, I::Vector{T}, unit::Symbol) where {N,T}
        if unit == :wavenumber
            ν = x
        elseif unit == :wavelength
            ν = wavelength_to_wavenumber.(x)
        else
            error("Specturm created with undefined units. Optinos for units are: :wavelength , :wavenumber")
        end

        if issorted(ν)
            return new{T}(ν, I)
        end

        # if they are not properly ordered, sort
        order = sortperm(ν)
        return new{T}(ν[order], I[order])
    end
end

Base.abs(s::Spectrum)  = Spectrum(wavenumbers(s), abs.(intensities(s)),  :wavenumber)
Base.abs2(s::Spectrum) = Spectrum(wavenumbers(s), abs2.(intensities(s)), :wavenumber)
Base.conj(s::Spectrum) = Spectrum(wavenumbers(s), conj.(intensities(s)), :wavenumber)
Base.sqrt(s::Spectrum) = Spectrum(wavenumbers(s), sqrt.(intensities(s)), :wavenumber)
Base.real(s::Spectrum) = Spectrum(wavenumbers(s), real.(intensities(s)), :wavenumber)

function +(s1::Spectrum, s2::Spectrum)
    if wavenumbers(s1) == wavenumbers(s2)
        return Spectrum(wavenumbers(s1), intensities(s1) .+ intensities(s2), :wavenumber) 
    else
        error("Spectra with unequal wavenumber arrays added")
    end
end

function *(s1::Spectrum, s2::Spectrum)
    if wavenumbers(s1) == wavenumbers(s2)
        return Spectrum(wavenumbers(s1), intensities(s1) .* intensities(s2), :wavenumber) 
    else
        error("Spectra with unequal wavenumber arrays multiplied")
    end
end


function wavenumbers(s::Spectrum)
    return s.ν
end

function wavelengths(s::Spectrum)
    return wavenumber_to_wavelength(s.ν)
end

function intensities(s::Spectrum; normalization = nothing)
    if isnothing(normalization)
        return s.I
    elseif normalization == :sum
        return s.I ./ sum(s.I)
    elseif normalization == :maximum
        return s.I ./ maximum(s.I)
    else
        error("normalization argument for intensities invalid. Options are: nothing, :maximum, :sum")
    end
end

Base.copy(s::Spectrum) = Spectrum(copy(s.ν), copy(s.I), :wavenumber)

function cut_spectral_range(spectrum::Spectrum, range::Tuple{Float64, Float64}, unit)
    if unit ==:wavenumber
        νmin, νmax = range
    elseif unit == :wavelength
        νmin, νmax = sort(wavelength_to_wavenumber.(range))
    else
        error("Unsupported unit for cutting spectral range. Options are :wavelength, :wavenumber")
    end

    i1 = argmin(abs.(νmin .- wavenumbers(spectrum)))
    i2 = argmin(abs.(νmax .- wavenumbers(spectrum)))
    return Spectrum(wavenumbers(spectrum)[i1:i2], intensities(spectrum)[i1:i2], :wavenumber)
end
