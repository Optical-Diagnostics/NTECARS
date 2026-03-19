import Base: +, -, *, /


"""
    Spectrum(x::Vector{N}, I::Vector{T}, unit::Symbol) where {N,T} -> Spectrum

Constructs a `Spectrum` from spectral positions `x` and intensities `I`.

# Constructor Arguments
- `x::Vector{N}`: Spectral positions at which the intensities are given.
  Units are specified by `unit`.
- `I::Vector{T}`: Intensities at spectral positions `x`.
- `unit::Symbol`: Units of `x`. Must be `:wavelength` (meters) or `:wavenumber` (cm⁻¹).

# Fields of returned type
- `ν::Vector{Float64}`: Spectral positions in wavenumbers (cm⁻¹), sorted from
  lowest to highest. Converted automatically if constructed with `:wavelength`.
- `I::Vector{T}`: Intensities corresponding to `ν`.

# Notes
- Spectral positions are always stored in wavenumbers internally, sorted from
  lowest to highest. If `x` is not sorted, the order is corrected automatically.
- Use `wavenumbers(spectrum)` and `intensities(spectrum)` to access the fields
  rather than accessing `ν` and `I` directly.

# Examples
```Julia
s = Spectrum([404e-9, 405e-9, 406e-9], [0.0, 1.0, 0.0], :wavelength)
s = Spectrum([20000.0, 21000.0], [0.5, 1.0], :wavenumber)
ν = wavenumbers(s)
λ = wavelengths(s)
I = intensities(s)
I = intensities(s, normalization = :maximum)
I = intensities(s, normalization = :sum)
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

"""
    cut_spectral_range(spectrum::Spectrum, range::Tuple{Float64, Float64}, unit::Symbol)

Extracts a sub-range of a spectrum by finding the closest grid points to the
specified range limits.

# Arguments
- `spectrum::Spectrum`: The spectrum to cut.
- `range::Tuple{Float64, Float64}`: The desired spectral range as `(min, max)`.
  For `:wavelength`, order does not matter as the tuple is sorted after conversion.
- `unit::Symbol`: Unit of `range`. Must be `:wavenumber` or `:wavelength` (in meters).

# Returns
- `Spectrum`: A new spectrum on a wavenumber grid, containing only the points
  within the specified range.

# Notes
- The cut is performed by nearest-neighbour lookup, so the returned range may
  differ slightly from the requested one if the grid spacing is coarse.
- The returned `Spectrum` is always in `:wavenumber` units, regardless of the
  input `unit`.

# Examples
```Julia
cut = cut_spectral_range(spectrum, (2200.0, 2400.0), :wavenumber)
cut = cut_spectral_range(spectrum, (4.0e-6, 4.5e-6), :wavelength)
```
"""
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
