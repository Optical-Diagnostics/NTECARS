"""
Bins a finely resolved spectrum defined at wavelengths ν with intensity I to the wavelengths ν_pixel
"""
function average_to_detector_pixels(S::Spectrum, ν_pixel)
    itp = Interpolations.LinearInterpolation(wavenumbers(S), intensities(S), extrapolation_bc=Line())
    S_pixel = Spectrum(ν_pixel, itp.(ν_pixel), :wavenumber)
    
    Δν = abs(ν_pixel[2] - ν_pixel[1])
    for i in eachindex(ν_pixel)
        ν_left  = ν_pixel[i] - 1/2 * Δν
        ν_right = ν_pixel[i] + 1/2 * Δν
        S_pixel.I[i] = mean(itp.(LinRange(ν_left,ν_right,10)))
    end
    
    return S_pixel
end

"""
    Samples the given spectrum at a given resolution via interpolation. Used for convolutions with
    the instrumental broadening profile.
"""
function resample_to_resolution(spec::Spectrum, Δν_ref)
    ν, I = wavenumbers(spec), intensities(spec)
    ν_min, ν_max = extrema(ν)

    # Generate symmetric sampling grid around 0 with step Δν_ref
    ν_sample = collect(ν_min:Δν_ref:ν_max)

    # Interpolate original spectrum onto new grid
    itp = interpolate((ν,), I, Gridded(Linear()))
    I_resampled = itp.(ν_sample)

    return Spectrum(ν_sample, I_resampled, :wavenumber)
end