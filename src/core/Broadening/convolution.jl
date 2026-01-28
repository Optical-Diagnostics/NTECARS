function convolute(spectrum::Spectrum, g::Spectrum)
    Δν_spectrum = abs(spectrum.ν[2] - spectrum.ν[1])

    kernel = resample_to_resolution(g, Δν_spectrum).I
    kernel ./= sum(kernel)

    I = intensities(spectrum)
    N = length(I)
    M = length(kernel)
    start = div(M, 2) + 1

    ν_convoluted = wavenumbers(spectrum)
    I_convoluted = DSP.conv(I, kernel)[start : start + N - 1] #trim of padded addition
    return Spectrum(ν_convoluted, I_convoluted, :wavenumber)
end