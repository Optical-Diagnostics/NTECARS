function wavelength_to_frequency(λ)
    c = 299_792_458
    return c ./ λ
end

function wavelength_to_angular_frequency(λ)
    return 2 * π * wavelength_to_frequency(λ)
end

function wavelength_to_wavenumber(λ)
    return 1 ./ λ .* 0.01
end

function wavenumber_to_wavelength(ν)
    return 1 ./ ν .* 0.01
end

function frequency_to_wavelength(ν)
    c = 299_792_458
    return c ./ ν
end

function angular_frequency_to_wavelength(ω)
    return frequency_to_wavelength.(ω ./ (2*π))
end

function angular_frequency_to_wavenumber(ω)
    return wavelength_to_wavenumber.(frequency_to_wavelength.(ω ./ (2*π)))
end

function frequency_to_wavenumber(ν)
    return 1/frequency_to_wavelength(ν) * 0.01
end

function Raman_shift(λ₁, λ₂)
    ν₁ = wavelength_to_wavenumber(λ₁)
    ν₂ = wavelength_to_wavenumber(λ₂)
    return (ν₁ .- ν₂)
end

function antistokes_wavenumber_from_wavelengths(λ_pu, λ_S, λ_pr)
    ν_pu = wavelength_to_wavenumber(λ_pu)
    ν_pr = wavelength_to_wavenumber(λ_pr)
    ν_S  = wavelength_to_wavenumber(λ_S) 

    return @. ν_pu .- ν_S .+ ν_pr
end

function antistokes_wavelength(λ_pu, λ_S, λ_pr)
    ν_pu = wavelength_to_wavenumber(λ_pu)
    ν_pr = wavelength_to_wavenumber(λ_pr)
    ν_S  = wavelength_to_wavenumber(λ_S)

    return wavenumber_to_wavelength(ν_pu .- ν_S .+ ν_pr)
end

function anti_Stokes_to_Stokes_wavelength(λ_aS, λ_pu, λ_pr)
    ν_pu = wavelength_to_wavenumber(λ_pu)
    ν_pr = wavelength_to_wavenumber(λ_pr)
    ν_aS = wavelength_to_wavenumber(λ_aS)
    ν_S  = @. ν_pu - ν_aS + ν_pr
    return wavenumber_to_wavelength.(ν_S)
end


