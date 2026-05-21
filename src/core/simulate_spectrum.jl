function simulate_spectrum(sim::CARSSimulator, ν_output = nothing)
    χ² = total_susceptibility(sim)
    
    if sim.lasers.stokes_profile isa Spectrum
        convolute_with_stokes_profile!(χ², sim.lasers)
    end
    
    if sim.instrument.profile isa Spectrum
        χ² = convolute(χ², sim.instrument.profile)
    end

    χ².ν = wavelength_to_wavenumber(wavelengths(χ²) .+ sim.wavelength_shift)

    if !isnothing(ν_output)
        χ² = average_to_detector_pixels(χ², ν_output)
    end

    χ².I = intensities(χ²) .* 1e100 .* sim.scaling_factor
    χ².I .+= sim.vertical_shift
    return χ²
end

function convolute_with_stokes_profile!(χ²::Spectrum, L::LaserConfiguration)
    stokes_interpol = LinearInterpolation(wavenumbers(L.stokes_profile), intensities(L.stokes_profile), extrapolation_bc = 0.0)
    ν_S = L.ν_1 + L.ν_2 .- wavenumbers(χ²)
    χ².I .*= stokes_interpol.(ν_S)
end

"""
    raman_shifts(; chi2, lasers) -> (ramanshift_ν1_νS, ramanshift_ν2_νS)

Computes the Raman shifts for the CARS-spectrum `chi2` relative to each
of the two laser wavenumbers.

# Arguments
- `chi2::Spectrum`: The CARS spectrum..
- `lasers::LaserConfiguration`: Laser configuration providing `ν_1` and `ν_2`
  in cm⁻¹.

# Returns
- `ramanshift_ν1_νS`: Raman shifts relative to the first laser,  `ν_1 - ν_S`, in cm⁻¹.
- `ramanshift_ν2_νS`: Raman shifts relative to the second laser, `ν_2 - ν_S`, in cm⁻¹.
"""
function raman_shifts(;chi2::Spectrum, lasers::LaserConfiguration)
    ramanshift_ν1_νS = wavenumbers(chi2) .- lasers.ν_2
    ramanshift_ν2_νS = wavenumbers(chi2) .- lasers.ν_1
    return ramanshift_ν1_νS, ramanshift_ν2_νS
end 
