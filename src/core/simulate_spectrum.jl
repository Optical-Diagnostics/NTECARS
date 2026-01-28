function simulate_spectrum(sim::CARSSimulator, ν_output = nothing)
    χ² = total_susceptibility(sim)
    
    if sim.lasers.stokes_profile isa Spectrum
        convolute_with_stokes_profile!(χ², sim.lasers)
    end
    
    if sim.instrument.profile isa Spectrum
        χ² = convolute(χ², sim.instrument.profile)
    end

    if !isnothing(ν_output)
        χ² = average_to_detector_pixels(χ², ν_output)
    end

    χ².I = intensities(χ², normalization = :maximum)
    χ².I .+= sim.vertical_shift
    return χ²
end

function convolute_with_stokes_profile!(χ²::Spectrum, L::LaserConfiguration)
    stokes_interpol = LinearInterpolation(wavenumbers(L.stokes_profile), intensities(L.stokes_profile), extrapolation_bc = 0.0)
    ν_S = L.ν_1 + L.ν_2 .- wavenumbers(χ²)
    χ².I .*= stokes_interpol.(ν_S)
end
"""


"""
function raman_shifts(;chi2::Spectrum, lasers::LaserConfiguration)
    raman_shift_1 = wavenumbers(chi2) .- lasers.ν_2
    raman_shift_2 = wavenumbers(chi2) .- lasers.ν_1
    return raman_shift_1, raman_shift_2
end 
