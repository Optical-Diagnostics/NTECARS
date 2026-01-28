const LINEWIDTH_DATA = TOML.parsefile(CO2_LINEWIDTH_PARAMETERS_PATH)

function total_linewidth(m, T, p, molar_fractions_dict)
    linewidth = 0
    for species_name in keys(molar_fractions_dict)
        broadening_coefficients  = LINEWIDTH_DATA["broadening_coefficients"][species_name]
        temperature_coefficients = LINEWIDTH_DATA["temperature_coefficients"][species_name]
        molar_fraction           = molar_fractions_dict[species_name]
        
        linewidth += molar_fraction * parameterized_linewidth(m, T, p, broadening_coefficients, temperature_coefficients)
    end
    return linewidth
end

"""
    parameterized_linewidth(J_initial, T, p, broadening_coefficients, temperature_coefficients)

Calculates the HWHM at a pressure p and temperature T for a set of broadening and temperature coefficients from
L. Rosenmann, et al. Appl. Opt. 27, 3902-3907 (1988) 
table III

m is equal to J + 1 for R-branch lines, -J for P-branch lines, and J for Q-branch lines.
"""
function parameterized_linewidth(m, T, p, broadening_coefficients, temperature_coefficients)
    a0, a1, a2 = broadening_coefficients
    b0, b1, b2 = temperature_coefficients

    γ_300K = (a0 + a1*m + a2*m^2) * 1e-3
    N      = b0 + b1*m + b2*m^2
    pressure_scaling = p/100_000

    return pressure_scaling * γ_300K * (300/T)^N
end

