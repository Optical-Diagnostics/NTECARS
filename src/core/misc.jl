function molar_fractions_dictionary(species::Array{T}) where {T<:CARSSpecies}
    Dict([S.name => S.molar_fraction for S in species])
end

function linewidth_extrema_of_spectrum(;conditions::GasConditions, species::Array{T}) where {T<:CARSSpecies}
    T_gas, p = conditions.T_gas, conditions.pressure
    molar_fractions_dict = molar_fractions_dictionary(species)

    # calculate all linewidths
    linewidths_of_species(S)     = [transition_linewidth(S, t, T_gas, p, molar_fractions_dict) for t in transitions(S)]
    linewidths                   = vcat([linewidths_of_species(S) for S in species]...)
    min_linewidth, max_linewidth = extrema(linewidths)

    return min_linewidth, max_linewidth
end

function normalize_by_stokes_intensity!(χ²::Spectrum, wavelength_1, wavelength_2, stokes_spectrum::Spectrum)
    ν_1 = wavelength_to_wavenumber(wavelength_1)
    ν_2 = wavelength_to_wavenumber(wavelength_2)
    inverse_stokes_interpol = LinearInterpolation(wavenumbers(stokes_spectrum), 1 ./ intensities(stokes_spectrum, normalization = :maximum), extrapolation_bc = 0.0)

    ν_aS = wavenumbers(χ²)
    ν_S = ν_1 + ν_2 .- ν_aS
    χ².I = χ².I .* inverse_stokes_interpol.(ν_S)
end


