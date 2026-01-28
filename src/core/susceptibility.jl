"""
    total_susceptibility(sim::CARSSimulator)

test2
"""
function total_susceptibility(sim::CARSSimulator)
    molar_fractions_dict = Dict([species.name => molar_fraction(species) for species in sim.species])

    ν_aS = sim.grid.ν_aS
    N    = length(ν_aS)
    χ1   = Spectrum(ν_aS, zeros(Complex, N), :wavenumber)
    χ2   = Spectrum(ν_aS, zeros(Complex, N), :wavenumber)
    ν_raman_1 = ν_aS .- sim.lasers.ν_1
    ν_raman_2 = ν_aS .- sim.lasers.ν_2

    for species in sim.species
        χ1.I .+= molar_fraction(species) .* (χ_CARS(
            ν_raman   = ν_raman_1,
            species   = species,
            T_gas     = sim.conditions.T_gas,
            p         = sim.conditions.pressure,
            molar_fractions_dict = molar_fractions_dict,
            buffer = sim.grid.buffer
        ) .+ non_resonant_susceptibility(species))

        χ2.I .+= molar_fraction(species) .* (χ_CARS(
            ν_raman   = ν_raman_2,
            species   = species,
            T_gas     = sim.conditions.T_gas,
            p         = sim.conditions.pressure,
            molar_fractions_dict = molar_fractions_dict,
            buffer = sim.grid.buffer
        ) .+ non_resonant_susceptibility(species))
    end

    if sim.grid isa AdaptiveGrid
        χ1 = interpolate_to_uniform(χ1, sim.grid.Δν_fine)
        χ2 = interpolate_to_uniform(χ2, sim.grid.Δν_fine)
    end

    g1 = sim.lasers.profile_1
    g2 = sim.lasers.profile_2

    if g1 isa DeltaProfile && g2 isa DeltaProfile
        return abs2(χ1 + χ2)
    elseif g1 isa Spectrum && g2 isa Spectrum
        spectrum = convolute(abs2(χ1), g1) + convolute(abs2(χ2), g2) + convolute(conj(χ1), g1) * convolute(χ2, g2) + convolute(χ1, g1) * convolute(conj(χ2), g2)
        return abs(spectrum) # just a change of datatypes since the complex part is zero now
    elseif g1 isa DeltaProfile && g2 isa Spectrum
        spectrum = abs2(χ1) + convolute(abs2(χ2), g2) + conj(χ1) * convolute(χ2, g2) + χ1 * convolute(conj(χ2), g2)
        return abs(spectrum) # just a change of datatypes since the complex part is zero now
    elseif g1 isa Spectrum && g2 isa DeltaProfile
        spectrum = abs2(χ2) + convolute(abs2(χ1), g1) + conj(χ2) * convolute(χ1, g1) + χ2 * convolute(conj(χ1), g1)
        return abs(spectrum) # just a change of datatypes since the complex part is zero now
    end
end


function χ_CARS(;
    ν_raman::Vector{Float64},
    T_gas::Number,
    p::Number,
    species::CARSSpecies,
    molar_fractions_dict,
    buffer)

    max_νfi = ν_raman[end]
    min_νfi = ν_raman[1]

    N_spec = length(ν_raman)

    buffer.χ      .= 0.0
    buffer.χ_real .= 0.0
    buffer.χ_imag .= 0.0
    buffer.narrowing_function_real .= 1.0
    buffer.narrowing_function_imag .= 0.0
    buffer.σ_real .= 0.0
    buffer.σ_imag .= 0.0

    transitions = transitions_band_dictionary(species)
    for band in keys(transitions)
        for t in transitions[band]
            ν_fi = transition_line_position(t)
            if ν_fi > max_νfi || ν_fi < min_νfi
                continue
            end
            ρ_i = rovibrational_population(species, initial_state(t))
            ρ_f = rovibrational_population(species, final_state(t))
            Δρ   = ρ_i - ρ_f

            γ_fi = transition_linewidth(species, t, T_gas, p, molar_fractions_dict)

            a = 1.0 / (ε_0 * h * c * cm⁻¹_to_m⁻¹) * Δρ * transition_polarizability(species, t)

            lorentz_line_profile_imag!(buffer.σ_imag, ν_raman, ν_fi, γ_fi)
            lorentz_line_profile_real!(buffer.σ_real, ν_raman, ν_fi, γ_fi)

            #apply_narrowing!(buffer, ρ_i, γ_fi, a)
            @tturbo @. buffer.χ_real += a * buffer.σ_real
            @tturbo @. buffer.χ_imag += a * buffer.σ_imag

            if use_collisional_narrowing(species)
                # split into real and imaginary: 1 + i*∑ ρ * γ_fi * (σ_real + i * σ_imag) = 1 - ∑ρ*γ_fi*σ_imag +  i*∑ ρ * γ_fi * σ_real
                @tturbo buffer.narrowing_function_real .+= -ρ_i .* γ_fi .* buffer.σ_imag
                @tturbo buffer.narrowing_function_imag .+= ρ_i .* γ_fi .* buffer.σ_real
            end
        end

        complex_div_and_add!(
            buffer.χ,
            buffer.χ_real,
            buffer.χ_imag,
            buffer.narrowing_function_real,
            buffer.narrowing_function_imag
        )
        # slower version:
        #buffer.χ += Complex(buffer.χ_real, buffer.χ_imag) / Complex(buffer.narrowing_function_real, buffer.narrowing_function_imag)
        
        
        buffer.χ_real .= 0.0
        buffer.χ_imag .= 0.0
        buffer.narrowing_function_real .= 1.0
        buffer.narrowing_function_imag .= 0.0
    end

    return buffer.χ
end

function complex_div_and_add!(χ::Vector{Complex},
                      σr::Vector{Float64}, σi::Vector{Float64},
                      ηr::Vector{Float64}, ηi::Vector{Float64})

    # faster version of χ += (σr +i*σi) / (ηr +i*ηi)
    Threads.@threads for i in eachindex(χ)
        denom = ηr[i]^2 + ηi[i]^2
        χ[i] += (σr[i]*ηr[i] + σi[i]*ηi[i]) / denom + 1im*(σi[i]*ηr[i] - σr[i]*ηi[i]) / denom
    end
end


function lorentz_line_profile_real!(σ, ν_raman, ν_if, γ_fi)
    @tturbo @. σ = (ν_if - ν_raman) / ((ν_if - ν_raman)^2 + γ_fi^2)
end

function lorentz_line_profile_imag!(σ, ν_raman, ν_if, γ_fi)
    @tturbo @. σ = γ_fi / ((ν_if - ν_raman)^2 + γ_fi^2)
end