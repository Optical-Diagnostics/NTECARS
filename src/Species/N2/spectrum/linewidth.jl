struct MEGL_constants
    α::Float64
    β::Float64
    δ::Float64
    n::Float64
    T₀::Float64
end

# pre-load and calculate
const MEGL_DATA = TOML.parsefile(LINEWIDTH_PARAMETERS_PATH)
const energies_buffer= [energy(QuantumNumbers(0,i)) for i in 0:100] .- energy(QuantumNumbers(0,0))


function γ_up(i, j, Tgas, p, C::MEGL_constants)
    # upward: i → j, j > i
    @unpack α, β, δ, n, T₀ = C

    E_i = energies_buffer[i+1] * cm⁻¹_to_J 
    E_j = energies_buffer[j+1] * cm⁻¹_to_J 
    ΔE = (E_j - E_i)   # always positive

    factor = ((1 + 1.5*E_i/(k_B*Tgas*δ)) / (1 + 1.5*E_i/(k_B*Tgas)))^2
    return α*p*(T₀/Tgas)^n * factor * exp(-β * ΔE/(k_B*Tgas))
end

function γ_down(i, j, Tgas, p, C::MEGL_constants)
    # downward: j → i, i < j
    @unpack α, β, δ, n, T₀ = C

    E_i = energies_buffer[i+1] * cm⁻¹_to_J 
    E_j = energies_buffer[j+1] * cm⁻¹_to_J 
    ΔE = (E_j - E_i)  # positive upward gap

    return γ_up(i, j, Tgas, p, C) * (2*i+1)/(2*j+1) * exp(ΔE/(k_B*Tgas))
end

function MEGL_linewidth(J, Tgas, p, C::MEGL_constants, J_max=70)
    γ = 0.0

    # downward transitions: J → K < J
    for K in 0:1:J-1
        γ += γ_down(K, J, Tgas, p, C)
    end

    # upward transitions: J → K > J
    for K in J+1:1:J_max
        γ += γ_up(J, K, Tgas, p, C)
    end
    
    return γ
end

function linewidth_constribution(Ji, Jf, Tgas, p, C::MEGL_constants)
    γ_fi = 1/2 * (
        MEGL_linewidth(Ji, Tgas, p, C) + 
        MEGL_linewidth(Jf, Tgas, p, C)
    )
    return γ_fi / 2 # /2 due to summing all J not every second
end



function total_linewidth(Ji, Jf, Tgas, p, molar_fractions_dict)
    γ = 0
    for (species_name, molar_fraction) in molar_fractions_dict
        MEGL_parameters     = copy(MEGL_DATA["MEGL_constants"][species_name])
        MEGL_parameters[1]  = MEGL_parameters[1]/atm_to_Pa # convert alpha to K/Pa
        linewidth_parameters = MEGL_constants(MEGL_parameters...)
        
        γ += molar_fraction * linewidth_constribution(Ji, Jf, Tgas, p, linewidth_parameters)
    end
    return γ
end