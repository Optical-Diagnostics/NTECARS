struct CO2Transition
    initial::State
    final::State
    ν_fi::Float64
    αα::Float64

    function CO2Transition(initial_state, final_state)
        ν_fi = total_energy(final_state) - total_energy(initial_state)
        return new(initial_state, final_state, ν_fi, αα(initial_state, final_state))
    end
end

"""
    Generate all transitions from the given initial state to all possible upper states
    the result is a vector of vector where each vector contains the transitions of a band
"""
function transitions_from_vibrational_state(vib_qn_initial::VibrationalQuantumNumbers, J_max; iso_ID = :O16C12O16)
    transitions_dict::Dict{String, Vector{CO2Transition}} = Dict()
    
    # Q branch
    @unpack v₁, v₂, l₂, v₃, r = vib_qn_initial
    v₁_final = v₁ + 1
    
    for r_final in 1:v₁_final+1
        v₁′, v₂′, l₂′, v₃′, r′ =  v₁_final, v₂, l₂, v₃, r_final

        identifier = "$(v₁)$(v₂)$(l₂)$(v₃)$(r)->$(v₁′)$(v₂′)$(l₂′)$(v₃′)$(r′)"
        transitions_of_band = []

        for initial_rotational_state in CO2.allowed_rotational_states(VibrationalQuantumNumbers(v₁, v₂, l₂, v₃, r); J_max=J_max)
            J = initial_rotational_state.qn.J
            
            initial = State(QuantumNumbers(v₁,  v₂,  l₂,  v₃,  r,  J , CO2_parity(J, l₂, v₃)))
            final   = State(QuantumNumbers(v₁′, v₂′, l₂′, v₃′, r′, J,  CO2_parity(J, l₂′, v₃′)))
            push!(transitions_of_band, CO2Transition(initial, final))
        end
        transitions_dict[identifier] = transitions_of_band
    end
    return transitions_dict
end

function get_transition_database_CO2(v₁_max, v₂_max, v₃_max, J_max; iso_ID = :O16C12O16)
    # CO2 transitions
    vib_states_to_consider = vibrational_quantum_numbers_in_range(v₁_max, v₂_max, v₃_max)
    transitions = [transitions_from_vibrational_state(vib_qn, J_max) for vib_qn in vib_states_to_consider]
    transitions = merge(transitions...)
    return transitions
end

function print_transition(t::CO2Transition)
    println("----------------------------------------------------------")
    E = angular_frequency_to_wavenumber(t.ω_fi)
    println("energy: $(E)")
    println("initial state:")
    print_wavefunction(t.initial)
    println("final state:")
    print_wavefunction(t.final)
    println("polarizability = $(sqrt(αα(t)))")
    println("----------------------------------------------------------")
end