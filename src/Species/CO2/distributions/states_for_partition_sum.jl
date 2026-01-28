function all_pure_vibrational_states(v₁_max, l₂_max, v₃_max)
    # return a vector containing all (rotationless) vibrational states, 
    # which are used for the calculation of the partition sums
    quantum_numbers    = vibrational_quantum_numbers_in_range(v₁_max, l₂_max, v₃_max)
    vibrational_states = [VibrationalState(qn) for qn in quantum_numbers]
    
    return vibrational_states
end

function vibrational_quantum_numbers_in_range(v₁_max, l₂_max, v₃_max)
    #returns a vector of all quantum numbers within the quantum number limits (in the AFGL-notation).
    vib_quantum_numbers::Vector{VibrationalQuantumNumbers} = []
    for (v₁, l₂, v₃) in Iterators.product(0:v₁_max, 0:l₂_max, 0:v₃_max)
        v₂              = l₂ # states with l₂ < v₂ are part of the polyad iwth different r
        N_polyad_states = v₁ + 1
        for r in 1:N_polyad_states
            push!(vib_quantum_numbers, VibrationalQuantumNumbers(v₁, v₂, l₂, v₃, r))
        end
    end
    return vib_quantum_numbers
end

function allowed_rotational_states(vib_qn::VibrationalQuantumNumbers; J_max = 100)
    #returns a vector of all states with the given quantum numbers that result in an allowed (symmetric state)
    allowed_states::Vector{State} = []
    for J in vib_qn.l₂:J_max
        qn = QuantumNumbers(vib_qn.v₁, vib_qn.v₂, vib_qn.l₂, vib_qn.v₃, vib_qn.r, J, CO2_parity(J, vib_qn.l₂, vib_qn.v₃))
        if !isnothing(allowed_term(qn))
            # for degenerate l₂ only add one state but with degeneracy 2
            push!(allowed_states, State(qn))
        end
    end
    return allowed_states
end

function database_of_allowed_states(;v₁_max = 2, l₂_max=4, v₃_max= 2, J_max = 100)
    vibrational_qns = vibrational_quantum_numbers_in_range(v₁_max, l₂_max, v₃_max)

    states = []
    for vibrational_qn in vibrational_qns
        @unpack v₁, v₂, l₂, v₃, r = vibrational_qn
        for rotational_state in CO2.allowed_rotational_states(VibrationalQuantumNumbers(v₁, v₂, l₂, v₃, r); J_max=J_max)
            push!(states, rotational_state)
        end
    end
    return states
end
