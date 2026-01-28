
mutable struct N2Species  <: CARSSpecies
    molar_fraction ::Float64
    χ_non_resonant ::Float64
    transitions    ::Dict{String, Vector{N2Transition}}
    distribution   ::N2Distribution
    use_collisional_narrowing ::Bool
    name           ::String

    function N2Species(;
        molar_fraction::Float64,
        v_max::Int = 4, 
        J_max::Int = 50, 
        distribution::Union{D}, 
        chi_non_resonant = 4.42e-51,
        use_collisional_narrowing = true) where {D<:N2Distribution}

        transitions = N2.get_transition_database(v_max, J_max)
        return new(molar_fraction, chi_non_resonant, transitions, distribution, use_collisional_narrowing, "N2")
    end
end

function transition_polarizability(::N2Species , transition::N2Transition)
    transition.αα
end

function transition_linewidth(::N2Species, transition::N2Transition, T::Number, p::Number, molar_fractions_dict)
    N2.total_linewidth(transition.initial.qn.J, transition.final.qn.J, T, p, molar_fractions_dict)
end

function transition_line_position(transition::N2Transition)
    return transition.ν_fi
end

function transitions_band_dictionary(N2::N2Species)
    N2.transitions
end

function transitions(N2::N2Species)
    vcat(collect(values(N2.transitions))...)
end

function non_resonant_susceptibility(N2::N2Species)
    return N2.χ_non_resonant
end

function initial_state(t::N2Transition)
    return t.initial
end

function final_state(t::N2Transition)
    return t.final
end

function rovibrational_population(N2::N2Species, state)
    return N2.distribution(state)
end

function use_collisional_narrowing(N2::N2Species)
    return N2.use_collisional_narrowing
end

function molar_fraction(N2::N2Species)
    return N2.molar_fraction
end

#####################################################################################################################
#                          Functions that are not needed for the calculation of spectra
#####################################################################################################################
"""
    Generate vibrational and rovibrational populations for a distribution function. If any parameter has
    an uncertainty using measurements.jl, then the results will have the corresponding uncertainties
"""
function rovibrational_populations_database(N2_molecule::N2Species) 
    all_states_simulated = []
    for t in transitions(N2_molecule)
        push!(all_states_simulated, t.initial)
        push!(all_states_simulated, t.final)
    end
    all_states_simulated = unique(all_states_simulated)

    # group them by the same vibrational quantum numbers
    grouped_NvJ = Dict()
    for s in all_states_simulated
        v = s.qn.v
        if !(v in keys(grouped_NvJ))
            f = rovibrational_population(N2_molecule, s)
            grouped_NvJ[v] = [s,f]
        else
            f = rovibrational_population(N2_molecule, s)
            grouped_NvJ[v] = hcat(grouped_NvJ[v], [s,f])
        end
    end

    grouped_Nv = Dict()
    for v in keys(grouped_NvJ)
        f = sum(grouped_NvJ[v][2,:])
        E_vib = N2.energy(N2.QuantumNumbers(v,0))
        grouped_Nv[v] = (E_vib, f)
    end
    return grouped_NvJ, grouped_Nv
end