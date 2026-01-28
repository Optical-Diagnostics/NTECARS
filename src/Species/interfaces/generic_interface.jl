# These functions have to be defined by each species 

function transition_linewidth(species::CARSSpecies, transition, T::Number, p::Number, molar_fractions_dict)
    throw(MethodError(transition_linewidth, (species, transition, T, p)))
end

function transition_polarizability(species::CARSSpecies, transition)
    throw(MethodError(transition_polarizability, (species, transition)))
end

function non_resonant_susceptibility(species::CARSSpecies)
    throw(MethodError(non_resonant_susceptibility, (species,)))
end

function initial_state(transition)
    throw(MethodError(initial_state, (transition)))
end

function final_state(transition)
    throw(MethodError(final_state, (transition)))
end

function rovibrational_population(distribution, state)
    throw(MethodError(rovibrational_population, (distribution, state)))
end

function use_collisional_narrowing(species::CARSSpecies)
    throw(MethodError(use_collisional_narrowing, (species,)))
end

function transition_line_position(transition)
    throw(MethodError(transition_line_position, (transition)))
end

function transitions_band_dictionary(species::CARSSpecies)
    throw(MethodError(transitions, (species,)))
end

function transitions(species::CARSSpecies)
    throw(MethodError(transitions, (species,)))
end

function molar_fraction(species::CARSSpecies)
    throw(MethodError(molar_fraction, (species,)))
end



