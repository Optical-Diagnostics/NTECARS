"""
    CO2Species(; molar_fraction, distribution, v_max=(1,2,1), J_max=100,
                 chi_non_resonant=4.5e-51, use_collisional_narrowing=true) -> CO2Species

Generates the transition database and then constructs a `CO2Species` and returns it.

# Constructor Arguments
- `molar_fraction::Float64`: Molar fraction of CO₂ in the gas mixture.
- `distribution::CO2Distribution`: Rovibrational population distribution, e.g.
  [`CO2.MultiTemperatureDistribution`](@ref).
- `v_max::Tuple{Int64,Int64,Int64}`: Maximum vibrational quantum numbers
  `(v1_max, v2_max, v3_max)` used to truncate the transition database. Higher.
  values increase accuracy at high temperatures at the cost of computation time.
- `J_max::Int`: Maximum rotational quantum number in the transition database.
- `chi_non_resonant`: Non-resonant susceptibility in m²V²/m³.
- `use_collisional_narrowing::Bool`: Whether to apply collisional narrowing.

# Fields of returned type
- `molar_fraction::Float64`: Molar fraction of CO₂ in the gas mixture.
- `χ_non_resonant::Float64`: Non-resonant susceptibility in m²V²/m³.
- `transitions::Dict{String, Vector{CO2Transition}}`: Transition database generated
  from `v_max` and `J_max` at construction time.
- `distribution::CO2Distribution`: Rovibrational population distribution.
- `use_collisional_narrowing::Bool`: Whether collisional narrowing is applied.
- `name::String`: Species identifier, always `"CO2"`. Required internally for
  constructing molar fraction dictionaries.

# Notes
- The transition database is computed once at construction via
  `CO2.get_transition_database_CO2` and stored in `transitions`. Increasing
  `v_max` can cause significant comuptation times since large hamiltonians
  have to be solved.
- `name` is always set to `"CO2"` and cannot be changed.

# Examples
```Julia
co2 = CO2Species(
    molar_fraction = 0.1,
    distribution   = CO2.MultiTemperatureDistribution(T_12 = 300.0, T_3 = 1800.0, T_rot = 300.0),
    v_max          = (1, 2, 1),
)
```
"""
mutable struct CO2Species  <: CARSSpecies
    molar_fraction            ::Float64
    χ_non_resonant            ::Float64
    transitions               ::Dict{String, Vector{CO2Transition}}
    distribution              ::CO2Distribution
    use_collisional_narrowing ::Bool
    name                      ::String #HAS TO BE CO2, is needed for creating dictionary with molar fractions

    function CO2Species(;
        molar_fraction::Float64,
        v_max::Tuple{Int64, Int64, Int64} = (1,2,1),
        J_max::Int = 100, 
        distribution::Union{D},
        chi_non_resonant = 4.5e-51,
        use_collisional_narrowing = true
        ) where {D<:CO2Distribution}

        transitions = CO2.get_transition_database_CO2(v_max..., J_max)
        return new(molar_fraction, chi_non_resonant, transitions, distribution, use_collisional_narrowing, "CO2")
    end
end

function transition_polarizability(::CO2Species , transition::CO2Transition)
    transition.αα
end

function transition_linewidth(::CO2Species, transition::CO2Transition, T::Number, p::Number, molar_fractions_dict)
    CO2.total_linewidth(transition.initial.qn.J, T, p, molar_fractions_dict)
end

function transition_line_position(transition::CO2Transition)
    return transition.ν_fi
end

function transitions_band_dictionary(CO2::CO2Species)
    CO2.transitions
end

function transitions(CO2::CO2Species)
    vcat(collect(values(CO2.transitions))...)
end

function non_resonant_susceptibility(CO2::CO2Species)
    return CO2.χ_non_resonant
end

function initial_state(t::CO2Transition)
    return t.initial
end

function final_state(t::CO2Transition)
    return t.final
end

function rovibrational_population(CO2::CO2Species, state)
    return CO2.distribution(state)
end

function use_collisional_narrowing(CO2::CO2Species)
    return CO2.use_collisional_narrowing
end

function molar_fraction(CO2::CO2Species)
    return CO2.molar_fraction
end

#####################################################################################################################
#                          Functions that are not needed for the calculation of spectra
#####################################################################################################################
"""
    Generate vibrational and rovibrational populations for a distribution function. If any parameter has
    an uncertainty using measurements.jl, then the results will have the corresponding uncertainties
"""
function rovibrational_populations_database(CO2_molecule::CO2Species) 
    # get quantumnumbers_of
    all_states_simulated = []
    for t in transitions(CO2_molecule)
        push!(all_states_simulated, t.initial)
        push!(all_states_simulated, t.final)
    end
    all_states_simulated = unique(all_states_simulated)

    # group them by the same vibrational quantum numbers
    grouped_NvJ = Dict()
    for s in all_states_simulated
        vib_qn = CO2.VibrationalQuantumNumbers(s.qn)
        if !(vib_qn in keys(grouped_NvJ))
            f = rovibrational_population(CO2_molecule, s)
            grouped_NvJ[vib_qn] = [s,f]
        else
            f = rovibrational_population(CO2_molecule, s)
            grouped_NvJ[vib_qn] = hcat(grouped_NvJ[vib_qn], [s,f])
        end
    end

    grouped_Nv = Dict()
    for vib_qn in keys(grouped_NvJ)
        f = sum(grouped_NvJ[vib_qn][2,:])
        E_vib = CO2.PartitionedEnergies(CO2.QuantumNumbers(vib_qn.v₁, vib_qn.v₂, vib_qn.l₂, vib_qn.v₃, vib_qn.r, 0, 1)).E_vib
        grouped_Nv[vib_qn] = (E_vib, f)
    end
    return grouped_NvJ, grouped_Nv
end