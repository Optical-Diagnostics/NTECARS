#include("states_for_partition_sum.jl")

# pre calculated states of CO2 for quicker calculation of the CO2 partition function.
USED_ISOS = [:O16C12O16]
v₁_max, l₂_max, v₃_max = 4,8,3#6,10,4
PRECALCULATED_VIBSTATES = Dict([iso => all_pure_vibrational_states( v₁_max, l₂_max, v₃_max) for iso in USED_ISOS])
PRECALCULATED_VIBSTATES_DICT = Dict([iso => Dict([state.qn => state for state in PRECALCULATED_VIBSTATES[iso]]) for iso in USED_ISOS])
