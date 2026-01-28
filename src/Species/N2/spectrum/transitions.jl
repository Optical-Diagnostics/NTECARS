struct N2Transition
    initial::State
    final  ::State
    ν_fi   ::Float64
    αα     ::Float64

    function N2Transition(initial::State, final::State)
        ν_fi = total_energy(final) - total_energy(initial)
        return new(initial, final, ν_fi, αα(initial, final) )
    end
end

"""
Get array with transitions connected to the state (depending on J it can be Q,O and P-branch)
"""
function transitions_from_state(initial::State)
    Q_branch = N2Transition(initial, State(v(initial)+1, J(initial)))
    S_branch = N2Transition(initial, State(v(initial)+1, J(initial)+2))

    transitions::Vector{N2Transition} = []

    if v(initial) == 0
        if J(initial) < 2
            transitions = [Q_branch, S_branch]
        else
            O_branch = N2Transition(initial, State(v(initial)+1, J(initial)-2))
            transitions =  [Q_branch, S_branch, O_branch]
        end
    else
        transitions = [Q_branch]
    end
    return transitions
end

function get_transition_database(v_max::Int, J_max::Int)
    transitions_dict::Dict{String, Vector{N2Transition}} = Dict()
    for v in 0:v_max
        band_identifier     = "$(v)->$(v+1)"
        transitions_of_band::Vector{N2Transition} = []
        for J in 0:J_max
            s = State(v,J)
            append!(transitions_of_band, transitions_from_state(s))
        end
        transitions_dict[band_identifier] = transitions_of_band
    end
    return transitions_dict
end

