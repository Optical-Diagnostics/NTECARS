struct MultiTemperatureDistribution <: N2Distribution
    T_vib::AbstractFloat             # vibrational temperature
    T_rot::AbstractFloat             # rotational temperature
    Q    ::AbstractFloat             # total partition sum (vib + rot)
    Q_rot::AbstractFloat            # rotational partition sum

    function MultiTemperatureDistribution(;T_vib::AbstractFloat, T_rot::AbstractFloat, vib_states = N2.vib_states(10))
        Q_rot = rot_partition_sum(T_rot)
        Q     = sum(Boltzmann_factor(s.E_vib, T_vib) * Q_rot for s in vib_states)
        return new(T_vib, T_rot, Q, Q_rot)
    end
end

function (df::MultiTemperatureDistribution)(state::State)
    f_vib = Boltzmann_factor(state.E_vib, df.T_vib)
    f_rot = Boltzmann_factor(state.E_rot, df.T_rot)
    return state.degen/df.Q * f_vib * f_rot
end

function vibrational_population(df::MultiTemperatureDistribution, v)
    rotational_ground_state = State(v,0)
    f_vib = Boltzmann_factor(rotational_ground_state.E_vib, df.T_vib)
    return df.Q_rot/df.Q * f_vib
end