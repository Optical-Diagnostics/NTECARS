"""
    N2.MultiTemperatureDistribution(; T_vib, T_rot, vib_states=N2.vib_states(10)) -> MultiTemperatureDistribution

Constructs a `MultiTemperatureDistribution` for N₂ assuming separate vibrational and
rotational temperatures.

# Constructor Arguments
- `T_vib::AbstractFloat`: Vibrational temperature in K.
- `T_rot::AbstractFloat`: Rotational temperature in K.
- `vib_states`: Vibrational states used to compute the partition sum.

# Fields of returned type
- `T_vib::AbstractFloat`: Vibrational temperature in K.
- `T_rot::AbstractFloat`: Rotational temperature in K.
- `Q::AbstractFloat`: Total partition sum (vibrational × rotational), computed at
  construction. Should not be set manually.
- `Q_rot::AbstractFloat`: Rotational partition sum, computed at construction.
  Should not be set manually.

# Notes
- For a system in full thermal equilibrium, set `T_vib = T_rot = T_gas`.

# Examples
```Julia
# Thermal equilibrium
dist = N2.MultiTemperatureDistribution(T_vib=1500.0, T_rot=1500.0)

# Non-equilibrium: vibrationally excited N₂
dist = N2.MultiTemperatureDistribution(T_vib=3000.0, T_rot=600.0)
```
"""
struct MultiTemperatureDistribution <: N2Distribution
    T_vib::AbstractFloat             # vibrational temperature
    T_rot::AbstractFloat             # rotational temperature
    Q    ::AbstractFloat             # total partition sum (vib + rot)
    Q_rot::AbstractFloat            # rotational partition sum

    function MultiTemperatureDistribution(;T_vib::AbstractFloat, T_rot::AbstractFloat, vib_states = N2.vib_states(16))
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