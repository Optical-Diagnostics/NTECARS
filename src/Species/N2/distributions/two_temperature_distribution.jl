"""
    N2.TwoTemperatureDistribution(; T_vib_cold, T_vib_hot, Rh, T_rot,
                                    vib_states=N2.vib_states(10)) -> TwoTemperatureDistribution

Constructs a `TwoTemperatureDistribution` for N₂ representing a bimodal vibrational
population consisting of a cold and a hot sub-population. 
(ref: J Kuhfeld et al 2021 J. Phys. D: Appl. Phys. 54 305205)

# Constructor Arguments
- `T_vib_cold::AbstractFloat`: Vibrational temperature of the cold population in K.
- `T_vib_hot::AbstractFloat`: Vibrational temperature of the hot population in K.
- `Rh::AbstractFloat`: Fraction of molecules in the hot vibrational population ∈ [0, 1].
  The cold population fraction is `1 - Rh`.
- `T_rot::AbstractFloat`: Rotational temperature in K, shared by both populations.
- `vib_states`: Vibrational states used to compute the partition sums. Defaults to
  `N2.vib_states(10)`, which includes the 10 lowest vibrational levels.

# Fields of returned type
- `T_vib_cold::AbstractFloat`: Vibrational temperature of the cold population in K.
- `T_vib_hot::AbstractFloat`: Vibrational temperature of the hot population in K.
- `Rh::AbstractFloat`: Fraction of molecules in the hot vibrational population.
- `T_rot::AbstractFloat`: Rotational temperature in K.
- `Q_cold::AbstractFloat`: Partition sum of the cold population, computed at construction.
- `Q_hot::AbstractFloat`: Partition sum of the hot population, computed at construction.
- `Q_rot::AbstractFloat`: Rotational partition sum, computed at construction.

# Examples
```Julia
# 10% of molecules in a hot vibrational population
dist = N2.TwoTemperatureDistribution(
    T_vib_cold = 500.0,
    T_vib_hot  = 5000.0,
    Rh         = 0.1,
    T_rot      = 600.0
)
```
"""
struct TwoTemperatureDistribution <: N2Distribution
    T_vib_cold::AbstractFloat         
    T_vib_hot::AbstractFloat
    Rh::AbstractFloat                # fraction in hot population
    T_rot::AbstractFloat             # rotational temperature
    Q_cold::AbstractFloat
    Q_hot::AbstractFloat
    Q_rot::AbstractFloat             # rotational partition sum

    function TwoTemperatureDistribution(;T_vib_cold::AbstractFloat, T_vib_hot::AbstractFloat,
                                        Rh::AbstractFloat, T_rot::AbstractFloat, vib_states = N2.vib_states(10))
        Q_rot = rot_partition_sum(T_rot)
        Q_cold = sum(Boltzmann_factor(s.E_vib, T_vib_cold) * Q_rot for s in vib_states)
        Q_hot  = sum(Boltzmann_factor(s.E_vib, T_vib_hot) * Q_rot for s in vib_states)
        return new(T_vib_cold, T_vib_hot, Rh, T_rot, Q_cold, Q_hot, Q_rot)
    end
end

function (df::TwoTemperatureDistribution)(state::State)
    f_vib = Boltzmann_factor(state.E_vib, df.T_vib)
    f_rot = Boltzmann_factor(state.E_rot, df.T_rot)
    return state.degen/df.Q * f_vib * f_rot
end

function fvib(df::TwoTemperatureDistribution, state::State)
    f_vib = (1-df.Rh)/df.Q_cold * Boltzmann_factor(state.E_vib, df.T_vib_cold) + 
            df.Rh/df.Q_hot * Boltzmann_factor(state.E_vib, df.T_vib_hot) 
    return f_vib
end

function vibrational_population(df::TwoTemperatureDistribution, v)
    f_vib = fvib(df, State(v,0))
    return df.Q_rot * f_vib
end