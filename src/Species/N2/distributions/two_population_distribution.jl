"""
    N2.TwoPopulationDistribution(; T_vib_cold, T_vib_hot, Rh, T_rot,
                                    vib_states=N2.vib_states(10)) -> TwoPopulationDistribution

Constructs a `TwoPopulationDistribution` for N₂ representing a bimodal vibrational
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
dist = N2.TwoPopulationDistribution(
    T_vib_cold = 500.0,
    T_vib_hot  = 5000.0,
    Rh         = 0.1,
    T_rot      = 600.0
)
```
"""
struct TwoPopulationDistribution <: N2Distribution
    T_vib_cold::AbstractFloat         
    T_vib_hot::AbstractFloat
    Rh::AbstractFloat                # fraction in hot population
    T_rot::AbstractFloat             # rotational temperature
    Q_cold::AbstractFloat
    Q_hot::AbstractFloat
    Q_rot::AbstractFloat             # rotational partition sum

    function TwoPopulationDistribution(;T_vib_cold::AbstractFloat, T_vib_hot::AbstractFloat,
                                        Rh::AbstractFloat, T_rot::AbstractFloat, vib_states = N2.vib_states(10))
        Q_rot = rot_partition_sum(T_rot)
        Q_cold = sum(Boltzmann_factor(s.E_vib, T_vib_cold) * Q_rot for s in vib_states)
        Q_hot  = sum(Boltzmann_factor(s.E_vib, T_vib_hot) * Q_rot for s in vib_states[2:end])
        return new(T_vib_cold, T_vib_hot, Rh, T_rot, Q_cold, Q_hot, Q_rot)
    end
end

function (df::TwoPopulationDistribution)(state::State)
    f_vib_cold = (1-df.Rh)/df.Q_cold * Boltzmann_factor(state.E_vib, df.T_vib_cold)
    f_vib_hot  = df.Rh/df.Q_hot * Boltzmann_factor(state.E_vib, df.T_vib_hot)
    f_rot      = Boltzmann_factor(state.E_rot, df.T_rot)
    return state.degen * f_rot * (f_vib_cold + f_vib_hot)
end

function vibrational_population(df::TwoPopulationDistribution, v)
    s = State(v,0)
    f_vib_cold = (1-df.Rh)/df.Q_cold * Boltzmann_factor(s.E_vib, df.T_vib_cold)
    f_vib_hot  = df.Rh/df.Q_hot * Boltzmann_factor(s.E_vib, df.T_vib_hot)
    return df.Q_rot * (f_vib_cold + f_vib_hot)
end