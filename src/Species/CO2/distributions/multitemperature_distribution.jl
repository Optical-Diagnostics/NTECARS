"""
    CO2.MultiTemperatureDistribution(; T_12, T_3, T_rot, iso_ID=:O16C12O16) -> df

Constructs a `MultiTemperatureDistribution` and returns it.

# Constructor Arguments
- `T_12::AbstractFloat`: Temperature of the v₁/v₂ vibrational modes in K.
- `T_3::AbstractFloat`: Temperature of the v₃ vibrational mode in K.
- `T_rot::AbstractFloat`: Rotational temperature in K.
- `iso_ID::Symbol`: Isotopologue identifier. Defaults to `:O16C12O16` which is the only one supported.

# Fields of returned type
- `T_12::AbstractFloat`: Temperature of the v₁/v₂ (symmetric stretch / bending)vibrational modes in K
- `T_3::AbstractFloat`: Temperature of the v₃ (asymmetric stretch) vibrationalmode in K.
- `T_rot::AbstractFloat`: Rotational temperature in K.
- `Q::AbstractFloat`: Partition sum, computed from the other fields at construction. Should not be set manually.
- `iso_ID::Symbol`: Isotopologue identifier. Defaults to `:O16C12O16` which is the only one supported.

# Notes
- For a system in full thermal equilibrium, set all three temperatures equal to the gas temperature.
- `Q` is computed via `partition_sum` and stored in the struct. It is
    recomputed automatically if you reconstruct the object with new temperatures.

# Examples
```Julia
    # Thermal equilibrium
    dist = CO2.MultiTemperatureDistribution(T_12=1500.0, T_3=1500.0, T_rot=1500.0)

    # Non-equilibrium: vibrationally excited ν₃ mode
    dist = CO2.MultiTemperatureDistribution(T_12=300.0, T_3=1800.0, T_rot=300.0)
```
"""
mutable struct MultiTemperatureDistribution <: CO2Distribution
    T_12   ::AbstractFloat
    T_3    ::AbstractFloat
    T_rot  ::AbstractFloat
    Q      ::AbstractFloat
    iso_ID ::Symbol

    function MultiTemperatureDistribution(;T_12, T_3, T_rot, iso_ID = :O16C12O16)
        df   = new(T_12 , T_3, T_rot, 1.0, iso_ID)
        df.Q = partition_sum(df)
        return df
    end
end

function fvib(df::MultiTemperatureDistribution, state::Union{State, VibrationalState})
    E_12, E_3 = state.energies.E_vib_12, state.energies.E_vib_3
    f_12 = boltzmann(E_12, df.T_12) 
    f_3  = boltzmann(E_3, df.T_3) 
    return f_12 * f_3
end

function frot(df::MultiTemperatureDistribution, state::State)
    boltzmann(state.energies.E_rot, df.T_rot)
end

