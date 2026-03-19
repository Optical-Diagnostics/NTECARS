"""
    GasConditions(; pressure, T_gas) -> GasConditions

Constructs a `GasConditions` the gas temperature and pressure used for linewidth 
calculations

# Constructor Arguments
- `pressure::AbstractFloat`: Gas pressure in Pa.
- `T_gas::AbstractFloat`: Translational temperature in K.

# Fields of returned type
- `pressure::AbstractFloat`: Gas pressure in Pa.
- `T_gas::AbstractFloat`: Translational temperature in K.

# Notes
- `T_gas` represents the translational temperature and is used for linewidth calculations. 
  It should match `T_rot` of the species distributions when assuming rotational equilibrium.

# Examples
```Julia
conditions = GasConditions(
    pressure = 15000.0,  # 150 mbar
    T_gas    = 600.0     # K
)
```
"""
mutable struct GasConditions
    pressure::AbstractFloat
    T_gas   ::AbstractFloat

    function GasConditions(;pressure::AbstractFloat, T_gas::AbstractFloat) 
        new(pressure, T_gas)
    end
end

