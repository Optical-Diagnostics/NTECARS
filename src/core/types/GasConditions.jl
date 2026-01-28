"""
    GasConditions(;pressure::AbstractFloat, T_gas::AbstractFloat) 

Contains the basis information of the pressure in Pa and translational temperature in K.

The pressure and translational temperature are used for calculating linewidths. When fitting,
it should be remembered to update `T_gas` accordingly.

# Examples
```Julia
conditions = GasConditions(
    pressure = 15000.0,
    T_gas    = 600.0   
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

