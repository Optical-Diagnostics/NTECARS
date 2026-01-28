"""
    CARSSimulator(; species, conditions, lasers, instrument, grid_type, vertical_shift)

struct that contains all information for the calculation of a CARS spectrum.

# Arguments
- `species::Vector{T} where {T<:CARSSpecies}`: Array of `CARSSpecies` such as `CO2Species` or `N2Species`
- `conditions::GasConditions`: bulk gas parameters used for calculating linewidths
- `lasers::LaserConfiguration`: information on laser wavelengths and profiles
- `instrument::InstrumentConfiguration`: inforamtino on instrumental broadening profile
- `grid_type ::Symbol`: Can be `:adaptive` or `:uniform`. the resolutions are automatically 
  determined from the linewidths of the given `CARSSpecies` at the temperature and pressure given in `GasConditions`
- `vertical_shift::Float64` A vertical offset that is added to the spectrum at the end. 
  The spectra are normalized to the maximum so ∈ [0,1]
"""
mutable struct CARSSimulator
    species   ::Vector{CARSSpecies}
    conditions::GasConditions
    lasers    ::LaserConfiguration
    instrument::InstrumentConfiguration
    grid      ::Union{UniformGrid, AdaptiveGrid}
    vertical_shift::Float64
end

function CARSSimulator(;
    species   ::Vector{T},
    conditions::GasConditions,
    lasers    ::LaserConfiguration,
    instrument::InstrumentConfiguration,
    grid_type ::Symbol = :adaptive,
    vertical_shift::Float64 = 0.0
    ) where {T<:CARSSpecies}

    if grid_type == :adaptive
        grid = AdaptiveGrid(species = species, lasers = lasers, conditions = conditions)
    elseif grid_type == :uniform
        grid = UniformGrid(species = species, lasers = lasers, conditions = conditions)
    end

    CARSSimulator(species, conditions, lasers, instrument, grid, vertical_shift)
end