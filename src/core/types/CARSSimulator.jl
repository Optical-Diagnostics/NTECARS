"""
    CARSSimulator(; species, conditions, lasers, instrument, grid_type=:adaptive,
                    vertical_shift=0.0, wavelength_shift=0.0) -> CARSSimulator

Constructs a `CARSSimulator` containing all information required to simulate a CARS spectrum.

# Constructor Arguments
- `species::Vector{T} where {T<:CARSSpecies}`: Vector of species to simulate.
- `conditions::GasConditions`: Bulk gas parameters used for calculating linewidths.
- `lasers::LaserConfiguration`: Laser wavelengths and spectral profiles.
- `instrument::InstrumentConfiguration`: Instrumental broadening profile.
- `grid_type::Symbol`: Can be `:adaptive` or `:uniform`.
- `vertical_shift::Float64`: Vertical offset added after the specturm is normalized to its maximum.
- `wavelength_shift::Float64`: Horizontal shift of the anti-Stokes wavelength axis in
  meters, i.e. the output is evaluated at `λ_aS + wavelength_shift`. Used to account
  for calibration offsets.

# Fields of returned type
- `species::Vector{CARSSpecies}`: Vector of species used in the simulation.
- `conditions::GasConditions`: Bulk gas conditions.
- `lasers::LaserConfiguration`: Laser configuration.
- `instrument::InstrumentConfiguration`: Instrumental broadening profile.
- `grid::Union{UniformGrid, AdaptiveGrid}`: Spectral grid.
- `vertical_shift::Float64`: Vertical offset of the normalised spectrum.
- `wavelength_shift::Float64`: Horizontal wavelength shift in meters.

# Notes
- The spectral grid is constructed once at initialisation. If `conditions`, `species`,
  or `lasers` are changed after construction, the grid is not automatically updated.
  Reconstruct the `CARSSimulator` if the grid needs to reflect the new parameters.
"""
mutable struct CARSSimulator
    species   ::Vector{CARSSpecies}
    conditions::GasConditions
    lasers    ::LaserConfiguration
    instrument::InstrumentConfiguration
    grid      ::Union{UniformGrid, AdaptiveGrid}
    vertical_shift::Float64
    wavelength_shift::Float64

    function CARSSimulator(;
        species   ::Vector{T},
        conditions::GasConditions,
        lasers    ::LaserConfiguration,
        instrument::InstrumentConfiguration,
        grid_type ::Symbol = :adaptive,
        vertical_shift::Float64 = 0.0,
        wavelength_shift::Float64 = 0.0
        ) where {T<:CARSSpecies}

        if grid_type == :adaptive
            grid = AdaptiveGrid(species = species, lasers = lasers, conditions = conditions)
        elseif grid_type == :uniform
            grid = UniformGrid(species = species, lasers = lasers, conditions = conditions)
        end

        new(species, conditions, lasers, instrument, grid, vertical_shift, wavelength_shift)
    end
end

