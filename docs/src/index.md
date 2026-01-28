# NETCARS: Non-Equilibrium CARS Spectra Fitting Code
The NTECARS.jl code is a flexible tool for the calculation and fitting CARS spectra under non-equilibrium conditions. A key strength of the code is the flexibility to model rovibrational distributions using either pre-implemented multi-temperature distribution functions, user-defined functions or to freely fit vibrational populations without the assumption of any distribution function. This flexibility enables the quantitative interpretation of CARS spectra under both equilibrium and strongly non-equilibrium conditions.

# Installation
```julia
using Pkg
Pkg.add("https://github.com/ChristianABusch/NTECARS")
```

# Minimal example for CARS spectrum
```@example
using NTECARS

conditions = GasConditions(
    pressure = 15000.0, # in Pa
    T_gas    = 600.0    # in K
)

lasers = LaserConfiguration(
    wavelength_1 = 532e-9,  # First laser wavelength in nm 
    wavelength_2 = 561e-9,  # Second laser wavelength in nm
    stokes_range = (603e-9, 611e-9),
)

instrument = InstrumentConfiguration(
    profile = Gaussian(0.5/2.35)
)

N2_species = N2Species(
    molar_fraction = 0.9, 
    distribution = N2.MultiTemperatureDistribution(T_vib = 2200.0, T_rot = 600.0), 
)

CO2_species = CO2Species(
    molar_fraction = 0.1, 
    distribution = CO2.MultiTemperatureDistribution(T_12 = 600.0, T_3 = 1800.0, T_rot = 600.0),
    v_max = (0,1,1)
)

sim  = CARSSimulator(
    species    = [CO2_species, N2_species], 
    conditions = conditions, 
    lasers     = lasers, 
    instrument = instrument,
)

synthetic_spectrum = simulate_spectrum(sim)
nothing # hide
```
