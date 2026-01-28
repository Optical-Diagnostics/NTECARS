# Fitting a CARS spectrum
In NTECARS, for fitting a spectrum, first a `CARSSimulator` struct has to be created for the given conditions. As an example, we will
create this object here, generate a synthetic spectrum and then fit it.

```@example fit
using NTECARS
conditions = GasConditions(
    pressure = 15000.0, # in Pa
    T_gas    = 600.0    # in K
)

lasers = LaserConfiguration(
    wavelength_1 = 532e-9,  # First laser wavelength in nm 
    wavelength_2 = 561e-9,  # Second laser wavelength in nm
    stokes_range = (600e-9, 610e-9),
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
nothing #hide
```

In NTECARS, to fit a set of parameters, an update function has to be defined which updates the `CARSSimulator` from an array of parameters
such as 
```@example fit
function update_function!(sim::CARSSimulator, param)
    T_12, T_3, T_N2vib, T_rot, CO2_frac = param
    
    # update translational temperature
    sim.conditions.T_gas = T_rot 
    
    # update CO2
    sim.species[1].molar_fraction = CO2_frac
    sim.species[1].distribution = CO2.MultiTemperatureDistribution(
        T_12 = T_12, T_3 = T_3, T_rot = T_rot)
    
    #update N2
    sim.species[2].molar_fraction = 1-CO2_frac
    sim.species[2].distribution = N2.MultiTemperatureDistribution(
        T_vib = T_N2vib, T_rot = T_rot)
end
nothing #hide
```
Fitting is then performed by giving the experimental data as a `Spectrum` to the `fit_spectrum()` function together with the
`CARSSimulator`, the update function and initial conditions and boundaries.
```@example fit
result = fit_spectrum(;
    spec_exp     = synthetic_spectrum,
    sim          = sim, 
    initial      = [500.0, 500.0, 500.0, 500.0, 0.1], 
    lower        = [0.0, 0.0, 0.0, 0.0, 0.0],
    upper        = [3000.0, 3000.0, 3000.0, 3000.0, 1.0],
    parameter_update_function! = update_function!
)
result.param
```

The fitted parameters, spectra and rovibrational populations together with their uncertainties can be saved using

```Julia
save_fit_results(folder_path, result)
```