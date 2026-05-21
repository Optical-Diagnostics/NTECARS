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

Now, the parameters that should be fitted have to be specified. In NTECARS, this is done by defining an update function that takes an array of parameters and updates the `CARSSimulator` from this array. This function is given to the solver. For this example, we will fit the fraction of CO2 in the gas `CO2_frac`, and the vibrational populations of CO2 and N2 using their respective `MultiTemperatureDistribution`. In total this corresponds to 5 fit parameters. The update function for this case is:
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
The NTECARS code then fits the model to the experiment. To match the arbitrary amplitude of the experimental data in the fitting process, the simulated spectrum is scaled by a constant that is automatically fitted together with the rest of the parameters. This approach was chosen over the initial approach used in the paper in which both experiment and simulation were normalized to their respective maximum, since the fitted scaling parameter leads to more stable results.

The fitted parameters, spectra and rovibrational populations together with their uncertainties can be saved using
```Julia
save_fit_results(folder_path, result)
```

## Weighted fitting
If the spectrum covers a large dynamic range with weak lines, a weighted least squares fitting has to be performed since the contribution of weak lines to the residual that is minimized is otherwise to small for the solver. For this, weights can be passed to the solver as well in the form a `Spectrum` that should be defined at the same points as the experimentally measured intensities. In the calculation fo the residual, the weights are applied as `w[i] * (I_exp[i]-I_sim[i])^2`.

For minimizing the relative intensity differences between the model and the experiment rather than absolute differences, the weight has to correspond to `1/I^2`. The weighted fit is then performed by calling:
```@example fit
weights = Spectrum(wavenumbers(synthetic_spectrum), 1 ./ intensities(synthetic_spectrum).^2, :wavenumber)

result = fit_spectrum(;
    spec_exp     = synthetic_spectrum,
    sim          = sim, 
    weights      = weights,
    initial      = [500.0, 500.0, 500.0, 500.0, 0.1], 
    lower        = [0.0, 0.0, 0.0, 0.0, 0.0],
    upper        = [3000.0, 3000.0, 3000.0, 3000.0, 1.0],
    parameter_update_function! = update_function!
)
nothing #hide
```

## Uncertainties
The NTECARS code can also calculate uncertainties for the fitted parameters and vibrational populations. If no uncertainties are provided in the form of weights, then the uncertainties are simply calculated from the Jacobian of the variance between the experimental data and the fit. If weights are supplied, then the weighted variance is used for the calculation of uncertainties. 

To calculate absolute uncertainties, the weights have to correspond to the inverse square of the uncertainty of the intensity e.g. w = 1/sigma^2. Furthermore, `fit_spectrum` has to be given the argument `absolute_sigma = true`. Setting `absolute_sigma = false` scales the uncertainties to yield a reduced chi-squared of unity. If the reduced chi-squared, which is also saved by `save_fit_results`, deviates from unity, the `absolute_sigma = false` option is likely more appropriate.
