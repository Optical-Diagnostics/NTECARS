
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

function update_function!(sim::CARSSimulator, param)
    println(param)
    T_12, T_3, T_N2vib, T_rot, CO2_frac = param
    
    # update translational temperature
    sim.conditions.T_gas = T_rot 
    
    # update CO2
    sim.species[1].molar_fraction = CO2_frac
    sim.species[1].distribution = CO2.MultiTemperatureDistribution(
        T_12  = T_12, 
        T_3   = T_3,
        T_rot = T_rot
        )
    
    #update N2
    sim.species[2].molar_fraction = 1-CO2_frac
    sim.species[2].distribution = N2.MultiTemperatureDistribution(
        T_vib = T_N2vib, 
        T_rot = T_rot
        )
end


@time result_free = fit_spectrum(;
    spec_exp     = synthetic_spectrum,
    sim          = sim, 
    initial      = [500.0, 500.0, 500.0, 500.0, 0.5], 
    lower        = [0.0, 0.0, 0.0, 0.0, 0.0],
    upper        = [3000.0, 3000.0, 3000.0, 3000.0, 1.0],
    parameter_update_function! = update_function!
)


using CairoMakie
with_theme(theme_latexfonts()) do
    fig = Figure(size = (800, 600), fontsize = 40, figure_padding = 30)
    ax = Axis(fig[1,1], 
        limits = ((494.0, 499), (-0.02, 1.1)),
        xlabel = L"$\mathrm{λ_{aS} \, (nm)}$",
        ylabel = L"$\mathrm{\sqrt{Intensity} \, (a.u.)}$",
        ygridvisible = false,
        xgridvisible = false,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(5),
    )
    lines!(ax, 1e9 .* wavelengths(synthetic_spectrum), sqrt.(intensities(synthetic_spectrum, normalization = :maximum)))
    fig
end

