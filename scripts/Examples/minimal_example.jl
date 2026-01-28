using NTECARS
using CairoMakie
CairoMakie.activate!()

###############################################################
#                           Setup
###############################################################
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

N2_species  = N2Species(
    molar_fraction  = 0.9, 
    distribution = N2.MultiTemperatureDistribution(T_vib = 2200.0, T_rot = 600.0), 
)

CO2_species = CO2Species(
    molar_fraction = 0.1, 
    distribution = CO2.MultiTemperatureDistribution(T_12 = 600.0, T_3 = 1800.0, T_rot = 600.0),
)

sim  = CARSSimulator(
    species    = [CO2_species, N2_species], 
    conditions = conditions, 
    lasers     = lasers, 
    instrument = instrument,
)

###############################################################
#                     Calculate spectrum
###############################################################
@time chi2 = simulate_spectrum(sim)

###############################################################
#                          Plotting
###############################################################
with_theme(theme_latexfonts()) do
    fig = Figure(size = (1000, 700), fontsize = 40, figure_padding = 30)
    ax = Axis(fig[1,1], 
        limits = ((494, 499), (-0.02, 1.1)),
        xlabel = L"$\mathrm{λ_{aS} \, (nm)}$",
        ylabel = L"$\mathrm{(CARS \, intensity)^{1/2}}$"
    )
    I_sqrt_normalized = sqrt.(chi2.I ./ maximum(chi2.I))
    lines!(ax, wavelengths(chi2) .* 1e9, I_sqrt_normalized, color = :black)
    fig
end
