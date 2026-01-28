# Getting Started with CARS-spectra Calculation
The NTECARS code works around the `CARSSimulator` object, which contains all data necessary to calculated a CARS spectrum.
This object is required both for the forward calculation of a CARS spectrum but also to fit. When fitting, the `CARSSimulator`
object is the starting point and gets iteratively modified.

## Spectra struct
All CARS spectra, laser profiles and instrumental profiles are `Spectrum` objects. For the laser profiles and instrumental profile,
the code implements constructors for different profiles like Gaussian, or Voigt profiles. Otherwise, the `Spectrum` has to be 
created from discrete data. A `Spectrum` can be constructed from either wavelengths (m) and intensities as an input:
```@setup ntecars
using NTECARS
using CairoMakie
```

```@example ntecars
λ  = [404e-9, 405e-9, 406e-9]
I  = [0.0, 1.0, 0.0]
s  = Spectrum(λ, I, :wavelength)
nothing #hide
```
or `Spectrum` can be constructed from wavenumbers (cm¹) and intensities as an input:
```@example ntecars
ν  = [20000, 21000, 22000]
I  = [0.0, 1.0, 0.0]
s  = Spectrum(ν, I, :wavenumber)
nothing #hide
```
When working with structs it is recommended to use the getter-functions which automatically handle unit conversions or normalizations
```@example ntecars
ν  = wavenumbers(s)
λ  = wavelengths(s)
I  = intensities(s)
I  = intensities(s, normalization = :maximum)
I  = intensities(s, normalization = :sum)
nothing #hide
```

## GasConditions
The bulk-gas parameters of pressure (Pa) and translational temperature (K) are used for the calculation and stored in `GasConditions`:
```@example ntecars
conditions = GasConditions(
    pressure = 15000.0, # in Pa
    T_gas    = 600.0    # in K
)
nothing #hide
```

## LaserConfiguration
Information about the wavelengths and profiles of the lasers is stored in `LaserConfiguration`. The minimum inputs are the central wavelengths of laser
1 and laser 2 and the spectral range of the Stokes-laser that determines the evaluated anti-Stokes frequencies. All inputs are in units of wavelength (m):

```@example ntecars
lasers = LaserConfiguration(
    wavelength_1    = 532e-9,
    wavelength_2    = 561e-9,
    stokes_range    = (600e-9, 610e-9)
)
nothing #hide
```
optional arguments are `profile_1`, `profile_2` which are `DeltaProfile()` by default and the `stokes_profile` which is `FlatProfile()` be default. All three
profiles can be passed a custom `Spectrum`. The narrow laser profiles `profile_1`, `profile_2` should be centerd around 0.0 and the `stokes_profile` should be 
at the correct absolute wavelength. For the convolutions, the profiles are interpolated, so they dont have to be defined at specific points. An example is:
```@example ntecars
lasers = LaserConfiguration(
    wavelength_1    = 532e-9,
    wavelength_2    = 561e-9,
    stokes_range    = (600e-9, 610e-9),
    profile_1       = Gaussian(0.2), # inputs are σ
    profile_2       = Voigt(0.1, 0.05), # inputs are σ, γ
    stokes_profile  = Spectrum([600e-9, 605e-9, 610e-0], [0.0, 1.0, 0.0], :wavelength)
)
nothing #hide
```

## Instrumental Broadening
just like in the LaserConfiguration, the instrumental profile can be a custom `Spectrum`, which should be centered around 0.0 or be passed pre-implemted profiles such as
`Gaussian` or `Voigt`. By default, teh instrumental profile is a `DeltaProfile()` e.g. it introduces no broadening.
```@example ntecars
instrument = InstrumentConfiguration(
    profile = Gaussian(0.1/2.35)
)
nothing #hide
```

## Molecular species
Each molecule is represented in the code as a struct. The minimum for their construction is a molar fraction and a distribution 
```@example ntecars
N2_species  = N2Species(
    molar_fraction = 0.9, 
    distribution = N2.MultiTemperatureDistribution(T_vib = 2200.0, T_rot = 600.0), 
)

CO2_species = CO2Species(
    molar_fraction = 0.1, 
    distribution = CO2.MultiTemperatureDistribution(T_12 = 600.0, T_3 = 1800.0, T_rot = 600.0),
)
nothing #hide
```
Here distribution functions can be set to different functions. For the implementation of custom distribution functions, please check the implementation of the current ones.
They need to be structs that are subtypes of the abstract types `CO2Distribution` and `N2.distribution` respecitvely, and be callable as `distribution(state)` where
`state` is a `CO2.State` or `N2.State` respectively. Especially for CO2, the calculation time of a spectrum rapidly increases with the maximum quantum numbers evaluated for a spectrum increasing. 
For this reason, the maximum quantum numbers considered for each spectrum can be set as an input. Furthermore if required, the non-resonant susceptibility can be changed. An example is:
```@example ntecars
N2_species  = N2Species(
    molar_fraction = 0.9, 
    distribution = N2.FreeVibrationalDistribution(delta_f_vib = [0.9, 0.09, 0.01], 
                        T_rot = 600.0, closure = :temperature_extrapolation), 
    v_max = 3,
    J_max = 50,
    chi_non_resonant = 4.2e-51
)

CO2_species = CO2Species(
    molar_fraction = 0.1, 
    distribution = CO2.MultiTemperatureDistribution(T_12 = 600.0, T_3 = 1800.0, T_rot = 600.0),
    v_max = (1,2,1),  # v1, v2, v3
    J_max = 80,
    chi_non_resonant = 4.5e-51
)
nothing #hide
```

## CARSSimulator
Finally, the created `GasConditions`, `LaserConfiguration`, `InstrumentConfiguration` and molecules such as `N2Species` and `CO2Species` are stored in the CARSSimulator struct
```@example ntecars
sim  = CARSSimulator(
    species    = [CO2_species, N2_species], 
    conditions = conditions, 
    lasers     = lasers, 
    instrument = instrument,
)
nothing #hide
```

## Calculation of spectrum
A spectrum is simply evaluated from the `CARSSimulator` using 
```@example ntecars
synthetic_spectrum = simulate_spectrum(sim)
nothing #hide
```
which returns a `Spectrum` containing the anti-Stokes wavenumbers and a intensity (not sqrt) of the anti-Stokes beam. 

```@example ntecars
with_theme(theme_latexfonts()) do
    fig = Figure(size = (1000, 700), fontsize = 40, figure_padding = 30)
    ax = Axis(fig[1,1], 
        limits = ((494, 499), (-0.02, 1.1)),
        xlabel = L"$\mathrm{λ_{aS} \, (nm)}$",
        ylabel = L"$\mathrm{(CARS \, intensity)^{1/2}}$"
    )
    lines!(ax, 1e9 .* wavelengths(synthetic_spectrum), sqrt.(intensities(synthetic_spectrum, normalization = :maximum)))
    fig
end
```
