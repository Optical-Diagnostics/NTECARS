"""
    InstrumentConfiguration(; profile=DeltaProfile()) -> InstrumentConfiguration

Constructs an `InstrumentConfiguration` storing the instrumental broadening profile
of the spectrometer or other measurement devices and returns it.

# Constructor Arguments
- `profile::Union{DeltaProfile, Spectrum}`: Instrumental broadening profile, centred
  around 0. Use `DeltaProfile()` if no instrumental broadening is present. For standard 
  profiles, convenience constructors are implemented: [`Gaussian`](@ref),
  [`Voigt`](@ref), [`PowerVoigt`](@ref). A spectrum from use data can be passed in the form
  of a [`Spectrum`](@ref).The profile does not have to be symmetric.

# Fields of returned type
- `profile::Union{DeltaProfile, Spectrum}`: The instrumental broadening profile.

# Examples
```Julia
# No instrumental broadening
instrument = InstrumentConfiguration()

# Gaussian broadening with FWHM = 0.2 cm⁻¹
instrument = InstrumentConfiguration(profile = Gaussian(0.2/2.355))

# Measured instrument profile from data
instrument = InstrumentConfiguration(
    profile = Spectrum(wavenumbers, intensities, :wavenumber))
```
"""
mutable struct InstrumentConfiguration
    profile::Union{DeltaProfile,Spectrum}

    function InstrumentConfiguration(;profile::Union{DeltaProfile,Spectrum} = DeltaProfile())
        new(profile)
    end
end