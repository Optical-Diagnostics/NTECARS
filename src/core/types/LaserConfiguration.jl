abstract type AbstractSpectralProfile end

struct DeltaProfile <: AbstractSpectralProfile end
struct FlatProfile <: AbstractSpectralProfile end

"""
    LaserConfiguration(; wavelength_1, wavelength_2, stokes_range,
                         profile_1=DeltaProfile(),
                         profile_2=DeltaProfile(),
                         stokes_profile=FlatProfile())

Stores the laser wavelengths and spectral profiles.

# Arguments
- `wavelength_1::Float64`: Central wavelength of the first laser (in meters).
- `wavelength_2::Float64`: Central wavelength of the second laser (in meters).
- `stokes_range::Tuple{Float64,Float64}`: Wavelength range over which the
  Stokes spectrum is defined (in meters). This determines the range of anti-Stokes 
  frequencies together with the central wavelengths.
- `profile_1::Union{DeltaProfile, Spectrum}`: Spectral profile of the first laser.
- `profile_2::Union{DeltaProfile, Spectrum}`: Spectral profile of the second laser.
- `stokes_profile::Spectrum`: Stokes spectral profile. The Stokes-spectrum has to be defined at
  the correct physical spectral range (not around 0!). A `FlatProfile` should be used if 
  the experimental data is already normalized by the stokes-proflie/non-resonant background.

# Notes
- `Spectrum` can be created from user data or Gaussian(σ), Voigt(σ, γ), PowerVoigt(σ, γ, n)
- Internally, wavelengths are converted to wavenumbers.
- The anti-Stokes limits are computed as `extrema(ν₁ + ν₂ − stokes_range)`.

# Examples
```Julia
lasers = LaserConfiguration(
    wavelength_1    = 532e-9,
    wavelength_2    = 561e-9,
    stokes_range    = (600e-9, 610e-9),
    profile_1       = Gaussian(0.2), # inputs are σ in cm^-1
    profile_2       = Voigt(0.1, 0.05), # inputs are σ, γ  in cm^-1
    stokes_profile  = Spectrum([600e-9, 605e-9, 610e-0], [0.0, 1.0, 0.0], :wavelength)
)
```
"""
mutable struct LaserConfiguration
    ν_1::Float64
    ν_2::Float64
    profile_1::Union{DeltaProfile, Spectrum} 
    profile_2::Union{DeltaProfile, Spectrum} 
    stokes_profile::Spectrum 
    ν_aS_limits::Tuple{Float64, Float64}
    
    function LaserConfiguration(;
        wavelength_1::Float64, 
        wavelength_2::Float64, 
        stokes_range::Tuple{Float64, Float64}, 
        profile_1::Union{DeltaProfile, Spectrum} = DeltaProfile(), 
        profile_2::Union{DeltaProfile, Spectrum} = DeltaProfile(),
        stokes_profile::Union{FlatProfile, Spectrum} = FlatProfile()
    )
        if stokes_profile isa FlatProfile
            stokes_profile = Spectrum([stokes_range...], [1.0, 1.0], :wavelength)
        end
        
        ν_1 = wavelength_to_wavenumber(wavelength_1)
        ν_2 = wavelength_to_wavenumber(wavelength_2)
        ν_S = wavenumbers(stokes_profile)
        ν_aS = ν_1 .+ ν_2 .- ν_S
        new(ν_1, ν_2, profile_1, profile_2, stokes_profile, extrema(ν_aS))
    end
end

