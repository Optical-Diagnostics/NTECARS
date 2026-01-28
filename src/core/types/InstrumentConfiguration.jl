"""
    InstrumentConfiguration(;profile::Union{DeltaProfile,Spectrum} = DeltaProfile())

Contains information about the broadening profile introduced by physical measurement devices (spectrometer). 

When no instrumental broadening is present, it should contain the default `DeltaProfile()`. Otherwise,
`profile` can be given profiles created by existing functions such as `Gaussian(σ)`, `Voigt(σ, γ)` or 
`PowerVoigt(σ, γ, n)`. Alternatively, a `Spectrum` containing a discrete spectrum can be used. 
The profiles should be centered around 0 but do not have to be symmetric.

# Examples
```Julia
InstrumentConfiguration(profile = Gaussian(0.2/2.35))
```
"""
struct InstrumentConfiguration
    profile::Union{DeltaProfile,Spectrum}

    function InstrumentConfiguration(;profile::Union{DeltaProfile,Spectrum} = DeltaProfile())
        new(profile)
    end
end