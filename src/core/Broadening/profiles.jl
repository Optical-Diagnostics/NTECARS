# The only point of these structs is a cleaner interface
"""
    Gaussian(σ, μ=0.0) -> Spectrum

Convenience constructor for a normalised Gaussian spectral profile.

# Arguments
- `σ::Float64`: Gaussian width (standard deviation, FWHM=2.355*σ) in cm⁻¹.
- `μ::Float64`: Centre offset in cm⁻¹. Defaults to `0.0`.

# Returns
- `Spectrum`: Sampled Gaussian profile centred at `μ`.

# Examples
```Julia
profile = Gaussian(0.2)          # σ = 0.2 cm⁻¹, centred at 0
profile = Gaussian(0.2, 1.0)     # σ = 0.2 cm⁻¹, centred at 1.0 cm⁻¹
```
"""
function Gaussian(σ, μ=0.0)
    power_voigt_spectrum(σ, 0, 1, μ)
end

"""
    Voigt(σ, γ, μ=0.0) -> Spectrum

Convenience constructor for a normalised Voigt spectral profile.

# Arguments
- `σ::Float64`: Gaussian width (standard deviation, FWHM=2.355*σ) in cm⁻¹.
- `γ::Float64`: Lorentzian half-width in cm⁻¹.
- `μ::Float64`: Centre offset in cm⁻¹. Defaults to `0.0`.

# Returns
- `Spectrum`: Sampled Voigt profile centred at `μ`.

# Examples
```Julia
profile = Voigt(0.1, 0.05)        # σ = 0.1 cm⁻¹, γ = 0.05 cm⁻¹
profile = Voigt(0.1, 0.05, 1.0)   # same but centred at 1.0 cm⁻¹
```
"""
function Voigt(σ, γ, μ=0.0)
    power_voigt_spectrum(σ, γ, 1.0, μ)
end

"""
    PowerVoigt(σ, γ, n, μ=0.0) -> Spectrum

Convenience constructor for a normalised power-Voigt spectral profile, i.e. a
Voigt profile raised to the power `n`. Setting `n = 1` recovers a standard Voigt.
Setting `γ = 0, n = 1` recovers a Gaussian.

# Arguments
- `σ::Float64`: Gaussian width (standard deviation, FWHM=2.355*σ) in cm⁻¹.
- `γ::Float64`: Lorentzian half-width in cm⁻¹.
- `n::Float64`: Exponent applied to the normalised Voigt profile.
- `μ::Float64`: Centre offset in cm⁻¹. Defaults to `0.0`.

# Returns
- `Spectrum`: Sampled power-Voigt profile centred at `μ`.

# Examples
```Julia
profile = PowerVoigt(0.1, 0.05, 2.0)       # squared Voigt profile
profile = PowerVoigt(0.1, 0.05, 2.0, 1.0)  # same but centred at 1.0 cm⁻¹
```
"""
function PowerVoigt(σ, γ, n, μ=0.0)
    power_voigt_spectrum(σ, γ, n, μ)
end

##################################################
#               Broanening profiles
##################################################
function power_voigt_spectrum(σ, γ, n, μ=0.0)
    threshold = 0.0001
    if voigt_super(1000.0, σ, γ, n) < threshold
        # find where profile decreases to less than 0.001% of the maximum
        ν_max = find_zero(x -> voigt_super(x, σ, γ, n) - 1e-5, (0.0, 1000.0))
    else
        ν_max = 100
    end

    # calculate spectrum in the determined range
    ν     = collect(LinRange(-ν_max, ν_max, 500))
    I     = voigt_super.(ν, σ, γ, n)
    return Spectrum(ν .+ μ, I, :wavenumber)
end

function voigt(x, σ, γ)
    z = (x + 1im*γ) / (√(2)*σ)
    Faddeeva(z) = erfcx(-1im*z)
    V = 1/(σ*2π) * real(Faddeeva(z))

    z0 = (0.0 + 1im*γ) / (√(2)*σ)
    V0 = 1/(σ*2π) * real(Faddeeva(z0))
    return V/V0
end

function voigt_super(x, σ, γ, n)
    return voigt(x, σ, γ) ^n
end
