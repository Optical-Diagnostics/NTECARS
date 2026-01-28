# The only point of these structs is a cleaner interface
function Gaussian(σ, μ=0.0)
    power_voigt_spectrum(σ, 0, 1, μ)
end

function Voigt(σ, γ, μ=0.0)
    power_voigt_spectrum(σ, γ, 1.0, μ)
end

function PowerVoigt(σ, γ, n, μ=0.0)
    power_voigt_spectrum(σ, γ, n, μ)
end

##################################################
#               Broanening profiles
##################################################
function power_voigt_spectrum(σ, γ, n, μ=0.0)
    # find where profile decreases to less than 1% of the maximum
    ν_max = find_zero(x -> voigt_super(x, σ, γ, n) - 0.0001, (0.0, 1000.0))

    # calculate spectrum in the determined range
    ν     = collect(LinRange(-ν_max, ν_max, 100))
    I     = voigt_super.(ν, σ, γ, n)
    return Spectrum(ν .+ μ, I .- minimum(I), :wavenumber)
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
