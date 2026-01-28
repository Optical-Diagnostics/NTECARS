module NTECARS
    using Dierckx, DSP, NonlinearSolve, Roots, Interpolations, LinearAlgebra
    using Measurements, Statistics, SpecialFunctions
    using CSV, DataFrames
    using LoopVectorization
    using Unitful, PhysicalConstants.CODATA2018
    using CairoMakie

    include(joinpath("core", "SpectralConverision", "spectral_conversion.jl"))

    include(joinpath("core", "constants.jl"))
    include(joinpath("core", "types", "abstract_types.jl"))

    include(joinpath("Species", "interfaces", "generic_interface.jl"))

    include(joinpath("Species", "N2", "N2.jl"))
    using .N2
    include(joinpath("Species", "interfaces", "N2_interface.jl"))

    include(joinpath("Species", "CO2", "CO2.jl"))
    using .CO2
    include(joinpath("Species", "interfaces", "CO2_interface.jl"))

    include(joinpath("core", "types", "Spectrum.jl"))
    include(joinpath("core", "types", "GasConditions.jl"))
    include(joinpath("core", "types", "LaserConfiguration.jl"))
    include(joinpath("core", "types", "InstrumentConfiguration.jl"))
    include(joinpath("core", "types", "Grids.jl"))
    include(joinpath("core", "types", "CARSSimulator.jl"))
    include(joinpath("core", "misc.jl"))

    include(joinpath("core/Broadening/profiles.jl"))
    include(joinpath("core/Broadening/convolution.jl"))
    include(joinpath("core/Broadening/resolution_change.jl"))

    include(joinpath("core", "susceptibility.jl"))
    include(joinpath("core", "simulate_spectrum.jl"))

    include(joinpath("fitting", "spectrum_fit.jl"))
    include(joinpath("fitting", "serialization", "identify_parameters_modifying_distributions.jl"))
    include(joinpath("fitting", "serialization", "save_results.jl"))

    include(joinpath("plotting", "spectrum_comparison.jl"))

    export GasConditions, LaserConfiguration, InstrumentConfiguration, Spectrum, N2, CO2, 
    CARSSimulator, simulate_spectrum, Spectrum, wavelengths, wavenumbers, intensities,
    gaussian, wavelength_to_wavenumber, wavenumber_to_wavelength, wavelength_to_angular_frequency,
    fit_spectrum, N2Species, CO2Species, raman_shifts, voigt_profile, voigt_super_profile, 
    gaussian_profile, voigt_asymmetric_profile, save_fit_results, spectra_compariso_with_residual,
    cut_spectral_range, total_susceptibility, DeltaProfile, FlatProfile, normalize_by_stokes_intensity!,
    Gaussian, Voigt, PowerVoigt
end



