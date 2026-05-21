"""
    save_spectrum(file_path, spec::Spectrum, sim::CARSSimulator)

Saves a spectrum to a CSV file with anti-Stokes wavelength, Raman shifts relative
to both lasers, and intensity columns.

# Arguments
- `file_path`: Path to the output CSV file, e.g. `"results/spectrum.csv"`.
- `spec::Spectrum`: The spectrum to save.
- `sim::CARSSimulator`: Simulator providing the laser configuration used to compute
  the Raman shifts.

# Output columns
- `anti-Stokes wavelength (nm)`: Anti-Stokes wavelength axis in nm.
- `Raman shift (cm⁻¹) [ω₁ - ωₛ]`: Raman shift relative to the first laser in cm⁻¹.
- `Raman shift (cm⁻¹) [ω₂ - ωₛ]`: Raman shift relative to the second laser in cm⁻¹.
- `intensity (a.u.)`: CARS intensity.
- `sqrt(intensity) (a.u.)`: Square root of the CARS intensity.

# Notes
- Prints the output path to the console upon successful completion.

# Examples
```Julia
save_spectrum("results/spectrum.csv", spec, sim)
```
"""
function save_spectrum(file_path, spec::Spectrum, sim::CARSSimulator)
    λ = wavelengths(spec)
    I = intensities(spec)
    ramanshift_ν1_νS, ramanshift_ν2_νS = raman_shifts(chi2 = spec, lasers = sim.lasers)

    # for calculating the sqrt of the CARS intensity
    I_sqrt = sqrt.(abs.(I ./ maximum(I)))

    field_names = [
        "anti-Stokes wavelength (nm)",
        "Raman shift (cm^-1) [omega_1-omega_S]",
        "Raman shift (cm^-1) [omega_2-omega_S]",
        "intensity (a.u.)",
        "sqrt(intensity) (a.u.)",
    ]
    field_values = [
        λ,
        ramanshift_ν1_νS,
        ramanshift_ν2_νS,
        I,
        I_sqrt,
    ]

    df = DataFrame(field_values, field_names)
    CSV.write(file_path, df)
    println("saved spectrum at: $(file_path)")
end