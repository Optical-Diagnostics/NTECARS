
"""
    save_fit_results(folderpath, result::FitResult,
                     parameter_labels=["fit param \$(i)" for i in eachindex(result.param)])

Saves the complete results of a `FitResult` to a folder structure of CSV files,
including the fitted and experimental spectra, fit parameters with uncertainties,
and rovibrational populations for each species.

# Arguments
- `folderpath`: Path to the output folder.
- `result::FitResult`: The fit result returned by [`fit_spectrum`](@ref).
- `parameter_labels`: Labels for the fit parameters used in the output file. 
  Defaults to `["fit param 1", "fit param 2", ...]`.

# Output files
- `fit_and_experiment_at_measurement_points.csv`: Fitted and experimental `√(I_CARS)`
  at the experimental wavenumber grid, including the residual. Contains anti-Stokes
  wavelength and Raman shift columns relative to both lasers.
- `fit_at_simulated_resolution.csv`: Fitted `√(I_CARS)` at the full simulation
  resolution, which is typically finer than the experimental grid.
- `fit_parameters.csv`: Fitted parameter values with absolute and relative
  uncertainties (%) for each parameter.
- `covariance_matrix.csv`: The covaraince matrix of the fit. If `absolute_chi=false`,
  the covariance matrix was rescaled with the redced chi squared.
- `reduced_chi2.csv`: The value of the reduced chi squared.
- `CO2_rovibrational_populations/`: Created if a `CO2Species` is present in
  `result.sim.species`. Contains one CSV per vibrational state with rotational
  populations `[CO2(v,J)]/[CO2]` and their uncertainties, plus a summary
  `vibrational_populations.csv`.
- `N2_rovibrational_populations/`: Created if an `N2Species` is present in
  `result.sim.species`. Contains one CSV per vibrational level with rotational
  populations `[N2(v,J)]/[N2]` and their uncertainties, plus a summary
  `vibrational_populations.csv`.

# Notes
- Rovibrational populations are computed with uncertainty propagation via the
  `Measurements.jl` package. Only the uncertainties of fit parameters and not
  the uncertainties of the polarizabilities and linewidths are considered.

# Examples
```Julia
save_fit_results("results/measurement_01", result)

# With custom parameter labels
save_fit_results("results/measurement_01", result,
    ["T_12 (K)", "T_3 (K)", "T_N2vib (K)", "T_rot (K)", "CO2 molar fraction"])
```
"""
function save_fit_results(folderpath, result::FitResult, parameter_labels = ["fit param $(i)" for i in eachindex(result.param)])
    # create folder structure for given folderpath if it doesnt exist
    mkpath(folderpath) 
    
    # for calculating the sqrt of the CARS intensity
    normalized_sqrt(I) = sqrt.(abs.(I ./ maximum(I)))

    #####################################################################################################
    #                    fitted and experimental CARS spectrum at measured wavelengths
    #####################################################################################################
    ν_raman_1, ν_raman_2 = raman_shifts(chi2 = result.experimental_spectrum, lasers = result.sim.lasers)
    field_names = [
        "anti-Stokes wavelength (nm)",
        "Raman shift (cm^-1) [omega_1-omega_S]",
        "Raman shift (cm^-1) [omega_2-omega_S]",
        "sqrt(I_CARS) [Experiment]",
        "sqrt(I_CARS) [Fit]",
        "residual"
    ]
    field_values = [
        wavelengths(result.experimental_spectrum) .* 1e9,
        ν_raman_1,
        ν_raman_2,
        normalized_sqrt(result.experimental_spectrum.I),
        normalized_sqrt(result.fitted_spectrum_at_measurement.I),
        normalized_sqrt(result.experimental_spectrum.I) .- normalized_sqrt(result.fitted_spectrum_at_measurement.I)
    ]

    df = DataFrame(field_values, field_names)
    save_path = joinpath(folderpath, "fit_and_experiment_at_measurement_points.csv")
    CSV.write(save_path, df)

    #####################################################################################################
    #                              Fitted spectrum at simulation resolution
    #####################################################################################################
    ν_raman_1, ν_raman_2 = raman_shifts(chi2 = result.fitted_spectrum, lasers = result.sim.lasers)
    field_names = [
        "anti-Stokes wavelength (nm)",
        "Raman shift (cm^-1) [omega_1-omega_S]",
        "Raman shift (cm^-1) [omega_2-omega_S]",
        "sqrt(I_CARS) [Fit]"
    ]
    field_values = [
        wavelengths(result.fitted_spectrum) .* 1e9,
        ν_raman_1,
        ν_raman_2,
        normalized_sqrt(result.fitted_spectrum.I)
    ]

    df = DataFrame(field_values, field_names)
    save_path = joinpath(folderpath, "fit_at_simulated_resolution.csv")
    CSV.write(save_path, df)

    #####################################################################################################
    #                                           Fit parameters
    #####################################################################################################
    field_names = [
        "Parameter",
        "Fitted values",
        "Absolute uncertainty",
        "relative_uncertainty (%)"
    ]

    field_values = [
        parameter_labels,
        result.param,
        result.uncertainties,
        result.uncertainties ./ result.param .* 100
    ]
    df = DataFrame(field_values, field_names)
    save_path = joinpath(folderpath, "fit_parameters.csv")
    CSV.write(save_path, df)

    #####################################################################################################
    #                            covariance matrix and reduced reduced_chi2
    ##################################################################################################### 
    df = DataFrame(result.pcov, :auto)
    save_path = joinpath(folderpath, "covariance_matrix.csv")
    CSV.write(save_path, df)

    df = DataFrame([[result.reduced_chi2]], ["reduced chi squared"])
    save_path = joinpath(folderpath, "reduced_chi2.csv")
    CSV.write(save_path, df)


    #####################################################################################################
    #                                     rovibrational_populations
    #####################################################################################################
    modifies_distribution = parameters_that_modify_distributions(result.sim, result.parameter_update_function!, result.param)
    param_with_error   = []
    for i in eachindex(result.param)
        if modifies_distribution[i]
            push!(param_with_error, result.param[i] ± result.uncertainties[i])
        else
            push!(param_with_error, result.param[i])
        end
    end
    
    result.parameter_update_function!(result.sim, param_with_error)
    for species in result.sim.species
        if typeof(species) == CO2Species
            CO2_results_folder = mkpath(joinpath(folderpath, "CO2_rovibrational_populations"))
            save_CO2_rovibrational_populations(CO2_results_folder, species)
        elseif typeof(species) == N2Species
            N2_results_folder = mkpath(joinpath(folderpath, "N2_rovibrational_populations"))
            save_N2_rovibrational_populations(N2_results_folder, species)
        end
    end
end

function save_CO2_rovibrational_populations(CO2_results_path, CO2species)
    NvJ, Nv = rovibrational_populations_database(CO2species)
    #####################################################################################################
    #                                     rovibrational_populations
    #####################################################################################################
    for vib_qn in keys(NvJ)
        J             = []
        E_total       = []
        E_rot         = []
        populations   = []
        uncertainties = []
        for i in 1:size(NvJ[vib_qn])[2]
            s,f = NvJ[vib_qn][:,i]
            push!(J, s.qn.J)
            push!(E_total, s.energies.E_total)
            push!(E_rot, s.energies.E_rot)
            push!(populations,   Measurements.value(f))
            push!(uncertainties, Measurements.uncertainty(f))
        end
        vib_qn_string = "$(vib_qn.v₁)$(vib_qn.v₂)$(vib_qn.l₂)$(vib_qn.v₃)$(vib_qn.r)"

        field_names = [
            "rotational quantum number J",
            "E_total($(vib_qn_string)J) (cm^-1))",
            "E_rot($(vib_qn_string)J)) (cm^-1)",
            "[CO2($(vib_qn_string)J)] / [CO2]",
            "absolute uncertainty of [CO2($(vib_qn_string)J)] / [CO2]",
            "relative uncertainty of [CO2($(vib_qn_string)J)] / [CO2]"
        ]
        field_values = [
            J,
            E_total,
            E_rot,
            populations,
            uncertainties,
            uncertainties ./ populations .* 100
        ]

        df = DataFrame(field_values, field_names)
        file_path = joinpath(CO2_results_path, vib_qn_string * "_J.csv")
        CSV.write(file_path, df)
    end

    #####################################################################################################
    #                                     vibrational_populations
    #####################################################################################################
    vib_qn_strings = []
    energies       = []
    populations    = []
    uncertainties  = []
    for vib_qn in keys(Nv)
        E_vib, f = Nv[vib_qn]
        push!(vib_qn_strings, "$(vib_qn.v₁)$(vib_qn.v₂)$(vib_qn.l₂)$(vib_qn.v₃)$(vib_qn.r)")
        push!(energies, E_vib)
        push!(populations,   Measurements.value(f))
        push!(uncertainties, Measurements.uncertainty(f))
    end

    order = sortperm(energies)
    field_names = [
        "v1,v2,l2,v3,r",
        "energy (cm^-1)",
        "[CO2(v1,v2,l2,v3,r)] / [CO2]",
        "absolute uncertainty of [CO2(v1,v2,l2,v3,r)] / [CO2]",
        "relative uncertainty of [CO2(v1,v2,l2,v3,r)] / [CO2]"
    ]

    field_values = [
        vib_qn_strings[order],
        energies[order],
        populations[order],
        uncertainties[order],
        uncertainties[order] ./ populations[order] .* 100
    ]

    df = DataFrame(field_values, field_names)
    file_path = joinpath(CO2_results_path, "vibrational_populations.csv")
    CSV.write(file_path, df)
end


function save_N2_rovibrational_populations(N2_results_path, N2species)
    NvJ, Nv = rovibrational_populations_database(N2species)
    #####################################################################################################
    #                                     rovibrational_populations
    #####################################################################################################
    for v in keys(NvJ)
        J             = []
        E_total       = []
        E_rot         = []
        populations   = []
        uncertainties = []
        for i in 1:size(NvJ[v])[2]
            s,f = NvJ[v][:,i]
            push!(J, s.qn.J)
            push!(E_total, s.E_tot)
            push!(E_rot, s.E_rot)
            push!(populations,   Measurements.value(f))
            push!(uncertainties, Measurements.uncertainty(f))
        end
        vib_qn_string = "$(v)"

        field_names = [
            "rotational quantum number J",
            "E_total(v=$(vib_qn_string),J) (cm^-1))",
            "E_rot(v=$(vib_qn_string),J)) (cm^-1)",
            "[N2(v=$(vib_qn_string),J)] / [N2]",
            "absolute uncertainty of [N2(v=$(vib_qn_string),J)] / [N2]",
            "relative uncertainty of [N2(v=$(vib_qn_string),J)] / [N2]"
        ]
        field_values = [
            J,
            E_total,
            E_rot,
            populations,
            uncertainties,
            uncertainties ./ populations .* 100
        ]

        df = DataFrame(field_values, field_names)
        file_path = joinpath(N2_results_path, "v=$(vib_qn_string)" * "_J.csv")
        CSV.write(file_path, df)
    end

    #####################################################################################################
    #                                     vibrational_populations
    #####################################################################################################
    vib_qn_strings = []
    energies       = []
    populations    = []
    uncertainties  = []
    for v in keys(Nv)
        E_vib, f = Nv[v]
        push!(vib_qn_strings, "$(v)"),
        push!(energies, E_vib)
        push!(populations, Measurements.value(f))
        push!(uncertainties, Measurements.uncertainty(f))
    end

    order = sortperm(energies)
    field_names = [
        "v",
        "energy (cm^-1)",
        "[N2(v)] / [N2]",
        "absolute uncertainty of [N2(v)] / [N2]",
        "relative uncertainty of [N2(v)] / [N2]"
    ]

    field_values = [
        vib_qn_strings[order],
        energies[order],
        populations[order],
        uncertainties[order],
        uncertainties[order] ./ populations[order] .* 100
    ]

    df = DataFrame(field_values, field_names)
    file_path = joinpath(N2_results_path, "vibrational_populations.csv")
    CSV.write(file_path, df)
end