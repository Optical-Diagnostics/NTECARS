
struct FitResult
    # results
    param::Vector{Float64}
    uncertainties::Vector{Float64}
    fitted_spectrum::Spectrum
    fitted_spectrum_at_measurement::Spectrum
    sim::CARSSimulator
    # for later serialization
    experimental_spectrum::Spectrum
    parameter_update_function!::Function
    intensity_eval_function::Function
end

"""
    fit_spectrum(; spec_exp, sim, parameter_update_function!, initial, lower, upper, 
        intensity_eval_function = x -> abs.(x ./ maximum(x)).^(1/2) , parameter_scaling_factor = initial,
        solver::Symbol = :LM, maxiters= 200)

Fits the the model to a measured spectrum.

# Arguments
- `spec_exp::Spectrum`: The measured spectrum that should be fitted
- `sim::CARSSimulator`: Should contain all required inputs such as linewidths
- `parameter_update_function!::Function`: A function with the signature `f(sim::CARSSimulator, param)`
  that take a `param` array of length `N` and updates `sim` using these parameters.
- `initial::Vector{Float64}`: Initial values of the fit parameters. Has to be of length `N`.
- `lower::Vector{Float64}`: Lower boundary of the fit parameters. Has to be of length `N`.
- `upper::Vector{Float64}`: upper boundary of the fit parameters. Has to be of length `N`.
- `intensity_eval_function::Function = x -> abs.(x ./ maximum(x)).^(1/2),` a function that is applied to 
  both the simulated and experimental spectrum before calculating the residuals. The default represents
  a fitting of the normalized sqaure-root of the CARS-intensity.
- `solver::Symbol`: Options are `:LM` for LevenbergMarquardt and `:IPOPT` for the IPOPT solver.
- `maxiters::Int64`: Maximum number of iterations

# return 
- `FitResult`: contains parameters, spectra, uncertainties etc...

# Example
Which parameters are fitted toegther is determined by the `parameter_update_function!`.

For a `sim::CARSSimulator` that contains the species `N2Species` as the first element in `sim.species`, 
the rotational and vibrational temperature can for example be fitted by defining the update function 
```Julia
function update_function!(sim::CARSSimulator, param)
    T_vib, T_rot = param
    sim.conditions.T_gas = T_rot
    sim.species[1].distribution = N2.MultiTemperatureDistribution(
        T_vib = T_N2vib, T_rot = T_rot)
end

result = fit_spectrum(;
    spec_exp     = experimental_spectrum,
    sim          = sim,
    initial      = [500.0, 500.0],
    lower        = [0.0, 0.0],
    upper        = [3000.0, 3000.0],
    parameter_update_function! = update_function!
)
``` 
"""
function fit_spectrum(;
    spec_exp::Spectrum,
    sim::CARSSimulator, 
    parameter_update_function!::Function, 
    initial::Vector{Float64}, 
    lower::Vector{Float64}, 
    upper::Vector{Float64}, 
    intensity_eval_function::Function = x -> abs.(x ./ maximum(x)).^(1/2),
    parameter_scaling_factors = initial,
    solver::Symbol = :IPOPT,
    maxiters::Int64 = 200
)
    # copy to not modify the users sim struct
    sim_fit = deepcopy(sim)
    auxiliary_parameters = []

    if solver == :LM
        #levenberg_parameter_transform(u) = @. (upper-lower)/2 * tanh(u) + (lower+upper)/2   
        #inv_levenberg_parameter_transform(u) = @. atanh((2*u - (upper + lower)) / (upper - lower))
        levenberg_parameter_transform(u) = @. lower + (upper-lower) * (atan(u)/π + 0.5)
        inv_levenberg_parameter_transform(u) = @. tan(π * ((u - lower) / (upper - lower) - 0.5))

        levenberg_residual(u, p) = residuals(
            sim_fit, 
            levenberg_parameter_transform(u), 
            spec_exp, 
            intensity_eval_function, 
            parameter_update_function!
        )

        prob = NonlinearLeastSquaresProblem(levenberg_residual, inv_levenberg_parameter_transform(initial))
        solver = LevenbergMarquardt(autodiff = ADTypes.AutoFiniteDiff())
        optimal_parameters = solve(prob, solver, abstol = 0.01, reltol = 0.100, maxiters = maxiters).u
        optimal_parameters = levenberg_parameter_transform(optimal_parameters)
    end

    if solver == :IPOPT
        normalize_parmeters(u) = @. u / parameter_scaling_factors
        denormalize_parmeters(u) = @. u * parameter_scaling_factors
        
        Ipopt_residual(u, p) = residuals(
            sim_fit, 
            denormalize_parmeters(u), 
            spec_exp, 
            intensity_eval_function, 
            parameter_update_function!
        )

        optf = SciMLBase.OptimizationFunction(
            (u,p) -> sum(abs2.(Ipopt_residual(u,p))), 
            ADTypes.AutoFiniteDiff()
        )

        prob = SciMLBase.OptimizationProblem(
            optf, 
            normalize_parmeters(initial), 
            auxiliary_parameters, 
            lb = normalize_parmeters(lower), 
            ub = normalize_parmeters(upper)
        )

        opt = IpoptOptimizer(
            hessian_approximation = "limited-memory", 
            acceptable_tol        = 1e-3,
            acceptable_iter       = 2, # number of iterations after solution is acceptable
            mu_strategy           = "monotone"
        )

        optimal_parameters = solve(prob, opt, maxiters = maxiters).u # scaling back to physical dimenions
        optimal_parameters = denormalize_parmeters(optimal_parameters)
    end


    param_uncertainties = uncertainties(
        sim_fit, 
        optimal_parameters, 
        spec_exp, 
        intensity_eval_function, 
        parameter_update_function!
    )     

    result = FitResult(
        optimal_parameters, 
        param_uncertainties, 
        simulate_spectrum(sim_fit),
        simulate_spectrum(sim_fit, wavenumbers(spec_exp)),
        sim_fit,
        spec_exp,
        parameter_update_function!,
        intensity_eval_function
    )

    return result
end


function residuals(sim, param, spec_exp, intensity_eval_function, update_function!)
    update_function!(sim, param)
    spec_sim = simulate_spectrum(sim, wavenumbers(spec_exp))
    I_sim = intensity_eval_function(spec_sim.I)
    I_exp = intensity_eval_function(spec_exp.I)

    r = I_exp .- I_sim
    return r
end

function uncertainties(sim_, optimal_params, spec_exp, intensity_eval_function, parameter_update_function!, use_svd = true)        
    sim = deepcopy(sim_)
    parameter_update_function!(sim, optimal_params)
    # estimate variance of data
    r  = residuals(sim, optimal_params, spec_exp, intensity_eval_function, parameter_update_function!)
    N  = length(wavenumbers(spec_exp))
    p  = length(optimal_params)
    σ² = sum(r .^ 2) / (N - p)

    # get local jacobian
    res(param) = residuals(sim, param, spec_exp, intensity_eval_function, parameter_update_function!)
    J          = FiniteDiff.finite_difference_jacobian(res, optimal_params)

    if use_svd
        # stable computation of inv(transpose(J) * J)
        F           = svd(J, full=false)
        U, s, V, Vt = F.U, F.S, F.V, F.Vt

        # Truncate values smaller than machine precision
        threshold  = eps(Float64) * max(size(J)...) * s[1]
        idx        = s .> threshold
        s          = s[idx]
        Vt         = Vt[idx, :]

        #J_reconstructed = U * (Diagonal(s) * Vt)
        pcov = σ² .* (V * Diagonal(1.0 ./ s.^2)) * Vt
    else
        pcov = σ² .* inv(transpose(J) * J)
    end

    uncertainties = sqrt.(diag(pcov))
    return uncertainties
end