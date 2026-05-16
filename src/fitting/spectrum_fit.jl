"""
    FitResult

Stores the complete result of a [`fit_spectrum`](@ref) call.

# Fields
- `param::Vector{Float64}`: Optimal fit parameters found by the solver.
- `uncertainties::Vector{Float64}`: Parameter uncertainties estimated from the
  Jacobian SVD at the optimal solution. Set to zeros if `compute_uncertainties = false`
  was passed to `fit_spectrum`.
- `fitted_spectrum::Spectrum`: Simulated spectrum evaluated on the simulator's
  internal wavenumber grid.
- `fitted_spectrum_at_measurement::Spectrum`: Simulated spectrum evaluated at the
  wavenumber positions of the experimental spectrum. Use this for direct
  comparison with `experimental_spectrum`.
- `sim::CARSSimulator`: Deep copy of the simulator after convergence, updated with
  the optimal parameters via `parameter_update_function!`.
- `experimental_spectrum::Spectrum`: The measured spectrum passed to `fit_spectrum`.
  Stored to allow serialization and later re-evaluation of the result.
- `parameter_update_function!::Function`: The update function used during fitting.
  Stored to allow re-simulation or re-fitting from the result.
- `weights::Function`: The weights used for weighted least-squares

# Notes
- `fitted_spectrum` and `fitted_spectrum_at_measurement` differ only in their
  wavenumber grids. For plotting against measured data, prefer
  `fitted_spectrum_at_measurement`.
- All fields needed to reproduce or serialize the fit are stored directly in the
  struct, so no reference to the original `sim` or functions is required after fitting.
"""
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
    weights::Union{Nothing, Spectrum}
    pcov::Matrix{Float64}
    reduced_chi2::Float64
end

"""
    fit_spectrum(; spec_exp, sim, parameter_update_function!, initial, lower, upper, 
        intensity_eval_function = x -> abs.(x ./ maximum(x)).^(1/2) , parameter_scaling_factor = initial,
        solver::Symbol = :LM, maxiters= 200, compute_uncertainties = true) -> FitResult

Fits the the model to a measured spectrum and returns a collection of results in a [`FitResult`](@ref).

# Arguments
- `spec_exp::Spectrum`: The measured spectrum that should be fitted
- `sim::CARSSimulator`: Should contain all required inputs such as linewidths
- `parameter_update_function!::Function`: A function with the signature `f(sim::CARSSimulator, param)`
  that take a `param` array of length `N` and updates `sim` using these parameters.
- `initial::Vector{Float64}`: Initial values of the fit parameters. Has to be of length `N`.
- `lower::Vector{Float64}`: Lower boundary of the fit parameters. Has to be of length `N`.
- `upper::Vector{Float64}`: upper boundary of the fit parameters. Has to be of length `N`.
- `weights::Union{Spectrum, Nothing}=nothing`: weights used for weighted least squares (formular: `residuals = weights .* (I_exp - I_sim).^2`).
  For absolute uncertainties, the weights represent w=1/σ^2 where σ is standard deviation of the measured intensity. 
- `absolute_sigma`::Bool`: If `true`, then the weights (sigmas) are used to give absolute uncertainties. Otherwise, 
  sigma is scaled to match the sample variance between the fitted and experimenal spectral.
- `solver::Symbol`: Options are `:LM` for LevenbergMarquardt and `:IPOPT` for the IPOPT solver.
- `maxiters::Int64`: Maximum number of iterations
- `compute_uncertainties`::Bool = true: Mostly an option to disable for cases is which the SVD throws an error. 
- `tolerance::Float64 = 1e-5`: Tolerance that is passed to the solver.

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
        T_vib = T_vib, T_rot = T_rot)
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
    parameter_scaling_factors = initial,
    solver::Symbol = :IPOPT,
    maxiters::Int64 = 200,
    compute_uncertainties = true,
    weights::Union{Spectrum, Nothing} = nothing,
    absolute_sigma = false,
    tolerance::Float64 = 1e-5,
)
    # copy to not modify the users sim struct
    sim_fit = deepcopy(sim)
    auxiliary_parameters = []

    if solver == :LM
        levenberg_parameter_transform(u) = @. lower + (upper-lower) * (atan(u)/π + 0.5)
        inv_levenberg_parameter_transform(u) = @. tan(π * ((u - lower) / (upper - lower) - 0.5))

        levenberg_residual(u, p) = residuals(
            sim_fit, 
            levenberg_parameter_transform(u), 
            spec_exp, 
            parameter_update_function!,
            weights,
            #transform = x -> sqrt(abs(x))
        )

        prob = NonlinearLeastSquaresProblem(levenberg_residual, inv_levenberg_parameter_transform(initial))
        solver = LevenbergMarquardt(autodiff = ADTypes.AutoFiniteDiff())
        optimal_parameters = solve(prob, solver, abstol = tolerance, reltol = tolerance, maxiters = maxiters).u
        optimal_parameters = levenberg_parameter_transform(optimal_parameters)
    end

    if solver == :IPOPT
        normalize_parmeters(u) = @. u / parameter_scaling_factors
        denormalize_parmeters(u) = @. u * parameter_scaling_factors
        
        Ipopt_residual(u, p) = residuals(
            sim_fit, 
            denormalize_parmeters(u), 
            spec_exp, 
            parameter_update_function!,
            weights,
            #transform = x -> sqrt(abs(x))
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
            acceptable_tol        = tolerance,
            acceptable_iter       = 4, # number of iterations after solution is acceptable
            mu_strategy           = "monotone"
        )

        optimal_parameters = solve(prob, opt, maxiters = maxiters).u # scaling back to physical dimenions
        optimal_parameters = denormalize_parmeters(optimal_parameters)
    end

    if compute_uncertainties
        param_uncertainties, pcov, χ² = uncertainties(
            sim_fit, 
            optimal_parameters, 
            spec_exp, 
            parameter_update_function!,
            weights,
            absolute_sigma,
            true
        )    
    else 
        Ni = length(optimal_parameters)
        param_uncertainties = zeros(Ni)
        pcov = zeros(Ni, Ni)
        χ² = 0
    end

    result = FitResult(
        optimal_parameters, 
        param_uncertainties, 
        simulate_spectrum(sim_fit),
        simulate_spectrum(sim_fit, wavenumbers(spec_exp)),
        sim_fit,
        spec_exp,
        parameter_update_function!,
        weights,
        pcov,
        χ²
    )

    return result
end


function residuals(sim, param, spec_exp::Spectrum, update_function!, weights::Union{Spectrum, Nothing}; transform::Function = x -> identity(x))
    # transform() is applied to help the solver consider peaks that may be too small by applying the sqrt()
    # to the spectra. Experimentally, I found that it works better to not apply the transform to the weights
    # when calculating the uncertainties, transform(x)=identity(x), so the ucnertainties are not
    # falsified
    
    # update sim to new parameters
    update_function!(sim, param)
    # get updated simulated spectrum
    spec_sim = simulate_spectrum(sim, wavenumbers(spec_exp))
    I_sim    = intensities(spec_sim, normalization = :maximum)

    # normalize exp data and weights which should be 1/sigma^2, together
    norm_exp = 1/maximum(intensities(spec_exp))
    I_exp    = norm_exp .* intensities(spec_exp)

    I_exp = transform.(I_exp)
    I_sim = transform.(I_sim)

    # apply transformation function. sqrt for better fitting and identity for the uncertainties
    r = I_exp .- I_sim

    if isnothing(weights)
        return r
    else
        σ = norm_exp .* 1 ./ sqrt.(abs.(intensities(weights)))
        #σ = transform.(σ)
        r_weighted = r ./ σ
        return r_weighted
    end
end

function uncertainties(sim_, optimal_params, spec_exp, parameter_update_function!, weights::Union{Spectrum, Nothing}, absolute_sigma = false, use_svd = true)        
    sim = deepcopy(sim_)
    parameter_update_function!(sim, optimal_params)
    # estimate variance of data
    r  = residuals(sim, optimal_params, spec_exp, parameter_update_function!, weights)
    N  = length(wavenumbers(spec_exp))
    p  = length(optimal_params)
    σ² = sum(r .^ 2) / (N - p) # reduced chi-squared

    # get local jacobian
    res(param) = residuals(sim, param, spec_exp, parameter_update_function!, weights)
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
        V          = V[:,idx]

        #J_reconstructed = U * (Diagonal(s) * Vt)
        pcov = (V * Diagonal(1.0 ./ s.^2)) * Vt
    else
        pcov = inv(transpose(J) * J)
    end

    if absolute_sigma == false || isnothing(weights)
        # if the weights given dont reflect absolute uncertainties
        pcov .*= σ²
    end

    uncertainties = sqrt.(diag(pcov))
    χ² = σ²
    return uncertainties, pcov, χ²
end