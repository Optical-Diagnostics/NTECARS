"""
    parameters_that_modify_distributions(sim_, parameter_update_function!, parameters)
    
Find out which parameters modify the distrution function of CO2 or N2
This is needed the know which paraeters can be extended with uncertainties
via Measurements.jl since this is only supported for the distribution functions right now.
"""
function parameters_that_modify_distributions(sim_, parameter_update_function!, parameters)
    # make sure that the sim struct has the given parameters
    sim = deepcopy(sim_)
    parameter_update_function!(sim, parameters)

    # array that track whether modified a parameter has changed any distribution function
    does_modify_distribution = fill(false, length(parameters))
    
    for i in eachindex(parameters)
        # modify the parameter slightly
        modified_params = deepcopy(parameters)
        modified_params[i] *= 0.95
        sim_copy = deepcopy(sim)
        parameter_update_function!(sim_copy, modified_params)

        for j in eachindex(sim_copy.species)
            if typeof(sim_copy.species[j]) == N2Species
                if !fieldequal(sim_copy.species[j].distribution, sim.species[j].distribution)
                    does_modify_distribution[i] = true
                end
            elseif typeof(sim_copy.species[j]) == CO2Species
                if !fieldequal(sim_copy.species[j].distribution, sim.species[j].distribution)
                    does_modify_distribution[i] = true
                end
            end
        end
    end
    return does_modify_distribution
end

function fieldequal(a, b)
    typeof(a) === typeof(b) || return false
    for name in fieldnames(typeof(a))
        getfield(a, name) == getfield(b, name) || return false
    end
    return true
end