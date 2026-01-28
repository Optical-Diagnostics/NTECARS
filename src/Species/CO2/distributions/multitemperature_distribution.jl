mutable struct MultiTemperatureDistribution <: CO2Distribution
    T_12   ::AbstractFloat
    T_3    ::AbstractFloat
    T_rot  ::AbstractFloat
    Q      ::AbstractFloat
    iso_ID ::Symbol

    function MultiTemperatureDistribution(;T_12, T_3, T_rot, iso_ID = :O16C12O16)
        df   = new(T_12 , T_3, T_rot, 1.0, iso_ID)
        df.Q = partition_sum(df)
        return df
    end
end

function fvib(df::MultiTemperatureDistribution, state::Union{State, VibrationalState})
    E_12, E_3 = state.energies.E_vib_12, state.energies.E_vib_3
    f_12 = boltzmann(E_12, df.T_12) 
    f_3  = boltzmann(E_3, df.T_3) 
    return f_12 * f_3
end

function frot(df::MultiTemperatureDistribution, state::State)
    boltzmann(state.energies.E_rot, df.T_rot)
end

