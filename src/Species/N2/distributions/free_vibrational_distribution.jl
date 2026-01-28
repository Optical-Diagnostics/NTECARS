struct FreeVibrationalDistribution <: N2Distribution
    f_vib::Vector{AbstractFloat}      # relative populations, f_vib[1] == 1.0 for v=0
    T_rot::AbstractFloat              # rotational temperature
    Q    ::AbstractFloat              # total partition sum
    Q_rot::AbstractFloat              # rotational partition sum (cached)
end

function FreeVibrationalDistribution(;
    f_vib      ::Union{Nothing, Vector{T}} = nothing,
    delta_f_vib::Union{Nothing, Vector{T}} = nothing,
    T_rot      ::AbstractFloat,
    closure    ::Union{Nothing, Symbol} = :zero
) where {T <: AbstractFloat}

    conflicting_input = !isnothing(f_vib)  && !isnothing(delta_f_vib)
    if conflicting_input
        derror("Define only the absolute poplations or the population differences")
    end

    from_absolute     = !isnothing(f_vib) && isnothing(delta_f_vib)
    from_differences  = isnothing(f_vib) && !isnothing(delta_f_vib)
    if from_absolute
        return free_distribution_from_absolute(f_vib = f_vib, T_rot = T_rot)
    elseif from_differences
        return free_distribution_from_differences(;delta_f_vib = delta_f_vib, T_rot = T_rot, closure = closure)
    end
end

function (df::FreeVibrationalDistribution)(state::State)
    Nv = length(df.f_vib)
    if v(state)+1 > Nv
        return 0.0
    end
    f_vib = df.f_vib[v(state)+1]
    f_rot = Boltzmann_factor(state.E_rot, df.T_rot)
    return state.degen/df.Q * f_vib * f_rot
end



function vibrational_population(df::FreeVibrationalDistribution, v::Int64)
    Q_rot = rot_partition_sum(df.T_rot)
    return  Q_rot/df.Q * df.f_vib[v+1]
end

function free_distribution_from_absolute(;f_vib::Vector{T}, T_rot::AbstractFloat, vib_states = N2.vib_states(length(f_vib))) where {T <: AbstractFloat}
    Q_rot = rot_partition_sum(T_rot) # rotational partition sum
    Q     = sum(f_vib .* Q_rot)      # Total partition sum over vibrational states 
    return FreeVibrationalDistribution(f_vib, T_rot, Q, Q_rot)
end

function free_distribution_from_differences(;delta_f_vib::Vector{T}, T_rot::AbstractFloat, closure::Symbol) where {T <: AbstractFloat}
    f_vib = vibrational_populations_from_differences(delta_f_vib, closure = closure)
    return free_distribution_from_absolute(f_vib = f_vib, T_rot = T_rot)
end

function vibrational_populations_from_differences(Δf_vib::Vector{T}; closure = :temperature_extrapolation) where {T <: AbstractFloat}
    N = length(Δf_vib)
    if N == 1
        return [1.0]
    end

    # highest intial and final state of highest transition
    vi_max = N-1
    vf_max = N
        
    if closure == :temperature_extrapolation && Δf_vib[N] < Δf_vib[N-1] && Δf_vib[N] > 0
        Ev(v) = State(v,0).E_vib

        # get temperature connecting the 3 states of the last 2 differences
        ΔE = Ev(vi_max) - Ev(vi_max-1)
        T̃  = -ΔE / (1/c₂ * log(Δf_vib[N]/Δf_vib[N-1]))

        # obtain population of final state of the highest transition
        ΔE           = Ev(vi_max+1) - Ev(vi_max)
        f_vib_vi_max = Δf_vib[N] / (1-exp(-c₂ * ΔE/T̃))
        f_vib_vf_max = f_vib_vi_max * exp(-c₂ * ΔE/T̃)
    elseif closure == :temperature_extrapolation && Δf_vib[N] >= Δf_vib[N-1]
        println("Warning: Δf_vib not decreasing monotonically, closure set to :zero")
        f_vib_vi_max = Δf_vib[N]
        f_vib_vf_max = 0.0
    elseif closure == :zero # set upper state of highest transition to 0
        f_vib_vi_max = Δf_vib[N]
        f_vib_vf_max = 0.0
    else
        f_vib_vi_max = Δf_vib[N]
        f_vib_vf_max = 0.0
    end

    # get all other populations in absolute values
    f_vib = [f_vib_vf_max, f_vib_vi_max]
    for i in N-1:-1:1
        # f(v) = [fv(v+1)-f(v)] + f(v+1)
        push!(f_vib, f_vib[end]+Δf_vib[i])
    end

    # reverse to have population in order ascending with quantum number v
    f_vib = reverse(f_vib)
    # normalize such that f(v=0) = 1
    f_vib ./= sum(f_vib)
    return f_vib
end