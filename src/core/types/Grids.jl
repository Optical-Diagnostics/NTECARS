mutable struct CARSBuffers
    χ      ::Vector{Complex}
    χ_real ::Vector{Float64}
    χ_imag ::Vector{Float64}
    narrowing_function_real::Vector{Float64}
    narrowing_function_imag::Vector{Float64}
    σ_real::Vector{Float64}
    σ_imag::Vector{Float64}

    function CARSBuffers(N)
        new(zeros(Complex, N), zeros(Float64, N), zeros(Float64, N), zeros(Float64, N), 
        zeros(Float64, N), zeros(Float64, N), zeros(Float64, N))
    end
end

struct AdaptiveGrid
    ν_aS         ::Vector{Float64}
    Δν_fine      ::Float64   # in nm
    Δν_coarse    ::Float64 # in nm
    window_width ::Float64    # in nm
    buffer       ::CARSBuffers
end

function AdaptiveGrid(;species::Vector{S}, lasers::LaserConfiguration, conditions::GasConditions) where {S<:CARSSpecies}
    γ_min, γ_max = linewidth_extrema_of_spectrum(conditions = conditions, species = species)
    
    Δν_fine      = γ_min
    Δν_coarse    = 20 * γ_min
    window_width = 20 * γ_max

    transitions_list = vcat([transitions(spec) for spec in species]...)

    ν_aS = adaptive_grid(
        transitions  = transitions_list, 
        ν_aS_limits  = lasers.ν_aS_limits, 
        ν_1          = lasers.ν_1, 
        ν_2          = lasers.ν_2, 
        window_width = window_width, 
        Δν_fine      = Δν_fine, 
        Δν_coarse    = Δν_coarse
    )
    return AdaptiveGrid(ν_aS, Δν_fine, Δν_coarse, window_width, CARSBuffers(length(ν_aS)))
end


struct UniformGrid
    ν_aS  ::Vector{Float64}
    Δν    ::Float64        # in nm
    buffer::CARSBuffers
end

function UniformGrid(;species::Vector{S}, lasers::LaserConfiguration, conditions::GasConditions) where {S<:CARSSpecies}
    γ_min, γ_max = linewidth_extrema_of_spectrum(conditions = conditions, species = species)
    ν_aS = collect(lasers.ν_aS_limits[1]:γ_min:lasers.ν_aS_limits[2])
    return UniformGrid(ν_aS, γ_min, CARSBuffers(length(ν_aS)))
end


##########################################################################
#                  Logic for construction and interpolation
##########################################################################

function interpolate_to_uniform(spec::Spectrum, Δν = 0.01)
    interpol  = linear_interpolation(wavenumbers(spec), intensities(spec), extrapolation_bc=Line())
    ν_uniform = collect(spec.ν[1]:Δν:spec.ν[end])
    I_uniform = interpol.(ν_uniform)

    spec_uniform = Spectrum(ν_uniform, I_uniform, :wavenumber)
    return spec_uniform
end

function adaptive_grid(;transitions, ν_aS_limits, ν_1, ν_2, window_width = 0.002e-9, Δν_fine = 0.01, Δν_coarse = 0.1)
    ν_if           = transition_line_position.(transitions)
    line_positions = vcat(ν_if .+ ν_1, ν_if .+ ν_2)
    line_intervals = intervals_around_points(line_positions, window_width)
    
    ν_aS = sampling_points_from_intervals(line_intervals, ν_aS_limits, Δν_fine, Δν_coarse)
    return unique(ν_aS)
end

"""
    creates an array of intervals that consist of (start_point, end_point) for all points
    in the positions array array with a width of window_width. Overlapping intervals ar merged
"""
function intervals_around_points(positions::Vector{T}, window_width::T) where {T<:Number}
    # create an interval around every line
    intervals = [(x-window_width/2, x+window_width/2) for x in sort(positions)]

    # if there are overlapping intervals, combine them
    N_intervals = length(intervals)
    for i in N_intervals:-1:2
        x_start = intervals[i][1]
        x_end   = intervals[i-1][2]
        
        intervals_overlap = x_end >= x_start 
        if intervals_overlap 
            intervals[i-1] = (intervals[i-1][1], intervals[i][2])
            deleteat!(intervals, i)
        end
    end
    return intervals
end


"""
    creates an array of sampling points in the range specified by "limits".
    In the specified intervals, the points are sampled with the fine resolution
    Δx_fine and inbetween them with Δx_coarse
"""
function sampling_points_from_intervals(intervals, limits, Δx_fine, Δx_coarse)
    sampling_points = []
    push!(sampling_points,limits[1])

    sample_fine(start, stop)   = collect(start:Δx_fine:stop)
    sample_coarse(start, stop) = collect(start:Δx_coarse:stop)

    for I in intervals
        append!(sampling_points, sample_coarse(sampling_points[end], I[1]))
        append!(sampling_points, sample_fine(I[1], I[2]))
    end

    append!(sampling_points, sample_coarse(intervals[end][2], limits[2]))
    push!(sampling_points, limits[end])

    # cut off points that are out of the limits
    sampling_points = filter(x -> (x >= limits[1]) && (x <= limits[2]), sampling_points)
    return unique(sampling_points)
end



