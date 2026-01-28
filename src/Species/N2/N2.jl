module N2
    using UnPack
    using PhysicalConstants.CODATA2018
    using Unitful
    using TOML
    
    include("constants.jl")
    include("file_paths.jl")

    include("states/states.jl")
    include("states/energies.jl")

    include("distributions/general.jl")
    include("distributions/multitemperature_distribution.jl")
    include("distributions/free_vibrational_distribution.jl")
    include("distributions/two_temperature_distribution.jl")

    include("spectrum/transitions.jl")
    include("spectrum/polarizability.jl")
    include("spectrum/linewidth.jl")

    export N2Distribution, N2Transition
end