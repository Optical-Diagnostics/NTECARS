

module CO2
    using TOML
    using Parameters
    using LinearAlgebra
    using StaticArrays
    using PhysicalConstants.CODATA2018
    using Unitful

    include("filepaths.jl")
    include("constants.jl")

    include("states/quantum_numbers.jl")

    # effective Hamiltonian
    include("effective_hamiltonian/load_energys_constants.jl")
    include("effective_hamiltonian/basis_and_transformation.jl")
    include("effective_hamiltonian/hamiltonian_construction.jl")
    include("effective_hamiltonian/hamiltonian_eigenstates.jl")

    #states

    include("states/metadata.jl")
    include("states/states.jl")
    include("states/symmetry.jl")

    #distributions
    include("distributions/general.jl")
    include("distributions/multitemperature_distribution.jl")
    include("distributions/states_for_partition_sum.jl")
    include("distributions/pre_calculations.jl")

    # Spectra
    include("spectrum/transitions.jl")
    include("spectrum/linewidth.jl")
    include("spectrum/polarizability.jl")

    export CO2Distribution, CO2Transition
end