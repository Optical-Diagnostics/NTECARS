struct HamiltonianEigenstate
    qn::QuantumNumbers
    E::Float64
    basis::Vector{WangQuantumNumbers}
    basis_coeff::Vector{Float64}
end

function HamiltonianEigenstate(qn::QuantumNumbers;)
    qn = QuantumNumbers(qn.v₁, qn.v₂, qn.l₂, qn.v₃, qn.r, qn.J, CO2_parity(qn.J, qn.l₂, qn.v₃))

    if haskey(_eigenstate_cache, qn)
        return _eigenstate_cache[qn]
    end
    
    assigned_eigenstates = all_eigenstates_of_hamiltonian(qn)
    if qn in keys(assigned_eigenstates)
        add_eigenstates_to_cache!(assigned_eigenstates)
        return assigned_eigenstates[qn]
    else
        assigned_eigenstates = all_eigenstates_of_hamiltonian(qn, only_diagonal=true)
        return assigned_eigenstates[qn]
    end
end



"""
Constructs the ffective Hamiltonian to which the state with the given QuantumNumbers qn belongs and
finds all of its Eigenstates.
"""
function all_eigenstates_of_hamiltonian(qn::QuantumNumbers; only_diagonal = false)
    P     = polyad_quantum_number(qn)
    wang_basis, HO_basis = polyad_basis_sets(P, qn.J)
    U     = wang_type_tranformation_matrix(wang_basis, HO_basis)
    H     = hamiltonian(HO_basis, only_diagonal = only_diagonal)

    H_eff = U' * H * U

    basis_e, Heff_e, basis_f, Heff_f = seperate_to_parities(H_eff, wang_basis)
    eigenstates::Dict{QuantumNumbers, HamiltonianEigenstate} = Dict()
    if !isnothing(basis_e)
        merge!(eigenstates, eigenstates_of_hamiltonian_block(Heff_e, basis_e, P))
    end
    if !isnothing(basis_f)
        merge!(eigenstates, eigenstates_of_hamiltonian_block(Heff_f, basis_f, P))
    end

    return eigenstates
end

function eigenstates_of_hamiltonian_block(Heff_block, wang_basis, P)
    eigenvalues, eigenvectors = eigen(Heff_block)
    eigenvectors .*= -1 # convention
    eigenvalues .-= groundstate_energy()
    eigenstates = assign_eigenstates(eigenvalues, eigenvectors, wang_basis, P)
    return eigenstates
end

function groundstate_energy()
    return hamiltonian_diagonal_term(HOQuantumNumbers(0,0,0,0,0))
end

function assign_eigenstates(eigenvalues, eigenvectors, basis, P)
    eigenvectors = [eigenvectors[:, i] for i in 1:length(eigenvalues)]

    # determine the quantum numbers with yet unkown ranking number
    quantum_numbers_without_ranking_index = []
    for (i, eigenvector) in enumerate(eigenvectors)
        leading_contributions = argmax(abs.(eigenvector))
        l₂ = basis[leading_contributions].l₂
        v₃ = basis[leading_contributions].v₃
        v₂ = l₂
        v₁ = (P - l₂ - 3*v₃)/2
        J  = basis[leading_contributions].J
        C  = basis[leading_contributions].C

        push!(quantum_numbers_without_ranking_index, [v₁, v₂, l₂, v₃, J, C])
    end

    # find all sets of equal quantum numbers to determine ranking numbers
    # store in dict that has qn array as keys and the indices where they appear as values
    sets_of_same_qn = Dict()
    for (i, qn) in enumerate(quantum_numbers_without_ranking_index)
        if qn in keys(sets_of_same_qn)
            push!(sets_of_same_qn[qn], i)
        else
            sets_of_same_qn[qn] = [i]
        end
    end
    
    eigenstates = Dict{QuantumNumbers, HamiltonianEigenstate}()
    for quantum_numbers in keys(sets_of_same_qn)
        indices_of_identical_qns = sets_of_same_qn[quantum_numbers]
        indices_ranked_by_energy = indices_of_identical_qns[sortperm(eigenvalues[indices_of_identical_qns]; rev=true)]
        for (r, index) in enumerate(indices_ranked_by_energy)
            v₁, v₂, l₂, v₃, J, C = Int64.(quantum_numbers)
            qn = QuantumNumbers(v₁, v₂, l₂, v₃, r, J, C)
            E  = eigenvalues[index]
            eigenstates[qn] = HamiltonianEigenstate(qn, E, basis, eigenvectors[index])
        end
    end

    return eigenstates
end

function seperate_to_parities(H_eff, wang_basis)
    N = length(wang_basis)
    
    if any(qn -> qn.C==1, wang_basis)
        e_range = findfirst(qn -> qn.C==1, wang_basis):findlast(qn -> qn.C==1, wang_basis)
        H_eff_e = H_eff[e_range, e_range]
        basis_e = wang_basis[e_range]
    else
        H_eff_e,basis_e= nothing, nothing
    end

    if any(qn -> qn.C==2, wang_basis)
        f_range = findfirst(qn -> qn.C==2, wang_basis):findlast(qn -> qn.C==2, wang_basis)
        H_eff_f = H_eff[f_range, f_range]
        basis_f = wang_basis[f_range]
    else
        H_eff_f, basis_f= nothing, nothing
    end
    return basis_e, H_eff_e, basis_f, H_eff_f
end

#######################################################################################################
#                                         Caching of Eigenstates                                         
#######################################################################################################
_eigenstate_cache = Dict{QuantumNumbers, HamiltonianEigenstate}()

function add_eigenstates_to_cache!(assigned_eigenstates::Dict{QuantumNumbers, HamiltonianEigenstate})
    global _eigenstate_cache
    for (qn, eigenstate) in assigned_eigenstates
        _eigenstate_cache[qn] = eigenstate
    end
end





#=
function seperate_to_parities(H_eff, wang_basis)
    N = length(wang_basis)
    
    if any(qn -> qn.C==1, wang_basis)
        e_range = findfirst(qn -> qn.C==1, wang_basis):findlast(qn -> qn.C==1, wang_basis)
        H_eff_e = H_eff[e_range, e_range]
        basis_e = wang_basis[e_range]
    else
        H_eff_e,basis_e= nothing, nothing
    end

    if any(qn -> qn.C==2, wang_basis)
        f_range = findfirst(qn -> qn.C==2, wang_basis):findlast(qn -> qn.C==2, wang_basis)
        H_eff_f = H_eff[f_range, f_range]
        basis_f = wang_basis[f_range]
    else
        H_eff_f, basis_f= nothing, nothing
    end
    return basis_e, H_eff_e, basis_f, H_eff_f
end

function assign_eigenstates2(eigenvalues, eigenvectors, basis, P, C)
    eigenvectors = [eigenvectors[:, i] for i in 1:length(eigenvalues)]

    # determine the quantum numbers with yet unkown ranking number
    quantum_numbers_without_ranking_index = []
    for (i, eigenvector) in enumerate(eigenvectors)
        leading_contributions = argmax(abs.(eigenvector))
        l₂ = basis[leading_contributions].l₂
        v₃ = basis[leading_contributions].v₃
        v₂ = l₂
        v₁ = (P - l₂ - 3*v₃)/2
        J  = basis[leading_contributions].J

        push!(quantum_numbers_without_ranking_index, [v₁, v₂, l₂, v₃, J])
    end

    # find all sets of equal quantum numbers to determine ranking numbers
    # store in dict that has qn array as keys and the indices where they appear as values
    sets_of_same_qn = Dict()
    for (i, qn) in enumerate(quantum_numbers_without_ranking_index)
        if qn in keys(sets_of_same_qn)
            push!(sets_of_same_qn[qn], i)
        else
            sets_of_same_qn[qn] = [i]
        end
    end
    
    eigenstates = Dict{QuantumNumbers, HamiltonianEigenstate}()
    for quantum_numbers in keys(sets_of_same_qn)
        indices_of_identical_qns = sets_of_same_qn[quantum_numbers]
        indices_ranked_by_energy = indices_of_identical_qns[sortperm(eigenvalues[indices_of_identical_qns]; rev=true)]
        for (r, index) in enumerate(indices_ranked_by_energy)
            v₁, v₂, l₂, v₃, J = Int64.(quantum_numbers)
            qn = QuantumNumbers(v₁, v₂, l₂, v₃, r, J, C)
            E  = eigenvalues[index]
            eigenstates[qn] = HamiltonianEigenstate(qn, E, basis, eigenvectors[index])
        end
    end

    return eigenstates
end
=#
