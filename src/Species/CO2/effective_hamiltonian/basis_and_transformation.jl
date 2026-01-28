struct HOQuantumNumbers <: AbstractQuantumNumbers
    # quantum numbers of harmonic oscillator, used for calc of effective hamiltonian
    v₁ ::Int16 # quantum number of symmetric stretch mode
    v₂ ::Int16 # quantum number of bending mode
    l₂ ::Int16 # effective rotation of bending modes
    v₃ ::Int16 # quantum number of asymmetric stretch mode
    J  ::Int16
end

struct WangQuantumNumbers <: AbstractQuantumNumbers
    # quantum numbers of harmonic oscillator, used for calc of effective hamiltonian
    v₁ ::Int16 # quantum number of symmetric stretch mode
    v₂ ::Int16 # quantum number of bending mode
    l₂ ::Int16 # effective rotation of bending modes
    v₃ ::Int16 # quantum number of asymmetric stretch mode
    J  ::Int16
    C  ::Int16 # parity, 1="e", 2="f"
end


"""
Requies that the state is specified in HITRAN notation with v₂ = l₂
"""
function polyad_quantum_number(qn::CO2.AbstractQuantumNumbers)
    v₁, v₂, v₃ = qn.v₁, qn.v₂, qn.v₃
    P = 2*v₁ + v₂ + 3*v₃  
end

function polyad_basis_sets(P, J)
    v₁_max = P ÷ 2 
    v₃_max = P ÷ 3

    # Get the basis vectors both in the HO-basis (l₂ = ±|l₂|) and wang basis HO-basis (|l₂|, "e"/"f")
    wang_basis = []
    HO_basis   = []

    parity(J,l₂,v₃) = (J+l₂+v₃) % 2 == 0 ? 1 : 2

    for v₃ in 0:v₃_max
        P′ = P - 3*v₃ 
        for v₁ in v₁_max:-1:0
            v₂ = P′ - 2*v₁
            for l₂ in v₂:-2:0
                if l₂ == 0
                    push!(wang_basis, WangQuantumNumbers(v₁, v₂, l₂, v₃, J, parity(J,l₂,v₃)))
                    push!(HO_basis, HOQuantumNumbers(v₁, v₂, l₂, v₃, J))
                else
                    push!(wang_basis, WangQuantumNumbers(v₁, v₂, l₂, v₃, J, 1))
                    push!(wang_basis, WangQuantumNumbers(v₁, v₂, l₂, v₃, J, 2))
                    push!(HO_basis, HOQuantumNumbers(v₁, v₂, l₂, v₃, J))
                    push!(HO_basis, HOQuantumNumbers(v₁, v₂, -l₂, v₃, J))
                end
            end
        end
    end

    # sort wang basis into e and f blocks
    idx_of_e = findall(qn -> qn.C == 1, wang_basis)
    idx_of_f = findall(qn -> qn.C == 2, wang_basis)
    wang_basis = vcat(wang_basis[idx_of_e], wang_basis[idx_of_f])

    return wang_basis, HO_basis 
end

function wang_type_tranformation_matrix(wang_basis, HO_basis)
    N = length(HO_basis)
    U = zeros(N,N)

    index_of_HOqn(target_qn) = findfirst(qn -> qn == target_qn, HO_basis)

    for (i,wqn) in enumerate(wang_basis)
        @unpack v₁, v₂, l₂, v₃, J, C = wqn
        if l₂ == 0
            j = index_of_HOqn(HOQuantumNumbers(v₁, v₂, l₂, v₃, J))
            U[j,i] = 1
        else
            j_plus  = index_of_HOqn(HOQuantumNumbers(v₁, v₂, l₂, v₃, J))
            j_minus = index_of_HOqn(HOQuantumNumbers(v₁, v₂, -l₂, v₃, J))
            U[j_plus,i]  = 1/√(2)
            U[j_minus,i] = 1/√(2) * (-1)^(C+1)
        end
    end
    return U
end