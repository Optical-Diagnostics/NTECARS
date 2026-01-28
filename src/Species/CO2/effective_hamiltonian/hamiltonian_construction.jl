const VIBROT_CONSTANTS::VibRotConstants = get_CO2_energy_constants(iso = :O16C12O16)

function hamiltonian(Ho_basis; only_diagonal=false)
    H = [hamiltonian_term(qn1, qn2, only_diagonal = only_diagonal) for qn1 in Ho_basis, qn2 in Ho_basis]
    return H
end

function hamiltonian_term(qn1, qn2; only_diagonal=true)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]

    is_diagonal          = Δqn == [ 0, 0, 0, 0, 0]
    is_l_doubling        = Δqn == [ 0, 0, 2, 0, 0] || Δqn == [ 0, 0, -2,0, 0]
    is_fermi_resonance_1 = Δqn == [-1, 2, 0, 0, 0] || Δqn == [ 1,-2, 0, 0, 0]
    is_fermi_resonance_2 = Δqn == [-2, 4, 0, 0, 0] || Δqn == [ 2,-4, 0, 0, 0]
    is_coriolis_1        = Δqn == [-1,-1, 1, 1, 0] || Δqn == [-1,-1,-1, 1, 0]
    is_coriolis_1_mirror = -Δqn==[-1,-1, 1, 1, 0]  || -Δqn==  [-1,-1,-1, 1, 0] # symmetric
    is_fermi_plus_l      = Δqn == [-1 ,2, 2, 0, 0] || Δqn == [-1, 2,-2, 0, 0]
    is_fermi_plus_l_mirror = -Δqn == [-1 ,2, 2, 0, 0] || -Δqn == [-1, 2,-2, 0, 0]

    
    if is_diagonal
        return hamiltonian_diagonal_term(qn1)
    end

    if only_diagonal
        return 0.0
    end

    if is_l_doubling
        return hamiltonian_l_doubling_term(qn1, qn2)
    elseif is_fermi_resonance_1
        return hamiltonian_Fermi_resonance_term(qn1, qn2)
    elseif is_fermi_resonance_2
        return hamiltonian_Fermi_resonance_term2(qn1, qn2)
    elseif is_coriolis_1
        return hamiltonian_coriolis_term(qn1, qn2)
    elseif is_coriolis_1_mirror
        return hamiltonian_coriolis_term(qn2, qn1)
    elseif is_fermi_plus_l
        return 0.0#hamiltonian_Fermi_plus_l_doubling_term(qn1, qn2)
    elseif is_fermi_plus_l_mirror
        return 0.0#hamiltonian_Fermi_plus_l_doubling_term(qn2, qn1)
    else
        return 0.0
    end
end

function hamiltonian_diagonal_term(qn::HOQuantumNumbers)
    @unpack v₁, v₂, l₂, v₃, J = qn
    @unpack ωi,xij,xll,yijk,yill,zijkm,zijll,zllll,Be,αi,γij,γll,εijk,εill,De,βi,ηij,ηll,He,δi,Ge,Fe,Fi,FJ,FJJ,Fij,FeIV,F1IV,F2IV,F3IV,FJIV,Le,LJ,Li,Lij,CJ,Ci = VIBROT_CONSTANTS

    gi = [1,2,1]
    vi = [v₁,v₂,v₃]

    vg_i    = @. vi + gi/2
    vg_ij   = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) for i in 1:3, j in 1:3]
    vg_ijk  = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) * (vi[k]+gi[k]/2) for i in 1:3, j in 1:3, k in 1:3]
    vg_ijkm = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) * (vi[k]+gi[k]/2) * (vi[m]+gi[m]/2) for i in 1:3, j in 1:3, k in 1:3, m in 1:3]

    E = sum(ωi .* vg_i) + sum(xij .* vg_ij) + xll*l₂^2 + sum(yijk .* vg_ijk) + sum(yill .* vg_i .* l₂^2) +
        sum(zijkm .* vg_ijkm) + sum(zijll .* vg_ij .* l₂^2) + zllll * l₂^4 +
        (Be - sum(αi .* vg_i) + sum(γij .* vg_ij) + γll * l₂^2 + sum(εijk .* vg_ijk) + sum(εill .* vg_i .* l₂^2)) * (J*(J+1) - l₂^2) -
        (De + sum(βi .* vg_i) + sum(ηij .* vg_ij) + ηll * l₂^2) * (J*(J+1) - l₂^2)^2 +
        (He + sum(δi .* vg_i)) * (J*(J+1) - l₂^2)^3 +
        Ge * (J*(J+1) - l₂^2)^4
    return E
end

function hamiltonian_Fermi_resonance_term(qn1::HOQuantumNumbers, qn2::HOQuantumNumbers)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]
    if Δqn == [-1,2,0,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
    elseif Δqn == [1,-2,0,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn2
    else
        error("Invalid combination of states for Fermi-resonance")
    end

    Fe, Fi, Fij = VIBROT_CONSTANTS.Fe, VIBROT_CONSTANTS.Fi, VIBROT_CONSTANTS.Fij
    FJ, FJJ = VIBROT_CONSTANTS.FJ , VIBROT_CONSTANTS.FJJ 

    gi = [1,2,1]
    vi = [v₁,v₂,v₃]
    Δv = [-1,2,0]

    vg_i    = @. vi + (gi + Δv)/2
    vg_ij   = [(vi[i]+(gi[i]+Δv[i])/2) * (vi[j]+(gi[j]+Δv[j])/2) for i in 1:3, j in 1:3]

    √(v₁*(v₂+l₂+2)*(v₂-l₂+2)) * (Fe + sum(Fi .* vg_i) + FJ*(J*(J+1) - l₂^2) + sum(vg_ij .* Fij) + FJJ * (J*(J+1) - l₂^2)^2)
end

function hamiltonian_Fermi_resonance_term2(qn1::HOQuantumNumbers, qn2::HOQuantumNumbers)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]
    if Δqn == [-2,4,0,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
    elseif Δqn == [2,-4,0,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn2
    else
        error("Invalid combination of states for Fermi-resonance")
    end

    FeIV = VIBROT_CONSTANTS.FeIV

    FeIV * √(v₁*(v₁-1)*(v₂+l₂+2)*(v₂+l₂+4)*(v₂-l₂+2)*(v₂-l₂+4))
end

function hamiltonian_l_doubling_term(qn1::HOQuantumNumbers, qn2::HOQuantumNumbers)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]
    if Δqn == [0,0,2,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = +1
    elseif Δqn == [0,0,-2,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = -1
    else
        error("Invalid combination of states for l-type doubling")
    end

    Le, LJ, Li, Lij = VIBROT_CONSTANTS.Le, VIBROT_CONSTANTS.LJ, VIBROT_CONSTANTS.Li, VIBROT_CONSTANTS.Lij 

    ±(x) = sign * x
    ∓(x) = -sign * x


    gi = [1,2,1]
    vi = [v₁,v₂,v₃]
    vg_i    = @. vi + gi/2
    vg_ij   = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) for i in 1:3, j in 1:3]

    term = √((v₂ + ±(l₂) + 2)*(v₂ + ∓(l₂))*(J*(J+1) - l₂*(l₂ + ±(1))) * (J*(J+1)-(l₂ + ±(1))*(l₂+ ±(2)))) * (Le+ sum(Li .* vg_i) + LJ*J*(J+1) + sum(Lij .* vg_ij))
    return term
end

function hamiltonian_Fermi_plus_l_doubling_term(qn1::HOQuantumNumbers, qn2::HOQuantumNumbers)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]
    if Δqn == [-1,2,2,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = +1
    elseif Δqn == [-1,2,-2,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = -1
    else
        error("Invalid combination of states for Fermi + l-type doubling")
    end

    ±(x) = sign * x
    ∓(x) = -sign * x

    return √(v₁*(v₂ + ±(l₂)+2)*(v₂+ ±(l₂)+4)*(J*(J+1) -l₂*(l₂+ ±(1)))*(J*(J+1) - (l₂+ ±(1)) * (l₂+ ±(2)))) * (-9.31e-6)
end


function hamiltonian_coriolis_term(qn1::HOQuantumNumbers, qn2::HOQuantumNumbers)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]

    if Δqn == [-1,-1,1,1,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = 1
    elseif Δqn == [-1,-1,-1,1,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = -1
    else
        error("Invalid combination of states for coriolis-interaction")
    end

    if J <= abs(l₂)
        return 0.0
    end

    Ce, CJ, Ci = VIBROT_CONSTANTS.Ce, VIBROT_CONSTANTS.CJ, VIBROT_CONSTANTS.Ci 

    ±(x)  = sign * x
    ∓(x)  = -sign * x

    gi = [1,2,1]
    vi = [v₁,v₂,v₃]
    Δv = [-1,-1,1]

    vg_i = @. vi + (gi + Δv)/2

    √(v₁*(v₂ + ∓(l₂))*(v₃+1)*(J*(J+1)-l₂*(l₂ + ±(1)))) * (Ce +sum(Ci .* vg_i) + CJ * J*(J+1))
end


function rotational_energy_from_database(qn::QuantumNumbers)
    Be, De, He = rotational_constants(qn)
    J = qn.J
    return Be * (J * (J + 1)) - De * (J * (J + 1))^2 + He * (J * (J + 1))^3
end