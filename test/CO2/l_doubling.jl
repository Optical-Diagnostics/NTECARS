function hamiltonian_l_doubling_term(qn1::CO2.HarmonicOscillatorQuantumNumbers, qn2::CO2.HarmonicOscillatorQuantumNumbers)
    Δqn = [qn2.v₁-qn1.v₁, qn2.v₂-qn1.v₂, qn2.l₂-qn1.l₂, qn2.v₃-qn1.v₃, qn2.J-qn1.J]
    println(Δqn)
    if Δqn == [0,0,2,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = +1
    elseif Δqn == [0,0,-2,0,0]
        @unpack v₁, v₂, l₂, v₃, J = qn1
        sign = -1
    else
        error("Invalid combination of states for l-type doubling")
    end

    #@unpack ωi,xij,xll,yijk,yill,zijkm,zijll,zllll,Be,αi,γij,γll,εijk,εill,De,βi,ηij,ηll,He,δi,Ge,Fe,Fi,FJ,FJJ,Fij,FeIV,F1IV,F2IV,F3IV,FJIV,Le,LJ,Li,Lij,Ce,CJ,Ci = VIBROT_CONSTANTS

    ±(x) = sign * x
    ∓(x) = -sign * x


    gi = [1,2,1]
    vi = [v₁,v₂,v₃]
    vg_i    = @. vi + gi/2
    vg_ij   = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) for i in 1:3, j in 1:3]

    term = √((v₂ + ±(l₂) + 2)*(v₂ + ∓(l₂))*(J*(J+1) - l₂*(l₂ + ±(1))) * (J*(J+1)-(l₂ + ±(1))*(l₂+ ±(2)))) -0.49e-1#* (Le + sum(Li .* vg_i) + LJ*J*(J+1) + sum(Lij .* vg_ij)) #1.4
    return term
end

qn1 = CO2.HarmonicOscillatorQuantumNumbers(0,1,-1,0,1)
qn2 = CO2.HarmonicOscillatorQuantumNumbers(0,1, 1,0,1)
hamiltonian_l_doubling_term(qn1,qn2)