# values from Lemus et al (2014) (nonlinear)
const a₁′  = 12.43e-42 # C*m^2/V
const a₂′′ = -2.83e-42 #C*m^2/V

function αα(t::CO2Transition)
    a = 0.0
    for (sf,cf) in zip(t.final.polyad_basis, t.final.polyad_coeff)
        for (si,ci) in zip(t.initial.polyad_basis, t.initial.polyad_coeff)
            a += (cf*ci)*transition_moment(sf, si)
        end
    end
    return abs2(a)
end

function αα(initial::State, final::State)
    a = 0.0
    for (sf,cf) in zip(final.polyad_basis, final.polyad_coeff)
        for (si,ci) in zip(initial.polyad_basis, initial.polyad_coeff)
            a += (cf*ci)*transition_moment(sf, si)
        end
    end
    return abs2(a)
end

function transition_moment(f::WangQuantumNumbers, i::WangQuantumNumbers)
    # scaling with quantum numbers can be found in J. Chem. Phys. 79, 5758–5768 (1983)
    if f == WangQuantumNumbers(i.v₁+1, i.v₂, i.l₂, i.v₃, i.J, i.C)
        return sqrt(1/2*(i.v₁+1))*a₁′
    elseif f == WangQuantumNumbers(i.v₁, i.v₂+2, i.l₂, i.v₃, i.J, i.C)
        return sqrt((i.v₂+2)^2 - i.l₂^2) *a₂′′/ 4 #/2 one two comes from the scaling and from taylor expansion in eq 14a
    else
        return zero(typeof(a₁′))
    end
end