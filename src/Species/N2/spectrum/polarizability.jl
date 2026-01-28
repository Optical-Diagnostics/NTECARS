# eq. 6.6.32-34 in The Raman Effect A Unified Treatment of the Theory of Raman Scattering by Molecules - Derek A. Long
@inline bJJ(J::Int)   = J*(J+1)/(2*J-1)/(2*J+3)
@inline bJJ_S(J::Int) = 3/2*(J+1)*(J+2)/(2J+3)/(2J+1)
@inline bJJ_O(J::Int) = 3/2*J*(J-1)/(2J+1)/(2J-1)

# from Appendix B of: Role of intramolecular interactions in Raman spectra of N 2 and O 2 molecules - Buldakov
angstrom3_to_cm3 = 1e-24
cgs_to_si        = 1/8.98e15
@inline a²(v::Int; dadQ = 8.52e-42 )  = (v+1)/2*(2*1.99/2358)*(1.871+0.0105*v)^2 * (angstrom3_to_cm3 * cgs_to_si)^2 # or # Table 5.2(g)) in The Raman Effect A Unified Treatment of the Theory of Raman Scattering by Molecules - Derek A. Long
@inline γ²(v::Int; dgammadQ=1.02e-41) = (v+1)/2*(2*1.99/2358)*(2.25+0.019*v)^2 * (angstrom3_to_cm3 * cgs_to_si)^2

function αα(initial::State, final::State; use_Herman_Wallis = true)
    is_Q_branch = v(final) == v(initial) + 1 && J(final) == J(initial)
    is_O_branch = v(final) == v(initial) + 1 && J(final) == J(initial)-2
    is_S_branch = v(final) == v(initial) + 1 && J(final) == J(initial)+2

    if is_Q_branch
        α²_term = a²(v(initial))
        γ²_term = bJJ(J(initial))*4/45*γ²(v(initial))
        
        # Herman-Wallis factor
        m  = J(initial) * (J(initial) + 1)
        Fα = 1 + (1.1e-5  - 0.61e-7*v(initial))*m
        Fγ = 1 + (0.14e-4 - 0.12e-6*v(initial))*m
        return Fα^2 * α²_term + Fγ^2 * γ²_term
    elseif is_O_branch
        γ²_term = bJJ_O(J(initial))*4/45*γ²(v(initial))

        # Herman-Wallis factor
        m  = 1-2*J(initial)
        Fγ = 1 + 1.04e-5 - 0.09e-6*v(initial) - (2.2e-3 + 0.37e-4*v(initial))*m +(0.47e-5 + 0.11e-7*v(initial))*m
        return Fγ^2 * γ²_term
    elseif is_S_branch
        γ²_term = bJJ_S(J(initial))*4/45*γ²(v(initial))

        # Herman-Wallis factor
        m  = 2*J(initial)+3
        Fγ = 1 + 1.04e-5 - 0.09e-6*v(initial) - (2.2e-3 + 0.37e-4*v(initial))*m +(0.47e-5 + 0.11e-7*v(initial))*m
        return Fγ^2 * γ²_term
    else
        return 0.0
    end
end

function αα(t::N2Transition)
    αα(t.initial, t.final)
end