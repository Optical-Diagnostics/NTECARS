@enum gerade_ungerade_symmetry begin
    gerade   = 1
    ungerade = -1
end

@enum positive_negative_symmetry begin
    positive       = 1
    p_m_degenerate = 0
    negative       = -1
end

struct VibrationalTerm 
    angular_momentum ::Int8
    g_u_symmetry     ::gerade_ungerade_symmetry 
    p_n_symmetry     ::positive_negative_symmetry
end

function allowed_term(qn::QuantumNumbers)
    # if a term exists for the given quantum numbers, that leads to asymmetric state, return it.
    # Otherwisee return nothing.
    for term in possible_vibrational_terms(VibrationalQuantumNumbers(qn.v₁, qn.v₂, qn.l₂, qn.v₃, qn.r))
        if is_symmetric(term, qn.J)
            return term
        end
    end
    return nothing
end

function possible_vibrational_terms(qn::VibrationalQuantumNumbers, return_pm_degenerate_state = true)
    # returns possible terms for the given quantum numbers, which can then be used to check for the symmetry 
    # of the states.
    @unpack v₁, v₂, l₂, v₃, r = qn

    v₁_term = VibrationalTerm(0, gerade, positive)
    v₂_term = VibrationalTerm(l₂, gerade_ungerade_symmetry((-1)^v₂), l₂>0 ? p_m_degenerate : positive)
    v₃_term = VibrationalTerm(0, gerade_ungerade_symmetry((-1)^v₃), positive) 

    vib_term = v₁_term * v₂_term * v₃_term # symmetry term multiplication is defined below

    if return_pm_degenerate_state
        [vib_term] 
    else
        if vib_term.p_n_symmetry == p_m_degenerate
            # with l₂ > 0 two degenrate bending modes exist
            positive_term = VibrationalTerm(vib_term.angular_momentum, vib_term.g_u_symmetry, positive)
            negative_term = VibrationalTerm(vib_term.angular_momentum, vib_term.g_u_symmetry, negative)
            return [positive_term, negative_term]
        else
            return [vib_term] 
        end
    end
end

function is_symmetric(term::VibrationalTerm, J)
    J_is_positive = (-1)^J==1

    if term.p_n_symmetry == positive
        if J_is_positive
            return term.g_u_symmetry == gerade ? true : false
        else
            return term.g_u_symmetry == gerade ? false : true
        end
    elseif term.p_n_symmetry == negative
        if J_is_positive
            return term.g_u_symmetry == gerade ? false : true
        else
            return term.g_u_symmetry == gerade ? true : false
        end
    elseif term.p_n_symmetry == p_m_degenerate
        return true
    end
end


Base.copy(t::VibrationalTerm) = VibrationalTerm(t.angular_momentum, t.g_u_symmetry, t.p_n_symmetry)

function Base.:*(term1::VibrationalTerm, term2::VibrationalTerm)
    l₂           = term1.angular_momentum + term2.angular_momentum
    g_u_symmetry = term1.g_u_symmetry * term2.g_u_symmetry
    p_m_symmetry = l₂>0 ? p_m_degenerate : term1.p_n_symmetry * term2.p_n_symmetry
    return VibrationalTerm(l₂, g_u_symmetry, p_m_symmetry)
end

function Base.:^(term::VibrationalTerm, n::Int64)
    t = copy(term)
    if n>1
        for i in 2:n
            t = t*term
        end
    end
    return t
end

function Base.:*(gu1::gerade_ungerade_symmetry, gu2::gerade_ungerade_symmetry)
    gerade_ungerade_symmetry(Int(gu1) * Int(gu2))
end

function Base.:*(pm1::positive_negative_symmetry, pm2::positive_negative_symmetry)
    positive_negative_symmetry(Int(pm1) * Int(pm2))
end