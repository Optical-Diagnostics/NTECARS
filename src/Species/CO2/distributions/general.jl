abstract type CO2Distribution end

function (df::CO2Distribution)(state::State)
    #fraction of CO2 molecules in the given state
    return state.g/df.Q * fvib(df, state) * frot(df, state)
end

function (df::CO2Distribution)(state::VibrationalState)
    #fraction of CO2 molecules in the given vibrational state
    f_vib = fvib(df, state)
    Q_rot = rotational_partition_sum_McDowell(df.T_rot, state.qn.l₂, df.iso_ID)
   
    return state.gᵢ/df.Q * f_vib * Q_rot
end

function partition_sum(df::CO2Distribution)
    Q = 0.0 # partition sum
    C = get_CO2_energy_constants() # load once because it is a bit slow
    for vibstate in PRECALCULATED_VIBSTATES[df.iso_ID]
        f_vib = fvib(df, vibstate)
        Qᵣ    = rotational_partition_sum_McDowell(df.T_rot, vibstate.qn, df.iso_ID, C)
        Q    += vibstate.gᵢ * Qᵣ * f_vib
    end
    Q
end



function rotational_partition_sum_McDowell(T_rot, qn, iso_ID, C::VibRotConstants = get_CO2_energy_constants())
    v₁, v₂, v₃, l₂ = qn.v₁, qn.v₂, qn.v₃, qn.l₂
    
    # calculate rotational constants of diagonal effective hamiltonian
    gi      = [1,2,1]
    vi      = [v₁,v₂,v₃]
    vg_i    = @. vi + gi/2
    vg_ij   = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) for i in 1:3, j in 1:3]
    vg_ijk  = [(vi[i]+gi[i]/2) * (vi[j]+gi[j]/2) * (vi[k]+gi[k]/2) for i in 1:3, j in 1:3, k in 1:3]

    B = C.Be - sum(C.αi .* vg_i) + sum(C.γij .* vg_ij) + C.γll * l₂^2 + sum(C.εijk .* vg_ijk) + sum(C.εill .* vg_i .* l₂^2)
    D = C.De + sum(C.βi .* vg_i) + sum(C.ηij .* vg_ij) + C.ηll * l₂^2
    H = C.He + sum(C.δi .* vg_i)
    
    # The partition sum of all rotational levels is Qᵣ = ∑_{J ∈ 0:1:∞} (2J+1) * exp(-Erot/kT) = exp(β/3)/β * (1 + β^2/90 + 8*β^3/2835) 
    β  = c₂ * B / T_rot
    Qᵣ = exp(β/3)/β * (1 + β^2/90 + 8*β^3/2835 + β^4*5/4536 + β^5*148/280655 )

    # The additional rotational constans D, H can be included by multiplying the partiotion sum with a correction factor
    d = D / B
    h = H / B
    centrifugal_correction = 1 + 2*d*(3-β)/(3*β) + 6*(2*d^2 - h) / β^2 + 120*d*(d^2-h)/β^3
    Qᵣ *= centrifugal_correction

    if l₂ > 0
        # degenerate bending mode, For each J there is a + and a - rotational level
        Qᵣ *= 2
    end
    
    if iso_ID == :O16C12O16 || iso_ID == :O16C13O16
        # For isotopologues that are linear and symmetric, either even or odd J result in a symmetric state
        # the sums over only odd or even J are approximately equal and half the sum over all J.
        Qᵣ /= 2
    end
    return Qᵣ
end

function boltzmann(E, T)
    exp(-c₂ * E/T)
end
