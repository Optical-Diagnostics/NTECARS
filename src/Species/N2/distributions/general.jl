abstract type N2Distribution end

function Boltzmann_factor(E_cm⁻1,T)
    exp(-c₂ * E_cm⁻1/T)
end

function rot_partition_sum(T_rot)
    odd_Js  = 3/2 * rotational_partition_sum_McDowell(T_rot) # g(J) = 3(2J + 1) when J is odd
    even_Js = 6/2 * rotational_partition_sum_McDowell(T_rot) # g(J) = 6(2J + 1) when J is even
    return odd_Js + even_Js
end

function rotational_partition_sum_McDowell(T_rot, B=1.99, D=5.76e-6, H=0.0)
    # Formulas from: "Rotational partition functions for linear molecules", McDowell 1988.
    # The partition sum of all rotational levels is 
    # Qᵣ = ∑_{J ∈ 0:1:∞} (2J+1) * exp(-Erot/kT) = exp(β/3)/β * (1 + β^2/90 + 8*β^3/2835) 
    β = c₂ * B / T_rot
    d = D / B
    h = H / B
    centrifugal_correction = 1 + 2*d*(3-β)/(3*β) + 6*(2*d^2 - h) / β^2 + 120*d*(d^2-h)/β^3
    Qᵣ = exp(β/3)/β * (1 + β^2/90 + 8*β^3/2835) * centrifugal_correction

    return Qᵣ
end



