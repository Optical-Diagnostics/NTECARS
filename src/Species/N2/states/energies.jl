# Dunham coefficients for calculating energy levels, load from file
DUNHAM_COEFFIECIENTS = reduce(vcat, permutedims.(TOML.parsefile(ROVIBRATIONAL_CONSTANTS_PATH)["Dunham_coefficients"]["Yij"]))

function energy(s::QuantumNumbers)
    energy = energy_from_dunham_coefficients(v(s),J(s))
    groundstate_energy = energy_from_dunham_coefficients(0,0)
    return energy - groundstate_energy
end

function energy_from_dunham_coefficients(v,J)
    energy = 0.0
    for k in 1:5
        for j in 1:4
            energy += DUNHAM_COEFFIECIENTS[k,j]*((v+0.5)^(k-1)*J^(j-1)*(J+1)^(j-1))
        end
    end
    return energy
end
