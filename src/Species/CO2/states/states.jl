######################################################################################################
#                                            Structs                                           
######################################################################################################
struct PartitionedEnergies
    E_total::Float32
    E_vib::Float32
    E_rot::Float32
    E_vib_12::Float32
    E_vib_3::Float32    
end

struct VibrationalState
    qn       ::VibrationalQuantumNumbers      # quantum numbers of the state
    energies ::PartitionedEnergies            # mode and total energies of the state
    gᵢ       ::Int64                          # state independant degeneracy
    polyad_basis::Vector{WangQuantumNumbers}
    polyad_coeff::Vector{Float64}
end

struct State
    qn       ::QuantumNumbers                 # quantum numbers of the state
    energies ::PartitionedEnergies            # mode and total energies of the state
    g        ::Int64                         # degeneracy of the state
    polyad_basis::Vector{WangQuantumNumbers}
    polyad_coeff::Vector{Float64}
end


total_energy(s::State) = s.energies.E_total

######################################################################################################
#                                            Constructor                                           
######################################################################################################
_CO2_metadata = MoleculeMetadata(:O16C12O16)
function State(qn::QuantumNumbers, mo_metadata::MoleculeMetadata = _CO2_metadata)
    # return State if the provided combination of term and quantum numbers lead to a symmetric state
    g_i = mo_metadata.gᵢ
    
    g_s = 1 # assumes that only allowed states are created.
    if mo_metadata.iso_ID == :O16C12O18 && qn.l₂ > 0
        g_s = 2
    end
    g  = g_i * g_s * (2*qn.J + 1)

    eigenstate                 = HamiltonianEigenstate(qn)
    polyad_basis, polyad_coeff = eigenstate.basis, eigenstate.basis_coeff
    partitioned_energies       = PartitionedEnergies(qn)

    return State(qn, partitioned_energies, g, polyad_basis, polyad_coeff)
end

function VibrationalState(vib_qn::VibrationalQuantumNumbers, mo_metadata::MoleculeMetadata = _CO2_metadata)
    # return State if the provided combination of term and quantum numbers lead to a symmetric state
    gᵢ = mo_metadata.gᵢ

    # get Quantum number with J=0
    @unpack v₁, v₂, l₂, v₃, r = vib_qn
    qn = QuantumNumbers(v₁, v₂, l₂, v₃, r, 0,1)
    eigenstate                 = HamiltonianEigenstate(qn)
    polyad_basis, polyad_coeff = eigenstate.basis, eigenstate.basis_coeff
    partitioned_energies       = PartitionedEnergies(qn)

    return VibrationalState(vib_qn, partitioned_energies, gᵢ, polyad_basis, polyad_coeff)
end

function PartitionedEnergies(qn::QuantumNumbers; use_Rothman = false, ignore_symmetry_parity = false)
    if use_Rothman
        E_vib    = HamiltonianEigenstate(QuantumNumbers(qn.v₁, qn.v₂, qn.l₂, qn.v₃, qn.r, 0, qn.C)).E
        E_rot    = rotational_energy_from_database(qn)
        E_total  = E_vib + E_rot
    else
        E_total  = HamiltonianEigenstate(qn).E
        E_vib    = HamiltonianEigenstate(QuantumNumbers(qn.v₁, qn.v₂, qn.l₂, qn.v₃, qn.r, 0, qn.C)).E
        E_rot    = E_total - E_vib 
    end

    E_vib_12 = HamiltonianEigenstate(QuantumNumbers(qn.v₁, qn.v₂, qn.l₂, 0, qn.r, 0, qn.C)).E
    E_vib_3  = E_vib - E_vib_12 
    return PartitionedEnergies(E_total, E_vib, E_rot, E_vib_12, E_vib_3)
end

######################################################################################################
#                                         Utilis for checks                                           
######################################################################################################

function print_wavefunction(s::State)
    qn = s.qn
    qn_string = "($(qn.v₁)$(qn.v₂)$(qn.l₂)$(qn.v₃)$(qn.r),$(qn.J))"

    wavefunction = "$(qn_string) =  "
    for i in eachindex(s.polyad_basis)
        qn = s.polyad_basis[i]
        qn_string = "($(qn.v₁)$(qn.v₂)$(qn.l₂)$(qn.v₃),$(qn.J))"
        wavefunction *= "$(round(s.polyad_coeff[i], digits = 3)) × $(qn_string)"
        if i != length(s.polyad_basis)
            wavefunction *= " + "
        end
    end
    println(wavefunction)

end