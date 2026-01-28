######################################################################################################
#                                          Quantum numbers                                           
######################################################################################################

abstract type AbstractQuantumNumbers end

struct HarmonicOscillatorQuantumNumbers <: AbstractQuantumNumbers
    # quantum numbers of harmonic oscillator, used for calc of effective hamiltonian
    v₁ ::Int16 # quantum number of symmetric stretch mode
    v₂ ::Int16 # quantum number of bending mode
    l₂ ::Int16 # effective rotation of bending modes
    v₃ ::Int16 # quantum number of asymmetric stretch mode
    J  ::Int16 
end

struct QuantumNumbers <: AbstractQuantumNumbers
    # full set of quantum numbers describing a state of CO2
    v₁ ::Int16 # quantum number of symmetric stretch mode
    v₂ ::Int16 # quantum number of bending mode
    l₂ ::Int16 # effective rotation of bending modes
    v₃ ::Int16 # quantum number of asymmetric stretch mode
    r  ::Int16 # ranking number in polyad. r==1 has highest energy.
    J  ::Int16 # rotational quantum number
    C  ::Int16 # paraity: 1 for "e" and 2 for "f"
end

struct VibrationalQuantumNumbers <: AbstractQuantumNumbers
    # quantum numbers of a vibrational state, used for culculating partition sums
    v₁ ::Int16 # quantum number of symmetric stretch mode
    v₂ ::Int16 # quantum number of bending mode
    l₂ ::Int16 # effective rotation of bending modes
    v₃ ::Int16 # quantum number of asymmetric stretch mode
    r  ::Int16 # ranking number in polyad. r==1 has highest energy.
end

function VibrationalQuantumNumbers(qn::QuantumNumbers)
    return VibrationalQuantumNumbers(qn.v₁, qn.v₂, qn.l₂, qn.v₃, qn.r)
end

function CO2_parity(J, l₂, v₃)
    # parity that an allowd state has at these quantum numbers
    return (J+l₂+v₃) % 2 == 0 ? 1 : 2
end


Base.show(io::IO, qn::VibrationalQuantumNumbers)        = print(io, "($(qn.v₁),$(qn.v₂),$(qn.l₂),$(qn.v₃),$(qn.r))")
Base.show(io::IO, qn::QuantumNumbers)                   = print(io, "($(qn.v₁),$(qn.v₂),$(qn.l₂),$(qn.v₃),$(qn.r), J=$(qn.J), C=$(qn.C))")
Base.show(io::IO, qn::HarmonicOscillatorQuantumNumbers) = print(io, "($(qn.v₁),$(qn.v₂),$(qn.l₂),$(qn.v₃),$(qn.J))")
