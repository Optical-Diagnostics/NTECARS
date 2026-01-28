abstract type AbstractQuantumNumers end
abstract type AbstractState end

struct QuantumNumbers <: AbstractQuantumNumers
    v::UInt8
    J::UInt16
end

struct State
    qn   ::QuantumNumbers
    E_tot::Float64
    E_vib::Float64
    E_rot::Float64
    degen::Int64
end

v(s::QuantumNumbers)::Int = s.v
J(s::QuantumNumbers)::Int = s.J
v(s::State)::Int = s.qn.v
J(s::State)::Int = s.qn.J

total_energy(s::State) = s.E_tot

degeneracy(s::QuantumNumbers) = isodd(J(s)) ? 3*(2*J(s)+1) : 6*(2*J(s)+1)

function State(v, J)
    qn    = QuantumNumbers(v, J)
    degen = degeneracy(qn)
    E_tot = energy(qn)
    E_vib = energy(QuantumNumbers(v,0))
    E_rot = E_tot - E_vib
    s     = State(qn, E_tot, E_vib, E_rot, degen)
    return s
end

function states(v_max, J_max)
    states = []
    for v in 0:v_max, J in 0:J_max
        push!(states, State(v, J))
    end
    return states
end 

function vib_states(v_max)
    return states(v_max, 0)
end 

