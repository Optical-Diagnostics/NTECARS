using Revise
using NTECARS
using CairoMakie
using CSV
using DataFrames
CairoMakie.activate!()

J = 60
qn = CO2.QuantumNumbers(1,2,2,0,1,J)
p = CO2.polyad_quantum_number(qn)
basis = CO2.polyad_basis(p, J)
H_eff = CO2.effective_Hamiltonian(basis)

H_eff