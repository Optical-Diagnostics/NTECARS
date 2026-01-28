using Revise
using NTECARS


qn = CO2.QuantumNumbers(0,0,0,0,1,1)
ψ  = CO2.HamiltonianEigenstate(qn)
CO2.print_wavefunction(ψ)