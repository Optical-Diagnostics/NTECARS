using Revise
using NTECARS

initial = CO2.State(CO2.QuantumNumbers(0,1,1,0,1,30, 1))
final   = CO2.State(CO2.QuantumNumbers(1,1,1,0,1,30, 1))

M = sqrt(CO2.αα(CO2.CO2Transition(initial, final)))