using Revise
#using NTECARS
using CairoMakie
using CSV
using DataFrames
CairoMakie.activate!()

qn_set = CO2.allowed_rotational_states(CO2.VibrationalQuantumNumbers(0,0,0,0,1))
Rothman_energies =  [CO2.PartitionedEnergies(s.qn, use_Rothman = true).E_total for s in qn_set]
Heff_energies    =  [CO2.PartitionedEnergies(s.qn, use_Rothman = false, ignore_symmetry_parity = true).E_total for s in qn_set]

fig = Figure(size = (600,750))
ax = Axis(fig[1,1])
scatter!(Rothman_energies)
scatter!(Heff_energies)
ax = Axis(fig[2,1])
scatter!(Rothman_energies .- Heff_energies)



qn_set = CO2.allowed_rotational_states(CO2.VibrationalQuantumNumbers(1,0,0,0,1))
Rothman_energies =  [CO2.PartitionedEnergies(s.qn, use_Rothman = true).E_total for s in qn_set]
Heff_energies    =  [CO2.PartitionedEnergies(s.qn, use_Rothman = false, ignore_symmetry_parity = true).E_total for s in qn_set]


ax = Axis(fig[3,1])
scatter!(Rothman_energies)
scatter!(Heff_energies)
ax = Axis(fig[4,1])
scatter!(Rothman_energies .- Heff_energies)


qn_set = CO2.allowed_rotational_states(CO2.VibrationalQuantumNumbers(0,2,2,0,1))
Rothman_energies =  [CO2.PartitionedEnergies(s.qn, use_Rothman = true).E_total for s in qn_set]
Heff_energies    =  [CO2.PartitionedEnergies(s.qn, use_Rothman = false, ignore_symmetry_parity = true).E_total for s in qn_set]


ax = Axis(fig[1,2])
scatter!(Rothman_energies)
scatter!(Heff_energies)
ax = Axis(fig[2,2])
scatter!(Rothman_energies .- Heff_energies)



qn_set = CO2.allowed_rotational_states(CO2.VibrationalQuantumNumbers(1,2,2,0,1))
Rothman_energies =  [CO2.PartitionedEnergies(s.qn, use_Rothman = true).E_total for s in qn_set]
Heff_energies    =  [CO2.PartitionedEnergies(s.qn, use_Rothman = false, ignore_symmetry_parity = true).E_total for s in qn_set]


ax = Axis(fig[3,2])
scatter!(Rothman_energies)
scatter!(Heff_energies)
ax = Axis(fig[4,2])
scatter!(Rothman_energies .- Heff_energies)

fig