using Revise
using NTECARS
using CairoMakie
using CSV
using DataFrames
CairoMakie.activate!()

df_HITRAN = CSV.read("data/validation/partition_sums/CO2.csv", delim=' ', ignorerepeated=true, DataFrame)

T_calc = df_HITRAN[1:10:end,1]
Q_calc = [CO2.MultiTemperatureDistribution(T_12 = T, T_3 = T, T_rot = T).Q for T in T_calc]

fig = Figure()
ax  = Axis(fig[1,1], limits = ((10, 10000),(10,1000000)), yscale = log10, xscale = log10, yminorticks = IntervalsBetween(9), yminorticksvisible = true,
      xminorticks = IntervalsBetween(9), xminorticksvisible = true)
lines!(ax, df_HITRAN[:,1], df_HITRAN[:,2], label = "HITRAN")
lines!(ax, T_calc, Q_calc, label = "This work")

ax2 = Axis(fig[2,1], limits = ((10,10000), nothing), xscale = log10, yscale = identity, yminorticks = IntervalsBetween(9), yminorticksvisible = true, 
yminorgridvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9), xminorticksvisible = true)
lines!(ax2, T_calc, abs.(df_HITRAN[1:10:end,2] .- Q_calc) ./ df_HITRAN[1:10:end,2] .* 100, label = "Relative error (%)")
fig