using Revise
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using NTECARS
using CairoMakie
using CSV
using DataFrames
CairoMakie.activate!()

df_HITRAN = CSV.read("data/validation/partition_sums/N2.csv", delim=' ', ignorerepeated=true, DataFrame)
Q_calc = [N2.MultiTemperatureDistribution(T_vib = T, T_rot = T).Q for T in Float64.(df_HITRAN[:,1])]

fig = Figure()
ax  = Axis(fig[1,1], limits = ((1, 10000),(1,50000)), yscale = log10, xscale = log10, yminorticks = IntervalsBetween(9), yminorticksvisible = true,
      xminorticks = IntervalsBetween(9), xminorticksvisible = true)
lines!(ax, df_HITRAN[:,1], df_HITRAN[:,2], label = "HITRAN")
lines!(ax, df_HITRAN[:,1], Q_calc, label = "This work")

ax2 = Axis(fig[2,1], limits = ((1,10000), nothing), xscale = log10)
lines!(ax2, df_HITRAN[:,1], abs.(df_HITRAN[:,2] .- Q_calc) ./ df_HITRAN[:,2] .* 100, label = "Relative error (%)")
fig