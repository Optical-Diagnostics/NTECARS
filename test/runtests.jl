using NTECARS
using Test
using CSV
using DataFrames

@testset "CO2" begin
    include("CO2.jl")
    include("N2.jl")
end
