
@testset "Multi-temperature partition sums" begin
    project_path = joinpath(@__DIR__, "..")
    file_path    = joinpath(project_path, "data", "validation", "partition_sums", "N2.csv")
    df_HITRAN    = CSV.read(file_path, delim=' ', ignorerepeated=true, DataFrame)
    T_HITRAN, Q_HITRAN = df_HITRAN[300:5:1500,1], df_HITRAN[300:5:1500,2]
    Q_calc = [N2.MultiTemperatureDistribution(T_vib = T, T_rot = T).Q for T in Float64.(T_HITRAN)]
    @test all(abs.((Q_HITRAN .- Q_calc)./Q_HITRAN) .< 0.01)
end