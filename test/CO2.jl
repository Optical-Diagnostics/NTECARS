
@testset "Parities" begin
    @test CO2.QuantumNumbers(0,0,0,0,1,0).C == 1
    @test CO2.QuantumNumbers(0,0,0,0,1,2).C == 1
    @test CO2.QuantumNumbers(0,0,0,0,1,4).C == 1

    @test CO2.QuantumNumbers(1,0,0,0,1,0).C == 1
    @test CO2.QuantumNumbers(1,0,0,0,1,2).C == 1
    @test CO2.QuantumNumbers(1,0,0,0,1,4).C == 1

    @test CO2.QuantumNumbers(0,1,1,0,1,1).C == 1
    @test CO2.QuantumNumbers(0,1,1,0,1,2).C == 2
    @test CO2.QuantumNumbers(0,1,1,0,1,3).C == 1
    @test CO2.QuantumNumbers(0,1,1,0,1,4).C == 2

    @test CO2.QuantumNumbers(0,1,1,1,1,0).C == 1
    @test CO2.QuantumNumbers(0,1,1,1,1,1).C == 2
    @test CO2.QuantumNumbers(0,1,1,1,1,2).C == 1
    @test CO2.QuantumNumbers(0,1,1,1,1,3).C == 2

    @test CO2.QuantumNumbers(0,0,0,0,1,0) == CO2.QuantumNumbers(0,0,0,0,1,0,1)
end

@testset "Energies" begin
    @test CO2.HamiltonianEigenstate(CO2.QuantumNumbers(0,0,0,0,1,0)).E == 0.0
    @test CO2.HamiltonianEigenstate(CO2.QuantumNumbers(0,0,0,0,1,100)).E > 3000
    @test 665 ≤ CO2.HamiltonianEigenstate(CO2.QuantumNumbers(0,1,1,0,1,0)).E  ≤ 670
    @test 1280 ≤ CO2.HamiltonianEigenstate(CO2.QuantumNumbers(1,0,0,0,2,0)).E  ≤ 1290
    @test 1385 ≤ CO2.HamiltonianEigenstate(CO2.QuantumNumbers(1,0,0,0,1,0)).E  ≤ 1390
    @test 2340 ≤ CO2.HamiltonianEigenstate(CO2.QuantumNumbers(0,0,0,1,1,0)).E  ≤ 2360
end

@testset "Wavefunctions" begin
    ψ = CO2.HamiltonianEigenstate(CO2.QuantumNumbers(0,0,0,0,1,0))
    @test ψ.basis == [NTECARS.CO2.WangQuantumNumbers(0, 0, 0, 0, 0, 1)]
    @test ψ.basis_coeff == [1.0] || ψ.basis_coeff == [-1.0]

    ψ = CO2.HamiltonianEigenstate(CO2.QuantumNumbers(1,0,0,0,1,0))
    @test NTECARS.CO2.WangQuantumNumbers(1, 0, 0, 0, 0, 1) in ψ.basis
    @test NTECARS.CO2.WangQuantumNumbers(0, 2, 2, 0, 0, 1) in ψ.basis
    @test NTECARS.CO2.WangQuantumNumbers(0, 2, 0, 0, 0, 1) in ψ.basis
    @test sum(ψ.basis_coeff .^ 2) ≈ 1
    @test abs(sum(ψ.basis_coeff)) < 1 # indirect check that [-0.67, 0.0, 0.73] or with inverse sign 
end

@testset "Multi-temperature partition sums" begin
    project_path = joinpath(@__DIR__, "..")
    file_path    = joinpath(project_path, "data", "validation", "partition_sums", "CO2.csv")
    df_HITRAN    = CSV.read(file_path, delim=' ', ignorerepeated=true, DataFrame)
    T_HITRAN, Q_HITRAN = df_HITRAN[300:5:1500,1], df_HITRAN[300:5:1500,2]
    Q_calc = [CO2.MultiTemperatureDistribution(T_12 = T, T_3 = T, T_rot = T).Q for T in T_HITRAN]
    @test all(abs.((Q_HITRAN .- Q_calc)./Q_HITRAN) .< 0.01)
end




