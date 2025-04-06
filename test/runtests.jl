using Con2Prim
using ForwardDiff
using Random
using StaticArrays
using Test

Random.seed!(0xb99747d6)
@testset "solve1d" begin
    f(x) = x^2 - 1
    xsol = 1.0
    for iter in 1:10
        xlo = rand(0.1:0.01:0.5)
        xhi = rand(1.5:0.01:2.5)
        @test solve1d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
    end
end

Random.seed!(0xadf8f6d0)
@testset "solve2d" begin
    f(x) = SVector(9 / 10 * x[1]^2 + 1 / 10 * x[2]^2 - 1, 1 / 10 * x[1]^2 + 9 / 10 * x[2]^2 - 1)
    xsol = SVector(1.0, 1.0)
    for iter in 1:10
        xlo = SVector(rand(0.1:0.01:0.5), rand(0.1:0.01:0.5))
        xhi = SVector(rand(1.5:0.01:2.5), rand(1.5:0.01:2.5))
        @test solve2d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
    end
end

Random.seed!(0x4809c1fe)
@testset "solve3d" begin
    f(x) = SVector(x[1]^2 - 1, x[2]^2 - 1, x[3]^2 - 1)
    xlo = SVector(0.0, 0.7, 0.3)
    xhi = SVector(1.3, 2.0, 1.6)
    xsol = SVector(1.0, 1.0, 1.0)
    @test solve3d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
end

@testset "con2prim" begin
    ρ = 1.0
    S = 1.0
    Ne = 0.0
    w = 0.1
    @show ρ, S, Ne, w

    D, De, M, T = prim2con(ρ, S, Ne, w)
    @show D, De, M, T

    prims = SVector(ρ, S, Ne, w)
    cons = prim2con(prims)
    dcons = ForwardDiff.jacobian(prim2con, prims)
    @show prims
    @show cons
    @show dcons

    atol = 1.0e-14

    ρ1, S1, Ne1, w1 = con2prim(D, De, M, T)
    @show ρ1, S1, Ne1, w1
    @test ρ1 ≈ ρ atol = atol
    @test S1 ≈ S atol = atol
    @test Ne1 ≈ Ne atol = atol
    @test w1 ≈ w atol = atol

    D1, De1, M1, T1 = prim2con(ρ1, S1, Ne1, w1)
    @show D1, De1, M1, T1
    @test D1 ≈ D atol = atol
    @test De1 ≈ De atol = atol
    @test M1 ≈ M atol = atol
    @test T1 ≈ T atol = atol
end
