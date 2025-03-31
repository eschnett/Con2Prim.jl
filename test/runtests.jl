using Con2Prim
using ForwardDiff
using StaticArrays
using Test

@testset "solve1d" begin
    f(x) = x^2 - 1
    xlo = 0.0
    xhi = 1.3
    xsol = 1.0
    @test solve1d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
end

@testset "solve2d" begin
    f(x) = SVector(x[1]^2 - 1, x[2]^2 - 1)
    xlo = SVector(0.0, 0.7)
    xhi = SVector(1.3, 2.0)
    xsol = SVector(1.0, 1.0)
    @test solve2d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
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
