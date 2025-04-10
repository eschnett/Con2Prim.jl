using Con2Prim
using ForwardDiff
using Random
using StaticArrays
using Test

Random.seed!(0xbe4b7c79)
@testset "boxes" begin
    bounds = Bounds{2,Float64}()
    add_point!(bounds, SVector(0, 0))
    add_point!(bounds, SVector(1, 2))
    add_point!(bounds, SVector(2, 1))
    @test all(length(bounds.bounds[d]) == 5 for d in 1:2)

    boxes = Boxes{2,Float64,Bool}(bounds)
    @test all(size(boxes.values) .== 4)

    set_past!(boxes, SVector(1, 1), true)
    @test boxes.values == Bool[1 1 0 0; 1 1 0 0; 0 0 0 0; 0 0 0 0]

    set_future!(boxes, SVector(1, 2), true)
    @test boxes.values == Bool[1 1 0 0; 1 1 0 0; 0 0 0 1; 0 0 0 1]
end

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
    f(x) = SVector(9/10 * x[1]^2 + 1/10 * x[2]^2 - 1, 1/10 * x[1]^2 + 9/10 * x[2]^2 - 1)
    xsol = SVector(1.0, 1.0)
    for iter in 1:10
        xlo = SVector(rand(0.1:0.01:0.5), rand(0.1:0.01:0.5))
        xhi = SVector(rand(1.5:0.01:2.5), rand(1.5:0.01:2.5))
        @test solve2d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
    end
end

Random.seed!(0x4809c1fe)
@testset "solve3d" begin
    f(x) = SVector(
        #TODO 97/100*x[1]^2 + 1/100*x[2]^2 + 2/100*x[3]^2 - 1,
        #TODO 1/100*x[1]^2 + 97/100*x[2]^2 + 2/100*x[3]^2 - 1,
        #TODO 1/100*x[1]^2 + 2/100*x[2]^2 + 97/100*x[3]^2 - 1,
        (x[1]-1) + 1/10*(x[2]-1) + 2/10*(x[3]-1),
        1/10 * (x[1]-1) + (x[2]-1) + 2/10*(x[3]-1),
        1/10 * (x[1]-1) + 2/10*(x[2]-1) + (x[3]-1),
    )
    xsol = SVector(1.0, 1.0, 1.0)
    for iter in 1:1 #TODO 10
        xlo = SVector(rand(0.1:0.01:0.5), rand(0.1:0.01:0.5), rand(0.1:0.01:0.5))
        xhi = SVector(rand(1.5:0.01:2.5), rand(1.5:0.01:2.5), rand(1.5:0.01:2.5))
        #TODO xlo = SVector(0.5, 0.5, 0.5)
        #TODO xhi = SVector(1.5, 1.5, 1.5)
        if !isapprox(solve3d(f, xlo, xhi).x, xsol; atol=1.0e-13)
            sol=solve3d(f, xlo, xhi).x
            @show xlo xhi sol f(sol)
        end
        @assert isapprox(solve3d(f, xlo, xhi).x, xsol; atol=1.0e-13)
        @test solve3d(f, xlo, xhi).x ≈ xsol atol = 1.0e-13
    end
end

#TODO @testset "con2prim" begin
#TODO     ρ = 1.0
#TODO     S = 1.0
#TODO     Ne = 0.0
#TODO     w = 0.1
#TODO     @show ρ, S, Ne, w
#TODO 
#TODO     D, De, M, T = prim2con(ρ, S, Ne, w)
#TODO     @show D, De, M, T
#TODO 
#TODO     prims = SVector(ρ, S, Ne, w)
#TODO     cons = prim2con(prims)
#TODO     dcons = ForwardDiff.jacobian(prim2con, prims)
#TODO     @show prims
#TODO     @show cons
#TODO     @show dcons
#TODO 
#TODO     atol = 1.0e-14
#TODO 
#TODO     ρ1, S1, Ne1, w1 = con2prim(D, De, M, T)
#TODO     @show ρ1, S1, Ne1, w1
#TODO     @test ρ1 ≈ ρ atol = atol
#TODO     @test S1 ≈ S atol = atol
#TODO     @test Ne1 ≈ Ne atol = atol
#TODO     @test w1 ≈ w atol = atol
#TODO 
#TODO     D1, De1, M1, T1 = prim2con(ρ1, S1, Ne1, w1)
#TODO     @show D1, De1, M1, T1
#TODO     @test D1 ≈ D atol = atol
#TODO     @test De1 ≈ De atol = atol
#TODO     @test M1 ≈ M atol = atol
#TODO     @test T1 ≈ T atol = atol
#TODO end
