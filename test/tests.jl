using Sorotoki
using Revise
using StaticArrays
using LinearAlgebra
using Test

## tests for Lie groups functions
#@testset "liegroup.jl" begin
    p₁ = SVector(1,2,3)
    p₂ = SVector(3,2,1)
    S = skew(p₁)
    @test S[1,2] == -S[2,1]
    @test S[1,3] == -S[3,1]
    @test S[2,3] == -S[3,2]
    @test S * p₁ == zeros(3)
    @test p₁ × p₂ == S * p₂

    η = SVector(1,2,3,4,5,6)
    Ξ = hat(η)
    @test det(Ξ) == 0.0
    adₓ = ad(η)
    @test adₓ * η == @SVector zeros(6)
#end

## testset for Shapes
#@testset "shapes.jl" begin
    #modes = [0,5,0,0,0,0]
    
    #shp = Shapes()

#end
