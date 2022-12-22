using Sorotoki
using Revise
using StaticArrays
using LinearAlgebra
using Test

## tests for Lie groups functions
@testset "liegroup.jl" begin
    v₁ = SVector(1,2,3)
    S = skew(v₁)
    @assert S * v₁ == zeros(3)
    @assert v₁ × v₁ == zeros(3)

    η = SVector(1,2,3,4,5,6)
    Ξ = hat(η)
    @assert det(Ξ) == 0.0
end

## testset for Shapes
#@testset "shapes.jl" begin
    modes = [0,5,0,0,0,0]
    
    #shp = Shapes()

#end
