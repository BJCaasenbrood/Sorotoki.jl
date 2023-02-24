using Sorotoki
using Revise
using Plots
# using StaticArrays
# using LinearAlgebra
using Test
using BenchmarkTools

#@testset "liegroup.jl" begin
    #theme(:solarized)
    shp = Shapes(Modes=[0 4 0 0 0 0], NNode = 250)
    shp = shp.rebuild()

    # do forward kinematics
    # q = 0.05 * (2 * (rand(shp.NDim) - 0.5 * ones(shp.NDim)))
    q = [0.05; 0.01; 0.01; 0.01]
    @time γ = shp.string(q);

    ## plot beam
    plot(γ[:,1],γ[:,3],linewidth=3,scale=:identity)

    L = sum(sqrt.(sum(diff(γ,dims=1).^2,dims=2)))
#end