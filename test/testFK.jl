using Sorotoki
using Revise
using Plots
using StaticArrays
#using SparseArrays
using LinearAlgebra
using Test
using BenchmarkTools

#@testset "liegroup.jl" begin
    #theme(:solarized)
    shp = Shapes(Modes=[0 3 2 0 0 0], NNode = 150)
    shp = shp.rebuild()

    ## plot beam
    q = rand(shp.NDim) * 0.075
    @time γ = shp.string(q);
    
    plot(γ[:,1],γ[:,2],γ[:,3],linewidth=3)

    L = sum(sqrt.(sum(diff(γ,dims=1).^2,dims=2)))
#end