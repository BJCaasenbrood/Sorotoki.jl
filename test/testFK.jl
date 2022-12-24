using Sorotoki
using Revise
using Plots
using StaticArrays
using LinearAlgebra
using Test

#@testset "liegroup.jl" begin
    theme(:solarized)
    shp = Shapes(Modes=[0 8 4 0 0 0], NNode = 150)
    @time plot(shp.Xₙ(),shp.θₑ,linewidth=5)
#end