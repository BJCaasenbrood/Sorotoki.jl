module Sorotoki

    using LinearAlgebra
    using StaticArrays
    #using SparseArrays
    #using ForwardDiff
    using Base

    # Sorotoki base packages
    include("packages.jl")
    include("./base/base.jl")
    include("./base/solver.jl")
    include("./base/polyspace.jl")
    export chebyshev
    export Id, skew, iskew, hat, SE3
    export ad, Ad, Ad⁻¹
    
    include("./base/tuplemath.jl")
    export ⊗, ⊕

    # Sorotoki packages for beam shapes
    include("./model/shapes.jl")
    export Shapes, forwardKinematicsCosseratODE

    # Sorotoki packages for (hyper-elastic) materials
    include("./material/material.jl")
    include("./material/neohookean.jl")
    export Material, NeoHookean
end
