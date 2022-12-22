module Sorotoki

    # Sorotoki base packages
    include("packages.jl")
    include("./base/base.jl")
    export skew, skewmat, hat, Id, ᵀ

    # Sorotoki packages for beam shapes
    include("./model/shapes.jl")
    export Shapes

    # Sorotoki packages for (hyper-elastic) materials
    include("./material/material.jl")
    include("./material/neohookean.jl")
    export Material, NeoHookean
end
