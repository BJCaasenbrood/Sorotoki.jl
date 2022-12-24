module Sorotoki
    # Sorotoki base packages
    include("packages.jl")
    include("./base/base.jl")
    export Id, skew, iskew, hat
    export ad, Ad, Ad⁻¹
    export chebyshev

    # Sorotoki packages for beam shapes
    include("./model/shapes.jl")
    export Shapes, forwardKinematicsCosseratODE

    # Sorotoki packages for (hyper-elastic) materials
    include("./material/material.jl")
    include("./material/neohookean.jl")
    export Material, NeoHookean
end
