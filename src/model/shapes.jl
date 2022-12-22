mutable struct Shapes
    NDim::Int16     # number of states/joints (dim(q) < Ƶ⁺)
    NDof::Int16     # number of degr. of freedom (NDof ≤ 6)
    Modes           # modes: [κ₁,κ₂,κ₃,ε₁,ε₂,ε₃]

    # zero-stress strain deformations 
    Γ₀::SVector{3}  # zero tangent vector
    K₀::SVector{3}  # zero curvature-torsion vector
    ξ₀::SVector{6}  # zero geom. strain vector

    # pre-eval matrices
    ξᵉ
    Θᵉ

    # material tensors 
    Ω₀  # reference volume Ω₀ ⊆ R³ (Sdf)
    L₀  # reference length L₀ > 0
    ρ₀  # volume density of volume ρ₀: Ω₀ → R₊
    Kᵪ  # stiffness tensor Kᵪ: Ω₀ × [0,T] → se*(3) × se*(3) (tangent)
    Mᵪ  # inertia tensor Mᵪ: Ω₀ × [0,T] → se*(3) × se*(3) 
    Rᵪ  # damping tensor Rᵪ: Ω₀ × [0,T] → se*(3) × se*(3) 
    
    function Shapes(;
        Modes = [0,1,0,0,0,0])
    
        this = new()

        print(typeof(Modes))
    
        if length(Modes) != 6
            msg = "Modes should be a ::Vector{6}, ::Tuple{6}, or ::Matrix{6}"
            throw(ArgumentError(msg))
        end
        
        # construct vector of intrinsic strains
        this.ξ₀ = @SVector [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    
        this.Modes = Modes
        this.NDof  = convert(Int16,sum(sign.(Modes)))
        this.NDim  = convert(Int16,sum(Modes))
    
        return this
    end
end