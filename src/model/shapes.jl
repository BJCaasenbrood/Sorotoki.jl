mutable struct Shapes
    NDim::Int16     # number of states/joints (dim(q) < Ƶ⁺)
    NDof::Int16     # number of degr. of freedom (NDof ≤ 6)
    NNode::Int16    # number of discrete points along g ∈ SE(3)
    Modes           # modes: [κ₁,κ₂,κ₃,ε₁,ε₂,ε₃]
    Basis           # list of (orthogonal) basis functions

    # domains and discrete samples
    X   # spatial domain 
    Xₙ  # normalize space domain 
    Δs  # spatial stepssize

    # zero-stress strain deformations 
    Γ₀  # zero tangent vector
    K₀  # zero curvature-torsion vector
    ξ₀  # zero geom. strain vector
    θ   # funtional basis (chebyshev 3rd-order by default)
    Bₐ  # matrix of active joints
    Iₐ  # matrix of modal number w.r.t basis {θ}ₖ

    # pre-eval matrices
    ξₑ  # evaluated strain vector
    θₑ  # evaluated POD matrix

    # material tensors 
    Ω₀  # reference volume Ω₀ ⊆ R³ (Sdf)
    L₀  # reference length L₀ > 0
    ρ₀  # volume density of volume ρ₀: Ω₀ → R₊
    Kᵪ  # stiffness tensor Kᵪ: Ω₀ × [0,T] → se*(3) × se*(3) (tangent)
    Mᵪ  # inertia tensor Mᵪ: Ω₀ × [0,T] → se*(3) × se*(3) 
    Rᵪ  # damping tensor Rᵪ: Ω₀ × [0,T] → se*(3) × se*(3) 

    rebuild::Function
    
    function Shapes(;
        Modes = [0,5,0,0,0,0],
        θ  = chebyshev(5),
        K₀ = SVector(0.,0.,0.),
        Γ₀ = SVector(1.,0.,0.),
        L₀ = 100.0,
        NNode = 100)
    
        this = new()
    
        @assert length(Modes) == 6 "Modes should be of size 6"
        @assert length(K₀) == 3 "Intrinsic curvature vector K₀ should be of size 3"
        @assert length(Γ₀) == 3 "Intrinsic linear vector Γ₀ should be of size 3"

        this.Modes = Modes
        this.NDof  = convert(Int16,sum(sign.(Modes)))
        this.NDim  = convert(Int16,sum(Modes))
        this.NNode = NNode
        this.K₀    = K₀
        this.Γ₀    = Γ₀
        this.ξ₀    = () -> vcat(this.K₀,this.Γ₀)
        this.L₀    = L₀

        this.θ = chebyshev(maximum(this.Modes))
        #------------------------------------------------------------------
        this.rebuild = function()
            #this.θ = 1
            this.Bₐ, this.Iₐ = buildActiveDofMatrix(this.Modes)
            this.X  = () -> range(0.0, this.L₀, this.NNode)
            this.Xₙ = () -> range(0.0, 1.0, this.NNode)
            this.Δs = this.L₀ / this.NNode 
            this.θₑ = this.θ(this.Xₙ())
            return this
        end
        #------------------------------------------------------------------

        this.rebuild()            

        return this
    end
end

function buildActiveDofMatrix(Modes)
    I6   = SMatrix{6,6}(Id(6.0));
    NDof = sum(Modes)
    Mmax = maximum(Modes)
    Xₐ   = []
    Xᵢ   = []
    # check which columns of I6 are needed
    # add to Xₐ for Modes[ii] > 0
    for ii = 1:6
        for jj = 1:Modes[ii]
            Xᵢ = push!(Xᵢ, ii)
            Xₐ = push!(Xₐ, jj)
        end
    end

    #print(Xₐ)
    return SMatrix{6,NDof}(Matrix(I6[:,Xᵢ])), Xₐ
end

function forwardKinematicsCosseratODE!(Z, shp::Shapes, k::Int64)

    γₖ, Rₖ, Jₖ = Z

    Adg = Ad(gₖ)
    Ξ   = shp.Bₐ * diagm(shp.θₑ[k, shp.Iₐ])
    ξ   = Ξ * q + ξₖ
    dξ  = Ξ * dq

    dγ = Rₖ * ξ[1:3]
    dR = Rₖ * skew(ξ[4:6])
    dJ = Adg * Ξ;

    Z = (dγ, dR, dJ)
end