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
    g₀  # initial base frame
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
    string::Function
    assemble::Function
    
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
        this.g₀    = SE3(Id(3),zeros(3,1))

        this.θ = chebyshev(maximum(this.Modes))
        #------------------------------------------------------------------
        this.rebuild = function()
            this.Bₐ, this.Iₐ = buildActiveDofMatrix(this.Modes)
            this.X  = (NStep) -> range(0.0, this.L₀, NStep)
            this.Xₙ = (NStep) -> range(0.0, 1.0, NStep)
            this.Δs = this.L₀ / this.NNode 
   
            # this.θₑ, ~ = gsogpoly(collect(this.θ(this.Xₙ())),
            #                    collect(this.X()))
            return this
        end
        #------------------------------------------------------------------
        this.string = function(q; ateach::Int64 = 1)

            γ = Matrix{Float64}(undef,this.NNode,3)
            #γ = []
            q = collect(q)  # ensure it is a column vector
            Z = (this.g₀[1:3,4], this.g₀[1:3,1:3], zeros(6,this.NDim))
            Q = (q, q)

            # set intermediate flow function for forward kinematics
            F = (Z,s) -> forwardKinematicsCosseratODE!(Z, Q, s, this.θ, this.Bₐ, this.Iₐ)
            S = this.Xₙ(this.NNode)

            @inbounds for i = 1:this.NNode 
                # leap-frog on tuple object
                Z_ = Z ⊕ (this.Δs ./ 2.0) ⊗ F(Z, S[i])
                Z  = Z ⊕ this.Δs ⊗ F(Z_, S[i] + this.Δs ./ 2.0)

                # forward euler on tuple object
                #Z = Z ⊕ this.Δs ⊗ F(Z, S[kk])

                if i == 1 || i == this.NNode || i % ateach == 0   
                    #push!(γ, Z[1])
                    γ[i,:] = Z[1]
                end
            end

            #return collect(reduce(hcat,γ)')
            return γ
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

    return SMatrix{6,NDof}(Matrix(I6[:,Xᵢ])), Xₐ
end

function forwardKinematicsCosseratODE!(Z, Q, s, PolyBasis, activeDof, structDof)
    γₖ, Rₖ = Z
    q, ~  = Q

    # compute the geometric strain vector
    Ξ  = activeDof * diagm(PolyBasis(s)[structDof]) 
    ξ  = Ξ * q 
    ξ[4] += 1.0

    # matrix differential equation for FK
    dγ = Rₖ * ξ[4:6]
    dR = Rₖ * skew(ξ[1:3])
    dJ = Ad( SE3(Rₖ, γₖ) ) * Ξ;

    # compose Tuple of state update
    return Z = (dγ, dR, dJ) 
end

function forwardDynamicsCosseratODE!(Z, Q, s, PolyBasis, activeDof, structDof, BeamTensor)
    
    q, dq = Q
    γₖ, Rₖ, Jₖ, dJdtₖ = Z 
    Mtt, Ktt = BeamTensor

    Ξ = activeDof * diagm(PolyBasis(s)[structDof]) 
    ξ = Ξ * q
    ξ[4] =+ 1.0

    A  = Ad( SE3(Rₖ, γₖ) )
    Ai = Ad⁻¹( SE3(Rₖ, γₖ) )

    # matrix differential equation for FK
    dγ   = Rₖ * ξ[4:6]
    dR   = Rₖ * skew(ξ[1:3])

    Jg    = Ai * Jₖ
    dJgdt = Ai * dJdtₖ
    η     = Jg * dq 
    adη   = ad(η)
    #adξ   = ad(ξ)

    dJ   = A * Ξ
    dJdt = A * adη * Ξ

    # body acceleration vector
    ag = [0,0,0,0,0,0]

    # matrix differential equation for EL
    dM =  Jg' * Mtt * Jg
    dC =  Jg' * ((Mtt * adη - adη' * Mtt) * Jg + Mtt * dJgdt)
    dG = -Jg' * (Ai * Mtt * ag)

    # tensor or energy integration
    dK  = (Ξ)' * Ktt * (Ξ) 
    dUe =  0.5 .* η' * Mtt * η
    dUg = Mtt[4,4] * γₖ * ag[4:6]

    # compose Tuple of state matrices update
    return Z = (dγ, dR, dJ, dJdt, dM, dC, dG) 
end