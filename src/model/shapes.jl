mutable struct Shapes
    NDim::Int8
    NJoint::Int8
    
    # intrinsic strain
    Γ₀::SVector{3}
    K₀::SVector{3}
    ξ₀

    # constructor
    function Shapes()
        this = new()

        this.K₀ = @SVector [0.0, 0.0, 0.0]
        this.Γ₀ = @SVector [1.0, 0.0, 0.0]
        this.ξ₀ = vcat(this.K₀,this.Γ₀)

        return this
    end
end