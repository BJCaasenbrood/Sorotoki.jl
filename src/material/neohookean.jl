mutable struct NeoHookean <: Material     # subtype of Material
    E       # linear Young's modulus (MPa)
    ν       # Poisson ratio -1 ≤ ν₀ < .5 (-)
    ρ       # volumetric density (kg/mm³)
    λ       # Lame constant Eν / (1+ν)(1-2ν)
    μ       # Lame constant E / 2(1+ν)

    getModulus::Function
    setModulus::Function

    function NeoHookean(; 
        E = 1.0,
        ν = 0.49,
        ρ = 1200E-12)

        # new struct 
        this = new()
        
        # assign material properties
        this.E = E
        this.ν = ν
        this.ρ = ρ

        # compute Lame constants
        this.λ, this.μ = lameConstants(E=this.E, ν=this.ν)
        #-----------------------------------------------------------------
        this.getModulus = function()
            return this.E
        end
        #------------------------------------------------------------------
        this.setModulus = function(x)
            this.E = x
            this.λ, this.μ = lameConstants(E=this.E, ν=this.ν)
            return this
        end
        #------------------------------------------------------------------
        return this
    end
end