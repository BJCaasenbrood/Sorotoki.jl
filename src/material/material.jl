abstract type Material end

"""
Returns Lame constants (λ,μ) from the (linear) material parameters E the Youngs' modulus and ν the Poisson ratio:
    λ := Eν/(1+ν)(1-2ν)
    μ := E/2(1+ν)
"""
function lameConstants(;E = 1.0, ν=0.499)
    @assert (-1 ≤ ν < 0.5) "Poisson ratio should be bounded by -1.0 ≤ ν₀ < 0.5"

    return (E*ν/((1.0+ν)*(1.0-2.0*ν)), E/(2.0*(1.0+ν)))
end