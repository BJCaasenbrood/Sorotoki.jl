# abstract type SO <: SMatrix{3,3, Float64, 9} end
# abstract type SO{N} <: SMatrix end?
# *(ᵀ)(x) = transpose(x)

function Id(N)
    if typeof(N) == Int64
        return Matrix(1I, N, N) 
    elseif typeof(N) == Float64
        return Matrix(1.0I, convert(Int64,N), convert(Int64,N))
    else
        return Matrix(1I, convert(Int64,N), convert(Int64,N)) 
    end
end

"""
The skewsymmetric operator ×(⋅): Rⁿ → so(n) that maps Rⁿ to the Lie algebra of SO(n) with n = {1,3}. Given A = skew(X), a property is that A + Aᵀ = 0. Also, any cross product between two vectors X,Y ∈ R³ can written as X × X = skew(X) * Y. 

Examples
≡≡≡≡≡≡≡≡≡≡

    X = [1,2,3]       % vector in R³
    Y = [3,2,1]       % vector in R³
    
    A = skew(X)       % skew-symmetric matrix, A = -Aᵀ ∈ so(3)
    b = skew(X) * Y   % cross product between x and y, b ∈ R³
"""
function skew(x::AbstractVector)
    N = length(x)
    if N == 1
        return @SMatrix[0 -x[1]; 0 x[1]]
    elseif N == 3
        return @SMatrix[0.0 -x[3] x[2]; x[3] 0.0 -x[1]; -x[2] x[1] 0.0]
    else
        msg = "Skew-symmetric matrix for N=$N is not feasible."
        throw(ArgumentError(msg))
    end
end

isskew(A::AbstractMatrix) = A' == -A

"""
The hat operator ∧(⋅): R⁶ → se(3) that maps R⁶ to the Lie algebra of SE(3) -- the rigid-bdy transformation group on R³. 
Usage:

    X = [0,0,0,0,1,0]       % vector in R⁶
    Y = hat(X)

    Y = 
"""
function hat(x::AbstractVector)
    if length(x) == 3
        return skew(x)
    elseif length(x) == 6
        return @SMatrix[0.0 -x[3] x[2] x[4];
                        x[3] 0.0 -x[1] x[5];
                        -x[2] x[1] .0 x[6];
                        0.0 0.0 0.0 0.0]
    end
end

"""
The adjoint action ad(⋅): se(3) → se(3) from one algebra onto another. 
Usage:
"""
function ad(x::AbstractVector)
    @assert length(x) == 6
    W = skew(x[1:3])
    U = skew(x[4:6])
    Z1 = hcat(W,abs.(0.0 * W))
    Z2 = hcat(U, W)
    return SMatrix{6,6}(vcat(Z1,Z2))
end

"""
The adjoint action ad(⋅): se(3) → se(3) from one algebra onto another. 
Usage:
"""
function Ad(X::AbstractMatrix)
    @assert size(X) == (4,4)
    R = X[1:3,1:3]
    S = skew(X[1:3,4])
    Z1 = hcat(R, 0. * R)
    Z2 = hcat(S * R, R)
    return SMatrix{6,6}(vcat(Z1,Z2))
end

"""
The Adjoint action Ad⁻¹G(⋅): of a group G onto se(3)
Usage:
"""

function Ad⁻¹(X::AbstractMatrix)
    @assert size(X) == (4,4)
    Rᵀ = (X[1:3,1:3])'
    Sᵀ = -skew(X[1:3,4])
    Z1 = hcat(Rᵀ, 0. * Rᵀ)
    Z2 = hcat(Rᵀ * Sᵀ, Rᵀ)
    return SMatrix{6,6}(vcat(Z1,Z2))
end