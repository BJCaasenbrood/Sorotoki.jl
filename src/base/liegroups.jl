function Id(N)
    @assert typeof(N) == Int64
    return SMatrix{3, 3, Float64, 9}([1 0 0; 0 1 0; 0 0 1])
end

"""
The skewsymmetric operator ×(⋅): R³ → so(3) that maps R³ to the Lie algebra of SO(3) -- the rotation group on R³. Given A = skew(X), a property is that A + Aᵀ = 0. Also, any cross product between two vectors X,Y ∈ R³ can written as X × X = skew(X) * Y. 

Usage:

    X = [1,2,3]       % vector in R³
    Y = [3,2,1]       % vector in R³
    
    A = skew(X)       % skew-symmetric matrix, A = -Aᵀ ∈ so(3)
    b = skew(X) * Y   % cross product between x and y, b ∈ R³
"""
function skew(x)
    @assert length(x) == 3
    return @SMatrix[0.0 -x[3] x[2]; x[3] 0.0 -x[1]; -x[2] x[1] 0.0]
end

(×)(x,y) = skew(x) * y
(×)(x) = skew(x)

"""
The hat operator ∧(⋅): R⁶ → se(3) that maps R⁶ to the Lie algebra of SE(3) -- the rigid-bdy transformation group on R³. 
Usage:

    X = [0,0,0,0,1,0]       % vector in R⁶
    Y = hat(X)

    Y = 
"""
function hat(x)
    if length(x) == 3
        return skew(x)
    elseif length(x) == 6
        return @SMatrix[0.0 -x[3] x[2] x[4];
                        x[3] 0.0 -x[1] x[5];
                        -x[2] x[1] .0 x[6];
                        0.0 0.0 0.0 0.0]
    end
end