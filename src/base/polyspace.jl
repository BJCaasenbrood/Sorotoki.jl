function chebyshev(k::Int64)
    Nts = 0.0:convert(Float64,k-1)
    θ = (x) -> cos.(Nts .* acos.( (2*clamp.(x,0,1) - 1 )));
    #return x -> SMatrix{length(x),k}(reduce(vcat,transpose.(θ.(x))))
    #return x -> SMatrix{length(x),k}(permutedims(hcat(θ.(x)...)))
    return x -> reduce(hcat,θ.(x))'
end

function gsogpoly(Y, x)
    d, n = size(Y)

    if x === nothing
        x = collect(range(0.0, 1.0, d))
    end
    m = min(d, n)
    R = Matrix{Float64}(I, m, n)
    Q = zeros(d, m)

    for ii in 1:m
        R = zeros(d)
        for jj in 1:ii-1
            D = sum(Q[:, jj] .* Q[:, jj])
            R .+= (sum(Y[:, ii] .* Q[:, jj]) / D) * Q[:, jj]
        end
        Q[:, ii] = Y[:, ii] - R
        Q[:, ii] /= sqrt(sum(Q[:, ii] .* Q[:, ii]))
    end
    
    return Q, R
end

function simpsons_rule(Y::Vector{T}, X::Vector{T}) where T
    @assert length(Y) == length(X) "Y and X must have the same length"
    Z , tmp = 0., 0.

    for i in 1:2:length(X)-2
        Z = (X[i+2] - X[i]) * (Y[i] + 4Y[i+1] + Y[i+2]) / 6
        Z += tmp
    end
    
    return Z
end
