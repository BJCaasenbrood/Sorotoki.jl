#

# function leapfrogMDE(f::Function, Z::Matrix{T}, t::T, dt::T) where {T <: AbstractFloat}
#     # Allocate intermediate arrays
#     Z_half = similar(Z)
#     Z_next = similar(Z)
    
#     Z_half .= Z + (dt / 2.0) * f(Z, t)
#     Z_next .= Z + dt * f(Z_half, t + dt / 2.0)
    
#     # Update the current state
#     Z .= Z_next
    
#     return Z
# end

function leapfrog(f::Function, Z::Vector{T}, t::T, dt::T) where {T <: AbstractFloat}
    # Allocate intermediate arrays
    Z_half = similar(Z)
    Z_next = similar(Z)
    
    Z_half .= Z + (dt / 2.0) * f(Z, t)
    Z_next .= Z + dt * f(Z_half, t + dt / 2.0)
    
    # Update the current state
    Z .= Z_next
    
    return Z
end

# function leapfrogAD(f::Function, Z::Vector{T}, t::T, dt::T) where {T <: AbstractFloat}
#     # Allocate intermediate arrays
#     Z11 = similar(Z)
#     Z12 = similar(Z)
    
#     # Compute the intermediate state
#     Z11 .= Z + (dt / 2.0) * f(Z, t)
    
#     # Define the residual function
#     function residual(Z_next)
#         return Z_next - (Z + dt * f(Z11, t + dt / 2.0))
#     end
    
#     # Solve NL eq using implicit scheme: Gauss-Newton :gn, 
#     # Levenberg-Marquardt :lm and the Trust Region :tr
#     Z12 = nlsolve(residual!, Z12, autodiff = :forward, 
#         ftol=1e-6, xtol=1e-6, maxiter=100, method=:gn)
    
#     Z .= Z12
    
#     return Z
# end
