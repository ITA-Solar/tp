#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Schemes.jl
#
#-------------------------------------------------------------------------------
# Module containing the various schemes for solving differential equations
#-------------------------------------------------------------------------------




#---------------------------------------#
# Integration of differential equations #
#---------------------------------------#---------------------------------------
function euler(pos::Vector{T} where {T<:Real},
               vel::Vector{T} where {T<:Real}, 
               acc::Vector{T} where {T<:Real},
               dt ::Real
               )
    nextPos = @. pos + vel * dt
    nextVel = @. vel + acc * dt
    return nextPos, nextVel
end # function euler
function euler(vel::Real,
               acc::Real,
               dt ::Real,
               )
    nextVel = @. vel + acc * dt
    return nextVel
end # function euler


function eulerCromer(pos::Vector{T} where {T<:Real},
                     vel::Vector{T} where {T<:Real}, 
                     acc::Vector{T} where {T<:Real},
                     dt ::Real
                     )
    nextVel = @. vel + acc * dt
    nextPos = @. pos + nextVel * dt
    return nextPos, nextVel
end # function eulerCromer


function positionHalfStep(pos::Vector{T} where {T<:Real},
                          vel::Vector{T} where {T<:Real},
                          dt ::Real
                          )
    return pos + 0.5dt*vel
end # function positionHalfStep


function euler( # rk1
    yn     ::Vector{T} where {T<:Real}, # 'y' at time step 'n'
    h      ::Real,         # The time step
    f      ::Function,        # the time derivative of 'y'
    args...                   # Variable number of arguments to pass to 'f'.
    )
    k1 = f(yn, args...)
    return yn + h*k1
end # function euler
#|
function euler( # rk1
    yn     ::Vector{T} where {T<:Real}, # 'y' at time step 'n'
    h      ::Real,         # The time step
    f      ::Vector{T} where {T<:Real}, # The  time derivative of 'y'
    )
    return yn + h*f
end # function euler


function eulerCromer( # Semi-implicit
    yn     ::Vector{T} where {T<:Real}, # 'y' at time step 'n'
    h      ::Real,         # The time step
    f      ::Function,        # the time derivative of 'y'
    args...                   # Variable number of arguments to pass to 'f'.
    )
    numdims = trunc(Int64, length(yn)/2)
    # Advance velocity using Euler
    svNext = euler(yn, h, f, args...)
    velNext = svNext[numdims+1:2numdims]
    # Advance position using updated velocity
    posNext = euler(yn[1:numdims], h, velNext)
    return [posNext; velNext]
end # function eulerCromer

function rk4(
    yn     ::Vector{T} where {T<:Real}, # 'y' at time step 'n'
    h      ::Real,         # The time step
    f      ::Function,        # the time derivative of 'y'
    args...                   # Variable number of arguments to pass to 'f'.
    )
    k1 = f(yn         , args...)
    k2 = f(yn + h*k1/2, args...)
    k3 = f(yn + h*k2/2, args...)
    k4 = f(yn + h*k3  , args...)
    return yn + h/6*(k1 + 2k2 + 2k3 + k4)
end # function rk4

"""
    vay(pos, vel, specie, bField, eField, dt, scheme)

Implementation of the Vay pusher. For integrating the relativistic Lorentz
eqution. Adapted from J.-L. Vay (2008). This is only the second step of the
implementation, where the relativistic velocity is advanced a time step. The
first step involves advancing the position half a time step using
`positionHalfStep`, and in the last step the position is advanced its final 
half using the same function but with the new relativistic velocity.
"""
function vay(vel   ::Vector{T} where {T<:Real},
             bField::Vector{T} where {T<:Real},
             eField::Vector{T} where {T<:Real},
             mass  ::Real,
             charge::Real,
             dt    ::Real
             )
    # Some factors which use is repeated
    factor1 = 0.5charge*dt/mass

    v = norm(vel)                  # Speed of the particle
    γ = √(1/(1 - v^2*cSqrdInv)) # The relativistic gamma-factor
    u = γ * vel                   # The relativistic velocity

    # Compute half-step in velocity
    uHalf = u + factor1*(eField + vel × bField)

    # Compute auxiliary quantities used to get the full step in velcity
    # See Vay 2008 for details
    uPrime = uHalf + factor1*eField
    τ = factor1*bField
    uPrimeNorm = norm(uPrime)
    τNorm = norm(τ)
    uStar = uPrime ⋅ τ/c
    γPrime = √(1 + uPrimeNorm^2*cSqrdInv)
    σ = γPrime^2 - τNorm^2
    # Compute the gamma-factor for full step
    γNext = √(0.5(σ + √(σ^2 + 4(τNorm^2 + uStar^2))))
    t = τ/γNext
    s = 1/(1 + norm(t)^2)
    # Finally, compute the full step in velocity
    uNext = s*(uPrime + (uPrime ⋅ t)t + uPrime × t)

    # Get the non-relativistic velocity and position
    velNext = uNext/γNext
    
    return velNext
end # function vay


"""
    boris(pos, vel, specie, bField, eField, dt, scheme)

Implementation of the Boris pusher. For integrating the relativistic Lorentz
eqution. Adapted form Ripperda et al. 2018. This is only the second step of the
implementation, where the relativistic velocity is advanced a time step. The
first step involves advancing the position half a time step using
`positionHalfStep`, and in the last step the position is advanced its final 
half using the same function but with the new relativistic velocity.
"""
function boris(vel   ::Vector{T} where {T<:Real},
               bField::Vector{T} where {T<:Real},
               eField::Vector{T} where {T<:Real},
               mass  ::Real,
               charge::Real,
               dt    ::Real
               )
    # Some factors which use is repeated
    factor1 = 0.5charge*dt/mass

    v = norm(vel)                  # Speed of the particle
    γ = √(1/(1 - v^2*cSqrdInv)) # The relativistic gamma-factor
    u = γ * vel                   # The relativistic velocity
    

    # First half of the electric field acceleration
    uMinus = u + factor1*eField
    uMinusNorm = norm(uMinus)
    
    # Compute auxiliary quantities
    γMinus = √(1 + uMinusNorm*cSqrdInv)
    t = factor1*bField/γMinus
    tNorm = norm(t)
    s = 2t/(1 + tNorm^2)

    # Rotation step
    uPlus = uMinus + (uMinus + (uMinus × t)) × s

    # Second half of electric field acceleration
    uNext = uPlus + factor1*eField
    uNextNorm = norm(uNext)
    γNext = √(1 + uNextNorm^2*cSqrdInv)

    # Get the non-relativistic velocity and position
    velNext = uNext/γNext
    
    return velNext
end # function boris
#-------------------------------------------------------------------------------


#-----------------#
# Differentiation #
#-----------------#-------------------------------------------------------------

"""
    dervateUpwind(
        field::Array{T, 3} where {T<:Real},
        dx   ::Vector{T} where {T<:Real},
        axis ::Tuple{Integer, Integer, Integer}
    )
Differentiates a 3D `field` with respect to a specified `axis` using the upwind
scheme. The grid size may be variable, hence given as the vector `dx`. End point
of result will be ill-calculated and the derivative will be defined at half grid
point higher than the input field.
"""
function derivateUpwind(
    field::Array{T, 3} where {T<:Real},
    xx   ::Vector{T} where {T<:Real},
    yy   ::Vector{T} where {T<:Real},
    zz   ::Vector{T} where {T<:Real},
    ;
    wfp::DataType=typeof(field[1])
    )
    ni, nj, nk = size(field)

    if length(xx) > 1
        df = circshift(field, (-1,0,0)) - field
        # Pad the Δx array with a copy of the first value at the end. I.e. assume
        # periodic boundary conditions
        dx = circshift(xx, -1) - xx
        ddx = df ./ dx
    else
        ddx = zeros(wfp, ni,nj,nk)
    end

    if length(yy) > 1
        df = circshift(field, (0,-1,0)) - field
        dy = circshift(yy, -1) - yy
        ddy = df ./ dy'
    else
        ddy = zeros(wfp, ni, nj, nk)
    end

    if length(zz) > 1
        df = circshift(field, (0,0,-1)) - field
        dz = circshift(zz, -1) - zz
        # I don't know of a fast method for the 3rd dimension
        ddz = zeros(wfp, ni, nj, nk)
        for k = 1:nk
            ddz[:,:,k] = df[:,:,k] / dz[k]
        end
    end
    
    return ddx, ddy, ddz
end #function derivateUpwind


"""
    derivateCentral(field, dx)
First and last grid point are ill calculated.
"""
function derivateCentral(field::Vector{T} where {T<:Real},
                         dx
                         )
    return (circshift(field, -1) - circshift(field, 1))/2dx
end # function derivateCentral
#|
function derivateCentral(field      ::Array{T, 3} where {T<:Real},
                         gridSpacing::Real,
                         axis       ::Tuple{Int64, Int64, Int64}
                         )
    ax1 = -1 .* axis
    return (circshift(field, ax1) - circshift(field, axis))/2gridSpacing
end # function derivateCentral


"""
    derivate4thOrder(field, dx)
First and last two grid point are ill calculated.
"""
function derivate4thOrder(field      ::Array{T, 3} where {T<:Real},
                          gridSpacing::Real,
                          axis       ::Tuple{Int64, Int64, Int64}
                          )
    ax1 = -2 .* axis
    ax2 = -1 .* axis
    ax3 =  2 .* axis
    return (-circshift(field, ax1) + 8.0 .* circshift(field, ax2) - 
        8.0 .* circshift(field, axis) + circshift(field, ax3))/12.0gridSpacing
end # function derivate4thOrder


"""
    ∇(field, dx, dy, dz, scheme)
The gradient operator. Calucaletes the gradient of a 3- or 2-dimensional scalar
field and the Jacobian of a 3D vector field. Requires grid spacing on all three 
axis and the numerical scheme as arguments. The scheme is given as a function 
type, e.g. Schemes.derivateCentral.
"""
function ∇(field::Array{T, 4} where {T<:Real},
           dx   ::Real,
           dy   ::Real,
           dz   ::Real,
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    #
    ∂fx∂x = scheme(fx, dx, (1,0,0))
    ∂fx∂y = scheme(fx, dy, (0,1,0))
    ∂fx∂z = scheme(fx, dz, (0,0,1))
    #
    ∂fy∂x = scheme(fy, dx, (1,0,0))
    ∂fy∂y = scheme(fy, dy, (0,1,0))
    ∂fy∂z = scheme(fy, dz, (0,0,1))
    #
    ∂fz∂x = scheme(fz, dx, (1,0,0))
    ∂fz∂y = scheme(fz, dy, (0,1,0))
    ∂fz∂z = scheme(fz, dz, (0,0,1))
    #
    _, nx, ny, nz = size(field)
    jacobian = zeros(wfp, 3, 3, nx, ny, nz)
    #
    jacobian[1, 1, :,:,:] = ∂fx∂x
    jacobian[1, 2, :,:,:] = ∂fx∂y
    jacobian[1, 3, :,:,:] = ∂fx∂z
    #
    jacobian[2, 1, :,:,:] = ∂fy∂x
    jacobian[2, 2, :,:,:] = ∂fy∂y
    jacobian[2, 3, :,:,:] = ∂fy∂z
    #
    jacobian[3, 1, :,:,:] = ∂fz∂x
    jacobian[3, 2, :,:,:] = ∂fz∂y
    jacobian[3, 3, :,:,:] = ∂fz∂z
    #
    return jacobian
end # function ∇ 
#|
function ∇(field::Array{T, 3} where {T<:Real},
           dx   ::Real,
           dy   ::Real,
           dz   ::Real,
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    ∂f∂x = scheme(field, dx, (1,0,0))
    ∂f∂y = scheme(field, dy, (0,1,0))
    ∂f∂z = scheme(field, dz, (0,0,1))
    nx, ny, nz = size(field)
    gradient = zeros(wfp, 3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient
end # function ∇ 
#|
function ∇(field::Array{T, 2} where {T<:Real},
           dx   ::Real,
           dy   ::Real,
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    dfdx = scheme(field, dx, (1,0,0))
    dfdy = scheme(field, dy, (0,1,0))
    nx, ny = size(field)
    gradient = zeros(wfp, 3, nx, ny)
    gradient[1, :, :] = ∂f∂x
    gradient[2, :, :] = ∂f∂y
    return gradient
end # function ∇ 
    
"""
Newer versions of the gradient, allowing for non-uniform structured grid.
"""
function ∇(field::Array{T, 4} where {T<:Real},
           xx   ::Vector{T} where {T<:Real},
           yy   ::Vector{T} where {T<:Real},
           zz   ::Vector{T} where {T<:Real},
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    #
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    #
    ∂fx∂x, ∂fx∂y, ∂fx∂z = scheme(fx, xx, yy, zz)
    ∂fy∂x, ∂fy∂y, ∂fy∂z = scheme(fy, xx, yy, zz)
    ∂fz∂x, ∂fz∂y, ∂fz∂z = scheme(fz, xx, yy, zz)
    #
    _, nx, ny, nz = size(field)
    jacobian = zeros(wfp, 3, 3, nx, ny, nz)
    #
    jacobian[1, 1, :,:,:] = ∂fx∂x
    jacobian[1, 2, :,:,:] = ∂fx∂y
    jacobian[1, 3, :,:,:] = ∂fx∂z
    #
    jacobian[2, 1, :,:,:] = ∂fy∂x
    jacobian[2, 2, :,:,:] = ∂fy∂y
    jacobian[2, 3, :,:,:] = ∂fy∂z
    #
    jacobian[3, 1, :,:,:] = ∂fz∂x
    jacobian[3, 2, :,:,:] = ∂fz∂y
    jacobian[3, 3, :,:,:] = ∂fz∂z
    #
    return jacobian
end # function ∇
#|
function ∇(field::Array{T, 3} where {T<:Real},
           xx   ::Vector{T} where {T<:Real},
           yy   ::Vector{T} where {T<:Real},
           zz   ::Vector{T} where {T<:Real},
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
        
    ∂f∂x, ∂f∂y, ∂f∂z = scheme(field, xx, yy, zz)

    nx, ny, nz = size(field)
    gradient = zeros(wfp, 3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient

end # function ∇


"""
    curlfield, gridsizes, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires uniform grid spacing on all three axis. 
arguments. The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(field      ::Array{T, 4} where {T<:Real},
              gridsizes  ::Tuple{Real, Real, Real},
              derivscheme::Function
              ;
              wfp::DataType=typeof(field[1])
              )
    dx, dy, dz = gridsizes
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    derivx = (1,0,0)
    derivy = (0,1,0)
    derivz = (0,0,1)
    ∂fx∂y = derivscheme(fx, dy, derivy)
    ∂fx∂z = derivscheme(fx, dz, derivz)
    ∂fy∂x = derivscheme(fy, dx, derivx)
    ∂fy∂z = derivscheme(fy, dz, derivz)
    ∂fz∂x = derivscheme(fz, dx, derivx)
    ∂fz∂y = derivscheme(fz, dy, derivy)
    _, nx, ny, nz = size(field)
    result = zeros(wfp, 3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # functin curl

"""
    curl(field, dx, dy, dz, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires grid spacing on all three axis (may be non-uniform)
The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(
    field      ::Array{T, 4} where {T<:Real},
    dx         ::Vector{T} where {T<:Real},
    dy         ::Vector{T} where {T<:Real},
    dz         ::Vector{T} where {T<:Real},
    derivscheme::Function
    ;
    wfp::DataType=typeof(field[1])
    )
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    derivx = (1,0,0)
    derivy = (0,1,0)
    derivz = (0,0,1)
    ∂fx∂y = derivscheme(fx, dy, derivy)
    ∂fx∂z = derivscheme(fx, dz, derivz)
    ∂fy∂x = derivscheme(fy, dx, derivx)
    ∂fy∂z = derivscheme(fy, dz, derivz)
    ∂fz∂x = derivscheme(fz, dx, derivx)
    ∂fz∂y = derivscheme(fz, dy, derivy)
    _, nx, ny, nz = size(field)
    result = zeros(wfp, 3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # function curl
#|
function curl(
    fx         ::Array{T, 3} where {T<:Real},
    fy         ::Array{T, 3} where {T<:Real},
    fz         ::Array{T, 3} where {T<:Real},
    xx         ::Vector{T} where {T<:Real},
    yy         ::Vector{T} where {T<:Real},
    zz         ::Vector{T} where {T<:Real},
    derivscheme::Function
    ;
    wfp::DataType=typeof(field[1])
    )
    ∂fx∂x, ∂fx∂y, ∂fx∂z = derivscheme(fx, xx, yy, zz)
    ∂fy∂x, ∂fy∂y, ∂fy∂z = derivscheme(fy, xx, yy, zz)
    ∂fz∂x, ∂fz∂y, ∂fz∂z = derivscheme(fz, xx, yy, zz)
    nx, ny, nz = size(fx)
    result = zeros(wfp, 3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # function curl


"""
    LinearAlgebra.cross(f, g)
Method for computing the cross product between two 3D vector-fields.
"""
function LinearAlgebra.cross(
    f::Array{T, 4} where {T<:Real},
    g::Array{T, 4} where {T<:Real}
    ;
    wfp::DataType=typeof(field[1])
    )
    _, ni, nj, nk = size(f)
    crossproduct = zeros(wfp, 3, ni, nj, nk)
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                crossproduct[:,i,j,k] = f[:,i,j,k] × g[:,i,j,k]
            end
        end
    end
    return crossproduct
end # function cross
