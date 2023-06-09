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


module Schemes

# Standard libraries
using LinearAlgebra
# Internal libraries
using WorkingPrecision: wpFloat, wpInt
using Constants:        c, cSqrdInv

export cross

#---------------------------------------#
# Integration of differential equations #
#---------------------------------------#---------------------------------------
function euler(pos::Vector{wpFloat},
               vel::Vector{wpFloat}, 
               acc::Vector{wpFloat},
               dt ::wpFloat
               )
    nextPos = @. pos + vel * dt
    nextVel = @. vel + acc * dt
    return nextPos, nextVel
end # function euler
function euler(vel::wpFloat,
               acc::wpFloat,
               dt ::wpFloat,
               )
    nextVel = @. vel + acc * dt
    return nextVel
end # function euler


function eulerCromer(pos::Vector{wpFloat},
                     vel::Vector{wpFloat}, 
                     acc::Vector{wpFloat},
                     dt ::wpFloat
                     )
    nextVel = @. vel + acc * dt
    nextPos = @. pos + nextVel * dt
    return nextPos, nextVel
end # function eulerCromer


function positionHalfStep(pos::Vector{wpFloat},
                          vel::Vector{wpFloat},
                          dt ::wpFloat
                          )
    return pos + 0.5dt*vel
end # function positionHalfStep


function euler( # rk1
    yn     ::Vector{wpFloat}, # 'y' at time step 'n'
    h      ::wpFloat,         # The time step
    f      ::Function,        # the time derivative of 'y'
    args...                   # Variable number of arguments to pass to 'f'.
    )
    k1 = f(yn, args...)
    return yn + h*k1
end # function euler
#|
function euler( # rk1
    yn     ::Vector{wpFloat}, # 'y' at time step 'n'
    h      ::wpFloat,         # The time step
    f      ::Vector{wpFloat}, # The  time derivative of 'y'
    )
    return yn + h*f
end # function euler


function eulerCromer( # Semi-implicit
    yn     ::Vector{wpFloat}, # 'y' at time step 'n'
    h      ::wpFloat,         # The time step
    f      ::Function,        # the time derivative of 'y'
    args...                   # Variable number of arguments to pass to 'f'.
    )
    numdims = trunc(wpInt, length(yn)/2)
    # Advance velocity using Euler
    svNext = euler(yn, h, f, args...)
    velNext = svNext[numdims+1:2numdims]
    # Advance position using updated velocity
    posNext = euler(yn[1:numdims], h, velNext)
    return [posNext; velNext]
end # function eulerCromer

function rk4(
    yn     ::Vector{wpFloat}, # 'y' at time step 'n'
    h      ::wpFloat,         # The time step
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
function vay(vel   ::Vector{wpFloat},
             bField::Vector{wpFloat},
             eField::Vector{wpFloat},
             mass  ::wpFloat,
             charge::wpFloat,
             dt    ::wpFloat
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
function boris(vel   ::Vector{wpFloat},
               bField::Vector{wpFloat},
               eField::Vector{wpFloat},
               mass  ::wpFloat,
               charge::wpFloat,
               dt    ::wpFloat
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
        field::Array{wpFloat, 3},
        dx   ::Vector{wpFloat},
        axis ::Tuple{wpInt, wpInt, wpInt}
    )
Differentiates a 3D `field` with respect to a specified `axis` using the upwind
scheme. The grid size may be variable, hence given as the vector `dx`. End point
of result will be ill-calculated and the derivative will be defined at half grid
point higher than the input field.
"""
function derivateUpwind(
    field::Array{wpFloat, 3},
    xx   ::Vector{wpFloat},
    yy   ::Vector{wpFloat},
    zz   ::Vector{wpFloat},
    )
    ni, nj, nk = size(field)

    if length(xx) > 1
        df = circshift(field, (-1,0,0)) - field
        # Pad the Δx array with a copy of the first value at the end. I.e. assume
        # periodic boundary conditions
        dx = circshift(xx, -1) - xx
        ddx = df ./ dx
    else
        ddx = zeros(wpFloat, ni,nj,nk)
    end

    if length(yy) > 1
        df = circshift(field, (0,-1,0)) - field
        dy = circshift(yy, -1) - yy
        ddy = df ./ dy'
    else
        ddy = zeros(wpFloat, ni, nj, nk)
    end

    if length(zz) > 1
        df = circshift(field, (0,0,-1)) - field
        dz = circshift(zz, -1) - zz
        # I don't know of a fast method for the 3rd dimension
        ddz = zeros(wpFloat, ni, nj, nk)
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
function derivateCentral(field::Vector{wpFloat},
                         dx
                         )
    return (circshift(field, -1) - circshift(field, 1))/2dx
end # function derivateCentral
#|
function derivateCentral(field      ::Array{wpFloat, 3},
                         gridSpacing::wpFloat,
                         axis       ::Tuple{wpInt, wpInt, wpInt}
                         )
    ax1 = -1 .* axis
    return (circshift(field, ax1) - circshift(field, axis))/2gridSpacing
end # function derivateCentral


"""
    derivate4thOrder(field, dx)
First and last two grid point are ill calculated.
"""
function derivate4thOrder(field      ::Array{wpFloat, 3},
                          gridSpacing::wpFloat,
                          axis       ::Tuple{wpInt, wpInt, wpInt}
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
function ∇(field::Array{wpFloat, 4},
           dx   ::wpFloat,
           dy   ::wpFloat,
           dz   ::wpFloat,
           scheme::Function
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
    jacobian = zeros(3, 3, nx, ny, nz)
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
function ∇(field::Array{wpFloat, 3},
           dx   ::wpFloat,
           dy   ::wpFloat,
           dz   ::wpFloat,
           scheme::Function
           )
    ∂f∂x = scheme(field, dx, (1,0,0))
    ∂f∂y = scheme(field, dy, (0,1,0))
    ∂f∂z = scheme(field, dz, (0,0,1))
    nx, ny, nz = size(field)
    gradient = zeros(3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient
end # function ∇ 
#|
function ∇(field::Array{wpFloat, 2},
           dx   ::wpFloat,
           dy   ::wpFloat,
           scheme::Function
           )
    dfdx = scheme(field, dx, (1,0,0))
    dfdy = scheme(field, dy, (0,1,0))
    nx, ny = size(field)
    gradient = zeros(3, nx, ny)
    gradient[1, :, :] = ∂f∂x
    gradient[2, :, :] = ∂f∂y
    return gradient
end # function ∇ 
    
"""
Newer versions of the gradient, allowing for non-uniform structured grid.
"""
function ∇(field::Array{wpFloat, 4},
           xx   ::Vector{wpFloat},
           yy   ::Vector{wpFloat},
           zz   ::Vector{wpFloat},
           scheme::Function
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
    jacobian = zeros(wpFloat, 3, 3, nx, ny, nz)
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
function ∇(field::Array{wpFloat, 3},
           xx   ::Vector{wpFloat},
           yy   ::Vector{wpFloat},
           zz   ::Vector{wpFloat},
           scheme::Function
           )
        
    ∂f∂x, ∂f∂y, ∂f∂z = scheme(field, xx, yy, zz)

    nx, ny, nz = size(field)
    gradient = zeros(wpFloat, 3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient

end # function ∇


"""
    ∇×(field, dx, dy, dz, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires grid spacing on all three axis and the numerical scheme as
arguments. The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(field      ::Array{wpFloat, 4},
              gridsizes  ::Tuple{wpFloat, wpFloat, wpFloat},
              derivscheme::Function
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
    result = zeros(3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # functin curl


end # module schemes
