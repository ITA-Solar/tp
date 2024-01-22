#-------------------------------------------------------------------------------
# Created 19.01.24. Code from back in 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#        ode_schemes.jl
#
#-------------------------------------------------------------------------------
# Contains various schemes for solving ordinary differential 
# equations
#-------------------------------------------------------------------------------


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
