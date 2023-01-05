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


function euler(pos::Vector{wpFloat},
               vel::Vector{wpFloat}, 
               acc::Vector{wpFloat},
               dt)
    nextPos = @. pos + vel * dt
    nextVel = @. vel + acc * dt
    return nextPos, nextVel
end # function euler


function eulerCromer(pos::Vector{wpFloat},
                     vel::Vector{wpFloat}, 
                     acc::Vector{wpFloat},
                     dt)
    nextVel = @. vel + acc * dt
    nextPos = @. pos + nextVel * dt
    return nextPos, nextVel
end # function eulerCromer


function positionHalfStep(pos::Vector{wpFloat},
                          vel::Vector{wpFloat},
                          dt
                          )
    return pos + 0.5dt*vel
end # function positionHalfStep


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

    # Compute various quantities used to get the full step in velcity
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
end # function vayAdvanceVel


end # module schemes
