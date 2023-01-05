#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Solvers.jl
#
#-------------------------------------------------------------------------------
# Module containing the various solvers.
#-------------------------------------------------------------------------------
module Solvers

# Standard library
using LinearAlgebra: ×, ⋅, norm
# Internal modules
using WorkingPrecision: wpInt, wpFloat
using Meshes
using Interpolations
using Particles: specieTable
using Constants: c, cSqrdInv

"""
    fullOrbit(pos, vel, specie, bField, eField, dt, scheme)

Solves the Lorentz equation of motion using an arbitrary numerical scheme
(defined by the argument `scheme`).
"""
function fullOrbit(pos         ::Vector{wpFloat},
                   vel         ::Vector{wpFloat},
                   specie      ::wpInt,
                   mesh        ::Mesh,
                   dt          ::wpFloat,
                   interpolator::Function
                   scheme      ::Function,
                   )
    # Extract particle mass and charge
    mass   = specieTable[specie, 1]
    charge = specieTable[specie, 2]
    bField, eField = Interpolations.grid(mesh,
                                         interpolator,
                                         pos)
    acc = charge/mass * (eField + vel × bField)
    newPos, newVel = scheme(pos, vel, acc, dt)
    return newPos, newVel
end # funcion fullOrbit


"""
    vay(pos, vel, specie, bField, eField, dt, scheme)

Implementation of the Vay pusher. For integrating the relativistic Lorentz
eqution. Adapted from J.-L. Vay (2008). This is only the second step of the
"""
function vay(pos         ::Vector{wpFloat},
             vel         ::Vector{wpFloat},
             specie      ::wpInt,
             mesh        ::Mesh,
             dt          ::wpFloat,
             interpolator::Function
             scheme      ::Function,
             )
    # Extract particle mass and charge
    mass   = specieTable[specie, 1]
    charge = specieTable[specie, 2]

    #
    # Step 1: Evaluate half-step in time for position
    #
    posHalf = pos + 0.5dt*vel
    # Interpolate fields to this location
    bField, efield = Interpolations.grid(mesh,
                                         interpolator,
                                         posHalf)

    # 
    # Step 2: Evaluate full time step in velocity, which is shceme-dependent.
    #
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

    #
    # Step 3: Evaluate second half of time step in position
    # 
    posNext = posHalf + 0.5dt*velNext

    return posNext, velNext
end

end # module Solvers
