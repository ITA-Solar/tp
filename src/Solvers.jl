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


function relFullOrbitExplLeapFrog(vel         ::Vector{wpFloat},
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
    posHalf = positionHalfStep(pos, vel, dt)
    # Interpolate fields to this location
    bField, efield = Interpolations.grid(mesh,
                                         interpolator,
                                         posHalf)

    # 
    # Step 2: Evaluate full time step in velocity, which is shceme-dependent.
    #
    velNext = scheme(vel, bField, eField, mass, charge, dt)
    
    #
    # Step 3: Evaluate second half of time step in position
    # 
    posNext =positionHalfStep(posHalf, velNext, dt)

    return posNext, velNext
end # function relFullOrbitExpLeapFrog

end # module Solvers
