#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Patches.jl
#
#-------------------------------------------------------------------------------
# Module containing the Patch struct and methods
#-------------------------------------------------------------------------------

module Patches

using WorkingPrecision: wpFloat, wpInt
using Meshes
using Particles
using Solvers
using Schemes
using Interpolations

#-------------#   
# Export      # 
#-------------#-----------------------------------------------------------------
export Patch
export run

#-------------#   
# Main struct # 
#-------------#-----------------------------------------------------------------
mutable struct Patch
    mesh        ::Mesh
    tp          ::ParticleSoA # The trace particles
    solver      ::Function
    scheme      ::Function
    interpolator::Function
    dt          ::wpFloat
    numSteps    ::wpInt
    numParticles::wpInt
end

#---------#
# Methods #
#---------#---------------------------------------------------------------------
function run!(patch::Patch)
    for i = 1:patch.numSteps # Over timesteps
        for j = 1:patch.numParticles # Over particles
            pos = patch.tp.pos[:, j, i]
            vel = patch.tp.vel[:, j, i]
            pos, vel = patch.solver(pos,
                                    vel,
                                    patch.tp.species[j],
                                    patch.mesh,
                                    patch.dt,
                                    patch.interpolator
                                    patch.scheme,
                                    )
            patch.tp.pos[:, j, i+1] = pos
            patch.tp.vel[:, j, i+1] = vel
       end # loop over particles (j)
    end # loop over timesteps (i)
end # function: run
#-------------------------------------------------------------------------------

end # module
