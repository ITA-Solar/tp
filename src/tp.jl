#-------------------------------------------------------------------------------
# Created 02.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 tp.jl
#
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

module tp

#-------------------------#   
#  use external libraries # 
#-------------------------#-----------------------------------------------------
using LinearAlgebra
using Random:           MersenneTwister
using Printf

using BifrostTools

"""
Abstract types: Ideas for next version
"""
#abstract type AbstractExperiment end
#abstract type AbstractPatch end
#abstract type AbstractParticle end
#abstract type AbstractProblem end
#abstract type AbstractProblemParameters end
abstract type AbstractMesh end
#abstract type AbstractSimulation end
#abstract type AbstractSolver end
#abstract type AbstractInitialConditions end

include("utils.jl")
include("constants.jl")
# ---
# Ideas for next version
#include("initial_conditions.jl")
#include("problems.jl")
#include("problem_parameters.jl")
#include("simulations.jl")
# --- 
include("solvers.jl")
include("meshes.jl")
include("tp_interpolations.jl")
include("particles.jl")
include("equations_of_motion.jl")
include("patches.jl")
include("experiments.jl")
	
# Exports
export Patch, run!
# meshes.jl
export Mesh
export PureMesh
export amplifyBfield! # Amplifies the magnetic field of the mesh by a factor
export amplifyEfield! # Amplifies the elctric field of the mesh by a factor
# particles.jl
export TraceParticle
export ParticleSoA # Particles represented as struct of arrays
export GCAParticleSoA
export specieTable # Maping specie to mass and charge
export getpos
export getvel
export reset!      # Resets particle positions to zero (except initial position)
export revive!     # Resets particle alive status to true
export setinitpos! # Sets the initial position of particles
export setinitvel! # Sets the initial velocity of particles
export kineticenergy # Computes the non-rel. kinetic energy at all time steps
export computeμ
# utils.jl
export randn
export rand
export initparticlesuniform
export initparticlesmaxwellian
export norm4
export createaxes, discretise!, mirroringfield, dipolefield
# solvers.jl
export cross
export euler, eulerCromer, rk4
export curl, normal3Donlyz, derivateCentral
# equations_of_motion.jl
export fullOrbit, relFullOrbitExplLeapFrog, vay, boris, GCA
# tp_interpolations.jl
export gridinterp
export locateCell
export trilinear
export trilinearGCA
export bilinear_xz
export bilinear_xzGCA
# experiments.jl
export Parameters
export Experiment
# Init-functions
export tp_init!
export tp_reinit_particles!
# Save/load functions
export tp_save
export tp_saveparams
export tp_savetp
export tp_savemesh
export tp_savebg
export tp_load
export tp_loadtp
export tp_loadtp!
export tp_loadmesh
export tp_loadbg
# Run functions
export tp_run!
# Set-functions
export tp_set_dt!
export tp_set_nsteps!
export tp_reset!


end
