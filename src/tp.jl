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
using DifferentialEquations
using Interpolations
using JLD
using Distributed

using BifrostTools

"""
Abstract types: Ideas for next version
"""
#abstract type AbstractExperiment end
abstract type AbstractPatch end
abstract type AbstractMesh end
abstract type AbstractParticle end
abstract type TraceParticle <: AbstractParticle end
abstract type AbstractProblemParameters end
abstract type AbstractInitialConditions end
#abstract type AbstractSimulation end
#abstract type AbstractSolver end

include("utils.jl")
include("constants.jl")
# ---
# Ideas for next version
#include("problems.jl")
#include("simulations.jl")
# --- 
include("callbacks.jl")
include("solvers.jl")
include("meshes.jl")
include("tp_interpolations.jl") 
include("particles.jl")
include("equations_of_motion.jl")
include("patches.jl")
include("initial_conditions.jl")
include("interpolations.jl")
include("problem_parameters.jl")
include("solve.jl")
include("experiments.jl")
	
# Exports
export FOStaticIC, LorentzForce, ODEParticle, 
        FOParams, GCAParams, GCAStaticIC
export Patch, run!, update, DEPatch
export compute_gradients, derivateUpwind, EMfield_itps, calc_GCA_IC_and_mu
export DE_init!

export VectorIC, FloatIC, AbstractInitialConditions, AbstractVariableIC
export PhaseSpace1DIC
export SingleCoefficientParams
export DuffingVanDerPolOscillatorDrift, DuffingVanDerPolOscillatorDiffusion
export OrnsteinUhlenbeckDrift, OrnsteinUhlenbeckDiffusion
export ConstantAdvection, ConstantDiffusion
export ProportionalAdvection, ProportionalDiffusion
export SDEParticle, euler_maruyama, milstein_central
export solve, solve_stoppingtime
export DiscreteTPCallback
export GCAPitchAngleScatteringParams 
export GCAPitchAngleScatteringIC
export GCAPitchAngleFriction, GCAPitchAngleDiffusion
export create_puremesh, get_fields_from_br, tp_get_initial_pos_and_vel!,
        EMfield_at_pos, calc_GCA_IC_and_mu, pitchangle 
export save_fast_PAS, tp_load_fast
export gca


# meshes.jl
export Mesh
export PureMesh
export amplifyBfield! # Amplifies the magnetic field of the mesh by a factor
export amplifyEfield! # Amplifies the elctric field of the mesh by a factor
# particles.jl
export get_problem
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
