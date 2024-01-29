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

#---------------#
#  Use packages #
#---------------#---------------------------------------------------------------
using LinearAlgebra
using Random

using DifferentialEquations
using ForwardDiff
using Interpolations

using Printf
using JLD

using Distributed

using BifrostTools
#-------------------------------------------------------------------------------

"""
Abstract types
"""
abstract type AbstractMesh end

include("utils.jl")
include("constants.jl")
include("mathematics.jl")
include("statistics.jl")
include("physics.jl")
include("equations_of_motion.jl")
include("solvers.jl")
include("callbacks.jl")
include("dataprocessing.jl")
include("numerical_methods/ode_schemes.jl")
include("numerical_methods/sde_schemes.jl")
include("numerical_methods/differentiations.jl")
include("numerical_methods/interpolations.jl")
include("io/save.jl")
include("io/load.jl")
include("io/bifrost_input.jl")


#---------#
# Exports #
#---------#---------------------------------------------------------------------
# physics.jl
export  get_guidingcentre,
        cosineof_pitchangle,
        kineticenergy
# equations_of_motion.jl
export  GCAPitchAngleFriction_lowmemory_2Dxz,
        GCAPitchAngleDiffusion_lowmemory_2Dxz,
        lorentzforce
# statistics.jl
export  maxwellianvelocitysample
# numerical_methods/sde_schemes
export  euler_maruyama,
        milstein_central
# solvers.jl
export  solve,
        solve_stoppingtime
# io/bifrost_input.jl
export  get_br_emfield_interpolator,
        get_br_emfield_numdensity_gastemp_interpolator,
        get_br_var_interpolator
# utils.jl
export dropdims


#-------------------------------------------------------------------------------
end
