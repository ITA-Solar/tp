#-------------------------------------------------------------------------------
# Created 12.02.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 doubleplasmoid.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Import external libraries
using LinearAlgebra

# Import internal modules from tp/src
using Constants:        k_B
using Meshes:           Mesh
using Patches:          Patch, run!
using Particles:        ParticleSoA, kineticenergy, specieTable
using Solvers          
using Schemes
using Interpolations_tp
using Utilities

"""
"""

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 10 # Number of particles to simulate    |
dt = 0.001        # Time step [s]                      |
tf = 3.0          # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(Int64, numparticles) # Specifies the species of the particles 

#...............................................
# INITIAL CONDITIONS
ic = "m" # u: uniform velocity)
         # m: Maxwell-Boltzmann distributed velocity
# Placement span
pos0 = [0.25, 0.0, 0.0]
posf = [0.75, 1.0, 0.0]
# Temperature and mass (for MB distribution of velocity
T = 1/(2k_B)
mass = specieTable[species[1], 1]
# Velocity span (for uniform velcity)
vel0 = [0.5, 0.0, 0.0]
velf = [0.5, 0.0, 0.0]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (0., 0., 0.)
# Upper bound of the three spatial axes
xif = (1., 1., 1.)
# Grid resolution of the axes
n = (10, 10, 2)
 
#...............................................
# MAGNETIC FIELD PARAMETERS
bamp     = 5.0            # Amplification factor
bconst   = [1., 0., 0]   # Constant term
bz       = 0.0             # z-component
λ, ϵ, B0 = (1.0, 1.0, 0.1) # Fadeev-equilibrium vector-potential parameters


#...............................................
# ELECTRIC FIELD PARAMETERS
# It will be static and homogeneous
Ex = 0.0
Ey = 0.0
Ez = 0.0

#...............................................
# SOLVER CONDITIONS
solver = Solvers.fullOrbit # Physics solver
scheme = Schemes.rk4       # Scheme for integrating the diff eqs.
interp = Interpolations_tp.trilinear # Interpolation scheme
pbc    = (true, true, true) # (x,y,z) Are mesh boundary conditions periodic?
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE MAGNETIC FIELD, AXES AND ELECTRIC FIELD
# Creating the vector potential also gives the axes of the experiment
axes, gridsizes, A = Utilities.fadeevEquilibrium(xi0, xif, n, λ, ϵ, B0)
# Derive the magnetic field from the curl of the vector-potential
B = Schemes.curl(A, gridsizes, Schemes.derivateCentral)
# Scale the field
@. B = bamp*B + bconst
@. B[3,:,:,:] = bz
# Create the electric field
E = zeros(Float64, size(B))
E[1, :, :, :] .= Ex
E[2, :, :, :] .= Ey
E[3, :, :, :] .= Ez

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numSteps = trunc(Int64, tf/dt)   # Number of timesteps in the simulation
println("Number of time steps = $numSteps.")

#-------------------------------------------------------------------------------
# MESH CREATION
# Create Mesh instance
xx, yy, zz = axes
mesh = Mesh(B, E, xx, yy, zz)

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
if ic == "u"
    pos, vel = Utilities.initparticlesuniform(
        numparticles,
        pos0,
        posf,
        vel0,
        velf
    )
elseif ic == "m"
    pos, vel = Utilities.initparticlesmaxwellianx(
        numparticles,
        pos0,
        posf,
        T,
        mass
    )
else
    println("Initial condition wrongly chosen.")
end
# Create ParticlesSoA-instance
particles = ParticleSoA(pos, vel, species, numSteps)

#-------------------------------------------------------------------------------
# CREATE PATCH
# Non-relativisitc Euler-Cromer
patch = Patch(mesh,
              particles,
              solver,
              scheme,
              interp,
              dt,
              numSteps,
              numparticles,
              pbc
              )

#...............................................................................
# RUN SIMULATION
@time run!(patch)
#...............................................................................
