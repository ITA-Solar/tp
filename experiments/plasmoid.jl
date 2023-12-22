#-------------------------------------------------------------------------------
# Created 13.01.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 plasmoid.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Import external libraries
#import Pkg; Pkg.add("Plots")
#ENV["PYTHON"]="/mn/alruba2/astro/local/mamba/envs/py310/python"
#ENV["PYTHON"]="python"
#Pkg.build("PyCall")
using LinearAlgebra

# Import internal modules from tp/src
using Constants:        k_B
using Meshes:           Mesh
using Patches:          Patch, run!
using Particles:        ParticleSoA, GCAParticleSoA, kineticenergy, specieTable
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
tf = 2.9 #n=100   # End time of simulation [s]         |
#tf = 2.3 #n=10    # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(Int64, numparticles) # Specifies the species of the particles 

#...............................................
# INITIAL CONDITIONS
ic = "d" # u: uniform velocity)
         # m: Maxwell-Boltzmann distributed velocity
         # d: Determined manually
# Placement span
pos0 = [0.10, 0.0, 0.0]
posf = [0.9, 0.4, 0.0]
# Temperature and mass (for MB distribution of velocity
T = 1/(2k_B)
mass = specieTable[species[1], 1]
# Velocity span (for uniform velcity)
vel0 = [0.5, 0.0, 0.0]
velf = [0.5, 0.0, 0.0]
# For manually determined positions and velocities
posd = [0.13  0.23  0.265 0.20 # At n = (100,100,2)
        0.03  0.05  0.12  0.04
        0.00  0.00  0.00  0.00]
#posd = [0.13  0.20  0.25  0.20  # At n = (10,10,2)
#        0.03  0.05  0.12  0.04
#        0.00  0.00  0.00  0.00]
veld = [0.50  0.50  0.50  1.00
        0.00  0.00  0.00  0.00
        0.00  0.00  0.00  0.00]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (0., 0., 0.)
# Upper bound of the three spatial axes
xif = (1., 1., 1.)
# Grid resolution of the axes
#n = (10, 10, 2)
n = (100, 100, 2)
 
#...............................................
# MAGNETIC FIELD PARAMETERS
bamp   = 10.0            # Amplification factor
bconst = [100., 0., 0]   # Constant term
bz     = 0.0             # z-component
μx = (xif[1] - xi0[1])/2 # 
μy = (xif[2] - xi0[2])/2 
σx = (xif[1] - xi0[1])/8
σy = (xif[2] - xi0[2])/8

#...............................................
# ELECTRIC FIELD PARAMETERS
# It will be static and homogeneous
Ex = 0.0
Ey = 0.0
Ez = 50.0

#...............................................
# SOLVER CONDITIONS
solver = Solvers.GCA # Physics solver
scheme = Schemes.rk4       # Scheme for integrating the diff eqs.
interp = Interpolations_tp.trilinearGCA # Interpolation scheme
pbc    = (true, true, true) # (x,y,z) Are mesh boundary conditions periodic?
#-------------------------------------------------------------------------------
seatpointsy = [0.1,0.17,0.23,0.28,0.32,0.34,0.35,0.36,0.37,0.375,0.38,0.39]
seatpointsy = [seatpointsy; [0.40,0.42,0.48, 0.55, 0.62, 0.71, 0.95]]
seatpoints = zeros(numdims, length(seatpointsy))
seatpoints[1,:] .= 0.5
seatpoints[2,:] .= seatpointsy


#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE MAGNETIC FIELD, AXES AND ELECTRIC FIELD
# Expectation and stdvalues. For creating vector potential with z-component
# normally distributed in x and y
μ  = (μx, μy)
σ  = (σx, σy)
# Creating the vector potential also gives the axes of the experiment
domainaxes, gridsizes, A = Utilities.normal3Donlyz(xi0, xif, n, μ, σ)
# Derive the magnetic field from the curl of the vector-potential
Bfield = Schemes.curl(A, gridsizes, Schemes.derivateCentral)
# Scale the field
@. Bfield = bamp*Bfield + bconst
@. Bfield[3,:,:,:] = bz
# Create the electric field
Efield = zeros(Float64, size(Bfield))
Efield[1, :, :, :] .= Ex
Efield[2, :, :, :] .= Ey
Efield[3, :, :, :] .= Ez

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numSteps = trunc(Int64, tf/dt)   # Number of timesteps in the simulation
println("Number of time steps = $numSteps.")

#-------------------------------------------------------------------------------
# MESH CREATION
# Create Mesh instance
xx, yy, zz = domainaxes
mesh = Mesh(Bfield, Efield, xx, yy, zz)

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
elseif ic == "d"
    pos, vel = posd, veld
    numparticles = size(posd)[2]
else
    println("Initial condition wrongly chosen.")
end
# Create ParticlesSoA-instance
if solver == Solvers.GCA
    R = pos
    vparal = zeros(numparticles)
    μ      = zeros(numparticles)
    for i = 1:numparticles
        m = specieTable[species[i], 1]
        (B⃗, E⃗, ∇B), _ = gridinterp(mesh,
                           interp,
                           R[:,i]
                           )
        B = norm(B⃗)
        b̂ = B⃗/B
        v = norm(vel[:,i])
        vparal[i] = vel[:,i] ⋅ b̂
        vperp = √(v^2 - vparal[i]^2)
        #vperp = Solvers.drift(b̂, E⃗, ∇B, B, μ, q)
        μ[i] = m*vperp^2/(2B)
    end
    particles = GCAParticleSoA(R, vparal, μ, species, numSteps)
else
    particles = ParticleSoA(pos, vel, species, numSteps)
end

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
