# Created 31.03.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 gradBdrift.jl
#
#-------------------------------------------------------------------------------
# This experiment tests the ∇B-drift term in the GCA against the full orbit
# solution in a static magnetic field where the gradient is perpendicular to the
# magnetic field direction. The parameters yield a Larmor radius that is much
# smaller than the characteristic length scale of the magnetic field
# gradient, such that the GCA is valid. The solutions is also compared to the
# approximate gradient drift presented in e.g. Chen (2016) or Aschwanden (2006).
#
# To be able to compare the GCA and full orbit, the simulation duration has to
# be a multiple of the gyration period.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 0.01         # Time step [s]                      |
# Use a final time equal to 10 gyrations for the full  |
# orbit, according to the parameters set (mass, charge,|
# vel0, pos0, B0, a)                                   |
tf = 10*1/1.7507  # End time of simulation [s]         |
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(Int64, numparticles) # Specifies the species of the particles 

#...............................................
# INITIAL CONDITIONS
pos0 = [0.5, 0.5, 0.5]
vel0 = [0.0, 0.1, 0.0]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (0., 0., 0.)
# Upper bound of the three spatial axes
xif = (1., 1., 1.)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 2)
 
#...............................................
# MAGNETIC FIELD PARAMETERS
a = 10.0 # gradient in magnetic field
B0 = 6.0 # Additional constant
function gradBfield(
    x::Float64,
    y::Float64,
    z::Float64
    )
    return [0.0, 0.0, a*y + B0]
end

#...............................................
# ELECTRIC FIELD PARAMETERS
# It will not exist
Ex = 0.0
Ey = 0.0

#...............................................
# SOLVER CONDITIONS
FOscheme = rk4  # Scheme for integrating the full-orbit eqs.gra
GCAscheme = euler # Scheme for integrating the GCA eqs.
FOinterp = trilinear # Interpolation scheme for full orbit
GCAinterp = trilinearGCA # Interpolation scheme for GCA
pbc    = (true, true, true) # (x,y,z) Are mesh boundary conditions periodic?


#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = createaxes(xi0, xif, ni)
Bfield = zeros(Float64, numdims, ni[1], ni[2], ni[3])
Efield = zeros(size(Bfield))
Efield[1,:,:,:] .= Ex
Efield[2,:,:,:] .= Ey
#Efield[1,:,1:30,:] .= 0.0
discretise!(Bfield, xx, yy, zz, gradBfield)

#-------------------------------------------------------------------------------
# MESH CREATION
# Create Mesh instance
mesh = Mesh(Bfield, Efield, xx, yy, zz)

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf/dt)   # Number of timesteps in the simulation
#println("Number of time steps = $numsteps.")

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
(B⃗, E⃗), _ = gridinterp(mesh,
                       trilinear,
                       pos0)
B = norm(B⃗)
b̂ = B⃗/B
v = norm(vel0)
vparal = [vel0 ⋅ b̂]
vperp = √(v^2 - vparal[1]^2)
mass = specieTable[species[1], 1]
μ = [mass*vperp^2/(2B)]

tpFO = ParticleSoA(pos0, vel0, species, numsteps)
tpGCA = GCAParticleSoA(pos0, vparal, μ, species, numsteps)

#-------------------------------------------------------------------------------
# CREATE PATCH
# Non-relativisitc Euler-Cromer
patchFO = Patch(mesh,
                tpFO,
                fullOrbit,
                FOscheme,
                FOinterp,
                dt,
                numsteps,
                numparticles,
                pbc
                )
patchGCA = Patch(mesh,
                 tpGCA,
                 GCA,
                 GCAscheme,
                 GCAinterp,
                 dt,
                 numsteps,
                 numparticles,
                 pbc
                 )

#...............................................................................
# RUN SIMULATION
run!(patchFO)
run!(patchGCA)
#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
charge = specieTable[species[1], 2]
# Initial gyrofrequency
B = norm(gradBfield(pos0[1],pos0[2],pos0[3]))
vperp = √(vel0[1]^2 + vel0[2]^2)
ff = charge*B/(mass*2π)
rL = mass*vperp/(charge*B)
L = a/B
vdrift = [-a*mass*vel0[2]^2/(2charge*B^2), 0.0, 0.0]
posf_anal = tf*vdrift + pos0
posf_GCA  = patchGCA.tp.R[:,1,end]
posf_FO   = patchFO.tp.pos[:,1,end]
@testset verbose = true "GCA: Euler" begin
    @test isapprox(posf_anal,posf_GCA, rtol=0.0001)
    @test isapprox(posf_FO,posf_GCA, rtol=0.001)
end # testset GCA: Euler
@testset verbose = true "Full orbit: RK4" begin
    @test isapprox(posf_anal,posf_FO, rtol=0.001)
end # testset GCA: Euler
    
