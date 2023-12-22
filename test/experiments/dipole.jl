# Created 14.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dipole.jl
#
#-------------------------------------------------------------------------------
# Runs a GCA and full orbit particle in a magnetic dipole. The setup is the same
# as in Ripperda et al., 2018. Initial positions is set such that the guiding
# centre of the particle is at x = [1.0, 0.0, 0.0]. See 'dipoleloop.jl' for the
# testsets that compares the result from the solvers with each other and the
# 'T_dipole' approximation. 
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 5.e-4        # Time step [s]                      |
tf = 100.0 #n=100   # End time of simulation [s]         |
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(Int64, numparticles) # Specifies the species of the particles 
mass = specieTable[species[1], 1]
charge = specieTable[species[1], 2]

#...............................................
# MAGNETIC FIELD PARAMETERS
#qMm = 20.0
#println("qMm: $qMm")
M = qMm*mass/charge

#...............................................
# INITIAL CONDITIONS
vel0 = [0.0, 1.0, 0.5]
# Initial position is hardcoded such that vperp is 1 and the guiding centre is
# at x⃗ = [1, 0, 0]. It is also assumed that the GCA is valid such that B is equal
# to B(x⃗) = 2M, directed in -ẑ. We may then compare with the Walt (1994)
# approximation of the azimuthal drift-period.
rL = mass/(charge*M) 
pos0 = [1.0 - rL, 0.0, 0.0]
R0 = [1.0, 0.0, 0.0]
pos0 = R0 - mass/(charge*M)*(vel0 × [0,0,-1.0])

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
a = 1.5
xi0 = (-a, -a, -a)
# Upper bound of the three spatial axes
xif = (a, a, a)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (256, 256, 256)

#...............................................
# ELECTRIC FIELD PARAMETERS
# It will not exist

#...............................................
# SOLVER CONDITIONS
FOscheme = rk4  # Scheme for integrating the full-orbit eqs.
GCAscheme = rk4 # Scheme for integrating the GCA eqs.
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
discretise!(Bfield, xx, yy, zz, dipolefield, M)

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
                       R0)
B = norm(B⃗)
b̂ = B⃗/B
v = norm(vel0)
vparal = [vel0 ⋅ b̂]
vperp = √(v^2 - vparal[1]^2)
mass = specieTable[species[1], 1]
μ = [mass*vperp^2/(2B)]
#μ = [0.0120]

tpFO = ParticleSoA(pos0, vel0, species, numsteps)
tpGCA = GCAParticleSoA(R0, vparal, μ, species, numsteps)

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
                 gca,
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
times = collect(range(0.0, step=dt, length=numsteps+1))
global ϕ_FO = zeros(size(times))
global ϕ_GCA = zeros(size(times))
@. ϕ_FO = atan(patchFO.tp.pos[2,1,:] ./ patchFO.tp.pos[1,1,:])
@. ϕ_GCA = atan(patchGCA.tp.R[2,1,:] ./ patchGCA.tp.R[1,1,:])

v = norm(vel0)
vparall = abs((B⃗ ⋅ vel0)/B)
vperp = √(v^2 - vparall^2)
α = atan(vperp/vparall)
R0 = 1.0

function T_dipole(q, m, M, R0, v, α)
    return @. 2π*q*M/(m*v^2*R0)*(1 - 1/3*sin(α)^0.62)
end # function T_dipole

T = T_dipole(charge, mass, M, R0, v, α)
global ϕ = 2π*times/T


#x = R0*cos(ϕ)
#y = R0*sin(ϕ)
#z = 0.0
#@testset verbose = true "GCA: Euler" begin
#    @test isapprox(rmse_GCA, 0.0, atol=0.001)
#end # testset GCA: Euler
#@testset verbose = true "Full orbit: RK4" begin
#    @test isapprox(rmse_FO, 0.0, atol=0.001)
#end # testset GCA: Euler
    
