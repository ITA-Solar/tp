# Created 02.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 mirroring.jl
#
#-------------------------------------------------------------------------------
# This experiments tests the mirroring term in the GCA solver and mirroring in
# the full orbit solver by comparing with an analytic expression of the motion
# parallel to the magnetic field. The analytic magnetic bottle field is obtained
# from Ripperda et al. 2018, and the motion of the particle is assume adiabatic
# such that the GCA holds.
#
# The GC is initialised in the centre of the bottle to avoid any other drifts.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 1.e-3        # Time step [s]                      |
tf = 30.0 #n=100   # End time of simulation [s]         |
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(Int64, numparticles) # Specifies the species of the particles 
mass = specieTable[species[1], 1]
charge = specieTable[species[1], 2]

#...............................................
# MAGNETIC FIELD PARAMETERS
B0 = 30.0 # Magnetic field strenth parameter
L = 0.4 # Mirroring length

#...............................................
# INITIAL CONDITIONS
vel0 = [0.0, 0.1, 0.1]
rL = mass*√(vel0[1]^2 + vel0[2]^2)/(charge*B0)
pos0 = [-rL, 0.0, 0.0]
R0 = [0.0, 0.0, 0.0]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
a = 2.1L
xi0 = (-a, -a, -a)
# Upper bound of the three spatial axes
xif = (a, a, a)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 100)
 
#...............................................
# ELECTRIC FIELD PARAMETERS
# It will not exist
Ex = 0.0
Ey = 0.0

#...............................................
# SOLVER CONDITIONS
FOscheme = rk4  # Scheme for integrating the full-orbit eqs.
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
discretise!(Bfield, xx, yy, zz, mirroringfield, B0, L)

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
Bmax = B0*norm(vel0)^2/(vel0[1]^2 + vel0[2]^2)
zmax = L*√(Bmax/B0 - 1)

function z(t, μ, mass, B0, L, A, ϕ)
    ω = √(2μ*B0/(mass*L^2))
    return @. A*sin(ω*t + ϕ)
end
ϕ = 0
A = zmax
times = collect(range(0.0, step=dt, length=numsteps+1))
z_anal = z(times, μ[1], mass, B0, L, A, ϕ)

z_GCA  = patchGCA.tp.R[3,1,:]
z_FO = patchFO.tp.pos[3,1,:]
rmse_GCA = √(sum((z_anal .- z_GCA).^2)/numsteps)
rmse_FO = √(sum((z_anal .- z_FO).^2)/numsteps)
@testset verbose = true "GCA: Euler" begin
    @test isapprox(rmse_GCA, 0.0, atol=0.001)
end # testset GCA: Euler
@testset verbose = true "Full orbit: RK4" begin
    @test isapprox(rmse_FO, 0.0, atol=0.001)
end # testset GCA: Euler
    
