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
tf = 100.0 #n=100   # End time of simulation [s]         |
tspan = (0.0, tf)
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(Int64, numparticles) # Specifies the species of the particles 
mass = 1.0
charge = 1.0

#...............................................
# MAGNETIC FIELD PARAMETERS
#qMm = 20.0
#println("qMm: $qMm")
M = qMm*mass/charge
itp_type = Gridded(Linear())
itp_bc = Flat()

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
# It will be zero

#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = createaxes(xi0, xif, ni)
Bfield = Array{Float64,  4}(undef, numdims, ni...)
Efield = zeros(Float64, size(Bfield))
discretise!(Bfield, xx, yy, zz, magneticdipolefield, M)
# Create interpolation objects
emfields = eachslice(vcat(Bfield, Efield), dims=(2,3,4))
emfields_itp = interpolate((xx, yy, zz), emfields, itp_type)
emfields_itp = extrapolate(emfields_itp, itp_bc)

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
B⃗ = emfields_itp(R0...)[1:3]
B = norm(B⃗)
b̂ = B⃗/B
v = norm(vel0)
vparal = vel0 ⋅ b̂
vperp = √(v^2 - vparal^2)
μ = mass*vperp^2/(2B)

fo_params = (charge, mass, emfields_itp)
gca_params = (charge, mass, μ, emfields_itp)

#-------------------------------------------------------------------------------
# CREATE PROBLEM
fo_prob = ODEProblem(lorentzforce!, [pos0; vel0], tspan, fo_params)
gca_prob = ODEProblem(
    guidingcentreapproximation!, [R0; vparal], tspan, gca_params
    )
#...............................................................................
# RUN SIMULATION
fo_sim = DifferentialEquations.solve(fo_prob)
gca_sim = DifferentialEquations.solve(gca_prob)
println("qMm = $qMm")
println("nof. FO-timesteps = $(length(fo_sim.t))")
println("nof. GCA-timesteps = $(length(gca_sim.t))")
#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
ϕ_FO = atan(fo_sim.u[end][2]/fo_sim.u[end][1])
ϕ_GCA = atan(gca_sim.u[end][2]/gca_sim.u[end][1])

v = norm(vel0)
vparall = abs((B⃗ ⋅ vel0)/B)
vperp = √(v^2 - vparall^2)
α = atan(vperp/vparall)
R0 = 1.0

function T_dipole(q, m, M, R0, v, α)
    return @. 2π*q*M/(m*v^2*R0)*(1 - 1/3*sin(α)^0.62)
end # function T_dipole

T = T_dipole(charge, mass, M, R0, v, α)
global ϕ = 2π*tf/T
