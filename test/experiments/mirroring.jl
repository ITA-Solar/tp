# Created 02.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
    mirroring.jl

This experiments tests the mirroring term in the GCA solver and mirroring in
the full orbit solver by comparing with an analytic expression of the motion
parallel to the magnetic field. The analytic magnetic bottle field is obtained
from Ripperda et al. 2018, and the motion of the particle is assume adiabatic
such that the GCA holds.

The GC is initialised in the centre of the bottle to avoid any other drifts.
"""

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 1.e-3        # Time step [s]                      |
tf = 30.0 #n=100   # End time of simulation [s]         |
tspan = (0, tf)
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
mass = 1
charge = 1

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
# RK4 for GCA? Euler-Cromer for full orbit?

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
discretise!(Bfield, xx, yy, zz, magneticmirrorfield, B0, L)
emfields = eachslice(vcat(Bfield, Efield), dims=(2,3,4))
emfields_itp = linear_interpolation((xx, yy, zz), emfields, 
    extrapolation_bc=Flat()
    )

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf/dt)   # Number of timesteps in the simulation
#println("Number of time steps = $numsteps.")

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
B⃗ = emfields_itp(pos0...)[1:3]
E⃗ = zeros(3)
B = norm(B⃗)
b̂ = B⃗/B
v = norm(vel0)
vparal = vel0 ⋅ b̂
vperp = √(v^2 - vparal[1]^2)
μ = mass*vperp^2/(2B)

#-------------------------------------------------------------------------------
# CREATE PROBLEM
prob_FO = ODEProblem(
    lorentzforce!, 
    [pos0; vel0],
    tspan,
    (charge, mass, emfields_itp),
    )
prob_GCA = ODEProblem(
    guidingcentreapproximation!,
    [R0; vparal],
    tspan,
    (charge, mass, μ, emfields_itp)
    )

#...............................................................................
# RUN SIMULATION
sol_FO = DifferentialEquations.solve(prob_FO)
sol_GCA = DifferentialEquations.solve(prob_GCA)

#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
Bmax = B0*norm(vel0)^2/(vel0[1]^2 + vel0[2]^2)
zmax = L*√(Bmax/B0 - 1)

function z(t, μ, mass, B0, L, A, ϕ)
    ω = √(2μ*B0/(mass*L^2))
    return @. A*sin(ω*t + ϕ)
end
ϕ = 0
A = zmax
times = collect(range(0.0, step=dt, length=numsteps+1))
z_anal_FO = z(sol_FO.t, μ, mass, B0, L, A, ϕ)
z_anal_GCA = z(sol_GCA.t, μ, mass, B0, L, A, ϕ)

z_GCA  = [u[3] for u in sol_GCA.u] 
z_FO  = [u[3] for u in sol_FO.u] 
rmse_GCA = √(sum((z_anal_GCA .- z_GCA).^2)/numsteps)
rmse_FO = √(sum((z_anal_FO .- z_FO).^2)/numsteps)
@testset verbose = true "GCA: Default alg." begin
    @test isapprox(rmse_GCA, 0.0, atol=0.001)
end # testset GCA: Euler
@testset verbose = true "Full orbit: Default alg." begin
    @test isapprox(rmse_FO, 0.0, atol=0.010)
end # testset GCA: Euler
    
