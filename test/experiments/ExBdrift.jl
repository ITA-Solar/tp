#-------------------------------------------------------------------------------
# Created 19.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 particleDrift.jl
#
#-------------------------------------------------------------------------------
# Experiment testing if the ExB-drift of particles in a static homogeneous 
# electromagnetic field follows the analytical solution.
#-------------------------------------------------------------------------------

"""
    testExBdrift.jl

    This experiment tests numerical particle paths produced by the tp-code 
against analytical solutions.
    An electron and proton in a static and homogenous electromagnetic field.
The electric field is normal to the magnetic field such that the particles will 
experience an ExB-drift. The resulting numerical paths are compared to the
analytic paths, which root mean square errors entails the test assertions. 

    In practice, this experiment also tests dependent type-constructors and
methods, such as:
    ParticleSoA(pos::Matrix, vel::Matrix, species::Vector, numSteps::Integer)
    Mesh(bField::Array{4},
         eField::Array{4},
         xx::Vector,
         yy::Vector,
         zz::Vector)
    run!(patch:Patch)
and the chosen numerical solver, scheme and interpolation chosen.
"""

#-------------------------------------------------------------------------------
# EXPERIMENT PARAMTERS
#
numDims = 3       # Number of spatial dimensions
numParticles = 2  # Number of particles to simulae
species = [1, 2]  # Specifies the species of the particles 
#   (electron, proton)
tf = 2*2π/(5.93096958e7) # s. End time of simulation
tspan = (0,tf)
#   multiple of electron (1eV) gyrofreq. at Bz=3.3e-4 G
numSteps = 1000   # Number of timesteps in the simulation
dt = tf/numSteps  # Size of timestep
# The electromagnetic field is static and homogeneous. Hence the mesh only 
# needs one cell, i.e. the boundary.
Bx = 0.0  # magnetic field component only in the positive z-direction
By = 0.0
Bz = 3.37213e-4 # G (so that electron larmor radius is 1 cm)
Ex = 0.0
Ey = 0.02Bz/tf  # Electric field component only in the positive y-direction
Ez = 0.0

# Particle conditions
qe = specieTable[1,2] # Electron charge
qp = specieTable[2,2] # Proton charge
me = specieTable[1,1] # Electron mass
mp = specieTable[2,1] # Proton mass
# Set initial position and velocities
vdriftx = Ey/Bz
# Electron initial velocity
vxe = 0.0#593_096.95848 # m/s (1eV electron
vye = 0.0
vze = 0.0
# Proton initial velocity
vxp = 13_841.122177 # m/s (1eV proton)
vyp = 0.0
vzp = 0.0
# Electron initial position
x0e = 0.0
y0e = 0.0
z0e = 0.0
# Proton initial position
x0p = 0.0
y0p = 0.0
z0p = 0.0

#-------------------------------------------------------------------------------
# ANALYTICAL SOLUTION 
times = LinRange(0, tf, numSteps + 1)

# Gyrofrequency
function ω(q, B, m)
    return q*B/m
end
# Larmor radius
function larmorRadius(ω, vperp)
    return vperp/ω
end
# x-position
function x(t, x0, ωt, rL, Ey, Bz)
    return @. rL*sin(ωt*t) + Ey/Bz*t
end
# y-positoin
function y(t, y0, ωt, rL, Ex, Bz)
    return @. -rL + rL*cos(ωt*t) - Ex/Bz*t
end
# z-position
function z(t, z0, vz0, ω, Ez, Bz)
    return @. z0 + 0.5ω/Bz*Ez*t^2 + vz0*t
end

# Define gyrofrequency and Larmor radius for both particles
ωe = ω(qe, Bz, me)
ωp = ω(qp, Bz, mp)
vperpe = vxe - vdriftx # Subtract drift velocity to find vperp
vperpp = vxp - vdriftx
rLe = larmorRadius(ωe, vperpe)
rLp = larmorRadius(ωp, vperpp)

# Compute time evolution of position of both particles
pose = zeros(Float64, numDims, numSteps + 1)
posp = zeros(Float64, numDims, numSteps + 1)
pose[1, :] .= x(times, x0e, ωe, rLe, Ey, Bz)
pose[2, :] .= y(times, y0e, ωe, rLe, Ex, Bz)
pose[3, :] .= z(times, z0e, vze,  ωe, Ez, Bz)
posp[1, :] .= x(times, x0p, ωp, rLp, Ey, Bz)
posp[2, :] .= y(times, y0p, ωp, rLp, Ex, Bz)
posp[3, :] .= z(times, z0p, vzp, ωp, Ez, Bz)


#-------------------------------------------------------------------------------
# NUMERICAL SOLUTION 
#-------------------------------------------------------------------------------
# MESH CREATION
B = zeros(Float64, numDims, 2, 2, 2)
E = zeros(Float64, numDims, 2, 2, 2)
B[3, :, :, :] .= Bz
E[1, :, :, :] .= Ex
E[2, :, :, :] .= Ey
E[3, :, :, :] .= Ez
# Create coordinates of cell
xx = [-0.1, 1.1]
yy = [-0.1, 1.1]
zz = [-0.1, 1.1]
# Create Mesh instance
mesh = Mesh(B, E, xx, yy, zz)
# Create mesh without fields (for new DifferentialEquations.jl implementation)
puremesh = PureMesh(xx, yy, zz)

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Full orbit particles
# Set initial position and velocities
pos = zeros(Float64, numDims, numParticles) # Position -> origin
vel = zeros(Float64, numDims, numParticles) # Velocity
vel[:, 1] = [vxe, vye, vze]  # Initial velocity electron
vel[:, 2] = [vxp, vyp, vzp]  # Initial velocity proton
pos[:, 1] = [x0e, y0e, z0e]  # Initial position electron
pos[:, 2] = [x0p, y0p, z0p]  # Initial position proton
particlesEC = ParticleSoA(pos, vel, species, numSteps)
particlesVay = ParticleSoA(pos, vel, species, numSteps)
#  Create DifferentialEquations.jl particle type
ic_FO = FO_IC(
    VectorIC(pos[1,:]), 
    VectorIC(pos[2,:]), 
    VectorIC(pos[3,:]), 
    VectorIC(vel[1,:]), 
    VectorIC(vel[2,:]), 
    VectorIC(vel[3,:])
    )
params_FO = FOParams(
    [-tp.e, tp.e],
    [tp.m_e, tp.m_p],
    puremesh,
    B,
    E,
    )
DEparticlesFO = ODEParticle(LorentzForce(), ic_FO, params_FO, numParticles)

# GCA particles
# Set initial position and velocities
R = zeros(Float64, numDims, numParticles) # Position -> origin
velGCA = zeros(Float64, numDims + 3, numParticles) # Velocity
R[:, 1] = [x0e, y0e, z0e]  # Initial position electron
R[:, 2] = [x0p, y0p, z0p]  # Initial position proton
vperpe = √(vxe^2 + vye^2)
vperpp = √(vxp^2 + vyp^2)
μₑ = me*vperpe^2/2Bz
μₚ = mp*vperpp^2/2Bz
vparal = zeros(numParticles) # Velocity parallell to the magnetic field
vparal[1] = vel[3,1] 
vparal[2] = vel[3,2] 
μ = [μₑ, μₑ]               # Magnetic moment of the particles

# Create ParticlesSoA-instance
particlesGCA = GCAParticleSoA(R, vparal, μ, species, numSteps)
#  Create DifferentialEquations.jl particle type
ic_GCA = GCA_IC(
    VectorIC(R[1,:]), 
    VectorIC(R[2,:]), 
    VectorIC(R[2,:]), 
    VectorIC(vparal)
    )
params_GCA = GCAParams(
    [-tp.e, tp.e],
    [tp.m_e, tp.m_p],
    [μₑ, μₚ],
    puremesh,
    B, E, mesh.∇B, mesh.∇b̂, mesh.∇ExB,
    )
DEparticlesGCA = ODEParticle(GCA(), ic_GCA, params_GCA, numParticles)


#-------------------------------------------------------------------------------
# CREATE PATCH
# Non-relativisitc Euler-Cromer
patchEC = Patch(mesh,
                particlesEC,
                fullOrbit,
                eulerCromer,
                trilinear,
                dt,
                numSteps,
                numParticles
                )
# Relativistic Vay pusher
patchVay = Patch(mesh,
                 particlesVay,
                 relFullOrbitExplLeapFrog,
                 vay,
                 trilinear,
                 dt,
                 numSteps,
                 numParticles)
# Relativistic Boris pusher
#patchBoris = Patch(mesh,
#              particles,
#              relFullOrbitExplLeapFrog,
#              boris,
#              trilinear,
#              dt,
#              numSteps,
#              numParticles)

# non-relativistic GCA
patchGCA = Patch(mesh,
                 particlesGCA,
                 gca,
                 rk4,
                 trilinearGCA,
                 dt,
                 numSteps,
                 numParticles)

# Using DifferentialEquations.jl
DEpatchFO = DEPatch(tspan, DEparticlesFO, puremesh)
DEpatchGCA = DEPatch(tspan, DEparticlesGCA, puremesh)

#-------------------------------------------------------------------------------
# RUN SIMULATION
run!(patchEC)
run!(patchVay)
#run!(patchBoris)
run!(patchGCA)
simFO = run!(DEpatchFO)
simGCA = run!(DEpatchGCA)

#-------------------------------------------------------------------------------
# CALCULATE DEVIATIONS FROM ANALYTICAL SOLUTION
# Calculate the root mean squared error between numerical and analytical 
#   trajectories for both electron and proton for all position components.

# Euler-Cromer
# Electron
numPose = patchEC.tp.pos[:, 1, :]
numPosp = patchEC.tp.pos[:, 2, :]
rmsErrexEC = √(sum((pose[1,2:numSteps] .-
    numPose[1,2:numSteps]).^2)/numSteps)
rmsErreyEC = √(sum((pose[2,2:numSteps] .- 
    numPose[2,2:numSteps]).^2)/numSteps)
rmsErrezEC = √(sum((pose[3,2:numSteps] .- 
    numPose[3,2:numSteps]).^2)/numSteps)
# Proton
rmsErrpxEC = √(sum((posp[1,2:numSteps] .- 
    numPosp[1,2:numSteps]).^2)/numSteps)
rmsErrpyEC = √(sum((posp[2,2:numSteps] .- 
    numPosp[2,2:numSteps]).^2)/numSteps)
rmsErrpzEC = √(sum((posp[3,2:numSteps] .- 
    numPosp[3,2:numSteps]).^2)/numSteps)

# Vay
# Electron
numPose = patchVay.tp.pos[:, 1, :]
numPosp = patchVay.tp.pos[:, 2, :]
rmsErrexVay = √(sum((pose[1,2:numSteps] .-
    numPose[1,2:numSteps]).^2)/numSteps)
rmsErreyVay = √(sum((pose[2,2:numSteps] .-
    numPose[2,2:numSteps]).^2)/numSteps)
rmsErrezVay = √(sum((pose[3,2:numSteps] .-
    numPose[3,2:numSteps]).^2)/numSteps)
# Proton
rmsErrpxVay = √(sum((posp[1,2:numSteps] .-
    numPosp[1,2:numSteps]).^2)/numSteps)
rmsErrpyVay = √(sum((posp[2,2:numSteps] .-
    numPosp[2,2:numSteps]).^2)/numSteps)
rmsErrpzVay = √(sum((posp[3,2:numSteps] .-
    numPosp[3,2:numSteps]).^2)/numSteps)

# GCA
# Electron
numPose = patchGCA.tp.R[:, 1, :]
numPosp = patchGCA.tp.R[:, 2, :]
rmsErrexGCA = √(sum((pose[1,2:numSteps] .-
    numPose[1,2:numSteps]).^2)/numSteps)
rmsErreyGCA = √(sum((pose[2,2:numSteps] .-
    numPose[2,2:numSteps]).^2)/numSteps)
rmsErrezGCA = √(sum((pose[3,2:numSteps] .-
    numPose[3,2:numSteps]).^2)/numSteps)
# Proton
rmsErrpxGCA = √(sum((posp[1,2:numSteps] .-
    numPosp[1,2:numSteps]).^2)/numSteps)
rmsErrpyGCA = √(sum((posp[2,2:numSteps] .-
    numPosp[2,2:numSteps]).^2)/numSteps)
rmsErrpzGCA = √(sum((posp[3,2:numSteps] .-
    numPosp[3,2:numSteps]).^2)/numSteps)

#-------------------------------------------------------------------------------
# PLOT RESULTS
if length(ARGS) >= 1 && ARGS[1] == "plot"
    using Plots
    println("in")
    pos = patchEC.tp.pos
    p1 = plot(pos[1, 1, :], pos[2, 1, :], label="Numerical",
              title="Electron", xlabel="x, m", ylabel="y, m")
    #p1 = plot!(pose[1, :], pose[2, :], label="Analytical", ls=:dash)
    p2 = plot(pos[1, 2, :], pos[2, 2, :], label="Numerical",
              title="Proton", xlabel="x, m", ylabel="y, m")
    #p2 = plot!(posp[1, :], posp[2, :], label="Analytical", ls=:dash)
    
    plot(p1, p2, layout=(2,1))
    gui()
end # if plot
#-------------------------------------------------------------------------------
# TEST RESULTS
@testset verbose = true "FO: Euler-Cromer" begin
    @test rmsErrexEC ≈  6.066588298648159e-5
    @test rmsErreyEC ≈ 5.5903789343654834e-5
    @test rmsErrezEC == 0.0
    @test rmsErrpxEC ≈  7.641789560048872e-14
    @test rmsErrpyEC ≈ 3.3694076090881336e-8
    @test rmsErrpzEC == 0.0
end # testset Euler-Cromer
#
@testset verbose = true "Vay" begin
    @test rmsErrexVay ≈  1.0902384513274832e-7
    @test rmsErreyVay ≈ 1.0818615735766856e-7
    @test rmsErrezVay == 0.0
    @test rmsErrpxVay ≈ 3.772326801578384e-14
    @test rmsErrpyVay ≈  2.802190198707962e-14
    @test rmsErrpzVay == 0.0
end # testset Vay
#
@testset verbose = true "GCA: RK4" begin
    @test isapprox(patchGCA.tp.R[1,1,end], 0.02, rtol=1e-13)
    @test isapprox(patchGCA.tp.R[2,1,end], 0.00, atol=1e-15)
    @test isapprox(patchGCA.tp.R[3,1,end], 0.00, atol=1e-15)
end # testset GCA: RK4

# Full orbit using DifferentialEquations.jl
@testset verbose = true "DE-FO" begin
    # Electron
    @test isapprox(simFO.u[1].u[end][1], pose[1,end], rtol=1e-3)
    @test isapprox(simFO.u[1].u[end][2], pose[2,end], atol=1e-5)
    @test isapprox(simFO.u[1].u[end][3], pose[3,end], rtol=1e-9)
    # Proton
    @test isapprox(simFO.u[2].u[end][1], posp[1,end], rtol=1e-9)
    @test isapprox(simFO.u[2].u[end][2], posp[2,end], rtol=1e-9)
    @test isapprox(simFO.u[2].u[end][3], posp[3,end], rtol=1e-9)
end

# GCA using DifferentialEquations.jl
@testset verbose = true "DE-GCA" begin
    # Electron
    @test isapprox(simGCA.u[1].u[end][1], 0.02, rtol=1e-13)
    @test isapprox(simGCA.u[1].u[end][2], 0.00, atol=1e-15)
    @test isapprox(simGCA.u[1].u[end][3], 0.00, rtol=1e-15)
end


