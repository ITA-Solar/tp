# Nullpoint simulation with 
#  Uniform distribution in position 
#  Maxwellian velocity distributions
using tp
using Random
using DifferentialEquations
using JLD2
using Interpolations
# Progress bar integration
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())


# SET EXPEIRMENT PARAMETERS
# -----------------------------------------------------------------------------
# Particle parameters
dtsnap = 1e-2 # in code units (hs)
tspan = (0.0, dtsnap*100)   # Simulation time-span
npart = convert(Int64, 2e4) # Number of particles
xbounds = (14.527344e6, 18.824219e6) # Bounds of of initial position in x
ybounds = (1e6, 1e6)                 # Bounds of of initial position in x
zbounds = (-7.78315e6,-3.8759463e6)  # Bounds of of initial position in x
charge = -tp.e              # Charge of particles
mass = tp.m_e               # Mass of particles
seed = 0                    # Random number generator seed
rng = MersenneTwister(seed) # Random number generator
wp_part=Float64             # Working precision of particles
units = "SI"               # Choice of units for Bifrost snapshot

# Bifrost-input parameters
br_expname="nullpoint"                        # Bifrost experiment name
br_expdir="/Users/eilifo/code/repos/data/8K"  # Bifrost experiment directory
br_snap=1490          # Bifrost experiment snapshot
wp_snap=Float32  # Working precision of Bifrost snapshot
itp_bc = (Flat(), Flat()) # Interpolation boundary conditions

# Parameters for saving data
tp_expname="nullpoint8K_30M_FF_1s_serial"    # Simulation experiment name
tp_expdir="/Users/eilifo/code/repos/data/8K" # Where to save results
 
# Other options
eom = gca_2Dxz! # Equation of motion


# DEFINE PARTICLE INITIAL CONDITIONS  AND ODE-PARAMETERS
# -----------------------------------------------------------------------------
# Background field interpolation object
println("Fetching Bifrost fields...")
@time fields_itp = get_br_emfield_vecof_interpolators(
    br_expname,
    br_snap,
    br_expdir
    ;
    itp_bc = itp_bc
    )
br_temp_itp = get_br_var_interpolator(
    br_expname,
    br_snap,
    br_expdir,
    "tg";
    units=units,
    itp_bc = itp_bc
    )

function prob_func(prob, _, _)
    x0 = rand(rng, wp_part, xbounds...)
    y0 = rand(rng, wp_part, ybounds...)
    z0 = rand(rng, wp_part, zbounds...)
    vel = maxwellianvelocitysample(rng, br_temp_itp, mass, x0, z0)
    fields_at_pos = [itp(x0, z0) for itp in fields_itp] # Has type Vector{Vector}
    R, vparal, magneticmoment = get_guidingcentre(
        [x0, y0, z0],
        vel,
        fields_at_pos[1:3],
        fields_at_pos[4:6],
        charge,
        mass
        )
    remake(prob, u0=[R; vparal],
        p=(charge=charge, mass=mass, μ=magneticmoment, itps=fields_itp)
    )
end
# DEFINE PROBLEMS
# -----------------------------------------------------------------------------
# Define the problem for the first particle
prob = ODEProblem(
    eom,
    zeros(5),
    tspan,
    )

# Define the ensamble problem for the whole lot
ensamble_prob = EnsembleProblem(
        prob, 
        prob_func=prob_func
        ;
        #safetycopy=false,
    )

# SOLVE PROBLEM
# -----------------------------------------------------------------------------
println("Running simulation...")
@time sim = DifferentialEquations.solve(
    ensamble_prob,
    EnsembleSerial(),
    ;
    trajectories=npart,
    progress=true,
    maxiters=10000,
    save_everystep=false,
    )

# SAVE THE RESULTS
# -----------------------------------------------------------------------------
println("Saving results...")
@save joinpath(tp_expdir, tp_expname) sim
