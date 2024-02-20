# Nullpoint simulation with 
#  Uniform distribution in position 
#  Maxwellian velocity distributions
using tp
using Random
using DifferentialEquations
using JLD2
using Interpolations
# Progress bar integration
#using Logging: global_logger
#using TerminalLoggers: TerminalLogger
#global_logger(TerminalLogger())


# SET EXPEIRMENT PARAMETERS
# -----------------------------------------------------------------------------
# Particle parameters
dtsnap = 1e-2 # in code units (hs)
tspan = (0.0, dtsnap*100)   # Simulation time-span
npart = convert(Int64, 3e6) # Number of particles
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
# Bifrost experiment directory
br_expdir="/mn/stornext/u3/eilifo/data@d9/bifrost/ohf/run0/8K-ext_damp"
br_snap=1490          # Bifrost experiment snapshot
wp_snap=Float32  # Working precision of Bifrost snapshot
itp_bc = (Flat(), Flat()) # Interpolation boundary conditions

# Parameters for saving data
#tp_expname="nullpoint8K_30M_FF_1s_serial.jld2"    # Simulation experiment name
tp_expname="test"    # Simulation experiment name
# Where to save the results
tp_expdir="/mn/stornext/u3/eilifo/data@d9/bifrost/ohf/run0/8K-ext_damp"
 
# Other options
eom = tp.gca_2Dxz! # Equation of motion
fieldfunc = tp.get_br_emfield_vecof_interpolators


# DEFINE PARTICLE INITIAL CONDITIONS  AND ODE-PARAMETERS
# -----------------------------------------------------------------------------
# Background field interpolation object
println("Fetching Bifrost fields...")
@time fields_itp = fieldfunc(
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

ninit = 0
lk = ReentrantLock()
function prob_func(prob, i, _)
    global ninit
    lock(lk)
    try
        ninit += 1
        print("Particle $ninit of $npart \r")
        flush(stdout)
    finally
        unlock(lk)
    end
    x0 = rand(rng, wp_part, xbounds...)
    y0 = rand(rng, wp_part, ybounds...)
    z0 = rand(rng, wp_part, zbounds...)
    vel = maxwellianvelocitysample(rng, br_temp_itp, mass, x0, z0)
    # Has type Vector{Vector}
    bfield_at_pos = [fields_itp[i](x0, z0) for i in 1:3]
    efield_at_pos = [fields_itp[i](x0, z0) for i in 4:6]
    R, vparal, magneticmoment = get_guidingcentre(
        [x0, y0, z0],
        vel,
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
        )
    remake(prob, u0=[R; vparal],
        p=(charge, mass, magneticmoment, fields_itp)
    )
end

# Define an output function - what to be saved from each particle solution
# The ODESolution type contains a lot of metadata, including alogorithm
# choice, number of algorithm switches, etc. In a large ensemble, this
# overhead will quickly consume a lot of memory.
#
# In this output function we save only the initial and final state of the
# particle, in addition to its timesteps. We set the required "rerun" return
# argument to false.
function output_func(sol, i)
    return ([first(sol), last(sol)], length(sol.t), last(sol.t)), false
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
ensemble_prob = EnsembleProblem(
        prob
        ;
        prob_func=prob_func,
        output_func=output_func,
        safetycopy=false,
    )

# SOLVE PROBLEM
# -----------------------------------------------------------------------------
# First define callbacks
cb = DiscreteCallback(tp.outside2dnullpointzoom, tp.killparticle!)

println("Running simulation...")
@time sim = DifferentialEquations.solve(
    ensemble_prob,
    EnsembleSerial(),
    #EnsembleThreads(),
    ;
    trajectories=npart,
    progress=false,
    maxiters=10000,
    #save_everystep=false,
    callback=cb
    )

# SAVE THE RESULTS
# -----------------------------------------------------------------------------
resultpath = joinpath(tp_expdir, tp_expname)
println("Saving results to $resultpath")
@save resultpath sim
