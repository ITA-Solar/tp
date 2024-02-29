# Nullpoint simulation with 
#  Uniform distribution in position 
#  Maxwellian velocity distributions
@everywhere using tp
@everywhere using Random
@everywhere using DifferentialEquations
@everywhere using JLD2
@everywhere using Interpolations
using Distributed
# Progress bar integration
#using Logging: global_logger
#using TerminalLoggers: TerminalLogger
#global_logger(TerminalLogger())


# SET EXPEIRMENT PARAMETERS
# -----------------------------------------------------------------------------
# Particle parameters
dtsnap = 1e-2 # in code units (hs)
tspan = (0.0, dtsnap*100)   # Simulation time-span
npart = convert(Int64, 3e7) # Number of particles
@everywhere xbounds = (14.527344e6, 18.824219e6) # Bounds of of initial position in x
@everywhere ybounds = (1e6, 1e6)                 # Bounds of of initial position in x
@everywhere zbounds = (-7.78315e6,-3.8759463e6)  # Bounds of of initial position in x
@everywhere charge = -tp.e              # Charge of particles
@everywhere mass = tp.m_e               # Mass of particles
@everywhere seed = 0                    # Random number generator seed
@everywhere rng = MersenneTwister(seed + myid()) # Random number generator
@everywhere wp_part=Float64             # Working precision of particles
units = "SI"               # Choice of units for Bifrost snapshot

# Bifrost-input parameters
br_expname="nullpoint"                        # Bifrost experiment name
# Bifrost experiment directory
br_expdir="/mn/stornext/u3/eilifo/data@d9/bifrost/ohf/run0/8K-ext_damp"
br_snap=1490          # Bifrost experiment snapshot
wp_snap=Float32  # Working precision of Bifrost snapshot
itp_bc = (Flat(), Flat()) # Interpolation boundary conditions

# Parameters for saving data
tp_expname="nullpoint8K_30M_FF_1s_16procs_lowmem.jld2"    # Simulation experiment name
#tp_expname="test.jld2"    # Simulation experiment name
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
    itp_bc = itp_bc,
    units=units,
    destagger=true,
    )
br_temp_itp = get_br_var_interpolator(
    br_expname,
    br_snap,
    br_expdir,
    "tg";
    units=units,
    itp_bc = itp_bc
    )
#remote_do.(()->(fields_itp, br_temp_itp), workers())
progresscount = 0
@everywhere updatecount() = global progresscount +=1;

@everywhere function prob_func(prob, i, _)
    remote_do(updatecount, 1)
    if myid() == 1
    end
    x0 = rand(rng, wp_part, xbounds...)
    y0 = rand(rng, wp_part, ybounds...)
    z0 = rand(rng, wp_part, zbounds...)
    vel = maxwellianvelocitysample(rng, prob.p[5], mass, x0, z0)
    # Has type Vector{Vector}
    bfield_at_pos = [prob.p[4][i](x0, z0) for i in 1:3]
    efield_at_pos = [prob.p[4][i](x0, z0) for i in 4:6]
    R, vparal, magneticmoment = get_guidingcentre(
        [x0, y0, z0],
        vel,
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
        )
    remake(prob, u0=[R; vparal],
    p=(charge, mass, magneticmoment, prob.p[4], prob.p[5])
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
@everywhere function output_func(sol, i)
    u0 = first(sol)
    uf = last(sol)
    B0_vec = [sol.prob.p[4][i](u0[1], u0[3]) for i in 1:3]
    E0_vec = [sol.prob.p[4][i](u0[1], u0[3]) for i in 4:6]
    Bf_vec = [sol.prob.p[4][i](uf[1], uf[3]) for i in 1:3]
    Ef_vec = [sol.prob.p[4][i](uf[1], uf[3]) for i in 4:6]
    return (
        u0,
        uf,
        sol.prob.p[3], # the magnetic moment
        B0_vec,
        E0_vec,
        Bf_vec,
        Ef_vec,
        length(sol.t), 
        last(sol.t)
        ), 
        false
end
# DEFINE PROBLEMS
# -----------------------------------------------------------------------------
# Define the problem for the first particle
prob = ODEProblem(
    eom,
    zeros(5),
    tspan,
    (charge, mass, nothing, fields_itp, br_temp_itp)
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
running = true
@async begin
    while running
        print("Particle $progresscount of $npart \r")
        sleep(0.1)
        flush(stdout)
    end
end
@time sim = DifferentialEquations.solve(
    ensemble_prob,
    EnsembleDistributed()
    #EnsembleSerial(),
    #EnsembleThreads(),
    ;
    trajectories=npart,
    progress=false,
    maxiters=10000,
    #save_everystep=false,
    callback=cb
    )
running = false

# SAVE THE RESULTS
# -----------------------------------------------------------------------------
resultpath = joinpath(tp_expdir, tp_expname)
println("Saving results to $resultpath...")
@save resultpath sim
println("success.")
