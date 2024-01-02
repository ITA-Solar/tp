# Nullpoint simulation with 
#  4092x4092x1 grid points (approx. 14km res)
#  Uniform distribution in position 
#  Maxwellian velocity distribution
using PyPlot
const plt = PyPlot
using tp
using LinearAlgebra

#__________________________SIMULATION TIME SPAN________________________________
# Let tspan be (0.0, dtsnap)
dt = 1e-3
dtsnap = 1e-2
dtsnapSI = tp.code2SI_t*dtsnap
tspan = (0.0, dtsnapSI)
display(tspan)
#__________________________NUMBER OF PARTICLES_________________________________
npart = convert(Int64, 1e4)


#__________________________SPATIAL SPAN OF PARTICLES___________________________
# Point in right exhaust
#posx = [16.98, 16.4]*code2SI_l
#posy = [1.0,1.0]*code2SI_l
#posz = [-5.71,-5.710]*code2SI_l
#
# Slicex4k and slicez4k
# Approximately the same as slice4k ends.
posx = [14.527344, 18.824219]*tp.code2SI_l
posy = [1.0,1.0]*tp.code2SI_l
posz = [-7.78315,-3.8759463]*tp.code2SI_l


#__________________________SET PARAMETERS HERE_________________________________
params = Parameters(
    patch_type="DE",
    tspan=tspan,
    charge = -tp.e,
    mass = tp.m_e,
    npart=npart,
    tp_expname="sde_proj_npart1e4_kd5e-11_8k_fixed_bc",
    #tp_expname="sde_proj_npart1e4_kd5e-8_8k",
    tp_expdir="/Users/eilifo/repos/data/sde",
    br_expname="nullpoint",
    br_expdir="/Users/eilifo/repos/data/8K",
    br_isnap=1490,
    solver=GCAPitchAngleFriction(),
    interp="bilinear_xz",
    wp_part=Float64,
    wp_snap=Float32,
    pos_distr="uniform",
    vel_distr="mb",
    posxbounds=posx,
    posybounds=posy,
    poszbounds=posz,
    bg_input="br",
    pbc=(false, true, false),
    SI_units=true,
    seed=convert(Int64, 0),
)

#exp = DE_init!(params)
# Get fields from Bifrost snap
mesh = create_puremesh(params)
bfield, efield, density, gas_temp = get_fields_from_br(params;
    get_bfield=true,
    get_efield=true,
    get_density=true,
    get_gas_temp=true,
    )
num_density = density ./ tp.m_p

# Create initial conditions for the particles
function init()
    pos, vel = tp_get_initial_pos_and_vel!(params, mesh)
    B_at_pos, E_at_pos = EMfield_at_pos(pos, bfield, efield, mesh)    
    R0, vparal0, mu0 = calc_GCA_IC_and_mu(
        pos, 
        vel, 
        params.charge, 
        params.mass,
        B_at_pos,
        E_at_pos,   
        )
    beta0 = pitchangle(params.mass, B_at_pos, vparal0, mu0)
    # Make half of the beta0-values negative
    beta0[1:convert(Int64, npart/2)] *= -1.0
    ic = InitialConditions(
        R0[1,:],
        R0[2,:],
        R0[3,:],
        vparal0,
        beta0,
        )
    # Create parameters for the EoM
    coulomb_logarithm = 5e-11 # 5e-11 = mfp = 6e9/5, 5e-8 = mfp 6e6/5
    gradB, gradb, gradExB = compute_gradients(bfield, efield,
        mesh.x, mesh.y, mesh.z, derivateUpwind
        )
    drift_params = GCAPitchAngleScatteringParams(
        params.mass, 
        params.charge,
        coulomb_logarithm,
        mesh,
        bfield,
        efield,
        gradB,
        gradb,
        gradExB,
        num_density,
        gas_temp,
        )
    
    # Create particle object
    particles = SDEParticle(
        GCAPitchAngleFriction(),
        GCAPitchAngleDiffusion(),
        ic,
        drift_params,
        drift_params,
        params.npart
        )
    return particles
end

#sol = tp_run!(exp);
# Callback
condition(u) = true
function affect!(u, i, n)
    if u[i+1, 5, n] < -1.0
        u[i+1, 5, n] = -1.0
    elseif u[i+1, 5, n] > 1.0
        u[i+1, 5, n] = 1
    end
    if u[i+1, 5, n] < -1.0 
        println("still to low")
    elseif u[i+1,5,n] > 1
        println("still to high")
    end
end
callback = DiscreteTPCallback(condition, affect!)

# SOLVE
function run()
    @time particles = init()
    @time u, t = solve(particles, dt, tspan, euler_maruyama; callback=callback,
        seed=params.seed)
    save_fast_PAS(u, params)
    return u, t
end
u, t = run()
function sde_plot()
    fig, ax = plt.subplots(3,2)
    for i = 1:params.npart
        ax[1,1].plot(u[:,1,i], u[:,3,i])
        ax[2,1].plot(t, u[:,1,i])
        ax[3,1].plot(t, u[:,2,i])
        ax[1,2].plot(t, u[:,3,i])
        ax[2,2].plot(t, u[:,4,i])
        ax[3,2].plot(t, u[:,5,i])
    end
    ax[1,1].set_title("z-plane")
    ax[2,1].set_title("x-pos")
    ax[3,1].set_title("y-pos")
    ax[1,2].set_title("z-pos")
    ax[2,2].set_title("vparal-pos")
    ax[3,2].set_title("beta-pos")
    fig.tight_layout()
end

# Save the results
#tp.tp_save_fast(sol.u, exp);
#nits = Array{Int64}(undef, npart)
#for i = 1:npart
#    nits[i] = length(sol.u[i].t)
#end
#avgits = sum(nits)/npart
#npart_over100 = length(findall(x -> x >100, nits))
