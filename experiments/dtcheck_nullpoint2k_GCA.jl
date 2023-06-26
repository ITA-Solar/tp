using TraceParticles
using LinearAlgebra

# Set simulation parameters
params = Parameters(
    dt=1e-8,
    nsteps=Int64(1e6),
    npart=5,
    tp_expname="nullpoint2k_dtrun",
    tp_expdir="/mn/stornext/u3/eilifo/code/repos/tp/data/dtrun/gca",
    solver="full-orbit",
    scheme="RK4",
    interp="bilinear_xz",
    br_expname="nullpoint2k",
    br_expdir="/mn/stornext/u3/eilifo/data/ohf/nullpoint2k",
    br_isnap=700,
    wp_part=Float64,
    wp_snap=Float32,
    pos_distr="point",
    posxbounds=[15.547e6, 15.739e6, 15.846e6, 15.947e6, 16.226e6],
    posybounds=[1e6, 1e6, 1e6, 1e6, 1e6],
    poszbounds=[-6.709e6, -6.519e6, -6.410e6, -6.314e6, -6.048e6],
    vel_distr="mb",
    pbc = (false, true, false),
    specie=[1,1,1,1,1],
)


tf = 0.01
ndt = 10
dts = 10 .^ (LinRange(-4,-10,ndt))
stepslist = round.(Integer, tf ./ dts)

exp = tp_init!(params)
for i = 1:ndt
    params.dt = dts[i]
    params.nsteps = stepslist[i]

    # Initialise experiment
    tp_reinit_particles!(exp)

    # Run experiment
    @time tp_run!(exp)

    basename = string(params.tp_expdir, "/", params.tp_expname)
    tp_savetp(exp, string(basename, "_$i.tp"))
    tp_saveparams(exp, string(basename, "_params_$i.jl"))
end
    
