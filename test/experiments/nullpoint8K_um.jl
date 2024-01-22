# Nullpoint simulation with 
#  8192x8192x1 grid points (approx. 7km res)
#  Uniform distribution in position 
#  Maxwellian velocity distributions
using tp


# Some usefull conversion factors
code2cgs_t = 1e2
code2cgs_l = 1e8
code2SI_t = code2cgs_t*tp.cgs2SI_t
code2SI_l = code2cgs_l*tp.cgs2SI_l

# Let tspan be (0.0, dtsnap)
dtsnap = 1e-2
dtsnapSI = code2SI_t*dtsnap
#
tspan = (0.0, dtsnapSI)
#

# choose position(s)
#
npart = convert(Int64, 3e7) # Number of particles
#
# Initial point positions
#
# Point in right exhaust
#posx = [16.98, 16.4]*code2SI_l
#posy = [1.0,1.0]*code2SI_l
#posz = [-5.71,-5.710]*code2SI_l
#
# Slicex4k and slicez4k
# Approximately the same as slice4k ends.
posx = [14.527344, 18.824219]*code2SI_l
posy = [1.0,1.0]*code2SI_l
posz = [-7.78315,-3.8759463]*code2SI_l
#

params = Parameters(
    patch_type="DE",
    tspan=tspan,
    charge = -tp.e,
    mass = tp.m_e,
    npart=npart,
    tp_expname="nullpoint8k_30M_um",
    tp_expdir="/mn/stornext/d9/data/eilifo/tp/nullpoint/8K-ext_damp",
    br_expname="nullpoint",
    br_expdir="/mn/stornext/d9/data/eilifo/bifrost/ohf/run0/8K-ext_damp",
    br_isnap=1490,
    solver=GCA(),
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
    seed=0,
)
 
exp = DE_init!(params);
sol = tp_run!(exp);

# Save the results
tp.tp_save_fast(sol.u, exp);
nits = Array{Int64}(undef, npart)
for i = 1:npart
    nits[i] = length(sol.u[i].t)
end
avgits = sum(nits)/npart
npart_over100 = length(findall(x -> x >100, nits))
npart_over1000 = length(findall(x -> x >1000, nits))
npart_over10000 = length(findall(x -> x >10000, nits))
npart_over90000 = length(findall(x -> x >90000, nits))
file = open(params.tp_expdir*"/stats.tp", "w+")
write(file, avgits)
write(file, npart_over100)
write(file, npart_over1000)
write(file, npart_over10000)
write(file, npart_over90000)
close(file)
file = open(params.tp_expdir*"/nits", "w+")
write(file, nits)
close(file)
