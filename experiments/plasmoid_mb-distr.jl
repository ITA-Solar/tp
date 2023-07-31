module PlasmoidParameters

using TraceParticles
using Schemes
using Utilities
using Constants: k_B

export plasmoidParameters

#...............................................
# INITIAL CONDITIONS
#  number of particles
npart = 1000
# Mass of particles
m = 1
# Average velocity of particles
v_avg = 0.5
# Temperature of particle ensamble
T = v_avg^2*π*m/(8k_B)
#...............................................
# SPATIAL PARAMETERS (x, y, z)
# Lower bounds of the three spatial axes
xi0 = (0., 0., 0.)
# Upper bound of the three spatial axes
xif = (1., 1., 1.)

#...............................................
# MAGNETIC FIELD PARAMETERS
bamp   = 10.0            # Amplification factor
bconst = [100., 0., 0]   # Constant term
bz     = 0.0             # z-component
μx = (xif[1] - xi0[1])/2 # 
μy = (xif[2] - xi0[2])/2 
σx = (xif[1] - xi0[1])/8
σy = (xif[2] - xi0[2])/8
μ  = (μx, μy)
σ  = (σx, σy)
time = 2.0

function plasmoidParameters(n, dt)
        # Creating the vector potential also gives the axes of the experiment
        domainaxes, gridsizes, A = Utilities.normal3Donlyz(xi0, xif, n, μ, σ)
        # Derive the magnetic field from the curl of the vector-potential
        Bfield = Schemes.curl(A, gridsizes, Schemes.derivateCentral)
        # Scale the field
        @. Bfield = bamp*Bfield + bconst
        @. Bfield[3,:,:,:] = bz
        # Extract the vector components
        Bx = Bfield[1,:,:,:]
        By = Bfield[2,:,:,:]
        Bz = Bfield[3,:,:,:]
        # ELECTRIC FIELD
        Ex = 0.0
        Ey = 0.0
        Ez = 50.0
        Efield = zeros(Float64, size(Bfield))
        Efield[1, :, :, :] .= Ex
        Efield[2, :, :, :] .= Ey
        Efield[3, :, :, :] .= Ez
        Ex = Efield[1,:,:,:]
        Ey = Efield[2,:,:,:]
        Ez = Efield[3,:,:,:]

        xx, yy, zz = domainaxes

        # Set simulation parameters
        params = Parameters(
        dt=dt,
        nsteps=round(Int, time/dt),
        npart=npart,
        tp_expname="plasmoid$(n[1])",
        tp_expdir="C:/Users/ixyva/data/test",
        solver="full-orbit",
        scheme="RK4",
        interp="trilinear",
        wp_part=Float32,
        wp_snap=Float64,
        pos_distr="uniform",
        vel_distr="mb-onetemp",
        posxbounds=[xi0[1], xif[1]],
        posybounds=[xi0[2], xif[2]],
        poszbounds=[xi0[3], xif[3]],
        x=xx,
        y=yy,
        z=zz,
        bx=Bx,
        by=By,
        bz=Bz,
        ex=Ex,
        ey=Ey,
        ez=Ez,
        nx=n[1],
        ny=n[2],
        nz=n[3],
        bg_input="manual",
        specie=4*ones(Int64, npart),
        pbc=(false, false, true),
        T=T,
        seed=0,
        )
        return params
end

println(abspath(PROGRAM_FILE) == @__FILE__)

using TPplots
import PyPlot
const plt = PyPlot
plt.pygui(true)

n = (100, 100, 2)
dt = 0.001

params = plasmoidParameters(n, dt)
# Initialise experiment
exp = tp_init!(params)
# Run experiment
tp_run!(exp)
nbins = 50
TPplots.plotenergydistr(PlasmoidParameters.exp.patch.tp, 1, 
        nbins, "Initial velocity distribution"; log=true)
TPplots.plotenergydistr(PlasmoidParameters.exp.patch.tp, params.nsteps+1, 
        nbins, "Final velocity distribution"; log=true)
plt.show()

end

