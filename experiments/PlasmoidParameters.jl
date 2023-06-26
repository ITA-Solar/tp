module PlasmoidParameters

using TraceParticles
using Schemes
using Utilities

export plasmoidParameters

#...............................................
# INITIAL CONDITIONS
# Placement span
pos0 = [0.10, 0.0, 0.0]
posf = [0.9, 0.4, 0.0]
# Velocity span (for uniform velcity)
vel0 = [0.5, 0.0, 0.0]
velf = [0.5, 0.0, 0.0]
# For manually determined positions and velocities
posd = [0.13  0.23  0.265 0.20 # At n = (100,100,2)
        0.03  0.05  0.12  0.04
        0.00  0.00  0.00  0.00]
#posd = [0.13  0.20  0.25  0.20  # At n = (10,10,2)
#        0.03  0.05  0.12  0.04
#        0.00  0.00  0.00  0.00]
veld = [0.50  0.50  0.50  1.00
        0.00  0.00  0.00  0.00
        0.00  0.00  0.00  0.00]
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
time = 2.9

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
        npart=4,
        tp_expname="plasmoid$(n[1])",
        tp_expdir="C:/Users/ixyva/data/test",
        solver="full-orbit",
        scheme="RK4",
        interp="trilinear",
        wp_part=Float32,
        wp_snap=Float32,
        pos_distr="point",
        vel_distr="point",
        posxbounds=posd[1,:],
        posybounds=posd[2,:],
        poszbounds=posd[3,:],
        velxbounds=veld[1,:],
        velybounds=veld[2,:],
        velzbounds=veld[3,:],
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
        specie=4*ones(Int64, 4),
        pbc=(false, false, true)
        )
        return params
end

println(abspath(PROGRAM_FILE) == @__FILE__)

if abspath(PROGRAM_FILE) == @__FILE__
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
        plot(exp.patch)
        plt.show()
end

end

