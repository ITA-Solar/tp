using tp

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
# Creating the vector potential also gives the axes of the experiment
n = (100, 100, 2)
domainaxes, gridsizes, A = normal3Donlyz(xi0, xif, n, μ, σ)
# Derive the magnetic field from the curl of the vector-potential
Bfield = curl(A, gridsizes, derivateCentral)
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


function plasmoidParameters(n, dt)
        # Set simulation parameters
        params = Parameters(
        dt=dt,
        nsteps=round(Int, time/dt),
        npart=4,
        tp_expname="plasmoid$(n[1])",
        tp_expdir="/mn/stornext/u3/eilifo/code/repos/tp",
        solver="full-orbit",
        scheme="RK4",
        interp="trilinear",
        wp_part=Float32,
        wp_snap=Float64,
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

function plasmoidParametersDE(n)
        # Set simulation parameters
        params = Parameters(
                patch_type="DE",
                tspan=(0.0, time),
                charge = 1.0,
                mass = 1.0,
                npart=4,
                tp_expname="plasmoid$(n[1])",
                tp_expdir="/mn/stornext/u3/eilifo/code/repos/tp",
                solver="full-orbit",
                wp_part=Float32,
                wp_snap=Float64,
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
# DE-run
q, m = (1.0, 1.)
npart = 4
pmesh = PureMesh(xx, yy, zz)
p  = FOParams(q, m, pmesh, Bfield, Efield)
foic = FOStaticIC(posd, veld)
part = ODEParticle(LorentzForce(), foic, p, 4)
tspan = (0.0, time)
DEpatch = DEPatch(tspan, part, pmesh)
println("Running DEpatch")
@time sim = update(DEpatch)

#  DE GCA-run
gradB, gradb, gradExB = compute_gradients(Bfield, Efield,
                                          pmesh.x,pmesh.y,pmesh.z,
                                          derivateUpwind
                                          )
B_itp, E_itp = EMfield_itps(pmesh, Bfield, Efield)
gcaB = Array{Float64}(undef, 3, npart)
gcaE = Array{Float64}(undef, 3, npart)
for i = 1:npart
        x,y,z = posd[:,i]
        gcaB[:,i] = [ B_itp[1](x,y,z),
                      B_itp[2](x,y,z),
                      B_itp[3](x,y,z) ]
        gcaE[:,i] = [ E_itp[1](x,y,z),
                      E_itp[2](x,y,z),
                      E_itp[3](x,y,z) ]
end     
R, vparal, mu = calc_GCA_IC_and_mu(posd, veld, q, m, gcaB, gcaE)
pGCA = GCAParams(q, m, mu, pmesh, Bfield, Efield, gradB, gradb, gradExB)
icGCA = GCAStaticIC(R, vparal)
gcapart = ODEParticle(GCA(), icGCA, pGCA, npart)
GCApatch = DEPatch(tspan, gcapart, pmesh)
println("Running GCApatch")
@time gcasim = update(GCApatch)

# Parameters DE FO run
FOparams = plasmoidParametersDE(n)
FOparams.solver = LorentzForce()
FOexp = DE_init!(FOparams)
println("Running FO exp")
@time sim2 = tp_run!(FOexp)

# Parameters DE GCA run
GCAparams = plasmoidParametersDE(n)
GCAparams.solver = GCA()
GCAexp = DE_init!(GCAparams)
println("Running GCA exp")
@time gcasim2 = tp_run!(GCAexp)

# Legacy-run
println("Running Legacy-patch")
dt = 0.001
params = plasmoidParameters(n, dt)
# Initialise experiment
exp = tp_init!(params)
# Run experiment
@time tp_run!(exp)

