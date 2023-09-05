#-------------------------------------------------------------------------------
# Created 05.06.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                TraceParticles.jl
#
#-------------------------------------------------------------------------------
# The tp (trace particles) package.
#-------------------------------------------------------------------------------

methodmap = Dict(
    "full-orbit"      => fullOrbit,
    "full-orbit-isf"  => fullOrbit_interstaticfield,
    "full-orbit-rel"  => relFullOrbitExplLeapFrog,
    "GCA"             => gca,
    "RK4"             => rk4,
    "trilinear"       => trilinear,
    "trilinear-GCA"   => trilinearGCA,
    "bilinear_xz"     => bilinear_xz,
    "bilinear_xz-GCA" => bilinear_xzGCA,
)

gcasolvers = (
              "GCA",
             )

#--------------------#
# Struct definitions #
#--------------------#----------------------------------------------------------
mutable struct Parameters
    #
    npart ::Integer
    nsteps::Integer
    dt    ::Real
    # Working precision on snap (mesh and fields) and part (trace particles)
    wp_snap::DataType
    wp_part::DataType
    # Solver, numerical scheme and interpolation scheme
    solver::String
    scheme::String
    interp::String
    # Location for saving simulation-data
    tp_expname::String
    tp_expdir ::String
    # Bifrost-snapshot
    br_expname::String
    br_expdir ::String
    br_isnap  ::Integer
    # Manual mesh and EM-fields
    x ::Vector{T} where {T<:Real}
    y ::Vector{T} where {T<:Real}
    z ::Vector{T} where {T<:Real}
    bx::Array{T, 3} where {T<:Real}
    by::Array{T, 3} where {T<:Real}
    bz::Array{T, 3} where {T<:Real}
    ex::Array{T, 3} where {T<:Real}
    ey::Array{T, 3} where {T<:Real}
    ez::Array{T, 3} where {T<:Real}
    # Number of grid points in the x, y, and z-direction
    nx::Integer
    ny::Integer
    nz::Integer
    # Particle-initialisation
    specie   ::Vector{T} where {T<:Integer}
    pos_distr::String # Particle position distribution "uniform": uniform
    vel_distr::String # Particle velocity distribution    ↑ or "mb": MB
    posxbounds::Vector{T} where {T<:Real} 
    posybounds::Vector{T} where {T<:Real} 
    poszbounds::Vector{T} where {T<:Real} 
    velxbounds::Vector{T} where {T<:Real} 
    velybounds::Vector{T} where {T<:Real} 
    velzbounds::Vector{T} where {T<:Real} 
    T         ::Real
    seed      ::Integer
    # Boundary conditions
    pbc::Tuple{Bool, Bool, Bool} # (x,y,z) If not periodic boundary conditions,
                                 # particles are
                                 # killed when crossing all boundaries.
    # "Private" parameters
    bg_input::String # "br": Bifrost snapshot
                     # "manual": Mesh and fields are defined manually in Julia.

    # Constructors
    #--------------------------------------------------------------------------
    # Bifrost-snap
    function Parameters(
        ;
        # Required parameters
        #-----------------------------------------|
        npart ::Integer,
        nsteps::Integer,
        dt    ::Real,
        # Working precision
        wp_snap::DataType,
        wp_part::DataType,
        #
        # Optionals
        #-----------------------------------------|
        solver::String="full-orbit",
        scheme::String="RK4",
        interp::String="trilinear",
        tp_expname=nothing,
        tp_expdir=nothing,
        # Either Bifrost-snapshot
        br_expname=nothing,
        br_expdir =nothing,
        br_isnap  =nothing,
        # or manual mesh and fields
        x =nothing,
        y =nothing,
        z =nothing,
        bx=nothing,
        by=nothing,
        bz=nothing,
        ex=nothing,
        ey=nothing,
        ez=nothing,
        # Either way, gridsize will be useful
        nx=nothing,
        ny=nothing,
        nz=nothing,
        # Particle-initialisation
        specie=ones(Integer, npart),
        pos_distr::String="uniform",
        vel_distr::String="point",
        # Default value of bounds depend on the mesh which may be given by
        # a Bifrost snapshot.
        posxbounds=nothing,
        posybounds=nothing,
        poszbounds=nothing,
        velxbounds=zeros(wp_part, 2),
        velybounds=zeros(wp_part, 2),
        velzbounds=zeros(wp_part, 2),
        T         =nothing,
        seed      =0,
        # Periodic boundary conditions
        pbc::Tuple{Bool, Bool, Bool}=(false, false, false),
        # "Private" parameters
        #-----------------------------------------|
        bg_input::String="br"
        ) 
        #
        println("tp.jl: Constructing Parameters...")
        #
        #  Create parameter instance with all parameters undefined
        params = new()

        #-----------------------------------------------------------------------
        # Check if mesh and EM-field paramters are present.
        # Bifrost snapshot is preferred by default.
        if all(i -> i != nothing, (br_expname, br_expdir, br_isnap))
            params.bg_input = "br"
            params.br_expname = br_expname
            params.br_expdir = br_expdir
            params.br_isnap = br_isnap
        elseif all(i -> i != nothing, (x, y, z, bx, by, bz, ex, ey, ez))
            params.bg_input = "manual"
            params.x = x
            params.y = y
            params.z = z
            params.nx = Int32(length(x))
            params.ny = Int32(length(y))
            params.nz = Int32(length(z))
            params.bx = bx
            params.by = by
            params.bz = bz
            params.ey = ex
            params.ex = ey
            params.ez = ez
        else
            error("Missing parameter(s): mesh or/and EM-field parameters")
        end 

        #-----------------------------------------------------------------------
        # Check if optional paramters without default are available
        if tp_expname != nothing
            params.tp_expname = tp_expname
        end
        if tp_expdir != nothing
            params.tp_expdir = tp_expdir
        end
        if posxbounds != nothing
            params.posxbounds = posxbounds
        end
        if posybounds != nothing
            params.posybounds = posybounds
        end
        if poszbounds != nothing
            params.poszbounds = poszbounds
        end
        if nx != nothing
            params.nx = nx
        end
        if ny != nothing
            params.ny = ny
        end
        if nz != nothing
            params.nz = nz
        end
        if T != nothing
            params.T = T
        end

        #-----------------------------------------------------------------------
        # Set required parameters
        params.npart = npart
        params.nsteps = nsteps
        params.dt = dt
        params.wp_snap = wp_snap
        params.wp_part = wp_part
        # Set optional parameters
        params.solver = solver
        params.scheme = scheme
        params.interp = interp
        params.specie = specie
        params.pos_distr = pos_distr
        params.vel_distr = vel_distr
        params.seed = seed
        params.pbc = pbc
        params.velxbounds = velxbounds
        params.velybounds = velybounds
        params.velzbounds = velzbounds
        # Set default value to optional parameters that are not passed and not
        # already defined above.

        # Construct methodmap

        # Return Parameters instance
        return params
    end # Constructor Parameters 
end # struct Parameters


mutable struct Experiment
    params::Parameters
    patch ::Patch
end


#------------------------------#
# Initialisation of Experiment #
#------------------------------#------------------------------------------------
function tp_init!(
    params::Parameters
    )
    # Check whether the required parameters are present. Since the struct is
    # mutable, it may have changed since constrcution.
    checkrequirements(params)


    println("tp.jl: Initialising experiment...")

    #---------------------------------------------------------------------------
    # Construct Experiment
    #---------------------------------------------------------------------------
    mesh = tp_createmesh!(params)
    particles = tp_createparticles!(params, mesh)
    patch = tp_createpatch(params, mesh, particles)

    exp = Experiment(params,
                     patch
                     )
    return exp
end # function tp_init


function tp_createmesh!(
    params::Parameters
    )
    #---------------------------------------------------------------------------
    # Construct mesh
    #---------------------------------------------------------------------------
    if params.bg_input == "br"
        mesh = Mesh(params.br_expname, params.br_expdir, params.br_isnap;
                    wfp=params.wp_snap)
        params.nx = length(mesh.xCoords)
        params.ny = length(mesh.yCoords)
        params.nz = length(mesh.zCoords)
    else params.bg_input == "manual"
        meshsize = (params.nx, params.ny, params.nz)
        bField = zeros(params.wp_snap, 3, meshsize...)
        eField = zeros(params.wp_snap, 3, meshsize...)
        bField[1,:,:,:] = params.bx
        bField[2,:,:,:] = params.by
        bField[3,:,:,:] = params.bz
        eField[1,:,:,:] = params.ex
        eField[2,:,:,:] = params.ey
        eField[3,:,:,:] = params.ez
        mesh = Mesh(bField, eField, params.x, params.y, params.z)
    end 
    return mesh
end


function tp_createparticles!(
    params::Parameters,
    mesh  ::Mesh,
    )
    # Define variables to avoid magic numbers
    numdims = 3
    #---------------------------------------------------------------------------
    # Construct particles
    #---------------------------------------------------------------------------
    #
    # First, create positions and velocities according to given parameters
    #
    # Go thorugh the various intial distributions
    # Positions
    if params.pos_distr == "uniform"
        # Set default bounds if not defined
        # Default bounds are the mesh domain boundaries.
        if !isdefined(params, :posxbounds)
            params.posxbounds = [mesh.xCoords[1], mesh.xCoords[end]]
        end
        if !isdefined(params, :posybounds)
            params.posybounds = [mesh.yCoords[1], mesh.yCoords[end]]
        end
        if !isdefined(params, :poszbounds)
            params.poszbounds = [mesh.zCoords[1], mesh.zCoords[end]]
        end
        pos = inituniform(params.npart,
                                    params.posxbounds,
                                    params.posybounds,
                                    params.poszbounds,
                                    params.wp_part
                                    ;
                                    params.seed
                                    )
    elseif params.pos_distr == "point"
        if !isdefined(params, :posxbounds)
            params.posxbounds = [(mesh.xCoords[1] - mesh.xCoords[end])/2.0]
        end
        if !isdefined(params, :posybounds)
            params.posybounds = [(mesh.yCoords[1] - mesh.yCoords[end])/2.0]
        end
        if !isdefined(params, :poszbounds)
            params.poszbounds = [(mesh.zCoords[1] - mesh.zCoords[end])/2.0]
        end
        pos = ones(params.wp_part, numdims, params.npart)
        # If particles are given by coordinates
        if (length(params.posxbounds) ==
            length(params.posybounds) ==
            length(params.poszbounds) == params.npart)
            pos[1,:] = params.posxbounds
            pos[2,:] = params.posybounds
            pos[3,:] = params.poszbounds
        else
            pos[1,:] *= params.posxbounds[1]
            pos[2,:] *= params.posybounds[1]
            pos[3,:] *= params.poszbounds[1]
        end
    else
        error("Invalid parameter value: pos_distr")
    end
    # IMPLEMENT e.g
    #
    # position distributed according to density
    #

    # Velocities
    if params.vel_distr == "mb"
        # Get temperature 
        if params.bg_input == "br"
            vel = ones(numdims, params.npart)
            # Temperature has code-units equal to Kelvin, no need for scaling.
            br_temp = dropdims(br_load_auxvariable(params.br_expname,
                                          [params.br_isnap],
                                          params.br_expdir,
                                          "tg",
                                          params.wp_snap
                                                   ),
                               dims=numdims + 1
                               )
            for i = 1:params.npart
                # Interpolate the temperature at the position of the particle
                x⃗ = pos[:,i]
                if params.interp == "bilinear_xz-GCA"
                    interp = methodmap["bilinear_xz"]
                elseif params.interp == "trilinear-GCA"
                    interp = methodmap["trilinear"]
                else
                    interp = methodmap[params.interp]
                end 
                t, _ = gridinterp(br_temp,
                               interp,
                               x⃗,
                               mesh.xCoords,
                               mesh.yCoords,
                               mesh.zCoords,
                               )
                # Standard deviation of velocity
                σ = √(k_B*t/specieTable[params.specie[i], 1])
                # Expectance value of particle velocity components is zero
                μ = 0.0
                # Draw velocity components from normal-distribution
                vel[:,i] = randn(μ, σ, (numdims); seed=params.seed+i)
            end
        else
            error(string("\"mb\" velocity distribution only available with ",
                         "Bifrost background. Use \"mb-onetemp\" to use one ",
                         "temperature for all particles"
                         )
                  )
        end
    elseif params.vel_distr == "mb-onetemp"
        if !isdefined(params, :T)
            error("Parameter not defined: T")
        end
        vel = ones(numdims, params.npart)
        for i = 1:params.npart
            # Standard deviation of velocity
            σ = √(k_B*params.T/specieTable[params.specie[i], 1])
            # Expectance value of particle velocity components is zero
            μ = 0.0
            # Draw velocity components from normal-distribution
            vel[:,i] = randn(μ, σ, (numdims); seed=params.seed+i)
        end
    elseif params.vel_distr == "uniform"
        vel = inituniform(params.npart,
                                    params.velxbounds,
                                    params.velybounds,
                                    params.velzbounds,
                                    params.wp_part
                                    ;
                                    seed=params.seed
                                    )
    elseif params.vel_distr == "point"
        vel = ones(params.wp_part, numdims, params.npart)
        if (length(params.velxbounds) ==
            length(params.velybounds) ==
            length(params.velzbounds) == params.npart)
            vel[1,:] = params.velxbounds
            vel[2,:] = params.velybounds
            vel[3,:] = params.velzbounds
        else
            vel = ones(numdims, params.npart)
            vel[2,:] *= params.velybounds[1]
            vel[3,:] *= params.velzbounds[1]
        end
    else
        error("Invalid parameter value: vel_distr")
    end

    #
    # Next, create the actual particle instances based on the solver (full orbit
    # or GCA).
    #
    # Full orbit
    if params.solver == "full-orbit"
        particles = ParticleSoA(pos,
                                vel,
                                params.specie,
                                params.nsteps)
    elseif params.solver == "full-orbit-isf"
        particles = ParticleSoA(pos,
                                vel,
                                params.specie,
                                params.nsteps)
    elseif params.solver == "full-orbit-relativistic"
        particles = ParticleSoA(pos,
                                vel,
                                params.specie,
                                params.nsteps)
    # GCA
    elseif params.solver == "GCA"
        vparal = zeros(params.wp_part, params.npart)
        magneticMoment = zeros(params.wp_part, params.npart)
        for i = 1:params.npart
            (B⃗, E⃗), _ = gridinterp(mesh,
                                   methodmap[params.interp],
                                   pos[:,i]
                                   )
            B = norm(B⃗)
            b̂ = B⃗/B
            v = norm(vel[:,i])
            vparal[i] = vel[:,i] ⋅ b̂
            vperp = √(v^2 - vparal[i]^2)
            mass = specieTable[params.specie[i], 1]
            magneticMoment[i] = mass*vperp^2/(2B)
        end
        particles = GCAParticleSoA(pos,
                                   vparal,
                                   magneticMoment,
                                   params.specie,
                                   params.nsteps
                                   )
    else
        error("Invalid solver-parameter: $params.solver")
    end
    return particles
end


function tp_createpatch(
    params   ::Parameters,
    mesh     ::Mesh,
    particles::TraceParticle
    )
    #---------------------------------------------------------------------------
    # Construct Patch
    #---------------------------------------------------------------------------
    patch = Patch(mesh,
                  particles,
                  methodmap[params.solver],
                  methodmap[params.scheme],
                  methodmap[params.interp],
                  params.dt,
                  params.nsteps,
                  params.npart,
                  params.pbc
                  )
    return patch
end


function tp_softinit!(
    exp::Experiment
    )
end


function tp_reinit_particles!(
    exp   ::Experiment,
    )
    particles = tp_createparticles!(exp.params, exp.patch.mesh)
    patch = tp_createpatch(exp.params, exp.patch.mesh, particles)
    exp.patch = patch
end


#----------------------------------------#
# Saving and loading Experiment or files #
#----------------------------------------#--------------------------------------
function tp_save(
    exp    ::Experiment
    ;
    expname=nothing,
    expdir=nothing,
    )
    #
    # Check whether path for saving files is defined.
    #
    if expname == nothing
        if !isdefined(exp.params, :tp_expname)
            error(string("Experiment name `tp_expname` needed for saving.",
                         " Specify it in `Parameters` or give as argument to",
                         " `tp_saveExp`."))
        else
            expname = exp.params.tp_expname
        end
    end
    if expdir == nothing
        if !isdefined(exp.params, :tp_expdir)
            error(string("Experiment directory `tp_expdir` needed for saving.",
                         " Specify it in `Parameters` or give as argument to",
                         " `tp_saveExp`."))
        else
            expdir = exp.params.tp_expdir
        end
    end

    # To avoid magic numbers
    numdims = 3

    println("tp.jl: Saving experiment...")
    #
    # Construct filenames
    #
    basename = string(expdir, "/", expname)
    tp_filename = string(basename, ".tp")
    mesh_filename = string(basename, ".mesh")
    bg_filename = string(basename, ".bg")
    params_filename = string(basename, "_params", ".jl")
    # SUGGESTION: Add some parameter-file? Not sure yet what to include and how
    # to make it. It would be nice with a params-file or script with same name
    # that could be run to produce the corresponding .tp file.

    #
    # Write particle positions, velocity and life-status
    #
    tp_savetp(exp, tp_filename)
    
    #
    # Write mesh (NB! here I seperate the fields from the mesh)
    #
    tp_savemesh(exp, mesh_filename)

    #
    # Write background fields
    #
    tp_savebg(exp, bg_filename)

    #
    # Write parameters
    #
    tp_saveparams(exp, params_filename)

    #
end # function tp_saveExp


function tp_saveparams(
    exp     ::Experiment,
    filename::String
    )
    paramsstring = "using tp\nparams = Parameters(\n"
    for p in fieldnames(Parameters)
        if isdefined(exp.params, p)
            if typeof(getfield(exp.params, p)) == String
                value = "\"$(getfield(exp.params, p))\""
            else
                value = "$(getfield(exp.params, p))"
            end
            spaces = " "^(11 - length("$p"))
            paramsstring = string(paramsstring,
                                  "\t$p", spaces, "= $value,\n")
        end
    end
    paramsstring = string(paramsstring, ")")
    #
    f = open(filename, "w+")
    write(f, paramsstring)
    close(f)
    #
    println("tp.jl: Wrote $(filename)")
end

function tp_savetp(
    exp::Experiment,
    filename::String,
    )
    f = open(filename, "w+")
    if usingGCA(exp.params)
        write(f, exp.patch.tp.R)
        write(f, exp.patch.tp.vparal)
    else
        write(f, exp.patch.tp.pos)
        write(f, exp.patch.tp.vel)
    end
    close(f)
    println("tp.jl: Wrote $filename")
end # tp_savetp


function tp_savemesh(
    exp::Experiment,
    filename::String,
    )
    f = open(filename, "w+")
    write(f, exp.patch.mesh.xCoords)
    write(f, exp.patch.mesh.yCoords)
    write(f, exp.patch.mesh.zCoords)
   # write(f, [exp.patch.mesh.xCoords[:];
   #           exp.patch.mesh.yCoords[:];
   #           exp.patch.mesh.zCoords[:];
   #           ]
   #       )
    close(f)
    println("tp.jl: Wrote $filename")
end # function tp_savemesh


function tp_savebg(
    exp::Experiment,
    filename::String,
    )
    f = open(filename, "w+")
    write(f, exp.patch.mesh.bField)
    write(f, exp.patch.mesh.eField)
    write(f, exp.patch.mesh.∇B)
    write(f, exp.patch.mesh.∇b̂)
    write(f, exp.patch.mesh.∇ExB)
    close(f)
    println("tp.jl: Wrote $filename")
end # function tp_savebg


function tp_load(
    params ::Parameters,
    ;
    expdir ::String=params.tp_expdir,
    expname::String=params.tp_expname,
    mesh_filename=nothing,
    bg_filename  =nothing,
    )
    #
    println("tp.jl: Loading experiment from file:")
    # Parse filenames
    basename = string(expdir, "/", expname)
    tp_filename = string(basename, ".tp")
    if mesh_filename == nothing
        mesh_filename = string(basename, ".mesh")
    end
    if bg_filename == nothing
        bg_filename = string(basename, ".bg")
    end

    #
    # Open tp-file
    #
    pos, vel, μ = tp_loadtp(params, tp_filename)
    #
    # Open mesh-file
    #
    x, y, z = tp_loadmesh(params, mesh_filename)
    #
    # Open bg-file
    #
    bField, eField, ∇B, ∇b̂, ∇ExB = tp_loadbg(params, bg_filename) 
    #
    # Create Mesh, ParticleSoA, Patch and Experiment
    #
    mesh = Mesh(bField, eField, ∇B, ∇b̂, ∇ExB, x, y, z)
    if usingGCA(params)
        particles = GCAParticleSoA(pos, vel, μ, params.specie)
    else
        particles = ParticleSoA(pos, vel, params.specie)
    end
    patch = Patch(mesh,
                  particles,
                  methodmap[params.solver],
                  methodmap[params.scheme],
                  methodmap[params.interp],
                  params.dt,
                  params.nsteps,
                  params.npart,
                  params.pbc
                  )
    exp = Experiment(params, patch);
    #
    return exp
end # function tp_load


function tp_loadtp(
    params::Parameters,
    filename::String
    )
    numdims = 3
    pos = Array{params.wp_part}(undef, numdims, params.npart, params.nsteps+1)
    #
    println("tp.jl: Loading particles...")
    f = open(filename)
    if usingGCA(params)
        vel = Array{params.wp_part}(undef, params.npart, params.nsteps+1)
        μ = Vector{params.wp_part}(undef, params.npart)
        read!(f, pos)
        read!(f, vel)
        read!(f, μ)
        close(f)
        return pos, vel, μ
    else
        vel = Array{params.wp_part}(undef, numdims, params.npart,
            params.nsteps+1)
        read!(f, pos)
        read!(f, vel)
        close(f)
        return pos, vel, nothing
    end 
end


function tp_loadtp!(
    exp::Experiment,
    filename::String
    )
    pos, vel, μ = tp_loadtp(exp.params, filename)
    if usingGCA(exp.params)
        exp.patch.tp.R .= pos
        exp.patch.tp.vparal .= vel
        exp.patch.tp.μ .= μ
    else
        exp.patch.tp.pos .= pos
        exp.patch.tp.vel .= vel
    end 
end
    
    
function tp_loadmesh(
    params::Parameters,
    filename::String
    )
    x = Vector{params.wp_snap}(undef, params.nx)
    y = Vector{params.wp_snap}(undef, params.ny)
    z = Vector{params.wp_snap}(undef, params.nz)
    #
    println("tp.jl: Loading mesh...")
    f = open(filename)
    read!(f, x)
    read!(f, y)
    read!(f, z)
    close(f)
    #
    return x, y, z
end
    

function tp_loadbg(
    params::Parameters,
    filename::String
    )
    meshsize = (params.nx, params.ny, params.nz)
    bField = zeros(params.wp_snap, 3, meshsize...)
    eField = zeros(params.wp_snap, 3, meshsize...)
    ∇B = zeros(params.wp_snap, 3, meshsize...)
    ∇b̂ = zeros(params.wp_snap, 3, 3, meshsize...)
    ∇ExB = zeros(params.wp_snap, 3, 3, meshsize...)
    #
    println("tp.jl: Loading fields...")
    f = open(filename)
    read!(f, bField)
    read!(f, eField)
    read!(f, ∇B)
    read!(f, ∇b̂)
    read!(f, ∇ExB)
    close(f)
    #
    return bField, eField, ∇B, ∇b̂, ∇ExB
end



#---------------------#
# Running Experiments #
#---------------------#---------------------------------------------------------
function tp_run!(
    exp::Experiment
    )
    statement = string("tp.jl: Running simulation:\n",
                       "\tnpart:  $(exp.params.npart)\n",
                       "\tnsteps: $(@sprintf("%.2e", exp.params.nsteps))\n",
                       "\tdt:     $(exp.params.dt)\n",
                       "\tNumber of iterations: ",
                       "$(@sprintf("%.2e",exp.params.npart*exp.params.nsteps))"
                       )
    println(statement)
    run!(exp.patch)
end # function tp_run


function tp_run(
    params::Parameters
    )
    exp = tp_init!(params)
    tp_run!(exp)
    return exp
end # function tp_run


#-----------------------------------------#
# Set and/or reset Parameters/Experiments #
#-----------------------------------------#-------------------------------------
function tp_set_dt!(
    exp::Experiment,
    dt ::Real,
    )
    exp.params.dt = dt
    exp.patch.dt = dt
end


function tp_set_nsteps!(
    exp::Experiment,
    nsteps::Integer,
    )
    exp.params.nsteps = nsteps
    tp_reinit_particles!(exp)
end


function tp_reset!(
    exp::Experiment,
    )
    revive!(exp.patch.tp)
    reset!(exp.patch.tp)
end


#-------------------#
# Utility functions #
#-------------------#-----------------------------------------------------------
function checkrequirements(
    params::Parameters
    )
    if params.npart == nothing
        error("Missing parameter: npart")
    elseif params.nsteps == nothing
        error("Missing parameter: nsteps")
    elseif params.dt == nothing
        error("Missing parameter: dt")
    elseif params.solver == nothing
        error("Missing parameter: solver")
    elseif params.scheme == nothing
        error("Missing parameter: scheme")
    elseif params.interp == nothing
        error("Missing parameter: interp")
    end
    if all(i -> isdefined(params, i), (:br_expname,
                                       :br_expdir,
                                       :br_isnap
                                       )
           )
    elseif all(i -> isdefined(params, i), (:x, :y, :z,
                                           :bx, :by, :bz,
                                           :ex, :ey, :ez
                                           )
               )
    else
        error("Missing parameter(s): mesh or/and EM-field")
    end
    if params.bg_input == nothing
        error("Missing paramter: bg_input")
    end
    if params.npart > length(params.specie)
        error("Number of particles outnumbers the length of specie-parameter.")
    end
end


function usingGCA(params::Parameters)
    if params.solver in gcasolvers
        return true
    else
        return false
    end
end

#-------------------------------------------------------------------------------
