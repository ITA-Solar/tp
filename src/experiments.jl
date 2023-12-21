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

gca_solvers = (
              "GCA",
              GCA,
              GCA()
             )
fo_solvers = (
    "full-orbit",
    "full-orbit-isf",
    "full-orbit-rel",
    LorentzForce,
    LorentzForce()
)

#--------------------#
# Struct definitions #
#--------------------#----------------------------------------------------------

mutable struct Parameters
    #
    patch_type::String
    npart ::Integer
    # Working precision on snap (mesh and fields) and part (trace particles)
    wp_snap::DataType
    wp_part::DataType
    # Required for standard patch
    nsteps::Integer
    dt    ::Real
    # Required for DE patch
    tspan ::Tuple{<:Real, <:Real}
    charge::Any
    mass  ::Any
    # Solver, numerical scheme and interpolation scheme
    solver::Any
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
    x ::Vector{<:Real}
    y ::Vector{<:Real}
    z ::Vector{<:Real}
    bx::Array{<:Real, 3}
    by::Array{<:Real, 3}
    bz::Array{<:Real, 3}
    ex::Array{<:Real, 3}
    ey::Array{<:Real, 3}
    ez::Array{<:Real, 3}
    # Number of grid points in the x, y, and z-direction
    nx::Integer
    ny::Integer
    nz::Integer
    # Particle-initialisation
    specie   ::Vector{<:Integer}
    pos_distr::String # Particle position distribution "uniform": uniform
    vel_distr::String # Particle velocity distribution    ↑ or "mb": MB
    posxbounds::Vector{<:Real}
    posybounds::Vector{<:Real}
    poszbounds::Vector{<:Real}
    velxbounds::Vector{<:Real}
    velybounds::Vector{<:Real}
    velzbounds::Vector{<:Real}
    T         ::Real
    seed      ::Integer
    # Boundary conditions
    pbc::Tuple{Bool, Bool, Bool} # (x,y,z) If not periodic boundary conditions,
                                 # particles are
                                 # killed when crossing all boundaries.
    # "Private" parameters
    bg_input::String # "br": Bifrost snapshot
                     # "manual": Mesh and fields are defined manually in Julia.
    SI_units::Bool

    # Constructors
    #--------------------------------------------------------------------------
    # Bifrost-snap
    function Parameters(
        ;
        # Required parameters
        #-----------------------------------------|
        npart ::Integer,
        # Working precision
        wp_snap::DataType,
        wp_part::DataType,
        #
        # Optionals
        #-----------------------------------------|
        # (required for standard patch)
        nsteps=nothing,
        dt    =nothing,
        # (required for DE patch)
        tspan = nothing,
        charge= nothing,
        mass  = nothing,
        #
        patch_type::String="STD",
        solver="full-orbit",
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
        SI_units::Bool=true,
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
        if tp_expname !== nothing
            params.tp_expname = tp_expname
        end
        if tp_expdir !== nothing
            params.tp_expdir = tp_expdir
        end
        if posxbounds !== nothing
            params.posxbounds = posxbounds
        end
        if posybounds !== nothing
            params.posybounds = posybounds
        end
        if poszbounds !== nothing
            params.poszbounds = poszbounds
        end
        if nx !== nothing
            params.nx = nx
        end
        if ny !== nothing
            params.ny = ny
        end
        if nz !== nothing
            params.nz = nz
        end
        if T !== nothing
            params.T = T
        end

        #-----------------------------------------------------------------------
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
        params.SI_units = SI_units
        # Set default value to optional parameters that are not passed and not
        # already defined above.
        
        # Set required parameters
        params.npart = npart
        if patch_type == "STD"
            params.patch_type = patch_type
            try
                params.dt = dt
                params.nsteps = nsteps
            catch
                error("Missing parameters: dt and/or nsteps")
            end
        #-----------------------------------------------------------------------
        # Additional parameters for DEPatch
        elseif patch_type == "DE"
            params.patch_type = patch_type
            try
                params.tspan = tspan
                params.charge = charge
                params.mass = mass
            catch
                error("Missing parameter(s). Check: tspan, charge, mass")
            end
        else
            error("Unrecognised Patch type")
        end

        # Return Parameters instance
        return params
    end # Constructor Parameters 
end # struct Parameters


mutable struct Experiment
    params::Parameters
    patch ::AbstractPatch
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------#
# Initialisation of Experiment #
#------------------------------#------------------------------------------------
function DE_init!(
    params::Parameters
    )
    # Check whether the required parameters are present. Since the struct is
    # mutable, it may have changed since constrcution.
    checkrequirements(params)


    println("tp.jl: Initialising experiment...")

    #---------------------------------------------------------------------------
    # Get mesh, fields, initial conditions and problem parameters from the 
    # experiment parameters
    mesh = create_puremesh(params)
    bfield, efield = get_fields(params) 
    println("tp.jl: Drawing initial posistions and velocities...")
    pos, vel = tp_get_initial_pos_and_vel!(params, mesh)
    if usingGCA(params)
        println("tp.jl: Getting GCA initial conditions...")
        ICs, mu = get_GCA_IC(params, pos, vel, bfield, efield, mesh)
        println("tp.jl: Computing field-gradients...")
        gradB, gradb, gradExB = compute_gradients(
                                    bfield, 
                                    efield,
                                    mesh.x, mesh.y, mesh.z,
                                    derivateUpwind
                                    )
        println("tp.jl: Getting GCA problem parameters...")
        prob_params = GCAParams(
            params.charge,
            params.mass,
            mu,
            mesh,
            bfield,
            efield,
            gradB,
            gradb,
            gradExB
        )
    else
        println("tp.jl: Getting FO initial conditions...")
        ICs = get_FO_IC(pos, vel) 
        println("tp.jl: Getting FO problem parameters...")
        prob_params = FOParams(params.charge, params.mass, mesh, bfield, efield)
    end

    #---------------------------------------------------------------------------
    # Construct instances of the ODEparticle (ensamble of particles), the 
    # the instances of the patch and the experiment.
    particle = ODEParticle(
                        params.solver, 
                        ICs, 
                        prob_params, 
                        params.npart
                        )
    patch = DEPatch(params.tspan, particle, mesh)
    exp = Experiment(params, patch)
    println("tp.jl: Initialisation of DE-experiment completed.")
    return exp
end


function tp_init!(
    params::Parameters
    )
    # Check whether the required parameters are present. Since the struct is
    # mutable, it may have changed since constrcution.
    checkrequirements(params)


    println("tp.jl: Initialising experiment...")

    #---------------------------------------------------------------------------
    # Construct Experiment
    mesh = tp_createmesh!(params)
    pos, vel = tp_get_initial_pos_and_vel!(params, mesh)
    particles = tp_createparticles!(params, pos, vel)
    patch = tp_createpatch(params, mesh, particles)

    exp = Experiment(params,
                     patch
                     )
    return exp
end # function tp_init


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function create_puremesh(params::Parameters)
    if params.bg_input == "br"
        return get_mesh_from_br(params)
    elseif params.bg_input == "manual"
        return PureMesh(params.x, params.y, params.z)
    else
        error("Invalid parameter: bg_input")
    end
end # function create_puremesh


function get_fields(params::Parameters)
    if params.bg_input == "br"
        return get_fields_from_br(params)
    elseif params.bg_input == "manual"
        B = [params.bx;;;; params.by;;;; params.bz]
        E = [params.ex;;;; params.ey;;;; params.ez]
        B = permutedims(B, (4, 1, 2, 3))
        E = permutedims(E, (4, 1, 2, 3))
        return B, E
    else
        error("Invalid parameter: bg_input")
    end
end

function EMfield_at_pos(
    pos::Matrix{<:Real},
    Bfield::Array{<:Real, 4},
    Efield::Array{<:Real, 4},
    pmesh ::PureMesh
    )
    npart = size(pos)[2]
    Bitp, Eitp = EMfield_itps(pmesh, Bfield, Efield)
    B_at_pos = Array{Float64}(undef, 3, npart)
    E_at_pos = Array{Float64}(undef, 3, npart)
    for i = 1:npart
        x,y,z = pos[:,i]
        B_at_pos[:,i] = Bitp(x,y,z)
        E_at_pos[:,i] = Eitp(x,y,z)
    end     
    return B_at_pos, E_at_pos
end

function get_GCA_IC(
    params::Parameters,
    pos::Matrix{<:Real},
    vel::Matrix{<:Real},
    Bfield::Array{<:Real, 4},
    Efield::Array{<:Real, 4},
    pmesh ::PureMesh
    )
    B_at_pos, E_at_pos = EMfield_at_pos(pos, Bfield, Efield, pmesh)    
    R, vparal, mu = calc_GCA_IC_and_mu(
        pos, 
        vel, 
        params.charge, 
        params.mass,
        B_at_pos,
        E_at_pos,
    )
    return GCA_IC(
        VectorIC(R[1,:]), 
        VectorIC(R[2,:]), 
        VectorIC(R[3,:]), 
        VectorIC(vparal)
        ), 
        mu
end

function get_FO_IC(
    pos::Matrix{<:Real},
    vel::Matrix{<:Real},
)
    return FO_IC(
        VectorIC(pos[1,:]),
        VectorIC(pos[2,:]),
        VectorIC(pos[3,:]),
        VectorIC(vel[1,:]),
        VectorIC(vel[2,:]),
        VectorIC(vel[3,:]),
        )
end

function get_mesh_from_br(
    params::Parameters
    )    
    println("tp.jl: Reading Bifrost mesh...")
    basename = string(
        params.br_expdir, "/", params.br_expname
        )
    idl_filename = string(basename, "_", params.br_isnap, ".idl")
    mesh_filename = string(basename, ".mesh")
    br_mesh = BifrostMesh(mesh_filename)
    br_params = br_read_params(idl_filename)
    # Scale axes
    code2cgs_l = params.wp_snap(br_params["u_l"])
    x = code2cgs_l * tp.cgs2SI_l * br_mesh.x
    y = code2cgs_l * tp.cgs2SI_l * br_mesh.y
    z = code2cgs_l * tp.cgs2SI_l * br_mesh.z
    if length(y) == 1
        y = [y[1],y[1]]
    end
    return PureMesh(x, y, z)
end

function get_fields_from_br(
    params::Parameters
    ;
    get_bfield  ::Bool=true, 
    get_efield  ::Bool=true,
    get_density ::Bool=false,
    get_gas_temp::Bool=false,
    )
    wfp = params.wp_snap
    ndims = 3

    file_basename = string(
        params.br_expdir, "/", params.br_expname
        )
    idl_filename = string(file_basename, "_", params.br_isnap, ".idl")
    snap_filename = string(file_basename, "_", params.br_isnap, ".snap")
    aux_filename = string(file_basename, "_", params.br_isnap, ".aux")

    println("tp.jl: Reading Bifrost fields...")
    println("           Reading parameters from $(basename(idl_filename))")
    br_params = br_read_params(idl_filename)
    println("           Loading snap       from $(basename(snap_filename))")
    br_snap = br_load_snapdata(snap_filename, br_params)
    println("           Loading auxiliares from $(basename(aux_filename))")
    br_aux = br_load_auxdata(aux_filename, br_params)

    pbc_x = Bool(br_params["periodic_x"])
    pbc_y = Bool(br_params["periodic_y"])
    pbc_z = Bool(br_params["periodic_z"])
    auxes = split(br_params["aux"])

    # Check which fields to get. This is not the best way to do it, because 
    # the amount of variables to return and which ones are unknown
    # What I should do is use br_load_snapvariable and have 1 function per
            # variable. Loading the parameters and mesh does not take 
            # that much time and you don't have to load all variables into
            # memory.
    if get_bfield
        # Magnetic field
        println("           Getting magnetic field from index...")
        bx_code = br_snap[:,:,:,6]
        by_code = br_snap[:,:,:,7]
        bz_code = br_snap[:,:,:,8]
    end
    if get_efield
        # Electric field
        ex_idx = findall(x -> x == "ex", auxes)
        ey_idx = findall(x -> x == "ey", auxes)
        ez_idx = findall(x -> x == "ez", auxes)
        e_idxs = [length(ex_idx), length(ey_idx), length(ez_idx)]
        if any(i -> i == 0, e_idxs)
            error("Did not find all electric field components in aux-file")
        end
        println("           Getting electric field from index...")
        ex_code = br_aux[:,:,:,ex_idx...]
        ey_code = br_aux[:,:,:,ey_idx...]
        ez_code = br_aux[:,:,:,ez_idx...]
    end
    if get_density
        println("           Getting density from index...")
        r_code = br_snap[:,:,:,1]
    end
    if get_gas_temp
        println("           Getting temperature from index...")
        tg_idx = findall(x -> x == "tg", auxes)
        tg_code = br_aux[:,:,:,tg_idx...]
    end


    #-----------------------------------------------------------------------
    # Destagger and scale variables to SI-units
    #
    # Vector quiantities are defined at the faces, while scalar quantities
    # are defined at cell centres. To get the bulk velcity at the centre
    # one thus have to interpolate the momentum to the cell centres.
    #
    # br_xup moves values half a grid in the x-direction.
    # params["u_u"] scales velocity from model/code-units to CGS units. 
    # params["u_b"] scales magnetic field from model/code-units to CGS units. 
    # cgs2SI_u scales velocity from CGS-units to SI-units
    # cgs2SI_b scales magnetic field from CGS-units to SI-units
    code2cgs_u = wfp(br_params["u_u"])
    code2cgs_b = wfp(br_params["u_B"])
    code2cgs_l = wfp(br_params["u_l"])
    code2cgs_e = code2cgs_u * code2cgs_b
    code2cgs_r = wfp(br_params["u_r"])

    if get_bfield
        # De-stagger and scale magnetic field
        println("           De-staggering: bx")
        bx = br_xup(bx_code, pbc_x)
        println("                          by")
        by = br_yup(by_code, pbc_y)
        println("                          bz")
        bz = br_zup(bz_code, pbc_z)
        bx *= code2cgs_b 
        by *= code2cgs_b 
        bz *= code2cgs_b 
        if params.SI_units
            bx *= wfp(cgs2SI_b)
            by *= wfp(cgs2SI_b)
            bz *= wfp(cgs2SI_b)
        end
    end
    if get_efield
        # De-stagger and scale electric field
        println("                          ex")
        ex = br_yup(br_zup(ex_code, pbc_z), pbc_y)
        println("                          ey")
        ey = br_zup(br_xup(ey_code, pbc_x), pbc_z)
        println("                          ez")
        ez = br_xup(br_yup(ez_code, pbc_y), pbc_x)
        ex *= code2cgs_e 
        ey *= code2cgs_e 
        ez *= code2cgs_e 
        if params.SI_units
            ex *= wfp(cgs2SI_e)
            ey *= wfp(cgs2SI_e)
            ez *= wfp(cgs2SI_e)
        end
    end
    if get_density
        r = r_code*code2cgs_r
        if params.SI_units
            r *= wfp(cgs2SI_r)
        end
    end
    if get_gas_temp
        # Nothing to be done here
    end

    # Scale to SI-units
    if params.SI_units
        println("            ** SCALING TO SI **")
    else
        println("            ** SCALING TO CGS **")
    end

    # Reshape to form necessary for tp (currently)
    meshsize = [size(bx)...]
    constant_dim_idx = nothing
    if any(i -> i == 1, meshsize)
        idxs = findall(x -> x == 1, meshsize)
        if length(idxs) > 1
            error("1D snaps is not yet supported.")
        end
        constant_dim_idx = idxs[1]
        meshsize[constant_dim_idx] = 2
    end 
    println("           Allocating arrays for reshaping...")
    if get_bfield
        bfield  = Array{wfp}(undef, ndims, meshsize...)
    end
    if get_efield
        efield  = Array{wfp}(undef, ndims, meshsize...)
    end
    if get_density
        density = Array{wfp}(undef, meshsize...)
    end
    if get_gas_temp
        gas_temp = Array{wfp}(undef, meshsize...)
    end
    println("           Reshaping...")
    if constant_dim_idx == 2
        println("               Bifrost simulation is 2D in the XZ-plane.")
        for i = 1:2
            if get_bfield
                bfield[1,:,i,:] = bx
                bfield[2,:,i,:] = by
                bfield[3,:,i,:] = bz
            end
            if get_efield
                efield[1,:,i,:] = ex
                efield[2,:,i,:] = ey
                efield[3,:,i,:] = ez
            end
            if get_density
                density[:,i,:] = r
            end
            if get_gas_temp
                gas_temp[:,i,:] = tg_code
            end
        end
    elseif constant_dim_idx == 1 && constant_dim_idx == 3
        error("2D Bifrost simualtions in XY and YZ planes not yet supported")
    else
        if get_bfield
            bfield[1,:,:,:] = bx
            bfield[2,:,:,:] = by
            bfield[3,:,:,:] = bz
        end
        if get_efield
            efield[1,:,:,:] = ex
            efield[2,:,:,:] = ey
            efield[3,:,:,:] = ez
        end
        if get_density
            density = r
        end
        if get_gas_temp
            gas_temp = tg_code
        end
    end
    
    # Temporary solution. All permutations not possible.
    # See comment above.
    if get_density
        return bfield, efield, density, gas_temp
    else
        return bfield, efield
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function tp_createmesh!(
    params::Parameters
    )
    #--------------------------------------------------------------------------
    # Construct mesh
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


function tp_get_initial_pos_and_vel!(
    params::Parameters,
    mesh  ::AbstractMesh,
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
    else
        error("Invalid parameter value: pos_distr")
    end

    pos = get_initial_positions(
                        params.npart,
                        params.pos_distr,
                        params.posxbounds,
                        params.posybounds,
                        params.poszbounds,
                        params.wp_part
                        ;
                        seed=params.seed
                        )
                        
    # Velocities
    stds = Array{params.wp_part}(undef, params.npart)
    if params.vel_distr == "mb"
        # Get temperature 
        if params.bg_input == "br"
            # Temperature has code-units equal to Kelvin, no need for scaling.
            br_temp = dropdims(br_load_auxvariable(params.br_expname,
                                          [params.br_isnap],
                                          params.br_expdir,
                                          "tg",
                                          params.wp_snap
                                                   ),
                               dims=numdims + 1
                               )
            println("           Finding temperatures at initial positions")
            @time for i = 1:params.npart
                # Interpolate the temperature at the position of the particle
                x = pos[:,i]
                if params.interp == "bilinear_xz-GCA"
                    interp = methodmap["bilinear_xz"]
                elseif params.interp == "trilinear-GCA"
                    interp = methodmap["trilinear"]
                else
                    interp = methodmap[params.interp]
                end 
                if params.patch_type == "STD"
                    xx = mesh.xCoords
                    yy = mesh.yCoords
                    zz = mesh.zCoords
                else
                    xx = mesh.x
                    yy = mesh.y
                    zz = mesh.z
                end 
                t, _ = gridinterp(br_temp, interp, x, xx, yy, zz) 
                # Standard deviation of velocity
                σ = √(k_B*t/specieTable[params.specie[i], 1])
                stds[i] = σ
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
        for i = 1:params.npart
            # Standard deviation of velocity
            σ = √(k_B*params.T/specieTable[params.specie[i], 1])
            stds[i] = σ
        end
    end
    vel = get_initial_velocities(
        params.npart,
        params.vel_distr,
        params.velxbounds,
        params.velybounds,
        params.velzbounds,
        params.wp_part
        ;
        seed=params.seed,
        stds=stds
    )
    return pos, vel
end

function tp_createparticles!(
    params::Parameters,
    pos::AbstractArray,
    vel::AbstractArray
)
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
function get_basename(
    params::Parameters
    ; 
    expname=nothing, 
    expdir=nothing
    )
    #
    # Check whether path for saving files is defined.
    #
    if expname === nothing
        if !isdefined(params, :tp_expname)
            error(string("Experiment name `tp_expname` needed for saving.",
                         " Specify it in `Parameters` or give as argument to",
                         " `tp_saveExp`."))
        else
            expname = params.tp_expname
        end
    end
    if expdir === nothing
        if !isdefined(params, :tp_expdir)
            error(string("Experiment directory `tp_expdir` needed for saving.",
                         " Specify it in `Parameters` or give as argument to",
                         " `tp_saveExp`."))
        else
            expdir = params.tp_expdir
        end
    end
    return string(expdir, "/", expname)
end

function tp_save(
    exp    ::Experiment
    ;
    expname=nothing,
    expdir=nothing,
    )

    basename = get_basename(exp.params; expname=expname, expdir=expdir)
    println("tp.jl: Saving experiment...")
    #
    # Construct filenames
    #
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
    if exp.params.patch_type == "DE"
        @warn "Use tp_save(::EnsambleSolution) to save particles of a DEPatch."
    else 
        tp_savetp(exp, tp_filename)
    end
    
    #
    # Write mesh (NB! here I seperate the fields from the mesh)
    #
    tp_savemesh(exp, mesh_filename)

    #
    # Write background fields
    #
    if exp.params.patch_type == "STD"
        tp_savebg(exp, bg_filename)
    end
    #
    # Write parameters
    #
    tp_saveparams(exp, params_filename)

    #
end # function tp_saveExp

function tp_save(
    sol::EnsembleSolution,
    params::Parameters,
    ;
    expname=nothing,
    expdir=nothing,
    )
    basename = get_basename(params; expname=expname, expdir=expdir)
    tp_filename = string(basename, ".tp")

    jldopen(tp_filename, "w") do file
        write(file, "npart", params.npart)
        for i = 1:params.npart
            write(file, "u$i", sol.u[i].u)
            write(file, "t$i", sol.u[i].t)
        end
    end
    println("tp.jl: Wrote $tp_filename")
end

function tp_save(
    u::Vector{ODESolution},
    params::Parameters,
    ;
    expname=nothing,
    expdir=nothing,
    )
    basename = get_basename(params; expname=expname, expdir=expdir)
    tp_filename = string(basename, ".tp")

    npart = length(u)
    jldopen(tp_filename, "w") do file
        write(file, "npart", npart)
        for i = 1:npart
            write(file, "u$i", u[i].u)
            write(file, "t$i", u[i].t)
        end
    end
    println("tp.jl: Wrote $tp_filename")
end

#-----------------------------------------------------
function getu0(
    u   ::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    ndof = length(u[1].u[1])
    RealT = typeof(u[1].u[1][1])
    u0 = Array{RealT}(undef, ndof, npart)
    for i = 1:npart
        u0[:,i] .= u[i].u[1]
    end
    return u0
end
   
function getuf(
    u   ::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    ndof = length(u[1].u[1])
    RealT = typeof(u[1].u[1][1])
    uf = Array{RealT}(undef, ndof, npart)
    for i = 1:npart
        uf[:,i] .= u[i].u[end]
    end
    return uf
end

function getenergy(
    mu   ::Vector{<:Real},
    mass ::Any,
    B_itp::Vector{AbstractInterpolation},
    u    ::Matrix{<:Real}=getu0(u),
    )
    R = u[1:3,:]
    vparal = u[4,:]
    vperp = getvperp(u[1:3,:], mu, mass, B_itp)

    nparts = length(mu)
    RealT = typeof(mu[1])
    v2 = Array{RealT}(undef, nparts)
    E_k = Array{RealT}(undef, nparts)
    
    @. v2 = vparal^2 + vperp^2
    @. E_k = 0.5*mass*v2
    
    return E_k
end 

function getvperp(
    R     ::Matrix{<:Real},
    mu    ::Vector{<:Real},
    mass  ::Any,
    B_itp::Vector{AbstractInterpolation},
    )
    ndims, npart = size(R)
    RealT = typeof(R).parameters[1]
    vperp = Array{RealT}(undef, npart)
    if typeof(mass) <: Real
        mass = ones(RealT, npart)*mass
    end
    for i = 1:npart
        B_vec = [ B_itp[1](R[:,i]...), B_itp[2](R[:,i]...), B_itp[3](R[:,i]...) ]
        B = norm(B_vec)
        vperp[i] = √(2mu[i]*B/mass[i])
    end # loop over i: particles
    return vperp
end # function getvperp


function tp_save_fast(
    u  ::Vector{<:ODESolution},
    exp::Experiment,
    ;
    expname=nothing,
    expdir=nothing,
    )
    basename = get_basename(exp.params; expname=expname, expdir=expdir)
    tp_filename = string(basename, ".tp")
    u0 = getu0(u)
    uf = getuf(u)
    Ek0 = getenergy(
        exp.patch.tp.p.mu,
        exp.patch.tp.p.m,
        exp.patch.tp.p.B,
        u0
    )
    Ekf = getenergy(
        exp.patch.tp.p.mu,
        exp.patch.tp.p.m,
        exp.patch.tp.p.B,
        uf
    )
    ndof, npart = size(u0)
    file = open(tp_filename, "w+")
    write(file, ndof)
    write(file, npart)
    write(file, u0)
    write(file, uf)
    write(file, Ek0)
    write(file, Ekf)
    close(file)
    println("tp.jl: Wrote $tp_filename")
end

#-----------------------------------------------------
    


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
    if exp.params.patch_type == "DE"
        write(f, exp.patch.mesh.x)
        write(f, exp.patch.mesh.y)
        write(f, exp.patch.mesh.z)
    else
        write(f, exp.patch.mesh.xCoords)
        write(f, exp.patch.mesh.yCoords)
        write(f, exp.patch.mesh.zCoords)
    end
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
    if mesh_filename === nothing
        mesh_filename = string(basename, ".mesh")
    end
    if bg_filename === nothing
        bg_filename = string(basename, ".bg")
    end

    #
    # Open tp-file
    #
    if params.patch_type == "DE"
        @warn "Load DEPatch results with tp_loadsol"
    else
        pos, vel, μ = tp_loadtp(params, tp_filename)
    end
    #
    # Open mesh-file
    #
    x, y, z = tp_loadmesh(params, mesh_filename)
    #
    # Open bg-file
    #
    if params.patch_type == "STD"
        bField, eField, ∇B, ∇b̂, ∇ExB = tp_loadbg(params, bg_filename) 
    end
    #
    # Create Mesh, ParticleSoA, Patch and Experiment
    #
    if params.patch_type == "STD"
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
        return exp
    else
        return x, y, z
    end
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


function tp_loadsol(
    filename::String,
    RealT   ::DataType
    )

    # Read size of EnsambleSolution
    sizes = Array{Int64}(undef, 3)
    f = open(filename)
    read!(f, sizes)
    npart, ndof, ndsteps = sizes

    u, t = jldopen(filename, "r") do file
        npart = read(file, "npart")
        u = Vector{Vector{Vector{RealT}}}(undef, npart)
        t = Vector{Vector{RealT}}(undef, npart)
        for i = 1:npart
            u[i] = read(file, "u$i")
            t[i] = read(file, "t$i")
        end
        return u, t
    end
    return u, t
end


function tp_load_fast(
    filename::String,
    RealT   ::DataType
    )
    file = open(filename)
    ndof = read(file, Int64)
    npart = read(file, Int64)
    u0 = Array{RealT}(undef, ndof, npart)
    uf = Array{RealT}(undef, ndof, npart)
    Ek0 = Array{RealT}(undef, npart)
    Ekf = Array{RealT}(undef, npart)
    read!(file, u0)
    read!(file, uf)
    read!(file, Ek0)
    read!(file, Ekf)
    close(file)
    return ndof, npart, u0, uf, Ek0, Ekf
end
    
#---------------------#
# Running Experiments #
#---------------------#---------------------------------------------------------
function tp_run!(
    exp::Experiment
    )
    if exp.params.patch_type == "STD"
        statement = string("tp.jl: Running simulation:\n",
            "\tnpart:  $(exp.params.npart)\n",
            "\tnsteps: $(@sprintf("%.2e", exp.params.nsteps))\n",
            "\tdt:     $(exp.params.dt)\n",
            "\tNumber of iterations: ",
            "$(@sprintf("%.2e",exp.params.npart*exp.params.nsteps))"
            )
    elseif exp.params.patch_type == "DE"
        statement = string("tp.jl: Running simulation:\n",
            "\tnpart:  $(exp.params.npart)\n",
            "\ttspan: $(@sprintf("(%.2e, %.2e)", exp.params.tspan...))\n",
            )
    end
    println(statement)
    res = run!(exp.patch)
    return res
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
    end
    if params.patch_type == "STD"
        if params.npart > length(params.specie)
            error("Number of particles outnumbers the length of species.")
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
    elseif params.patch_type == "DE"
        if params.npart < length(params.charge)
            error("Number of charges outnumbers the length of particles.")
        elseif params.tspan == nothing
            error("Missing parameter: tspan")
        elseif params.charge == nothing
            error("Missing parameter: charge")
        elseif params.mass == nothing
            error("Missing parameter: mass")
        end
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
end


function usingGCA(params::Parameters)
    if params.solver in gca_solvers
        return true
    else
        return false
    end
end

#-------------------------------------------------------------------------------
