#-------------------------------------------------------------------------------
# Created 02.12.22
# email: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Meshes.jl
#
#-------------------------------------------------------------------------------
# Module containing mesh-structs and methods
#-------------------------------------------------------------------------------




#-------------#
# Structs     #
#-------------#-----------------------------------------------------------------
struct Mesh
    bField ::Array{T, 4} where {T<:Real} # The magnetic field
    eField ::Array{T, 4} where {T<:Real} # The eletric field
    ∇B     ::Array{T, 4} where {T<:Real} # The gradient of the magnetic field
    ∇b̂     ::Array{T, 5} where {T<:Real} # The gradient of the magnetic field
    ∇ExB   ::Array{T, 5} where {T<:Real} # The gradient of the magnetic field
    xCoords::Vector{T} where {T<:Real} # The cartesian x-coordinates of the grid points
    yCoords::Vector{T} where {T<:Real} # The cartesian x-coordinates of the grid points
    zCoords::Vector{T} where {T<:Real} # The cartesian x-coordinates of the grid points
    domain ::Matrix{T} where {T<:Real} # Contains the extent of the numerical domain
    numdims::Integer

    # Constructors
    #--------------------------------------------------------------------------
    function Mesh(bField ::Array{T, 4} where {T<:Real},
                  eField ::Array{T, 4} where {T<:Real},
                  ∇B     ::Array{T, 4} where {T<:Real},
                  ∇b̂     ::Array{T, 5} where {T<:Real},
                  ∇ExB   ::Array{T, 5} where {T<:Real},
                  xCoords::Vector{T} where {T<:Real},
                  yCoords::Vector{T} where {T<:Real},
                  zCoords::Vector{T} where {T<:Real}
                  )
        domain = [xCoords[1] xCoords[end]
                  yCoords[1] yCoords[end]
                  zCoords[1] zCoords[end]]
        numdims = 3
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{T, 4} where {T<:Real},
                  eField ::Array{T, 4} where {T<:Real},
                  xCoords::Vector{T} where {T<:Real},
                  yCoords::Vector{T} where {T<:Real},
                  zCoords::Vector{T} where {T<:Real}
                  )
        domain = [xCoords[1] xCoords[length(xCoords)]
                  yCoords[1] yCoords[length(yCoords)]
                  zCoords[1] zCoords[length(zCoords)]]
        numdims = 3
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 xCoords,
                                 yCoords,
                                 zCoords,
                                 )
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{T, 4} where {T<:Real},
                  eField ::Array{T, 4} where {T<:Real}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        zCoords = collect(LinRange(0,1, size(bField)[4]))
        domain = [xCoords[1] xCoords[end]
                  yCoords[1] yCoords[end]
                  zCoords[1] zCoords[end]]
#        ∇B = compute∇B(bField, 
#                       xCoords,
#                       yCoords,
#                       zCoords)
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 xCoords,
                                 yCoords,
                                 zCoords)
        numdims = 3
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{T, 3} where {T<:Real},
                  eField ::Array{T, 3} where {T<:Real}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        dx = xCoords[2] - xCoords[1]
        dy = yCoords[2] - yCoords[1]
        #∇B = ∇(bField, dx, dy, derivate4thOrder) # Won't work. need norm3
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField, # Won't work. Need norm3 etc.
                                 xCoords,
                                 yCoords,
                                 zCoords)
        numdims = 2
        # Won't work. missing ∇b̂ and ∇ExB
        return new(bField, eField, ∇B, xCoords, yCoords, numdims)
    end # constructor 

    function Mesh(bField ::Array{T, 2} where {T<:Real},
                  eField ::Array{T, 2} where {T<:Real}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        dx = xCoords[2] - xCoords[1]
        ∇B = derivateCentral(bField, dx, (1,0,0)) # Won't work. need norm2
        numdims = 1
        # Won't work. missing ∇b̂ and ∇ExB
        return new(bField, eField, ∇B, xCoords, numdims)
    end # constructor 

    """
         Using Bifrost input
    """
    function Mesh(
        expname::String,
        expdir ::String,
        snap   ::Integer,
        ;
        SI_units::Bool=true,
        wfp      ::DataType=Float32 # Working float precision
        )

        numdims = 3
        #-----------------------------------------------------------------------
        # Load Bifrost-snap

        # Parse filenames
        snapstr = "_$snap"
        basename = string(expdir, "/", expname, "_$snap")
        idl_filename = string(basename, ".idl")
        snap_filename = string(basename, ".snap")
        aux_filename = string(basename, ".aux")
        mesh_filename = string(expdir, "/", expname, ".mesh")

        br_mesh = BifrostMesh(mesh_filename)
        br_params = br_read_params(idl_filename)
        br_snap = br_load_snapdata(snap_filename, br_params)
        br_aux = br_load_auxdata(aux_filename, br_params)

        auxvars = split(br_params["aux"])
        aux_avail = [0,                # η_total
                     0, 0, 0,  # ηx, ηy, ηz
                     0, 0, 0,  # Jx, Jy, Jz
                     0, 0, 0,  # Ex, Ey, Ez
                     ]
        for i = 1:length(auxvars)
            if auxvars[i] == "eta_total"
                aux_avail[1] = i
            elseif auxvars[i] == "etax"
                aux_avail[3] = i
            elseif auxvars[i] == "etay"
                aux_avail[4] = i
            elseif auxvars[i] == "etaz"
                aux_avail[2] = i
            elseif auxvars[i] == "ix"
                aux_avail[5] = i
            elseif auxvars[i] == "iy"
                aux_avail[6] = i
            elseif auxvars[i] == "iz"
                aux_avail[7] = i
            elseif auxvars[i] == "ex"
                aux_avail[8] = i
            elseif auxvars[i] == "ey"
                aux_avail[9] = i
            elseif auxvars[i] == "ez"
                aux_avail[10] = i
            end
        end

        #-----------------------------------------------------------------------
        # Allocate memory for simple variables
        meshsize = (br_mesh.mx, br_mesh.my, br_mesh.mz)

        # Density
        ρ_code  = br_snap[:,:,:,1]
        # Bulk momentum
        px_code = br_snap[:,:,:,2]
        py_code = br_snap[:,:,:,3]
        pz_code = br_snap[:,:,:,4]
        # Magnetic field
        bx_code = br_snap[:,:,:,6]
        by_code = br_snap[:,:,:,7]
        bz_code = br_snap[:,:,:,8]

        pbc_x = Bool(br_params["periodic_x"])
        pbc_y = Bool(br_params["periodic_y"])
        pbc_z = Bool(br_params["periodic_z"])

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

        # De-stagger and scale velocity
        ux = code2cgs_u * br_xup(px_code, pbc_x) ./ ρ_code
        uy = code2cgs_u * br_yup(py_code, pbc_y) ./ ρ_code
        uz = code2cgs_u * br_zup(pz_code, pbc_z) ./ ρ_code

        # De-stagger and scale magnetic field
        bx = code2cgs_b * br_xup(bx_code, pbc_x)
        by = code2cgs_b * br_xup(by_code, pbc_y)
        bz = code2cgs_b * br_xup(bz_code, pbc_z)
        
        # Scale axes
        x = code2cgs_l * br_mesh.x
        y = code2cgs_l * br_mesh.y
        z = code2cgs_l * br_mesh.z

        if SI_units
            ux *= wfp(Constants.cgs2SI_u)
            uy *= wfp(Constants.cgs2SI_u)
            uz *= wfp(Constants.cgs2SI_u)
            bx *= wfp(Constants.cgs2SI_b)
            by *= wfp(Constants.cgs2SI_b)
            bz *= wfp(Constants.cgs2SI_b)
            x  *= wfp(Constants.cgs2SI_l)
            y  *= wfp(Constants.cgs2SI_l)
            z  *= wfp(Constants.cgs2SI_l)
        end

        dx = diff(x)
        dy = diff(y)
        dz = diff(z)

        domain = [x[1] x[end]
                  y[1] y[end]
                  z[1] z[end]]

        #-----------------------------------------------------------------------
        # Compute the electric field and current density
        calcEfield = true
        resistiveMHD = true

        if (aux_avail[8] > 0) & (aux_avail[9] > 0) & (aux_avail[10] > 0)
            ex = br_aux[:,:,:,aux_avail[8]]
            ey = br_aux[:,:,:,aux_avail[9]]
            ez = br_aux[:,:,:,aux_avail[10]]
            calcEfield = false
        elseif  (aux_avail[2] > 0) & (aux_avail[3] > 0) & (aux_avail[4] > 0)
            # Calculate η_total? or J using components?
            println(string("Error: Computing E-field using η-components is not",
                           " implemented"))
            return
        elseif aux_avail[1] > 0
            η = br_aux[:,:,:,aux_avail[1]]
            # Do nothing
        else
            @warn string("No resisstance nor electric field in aux-vars.",
                         "Proceeding assuming ideal MHD")
            resistiveMHD = false
        end

        # Bring magnetic field to cell centres
        bField  = zeros(wfp, 3, meshsize...)
        bField[1,:,:,:] = bx
        bField[2,:,:,:] = by
        bField[3,:,:,:] = bz

        if calcEfield
            # This allocation is not needed if cross() allows components and
            # not only the vectors.
            bulkvel = zeros(wfp, 3, meshsize...)
            bulkvel[1,:,:,:] = ux
            bulkvel[2,:,:,:] = uy
            bulkvel[3,:,:,:] = uz
            if resistiveMHD
                # Calculate current density
                if (aux_avail[5] > 0) & (aux_avail[6] > 0) & (aux_avail[7] > 0)
                    # Current density available in aux-variables
                    J = zeros(wfp, 3, meshsize...)
                    J[1,:,:,:] = br_aux[:,:,:,aux_avail[5]]
                    J[2,:,:,:] = br_aux[:,:,:,aux_avail[6]]
                    J[3,:,:,:] = br_aux[:,:,:,aux_avail[7]]
                elseif aux_avail[1] > 0
                    # Calculate J using the curl of the magnetic field.
                    # Use magnetic field on cell faces and upwind scheme
                    J = Schemes.curl(bx, by, bz,
                                     br_mesh.x, br_mesh.y, br_mesh.z,
                                     Schemes.derivateUpwind
                                     )/Constants.μ_0
                end
                eField = J - cross(bulkvel, bField)
                eField[1,:,:,:] .*= η
                eField[2,:,:,:] .*= η
                eField[3,:,:,:] .*= η
            else
                eField = -cross(bulkvel, bField)
            end
        end

        # Scale electric field. No need to de-stagger. Should be cell centred
        if SI_units
            ex *= Constants.cgs2SI_e
            ey *= Constants.cgs2SI_e
            ez *= Constants.cgs2SI_e
        end
        eField  = zeros(wfp, 3, meshsize...)
        eField[1,:,:,:] = code2cgs_e * ex
        eField[2,:,:,:] = code2cgs_e * ey
        eField[3,:,:,:] = code2cgs_e * ez


        #-----------------------------------------------------------------------
        # Define/compute other parameters/variables.

        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 x,
                                 y,
                                 z,
                                 derivateUpwind,
                                 )

        #-------------------------------------------------------

        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   x, y, z,
                   domain, numdims)
    end # constructor Mesh

end # Struc

#-------------------#
# Utility functions #
#-------------------------------------------------------------------------------
function compute∇B(bField ::Array{T, 4} where {T<:Real},
                   xCoords::Vector{T} where {T<:Real},
                   yCoords::Vector{T} where {T<:Real},
                   zCoords::Vector{T} where {T<:Real}
                   )
    dx = xCoords[2] - xCoords[1]
    dy = yCoords[2] - yCoords[1]
    dz = zCoords[2] - zCoords[1]
    ∇B = ∇(norm4(bField), dx, dy, dz, derivate4thOrder) 
    return ∇B
end # function compute∇B


function compute∇s(bField ::Array{T, 4} where {T<:Real},
                   eField ::Array{T, 4} where {T<:Real},
                   xCoords::Vector{T} where {T<:Real},
                   yCoords::Vector{T} where {T<:Real},
                   zCoords::Vector{T} where {T<:Real}
                   )
    wfp = typeof(bField[1])
    _, nx, ny, nz = size(bField)
    dx = xCoords[2] - xCoords[1]
    dy = yCoords[2] - yCoords[1]
    dz = zCoords[2] - zCoords[1]
    BB = norm4(bField)
    b̂ = zeros(wfp, 3, nx, ny, nz)
    ExBdrift = zeros(wfp, 3, nx, ny, nz)
    for i = 1:nx
        for j= 1:ny
            for k = 1:nz
                B⃗ = bField[:,i,j,k]
                E⃗ = eField[:,i,j,k]
                B = BB[i,j,k]
                b̂[:,i,j,k]  .= B⃗ ./ B
                ExBdrift[:,i,j,k] .= (E⃗ × B⃗) ./ B^2
            end
        end
    end
    ∇B = ∇(BB, dx, dy, dz, derivate4thOrder) 
    ∇b̂ = ∇(b̂,  dx, dy, dz, derivate4thOrder)
    ∇ExBdrift= ∇(ExBdrift, dx, dy, dz, derivate4thOrder)
    return ∇B, ∇b̂, ∇ExBdrift
end # function compute∇s
#|
function compute∇s(bField ::Array{T, 4} where {T<:Real},
                   eField ::Array{T, 4} where {T<:Real},
                   xCoords::Vector{T} where {T<:Real},
                   yCoords::Vector{T} where {T<:Real},
                   zCoords::Vector{T} where {T<:Real},
                   scheme ::Function
                   )
    wfp = typeof(bField[1])
    _, nx, ny, nz = size(bField)
    BB = norm4(bField)
    b̂ = zeros(wfp, 3, nx, ny, nz)
    ExBdrift = zeros(wfp, 3, nx, ny, nz)
    for i = 1:nx
        for j= 1:ny
            for k = 1:nz
                B⃗ = bField[:,i,j,k]
                E⃗ = eField[:,i,j,k]
                B = BB[i,j,k]
                b̂[:,i,j,k]  .= B⃗ ./ B
                ExBdrift[:,i,j,k] .= (E⃗ × B⃗) ./ B^2
            end
        end
    end
    ∇B = ∇(BB, xCoords, yCoords, zCoords, scheme) 
    ∇b̂ = ∇(b̂,  xCoords, yCoords, zCoords, scheme)
    ∇ExBdrift= ∇(ExBdrift, xCoords, yCoords, zCoords, scheme)
    return ∇B, ∇b̂, ∇ExBdrift
end # function compute∇s
#------------------#
# Mesh set-methods #
#-------------------------------------------------------------------------------
function amplifyBfield!(mesh  ::Mesh,
                        factor::Real)
    @. mesh.bField *= factor
    mesh.∇B = compute∇B(mesh.bField, 
                        mesh.xCoords,
                        mesh.yCoords,
                        mesh.zCoords)
end # function amplifyBfield
                   

function amplifyEfield!(mesh  ::Mesh,
                        factor::Real)
    @. mesh.eField *= factor
end # function amplifyEfield

