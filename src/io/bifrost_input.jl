#-------------------------------------------------------------------------------
# Created 03.01.24 
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 io/bifrost_input.jl
#
#-------------------------------------------------------------------------------
# Contains functions for reading Bifrost snapshots and creating interpoaltion
# objects from the snapshot-fields. Depends on BifrostTools.jl.
#-------------------------------------------------------------------------------


function get_br_var_interpolator(
        expname ::String,
        snap    ::Integer,
        expdir  ::String,
        variable::String,
        ;
        itp_type=Gridded(Linear()),
        itp_bc=Flat(),
        kwargs...
    )
    var = get_var(expname, snap, expdir, variable; kwargs...)
    var = dropdims(var)
    br_axes = load_br_axes(expname, snap, expdir)
    br_axes = dropdims(br_axes)
    var_itp = linear_interpolation(
        br_axes, var, extrapolation_bc=itp_bc
        )
    #var_itp = interpolate(br_axes, var, itp_type)
    #var_itp = extrapolate(var_itp, itp_bc)
    return var_itp
end


function get_br_emfield_interpolator(
        expname::String,
        snap   ::Integer,
        expdir ::String,
        ;
        itp_type=Gridded(Linear()),
        itp_bc=Flat(),
        )
    bx = get_var(expname, snap, expdir, "bx"; units="SI", destagger=true)
    by = get_var(expname, snap, expdir, "by"; units="SI", destagger=true)
    bz = get_var(expname, snap, expdir, "bz"; units="SI", destagger=true)
    ex = get_var(expname, snap, expdir, "ex"; units="SI", destagger=true)
    ey = get_var(expname, snap, expdir, "ey"; units="SI", destagger=true)
    ez = get_var(expname, snap, expdir, "ez"; units="SI", destagger=true)
    #
    emfield = eachslice(stack([bx, by, bz, ex, ey, ez]), dims=(1,2,3))
   # nx,ny,nz = size(bx)
   # emfield = [[bx[i,j,k],by[i,j,k],bz[i,j,k],ex[i,j,k],ey[i,j,k],ez[i,j,k]] for i=1:nx, j = 1:ny, k = 1:nz]
    #
    # Make interpolation-object
    emfield = dropdims(emfield)
    br_axes = load_br_axes(expname, snap, expdir)
    br_axes = dropdims(br_axes)
    emfields_itp = linear_interpolation(
        br_axes, emfield, extrapolation_bc=itp_bc
        )
    #emfields_itp = interpolate(br_axes, emfield, itp_type)
    #emfields_itp = extrapolate(emfields_itp, itp_bc)
    #
    return emfields_itp
end


function get_br_emfield_numdensity_gastemp_interpolator(
        expname::String,
        snap   ::Integer,
        expdir ::String,
        ;
        itp_type=Gridded(Linear()),
        itp_bc=Flat(),
        )
    #
    bx = get_var(expname, snap, expdir, "bx"; units="SI", destagger=true)
    by = get_var(expname, snap, expdir, "by"; units="SI", destagger=true)
    bz = get_var(expname, snap, expdir, "bz"; units="SI", destagger=true)
    ex = get_var(expname, snap, expdir, "ex"; units="SI", destagger=true)
    ey = get_var(expname, snap, expdir, "ey"; units="SI", destagger=true)
    ez = get_var(expname, snap, expdir, "ez"; units="SI", destagger=true)
    rho= get_var(expname, snap, expdir, "r"; units="SI")
    tg = get_var(expname, snap, expdir, "tg"; units="SI", destagger=false)
    #
    # Convert mass density to number density
    rho ./= tp.m_p
    fields = eachslice(
        [bx;;;; by;;;; bz;;;; ex;;;; ey;;;; ez;;;; rho;;;; tg],
        dims=(1,2,3)
        )       
    #
    # Make interpolation-object
    fields = dropdims(fields)
    br_axes = load_br_axes(expname, snap, expdir)
    br_axes = dropdims(br_axes)
    fields_itp = interpolate(br_axes, fields, itp_type)
    fields_itp = extrapolate(fields_itp, itp_bc)
    #
    return fields_itp
end


function load_br_axes(
        expname::String,
        snap   ::Integer,
        expdir ::String,
        )
    basename = string(
        expdir, "/", expname
        )
    idl_filename = string(basename, "_", snap, ".idl")
    mesh_filename = string(basename, ".mesh")
    mesh = BifrostMesh(mesh_filename)
    params = read_params(idl_filename)
    # Scale br_axes
    code2cgs_l = parse(Float32, params["u_l"]) # exploit promotion
    x = code2cgs_l * tp.cgs2SI_l * mesh.x
    y = code2cgs_l * tp.cgs2SI_l * mesh.y
    z = code2cgs_l * tp.cgs2SI_l * mesh.z
    return (x, y, z)
end
