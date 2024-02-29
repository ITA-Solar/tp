#-------------------------------------------------------------------------------
# Created 29.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dataprocessing.jl
#
#-------------------------------------------------------------------------------
# Contains functions for processing results.
#-------------------------------------------------------------------------------
"""
    getinitialstate(
        u::Vector{<:ODESolution},
        )
Return the initial state (or initial condition) of a vector of `ODEsolution`'s.
"""
function getinitialstate(
    u::Vector{<:ODESolution},
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


"""
    getinitialstate(
        u::Array{<:Real, 3},:
        )
Return the initial state (or initial condition) of a solution stored as a
3D-array.
"""
function getinitialstate(
    u::Array{<:Real, 3},
    )
    return u[1,:,:]
end


"""
    getfinalstate(
        u::Vector{<:ODESolution},
        )
Return the final state of a vector of `ODEsolution`'s.
"""
function getfinalstate(
    u::Vector{<:ODESolution},
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


"""
    getfinalstate(
        u::Array{<:Real, 3},:
        )
Return the final state of a solution stored as a 3D-array.
"""
function getfinalstate(
    u::Array{<:Real, 3},
    )
    return u[end,:,:]
end


"""
    getfinaltime(
        u::Vector{<:ODESolution},
        )
Return the final times of a vector of `ODEsolution`'s.
"""
function getfinaltime(
    u::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    RealT = typeof(u[1].u[1][1])
    tf = Vector{RealT}(undef, npart)
    for i = 1:npart
        tf[i] = u.t[end]
    end
    return tf
end


function calculateenergy(
    solution::Vector{Tuple{Any, Any, Any, Any, Any}},
    expname ::String,
    snap    ::Integer,
    expdir  ::String,
    ;
    mass  ::Real=tp.m_e,
    charge::Real=tp.e,
    units ::String="SI",
    destagger::Bool=true,
    itp_bc=(Flat(), Flat())
    )
    itpvec = get_br_emfield_vecof_interpolators(
        expname, snap, expdir;
        units=units, destagger=destagger, itp_bc=itp_bc
        )
    nofparticles = length(solution)
    ekf = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    for j in eachindex(solution)
        u0 = solution[j][1]
        uf = solution[j][2]
        B0_vec = [itpvec[i](u0[1], u0[3]) for i in 1:3]
        E0_vec = [itpvec[i](u0[1], u0[3]) for i in 4:6]
        Bf_vec = [itpvec[i](uf[1], uf[3]) for i in 1:3]
        Ef_vec = [itpvec[i](uf[1], uf[3]) for i in 4:6]
        ek0[j] = kineticenergy(
            B0_vec, E0_vec, u0[1:3], u0[4], solution[j][3],
            mass, charge
            )
        ekf[j] = kineticenergy(
            Bf_vec, Ef_vec, uf[1:3], uf[4], solution[j][3],
            mass, charge
            )
    end
    return ek0, ekf
end


function calculateenergy(
    solution::Vector{Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Any}}
    ;
    mass  ::Real=tp.m_e,
    charge::Real=tp.e
    )
    nofparticles = length(solution)
    ek0 = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    ekf = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    for part in solution
        ek0 = kineticenergy(
            part[4],
            part[5],
            part[1][1:3],
            part[1][4],
            part[3],
            mass,
            charge
            )
        ekf = kineticenergy(
            part[6],
            part[7],
            part[2][1:3],
            part[2][4],
            part[3],
            mass,
            charge
            )
    end
    return ek0, ekf
end
