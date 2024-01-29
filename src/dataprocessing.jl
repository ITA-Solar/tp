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


