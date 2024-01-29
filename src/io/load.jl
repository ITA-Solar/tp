#-------------------------------------------------------------------------------
# Created 25.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 load.jl
#
#-------------------------------------------------------------------------------
# Contains functions for loading data, both generally and for data saved in 
# specific format defined by functions in save.jl
#-------------------------------------------------------------------------------
"""
    load(
        filename::String,
        ;
        kwargs...
        )
Load `filename`. Assumes the filextension is either of:

- .jld : for loading `ODESolution` objects.
- .fo  : for loading 3D solution arrays.
- .gca : for loading 3D solution arrays with magnetic moment.
- .lw  : for loading results stored using the `save_lightweight` function.
"""
function load(
    filename::String,
    ;
    kwargs...
    )
    file_ext = splitext(filename)[1]
    if file_ext == ".jld"
        load_jld(filename; kwargs...)
    elseif file_ext == ".fo"
        load_array(filename; kwargs...)
    elseif file_ext == ".gca"
        load_gca(filename; kwargs...)
    elseif file_ext == ".lw"
        load_lightweight(filename; kwargs...)
    end
end


"""
    load_jld(
        filename::String
        ;
        RealT   ::DataType
        )
Load `ODESolutions` using jld. Not sure if works. Maybe there already are some
functionality from `DifferentialEquations.jl`.
"""
function load_jld(
    filename::String
    ;
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


"""
    load_array(
        filename::String
        ;
        ndof::Integer,
        npart::Integer,
        nsteps::Integer,
        precision::DataType=Float64,
        )
For loading solutions stored in a 3D array with dimesnions `(nsteps + 1,
ndof, npart`). Double precision is assumed by default.
"""
function load_array(
    filename::String
    ;
    ndof::Integer,
    npart::Integer,
    nsteps::Integer,
    precision::DataType=Float64,
    )
    data = Array{precision, 3}(undef, nsteps+1, ndof, npart)
    f = open(filename)
    read!(f, data)
    close(f)
    return data
end


"""
    load_gca(
        filename::String
        ;
        ndof::Integer=5,
        npart::Integer,
        nsteps::Integer,
        precision::DataType=Float64,
        )
For loading GCA solutions stored in a 3D array with dimesnions `(nsteps + 1,
ndof, npart`), and their magnetic moments. Double precision is assumed by
default.
"""
function load_gca(
    filename::String
    ;
    ndof::Integer=5,
    npart::Integer,
    nsteps::Integer,
    precision::DataType=Float64,
    )
    data = Array{precision, 3}(undef, nsteps+1, ndof, npart)
    magnetic_moment = Vector{precision}(undef, npart)
    f = open(filename)
    read!(f, data)
    read!(f, magneticmoment)
    close(f)
    return data, magnetic_moment
end


"""
    load_lightweight(
        filename::String
        ;
        RealT   ::DataType
        )
Loads data stored in the lightweight format by `save_lightweight`.
"""
function load_lightweight(
    filename::String
    ;
    RealT   ::DataType
    )
    file = open(filename)
    ndof = read(file, Int64)
    npart = read(file, Int64)
    u0 = Array{RealT}(undef, ndof, npart)
    uf = Array{RealT}(undef, ndof, npart)
    tf = Vector{RealT}(undef, npart)
    read!(file, u0)
    read!(file, uf)
    read!(file, tf)
    close(file)
    return ndof, npart, u0, uf, tf
end
