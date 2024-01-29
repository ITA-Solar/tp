#-------------------------------------------------------------------------------
# Created 25.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 save.jl
#
#-------------------------------------------------------------------------------
# Contains functions for saving data.
#-------------------------------------------------------------------------------
"""
    save(
        expname::String,
        expdir ::String,
        args...
        )
Save the data slurped in `args` with `expname` at `expdir`.
"""
function save(
    expname::String,
    expdir ::String,
    args...
    )
    basename = joinpath(expdir, expname)
    filename = string(basename, ".tp")
    save(filename, args...)
end


"""
    save(
        filename::String,
        solution::EnsembleSolution,
        )
Save an `EnsambleSolution` using `JLD`.
"""
function save(
    filename::String,
    solution::EnsembleSolution,
    )
    save(filename, solution.u)
end

"""
    save(
        filename::String,
        solutionvector::Vector{<:ODESolution},
        )
Save a vector of `ODESolution`'s using `JLD`.
"""
function save(
    filename::String,
    solutionvector::Vector{<:ODESolution},
    )
    npart = length(solutionvector)
    jldopen(filename, "w") do file
        write(file, "npart", npart)
        for i = 1:npart
            write(file, "u$i", solutionvector[i].u)
            write(file, "t$i", solutionvector[i].t)
        end
    end
    println("tp.jl: Wrote $filename")
end


"""
    save(
        filename::String,
        solution::Array{<:Real, 3},
        )
Save the numerical solution which has the form of a 3D-array. E.g. the solution
to the Loretnz equation.
"""
function save(
    filename::String,
    solution::Array{<:Real, 3},
    )
    f = open(filename, "w+")
    write(f, solution)
    close(f)
    println("tp.jl: Wrote $filename")
end 


"""
    save(
        filename::String,
        solution::Array{<:Real, 3},
        magneticmoment::Vector{<:Real}
        )
Save the numerical solution of the guiding centre approximation - which has the
form of a 3D-array - and the magnetic moment - which is a parameter of the
equations.
"""
function save(
    filename::String,
    solution::Array{<:Real, 3},
    magneticmoment::Vector{<:Real}
    )
    f = open(filename, "w+")
    write(f, solution)
    write(f, magneticmoment)
    close(f)
    println("tp.jl: Wrote $filename")
end 


"""
    save_lightweight(
        filename::String,
        solutionvector::Vector{<:ODESolution},
        )
Save only the initial state, final state, and final time of the ODE-solutions.
"""
function save_lightweight(
    filename::String,
    solutionvector::Vector{<:ODESolution},
    )
    u0 = getinitialstate(solutionvector)
    uf = getfinalstate(solutionvector)
    tf = getfinaltime(solutionvector)
    ndof, npart = size(u0)
    file = open(filename, "w+")
    write(file, ndof)
    write(file, npart)
    write(file, u0)
    write(file, uf)
    write(file, tf)
    close(file)
    println("tp.jl: Wrote $filename")
end


"""
    save(
        filename::String,
        structure::DataType,
        )
Save a `DataType` of type `structure` in a runnable Julia-script.

This is a rather complicated way of saving a structure as a Julia-script, but
it looks nice. A simpler way would be to just write

    f = open(filename, "w+")
    write(f, string("instance = ", structure))
    close(f)

and that would also be a runnable script.
"""
function save(
    filename::String,
    structure::DataType,
    )
    writestring = "instance = $(typeof(structure))(\n"
    for field in fieldnames(structure)
        if isdefined(structure, field)
            if typeof(getfield(structure, field)) == String
                value = "\"$(getfield(structure, field))\""
            else
                value = "$(getfield(structure, field))"
            end
            spaces = " "^(11 - length("$field"))
            writestring = string(writestring,
                                  "\t$field", spaces, "= $value,\n")
        end
    end
    writestring = string(writestring, ")")
    #
    f = open(filename, "w+")
    write(f, writestring)
    close(f)
    #
    println("tp.jl: Wrote $(filename)")
end
