#-------------------------------------------------------------------------------
# Created 12.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             Interpolations_tp.jl
#
#-------------------------------------------------------------------------------
# Module containing interpolation methods.
#-------------------------------------------------------------------------------

module Interpolations_tp

using WorkingPrecision: wpFloat, wpInt
using Meshes

export gridinterp
export locateCell
export trilinear
export bilinear_xz


#-----------------------#
# Type of interpolation #
#-------------------------------------------------------------------------------
function gridinterp(
    tensorfield ::AbstractArray{wpFloat},
    interpolator::Function,
    pos         ::Vector{wpFloat},
    xx          ::Vector{wpFloat},
    yy          ::Vector{wpFloat},
    zz          ::Vector{wpFloat}
    )
    numDims = 3
    x, y, z = pos
    # Find which cell the position is in, and store the indices of its corner.
    i = locateCell(xx, x)
    j = locateCell(yy, y)
    k = locateCell(zz, z)
    interpolatedfield = interpolator(
        tensorfield,
        xx, yy, zz,
        (i,j,k),
        (x,y,z)
    )
    return interpolatedfield, i, j, k
end # function gridinterp
#|
function gridinterp(
    mesh        ::Mesh,
    interpolator::Function,
    pos         ::Vector{wpFloat}
    )
    numDims = 3
    x, y, z = pos
    # Find which cell the position is in, and store the indices of its corner.
    i = locateCell(mesh.xCoords, x)
    j = locateCell(mesh.yCoords, y)
    k = locateCell(mesh.zCoords, z)
    # Interpolate
    interpolatedfields = interpolator(mesh, (i,j,k), (x,y,z))
    return interpolatedfields, i, j, k
end # function gridinterp



#-----------#
# Utilities #
#-------------------------------------------------------------------------------
"""
    locateCell(coords, point)

Finds the position of `point` in the vector `coords` using binary search.
Returns its lower neighbour.
"""
function locateCell(coords::Vector{wpFloat},
                    point ::Number
                    )
    # Initial bounds
    low = 1
    high = length(coords)
    while  high > low + 1 # While we haven't found the cell
        mid  = floor(wpInt, (high + low)/2)
        if point > coords[mid]
            low = mid 
        elseif point < coords[mid]
            high = mid
	else
	    low = mid - 1
	    high = mid
        end # if
    end # while
    return low 
end # function locateCell


#---------------------#
# Interpolation steps #
#-------------------------------------------------------------------------------
function trilinear( # Arbitrary vector-field
    tensorfield ::AbstractArray{wpFloat},
    xx          ::Vector{wpFloat},
    yy          ::Vector{wpFloat},
    zz          ::Vector{wpFloat},
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
    )
    coefficients = trilinearcoefficients(xx,
                                         yy,
                                         zz,
                                         (i,j,k),
                                         (x,y,z))
    f = trilinearsum(tensorfield,
                     (i,j,k),
                     coefficients)
    return f
end # function trilinear
#|
function trilinear( # Method for passing the mesh-struct
    mesh   ::Mesh,
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
    )
    coefficients = trilinearcoefficients(mesh.xCoords,
                                         mesh.yCoords,
                                         mesh.zCoords,
                                         (i,j,k),
                                         (x,y,z))
    B = trilinearsum(mesh.bField,
                              (i,j,k),
                              coefficients)
    E = trilinearsum(mesh.eField,
                              (i,j,k),
                              coefficients)
    return B, E
end # function trilinear


function bilinear_xz( # Arbitrary vector-field
    tensorfield ::AbstractArray{wpFloat},
    xx          ::Vector{wpFloat},
    yy          ::Vector{wpFloat},
    zz          ::Vector{wpFloat},
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
    )
    coefficients = bilinearcoefficients(xx,
                                        zz,
                                        (i,j),
                                        (x,z))
    f = bilinearsum(dropdims(tensorfield, dims=ndims(tensorfield)-1),
                    (i,j),
                    coefficients)
    return f
end # function trilinear
#|
function bilinear_xz( 
    mesh   ::Mesh,
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
    )
    coefficients = bilinearcoefficients(mesh.xCoords,
                                        mesh.zCoords,
                                        (i,j),
                                        (x,z))
    B = bilinearsum(dropdims(mesh.bField, dims=3),
                    (i,j),
                    coefficients)
    E = bilinearsum(dropdims(mesh.eField, dims=3),
                    (i,j),
                    coefficients)
    return B, E
end # function trilinear


function trilinearGCA(
    mesh   ::Mesh,
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
    )
    coefficients = trilinearcoefficients(mesh.xCoords,
                                         mesh.yCoords,
                                         mesh.xCoords,
                                         (i,j,k),
                                         (x,y,z))
    B = trilinearsum(mesh.bField,
                              (i,j,k),
                              coefficients)
    E = trilinearsum(mesh.eField,
                              (i,j,k),
                              coefficients)
    ∇B = trilinearsum(mesh.∇B,
                              (i,j,k),
                              coefficients)
    ∇b̂ = trilinearsum(mesh.∇b̂,
                              (i,j,k),
                              coefficients)
    ∇ExB = trilinearsum(mesh.∇ExB,
                              (i,j,k),
                              coefficients)

    return B, E, ∇B, ∇b̂, ∇ExB
end # function trilinearGCA


function trilinearcoefficients(
    xCoords::Vector{wpFloat},
    yCoords::Vector{wpFloat},
    zCoords::Vector{wpFloat},
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
    )
    t = (x - xCoords[i])/(xCoords[i + 1] - xCoords[i])
    u = (y - yCoords[j])/(yCoords[j + 1] - yCoords[j])
    v = (z - zCoords[k])/(zCoords[k + 1] - zCoords[k])
    #
    c0 = (1 - t)*(1 - u)*(1 - v)
    c1 = t*(1 - u)*(1 - v)
    c2 = t*u*(1 - v)
    c3 = (1 - t)*u*(1 - v)
    c4 = (1 - t)*(1 - u)*v
    c5 = t*(1 - u)*v
    c6 = t*u*v
    c7 = (1 - t)*u*v
    #
    return c0, c1, c2, c3, c4, c5, c6, c7
end # function trinlinearcoefficients


function bilinearcoefficients(
    xCoords::Vector{wpFloat},
    yCoords::Vector{wpFloat},
    (i,j)::Tuple{Integer, Integer},
    (x,y)::Tuple{wpFloat, wpFloat}
    )
    t = (x - xCoords[i])/(xCoords[i + 1] - xCoords[i])
    u = (y - yCoords[j])/(yCoords[j + 1] - yCoords[j])
    #
    c0 = (1 - t)*(1 - u)
    c1 = t*(1 - u)
    c2 = t*u
    c3 = (1 - t)*u
    #
    return c0, c1, c2, c3
end # bilinearcoefficients

function trilinearsum(
    tensorfield::Array{wpFloat, 5},
    (i,j,k)    ::Tuple{Integer, Integer, Integer},
    c          ::NTuple{8, wpFloat}
    )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = tensorfield[:, :,   i,   j,   k]
    A1 = tensorfield[:, :, i+1,   j,   k]
    A2 = tensorfield[:, :, i+1, j+1,   k]
    A3 = tensorfield[:, :,   i, j+1,   k]
    A4 = tensorfield[:, :,   i,   j, k+1]
    A5 = tensorfield[:, :, i+1,   j, k+1]
    A6 = tensorfield[:, :, i+1, j+1, k+1]
    A7 = tensorfield[:, :,   i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearsum
#|
function trilinearsum(
    vectorfield::Array{wpFloat, 4},
    (i,j,k)    ::Tuple{Integer, Integer, Integer},
    c          ::NTuple{8, wpFloat}
    )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = vectorfield[:,   i,   j,   k]
    A1 = vectorfield[:, i+1,   j,   k]
    A2 = vectorfield[:, i+1, j+1,   k]
    A3 = vectorfield[:,   i, j+1,   k]
    A4 = vectorfield[:,   i,   j, k+1]
    A5 = vectorfield[:, i+1,   j, k+1]
    A6 = vectorfield[:, i+1, j+1, k+1]
    A7 = vectorfield[:,   i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearsum
#|
function trilinearsum(
    scalarfield::Array{wpFloat, 3},
    (i,j,k)    ::Tuple{Integer, Integer, Integer},
    c          ::NTuple{8, wpFloat}
    )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = scalarfield[  i,   j,   k]
    A1 = scalarfield[i+1,   j,   k]
    A2 = scalarfield[i+1, j+1,   k]
    A3 = scalarfield[  i, j+1,   k]
    A4 = scalarfield[  i,   j, k+1]
    A5 = scalarfield[i+1,   j, k+1]
    A6 = scalarfield[i+1, j+1, k+1]
    A7 = scalarfield[  i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearsum


function bilinearsum(
    tensorfield::Array{wpFloat, 4},
    (i,j)      ::Tuple{Integer, Integer},
    c          ::NTuple{4, wpFloat}
    )
    c0, c1, c2, c3 = c
    A0 = tensorfield[:, :,   i,   j]
    A1 = tensorfield[:, :, i+1,   j]
    A2 = tensorfield[:, :, i+1, j+1]
    A3 = tensorfield[:, :,   i, j+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 
    return A
end # bilinearsum
#|
function bilinearsum(
    vectorfield::Array{wpFloat, 3},
    (i,j)      ::Tuple{Integer, Integer},
    c          ::NTuple{4, wpFloat}
    )
    c0, c1, c2, c3 = c
    A0 = vectorfield[:,   i,   j]
    A1 = vectorfield[:, i+1,   j]
    A2 = vectorfield[:, i+1, j+1]
    A3 = vectorfield[:,   i, j+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 
    return A
end # bilinearsum
#|
function bilinearsum(
    scalarfield::Array{wpFloat, 2},
    (i,j)      ::Tuple{Integer, Integer},
    c          ::NTuple{4, wpFloat}
    )
    c0, c1, c2, c3 = c
    A0 = scalarfield[  i,   j]
    A1 = scalarfield[i+1,   j]
    A2 = scalarfield[i+1, j+1]
    A3 = scalarfield[  i, j+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 
    return A
end # bilinearsum

end # module Interpolations_tp
