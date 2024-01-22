# Created 05.01.24 by eilfso
# Author: e.s.oyre@astro.uio.no
#
#           mathematics.jl
#-------------------------------------------------------------------------------
# Contains calculus
#-------------------------------------------------------------------------------


"""
    skewsymmetric_matrix(vector::Vector)
Returns the skew-symmetric matrix of a 3D vector. The matrix is used to 
represent a cross product between two 3D vectors as a matrix product between 
the skew-symmetric matrix and the other vector. It is probably useful for 
other stuff too.
"""
function skewsymmetric_matrix(vector::Vector{Float64})
    return [       0.0 -vector[3]  vector[2]
             vector[3]        0.0 -vector[1]
            -vector[2]  vector[1]        0.0]
end
function skewsymmetric_matrix(vector::Vector{Float32})
    return [0.0f0 -vector[3] vector[2]
            vector[3] 0.0f0 -vector[1]
            -vector[2] vector[1] 0.0f0]
end


"""
    norm2(field, axis=1)
Calculates the p=2 norm of the vectors in a 1D vector field, i.e. the field
strength. The argument 'axis' determines whether the vector components are
stored in the first or second dimension of the array storing the field.
"""
function norm2(
    field::Matrix{T} where {T<:Real},
    axis ::Integer=1
    ;
    wfp  ::DataType=typeof(field[1])
    )
    dims = size(field)
    if axis == 1
        fieldstrength = zeros(wfp, dims[2])
        for i = 1:dims[2]
            fieldstrength[i] = norm(field[:, i])
        end # loop i
    elseif axis == 2
        fieldstrength = zeros(wfp, dims[1])
        for i = 1:dims[1]
            fieldstrength[i] = norm(field[i, :])
        end
    else
        println("Error: Your axes are wierd...")
    end # if
    return fieldstrength
end # function norm2


"""
    norm3(field)
Calculates the p=2 norm of the vectors in a 2D vector field, i.e. the field
strength. The function assumes the vector components are store in the first
dimension of the field array.
"""
function norm3(
    field::Array{T} where {T<:Real}
    ;
    wfp  ::DataType=typeof(field[1])
    ) 
    dims = size(field)
    fieldstrength = zeros(wfp, dims[2:3])
    for i = 1:dims[2]
        for j = 1:dims[3]
            fieldstrength[i,j] = norm(field[:, i, j])
        end # loop j
    end # loop i
    return fieldstrength
end # function norm2


"""
    norm4(field, axis)
Calculates the p=2 norm of the vectors in a 3D vector field, i.e. the field
strength. The argument 'axis' determines whether the vector components are
stored in the first or fourth dimension of the array storing the field.
"""
function norm4(
    field::Array{T} where {T<:Real},
    axis ::Integer=1
    ;
    wfp  ::DataType=typeof(field[1])
               )
    if axis == 1
        dims = size(field[1,:,:,:])
        fieldstrength = zeros(wfp, dims)
        for i = 1:dims[1]
            for j = 1:dims[2]
                for k = 1:dims[3]
                    fieldstrength[i,j,k] = √(field[1,i,j,k]^2 + 
                                             field[2,i,j,k]^2 +
                                             field[3,i,j,k]^2)
                end # loop k
            end # look j
        end # loop i
    elseif axis == 4
        dims = size(field[:,:,:,1])
        fieldstrength = zeros(wfp, dims)
        for i = 1:dims[1]
            for j = 1:dims[2]
                for k = 1:dims[3]
                    fieldstrength[i,j,k] = √(field[i,j,k,1]^2 + 
                                             field[i,j,k,2]^2 +
                                             field[i,j,k,3]^2)
                end # loop k
            end # look j
        end # loop i
    else
        println("Error: Yours axes are wierd...")
    end # if
    return fieldstrength
end # function norm4
#|
function norm4(
    fx::Array{T} where {T<:Real},
    fy::Array{T} where {T<:Real},
    fz::Array{T} where {T<:Real},
    ;
    wfp  ::DataType=typeof(fx[1])
               )
    dims = size(fx[:,:,:])
    fieldstrength = zeros(wfp, dims)
    for i = 1:dims[1]
        for j = 1:dims[2]
            for k = 1:dims[3]
                fieldstrength[i,j,k] = √(fx[i,j,k]^2 + 
                    fy[i,j,k]^2 +
                    fz[i,j,k]^2)
            end # loop k
        end # look j
    end # loop i
    return fieldstrength
end # function norm4
