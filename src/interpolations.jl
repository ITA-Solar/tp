
function tensor_interpolate(
    axes       ::Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}},
    tensorfield::Array{<:Real, 4},
    itp_type   ::Interpolations.InterpolationType,
    itp_bc     ::Interpolations.BoundaryCondition,
    ) 
    num_components, _, _, _ = size(tensorfield)
    itp = Array{AbstractInterpolation}(undef, num_components)
    for i = 1:num_components
        itp[i] = interpolate(axes, tensorfield[i,:,:,:], itp_type)
        itp[i] = extrapolate(itp[i], itp_bc)
    end 
    return itp
end

function tensor_interpolate(
    axes       ::Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}},
    tensorfield::Array{<:Real, 5},
    itp_type   ::Interpolations.InterpolationType,
    itp_bc     ::Interpolations.BoundaryCondition,
    ) 
    num_components_m, num_components_n, _, _, _ = size(tensorfield)
    num_components_i = size(tensorfield)[3]
    itp = Array{AbstractInterpolation}(undef, num_components_m, num_components_n)
    for n = 1:num_components_n
        for m = 1:num_components_m
            itp[m,n] = interpolate(axes, tensorfield[m,n,:,:,:], itp_type)
            itp[m,n] = extrapolate(itp[m,n], itp_bc)
        end
    end 
    return itp
end

function EMfield_itps(
    mesh   ::AbstractMesh,
    B      ::AbstractArray,
    E      ::AbstractArray,
    )
    itp_type = Gridded(Linear())
    itp_bc = Flat()
    axes = (mesh.x, mesh.y, mesh.z)
    B_itp = tensor_interpolate(axes, B, itp_type, itp_bc)
    E_itp = tensor_interpolate(axes, E, itp_type, itp_bc)
    return B_itp, E_itp
end
function EMfield_itps(
    mesh   ::AbstractMesh,
    B      ::AbstractArray,
    E      ::AbstractArray,
    gradB  ::AbstractArray,
    gradb  ::AbstractArray,
    gradExB::AbstractArray,
    )
    itp_type = Gridded(Linear())
    itp_bc = Flat()
    axes = (mesh.x, mesh.y, mesh.z)
    B_itp = tensor_interpolate(axes, B, itp_type, itp_bc)
    E_itp = tensor_interpolate(axes, E, itp_type, itp_bc)
    gradB_itp = tensor_interpolate(axes, gradB, itp_type, itp_bc)
    gradb_itp = tensor_interpolate(axes, gradb, itp_type, itp_bc)
    gradExB_itp = tensor_interpolate(axes, gradExB, itp_type, itp_bc)
    return B_itp, E_itp, gradB_itp, gradb_itp, gradExB_itp
end
