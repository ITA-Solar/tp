
struct InterpolateTensor
    itp::Array{AbstractInterpolation}

    function InterpolateTensor(
        axes       ::Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}},
        tensorfield::Array{<:Real, 3},
        itp_type   ::Interpolations.InterpolationType,
        itp_bc     ::Interpolations.BoundaryCondition,
        ) 
        num_components = 1
        itp = Array{AbstractInterpolation}(undef, num_components)
        itp[1] = interpolate(axes, tensorfield[:,:,:], itp_type)
        itp[1] = extrapolate(itp[1], itp_bc)
        new(itp)
    end
    function InterpolateTensor(
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
        new(itp)
    end
    function InterpolateTensor(
        axes       ::Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}},
        tensorfield::Array{<:Real, 5},
        itp_type   ::Interpolations.InterpolationType,
        itp_bc     ::Interpolations.BoundaryCondition,
        ) 
        num_components_m, num_components_n, _, _, _ = size(tensorfield)
        itp = Array{AbstractInterpolation}(undef, num_components_m, num_components_n)
        for n = 1:num_components_n
            for m = 1:num_components_m
                itp[m,n] = interpolate(axes, tensorfield[m,n,:,:,:], itp_type)
                itp[m,n] = extrapolate(itp[m,n], itp_bc)
            end
        end 
        new(itp)
    end
end
function (itp::InterpolateTensor)(x,y,z)
    field = Array{typeof(x)}(undef, size(itp.itp)...)
    for i = eachindex(itp.itp)
        field[i] = itp.itp[i](x, y, z)
    end
    return field
end
