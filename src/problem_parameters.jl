#-------------------------------------------------------------------------------
# Created 23.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                problem_parameters.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
struct FOParams{qType} <: AbstractProblemParameters
    q ::qType
    m ::qType
    #B ::Vector{AbstractInterpolation}
    #E ::Vector{AbstractInterpolation}
    bx::AbstractInterpolation
    by::AbstractInterpolation
    bz::AbstractInterpolation
    ex::AbstractInterpolation
    ey::AbstractInterpolation
    ez::AbstractInterpolation

    function FOParams(
        q   ::Any,
        m   ::Any,
        mesh::AbstractMesh,
        B   ::AbstractArray,
        E   ::AbstractArray,
        )
        itp_type = Gridded(Linear())
        bx_interp = interpolate((mesh.x, mesh.y, mesh.z),B[1,:,:,:], itp_type)
        by_interp = interpolate((mesh.x, mesh.y, mesh.z),B[2,:,:,:], itp_type)
        bz_interp = interpolate((mesh.x, mesh.y, mesh.z),B[3,:,:,:], itp_type)
        ex_interp = interpolate((mesh.x, mesh.y, mesh.z),E[1,:,:,:], itp_type)
        ey_interp = interpolate((mesh.x, mesh.y, mesh.z),E[2,:,:,:], itp_type)
        ez_interp = interpolate((mesh.x, mesh.y, mesh.z),E[3,:,:,:], itp_type)
        bx_interp = extrapolate(bx_interp, Periodic())
        by_interp = extrapolate(by_interp, Periodic())
        bz_interp = extrapolate(bz_interp, Periodic())
        ex_interp = extrapolate(ex_interp, Periodic())
        ey_interp = extrapolate(ey_interp, Periodic())
        ez_interp = extrapolate(ez_interp, Periodic())
        return new{typeof(q)}(q, m,
                    bx_interp, by_interp, bz_interp, 
                    ex_interp, ey_interp, ez_interp
                    )
    end 
end
function (p::FOParams{<:Real})(i::Int64=1)
    return (p.q, p.m, p.bx, p.by, p.bz, p.ex, p.ey, p.ez)
end
function (p::FOParams{<:Vector{<:Real}})(i::Int64=1)
    return (p.q[i], p.m[i], p.bx, p.by, p.bz, p.ex, p.ey, p.ez)
end


struct GCAParams{qType,RealT} <: AbstractProblemParameters
    q      ::qType
    m      ::qType
    mu     ::Vector{RealT}
    B      ::tensor_interpolate
    E      ::tensor_interpolate
    gradB  ::tensor_interpolate
    gradb  ::tensor_interpolate
    gradExB::tensor_interpolate

    function GCAParams(
        q      ::Any,
        m      ::Any,
        mu     ::Vector{<:Real},
        mesh   ::AbstractMesh,
        B      ::AbstractArray,
        E      ::AbstractArray,
        gradB  ::AbstractArray,
        gradb  ::AbstractArray,
        gradExB::AbstractArray,
        )
        new{typeof(q), typeof(mu).parameters[1]}(
            q, m, mu, 
            EMfield_itps(mesh, B, E, gradB, gradb, gradExB)...
            ) 
    end
end
function (p::GCAParams{<:Real,<:Real})(i::Int64=1)
    return (p.q, p.m, p.mu[i], p.B, p.E, p.gradB, p.gradb, p.gradExB)
end 
function (p::GCAParams{<:Vector{<:Real},<:Real})(i::Int64=1)
    return (p.q[i], p.m[i], p.mu[i], p.B, p.E, p.gradB, p.gradb, p.gradExB)
end 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For stochastic differential equations (SDEs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct SingleCoefficientParams{RealT} <: AbstractProblemParameters
    alpha::RealT
end
function (p::SingleCoefficientParams{<:Real})(_::Int64)
    return (p.alpha)
end

struct DoubleCoefficientParams{RealT} <: AbstractProblemParameters
    alpha::RealT
    beta ::RealT
end
function (p::DoubleCoefficientParams{<:Real})(_::Int64)
    return (p.alpha, p.beta)
end

struct NoParams
end
function (p::NoParams)(_::Int64)
    return ()
end

struct GCAPitchAngleScatteringParams{qType} <: AbstractProblemParameters
    q      ::qType
    m      ::qType
    lnΛ    ::qType
    B      ::Any
    E      ::Any
    gradB  ::Any
    gradb  ::Any
    gradExB::Any
    n      ::Any
    T      ::Any

    function GCAPitchAngleScatteringParams(
        q      ::Any,
        m      ::Any,
        lnΛ    ::Any,
        mesh   ::AbstractMesh,
        B      ::AbstractArray,
        E      ::AbstractArray,
        gradB  ::AbstractArray,
        gradb  ::AbstractArray,
        gradExB::AbstractArray,
        n      ::AbstractArray,
        T      ::AbstractArray,
        )
        new{typeof(q)}(
            q, m, lnΛ,
            EMfield_itps(mesh, B, E, gradB, gradb, gradExB, n, T)...
            )
    end
    function GCAPitchAngleScatteringParams(
        q      ::Any,
        m      ::Any,
        lnΛ    ::Any,
        B      ::Function,
        E      ::Function,
        gradB  ::Function,
        gradb  ::Function,
        gradExB::Function,
        n      ::Function,
        T      ::Function,
        )
        new{typeof(q)}(
            q, m, lnΛ, B, E, gradB, gradb, gradExB, n, T
            )
    end
    function GCAPitchAngleScatteringParams(
        q      ::Any,
        m      ::Any,
        lnΛ    ::Any,
        B      ::tensor_interpolate,
        E      ::tensor_interpolate,
        gradB  ::tensor_interpolate,
        gradb  ::tensor_interpolate,
        gradExB::tensor_interpolate,
        n      ::tensor_interpolate,
        T      ::tensor_interpolate,
        )
        new{typeof(q)}(
            q, m, lnΛ, B, E, gradB, gradb, gradExB, n, T
            )
    end
end
function (p::GCAPitchAngleScatteringParams{<:Real})(_::Int64=1)
    return (p.q, p.m, p.lnΛ, p.B, p.E, p.gradB, p.gradb, p.gradExB, p.n, p.T)
end
function (p::GCAPitchAngleScatteringParams{<:Vector{<:Real}})(i::Int64=1)
    return (p.q[i], p.m[i]. p.lnΛ[i], p.B, p.E, p.gradB, p.gradb, p.gradExB,
        p.n, p.T
        )
end
