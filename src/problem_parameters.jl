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
    B ::Any
    E ::Any

    function FOParams(
        q   ::Any,
        m   ::Any,
        mesh::AbstractMesh,
        B   ::AbstractArray,
        E   ::AbstractArray,
        )
        itp_type = Gridded(Linear())
        itp_bc = Periodic()
        axes = (mesh.x, mesh.y, mesh.z)
        B_itp = InterpolateTensor(axes, B, itp_type, itp_bc)
        E_itp = InterpolateTensor(axes, E, itp_type, itp_bc)
        return new{typeof(q)}(q, m, B_itp, E_itp)
    end 
end
function (p::FOParams{<:Real})(_::Int64=1)
    return (p.q, p.m, p.B, p.E)
end
function (p::FOParams{<:Vector{<:Real}})(i::Int64=1)
    return (p.q[i], p.m[i], p.B, p.E)
end


struct GCAParams{qType,RealT} <: AbstractProblemParameters
    q      ::qType
    m      ::qType
    mu     ::Vector{RealT}
    B      ::InterpolateTensor
    E      ::InterpolateTensor
    gradB  ::InterpolateTensor
    gradb  ::InterpolateTensor
    gradExB::InterpolateTensor

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
        itp_type = Gridded(Linear())
        itp_bc = Flat()
        axes = (mesh.x, mesh.y, mesh.z)
        #println("tp.jl: Getting GCA problem parameters...")
        #println("           Creating interpolation objects: B")
        B_itp = InterpolateTensor(axes, B, itp_type, itp_bc)
        #println("                                           E")
        E_itp = InterpolateTensor(axes, E, itp_type, itp_bc)
        #println("                                           grad B")
        gradB_itp = InterpolateTensor(axes, gradB, itp_type, itp_bc)
        #println("                                           grad b")
        gradb_itp = InterpolateTensor(axes, gradb, itp_type, itp_bc)
        #println("                                           grad ExB")
        gradExB_itp = InterpolateTensor(axes, gradExB, itp_type, itp_bc)
        new{typeof(q), typeof(mu).parameters[1]}(
            q, m, mu, 
            B_itp, E_itp, gradB_itp, gradb_itp, gradExB_itp
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

        itp_type = Gridded(Linear())
        itp_bc = Flat()
        axes = (mesh.x, mesh.y, mesh.z)
        #println("tp.jl: Getting GCA problem parameters...")
        #println("           Creating interpolation objects: B")
        B_itp = InterpolateTensor(axes, B, itp_type, itp_bc)
        #println("                                           E")
        E_itp = InterpolateTensor(axes, E, itp_type, itp_bc)
        #println("                                           grad B")
        gradB_itp = InterpolateTensor(axes, gradB, itp_type, itp_bc)
        #println("                                           grad b")
        gradb_itp = InterpolateTensor(axes, gradb, itp_type, itp_bc)
        #println("                                           grad ExB")
        gradExB_itp = InterpolateTensor(axes, gradExB, itp_type, itp_bc)
        n_itp = InterpolateTensor(axes, n, itp_type, itp_bc)
        T_itp = InterpolateTensor(axes, T, itp_type, itp_bc)
        new{typeof(q)}(
            q, m, lnΛ,
            B_itp, E_itp, gradB_itp, gradb_itp, gradExB_itp, n_itp, T_itp
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
