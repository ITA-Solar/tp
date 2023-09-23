#-------------------------------------------------------------------------------
# Created 23.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                initial_conditions.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

abstract type AbstractVariableIC <: AbstractInitialConditions end

struct FOStaticIC{RealT<:Real} <: AbstractInitialConditions
    pos::Matrix{RealT}
    vel::Matrix{RealT}
end
function (ic::FOStaticIC)(i::Integer=1)
    return [ic.pos[:,i]; ic.vel[:,i]]
end 

struct GCAStaticIC{RealT<:Real} <: AbstractInitialConditions
    R::Matrix{RealT}
    vparal::Vector{RealT}
end
function (ic::GCAStaticIC)(i::Integer=1)
    return [ic.R[:,i]; ic.vparal[i]]
end


struct FO_IC <: AbstractInitialConditions
    x ::AbstractVariableIC
    y ::AbstractVariableIC
    z ::AbstractVariableIC
    vx::AbstractVariableIC
    vy::AbstractVariableIC
    vz::AbstractVariableIC
end
function (ic::FO_IC)(i::Integer=1)
    return [ic.x(i), ic.y(i), ic.z(i), ic.vx(i), ic.vy(i), ic.vz(i)]
end 

struct GCA_IC <: AbstractInitialConditions
    rx    ::AbstractVariableIC
    ry    ::AbstractVariableIC
    rz    ::AbstractVariableIC
    vparal::AbstractVariableIC
end
function (ic::GCA_IC)(i::Integer=1)
    return [ic.rx(i), ic.ry(i), ic.rz(i), ic.vparal(i)]
end 

struct FloatIC{RealT<:Real} <: AbstractVariableIC
    float::RealT
end
function (ic::FloatIC)(i::Integer=1)
    return ic.float
end

struct VectorIC{RealT<:Real} <: AbstractVariableIC
    vec::Vector{RealT}
end
function (ic::VectorIC)(i::Integer=1)
    return ic.vec[i]
end

#-------------------------------------------------------------------------------
# On the fly generalisation of initial conditions
#
struct RandomUniformIC{RealT<:Real} <: AbstractVariableIC
    init_seed::Int64
    curr_seed::Int64
    lower_bound::RealT
    extent     ::RealT
end
function (ic::RandomUniformIC)(i::Integer=1)
    ic.curr_seed += 1
    r = rand(MersenneTwister(ic.curr_seed),
             typeof(ic).parameters[1],
             )
    return lower_bound + extent * r
end

#struct RandomMaxwellianIC{RealT<:Real} <: AbstractVariableIC
#    init_seed::Int64
#    curr_seed::Int64
#    lower_bound::
#    extent     ::RealT
#end


function calc_GCA_IC_and_mu(
    pos  ::Matrix{<:Real},
    vel  ::Matrix{<:Real},
    q    ::Any,
    m    ::Any,
    B_vec::Matrix{<:Real},
    E_vec::Matrix{<:Real},
    )
    ndims, npart = size(pos)
    RealT = typeof(vel).parameters[1]
    vparal = Array{RealT}(undef, npart)
    mu = Array{RealT}(undef, npart)
    R = Array{RealT}(undef, ndims, npart)
    if typeof(q) == typeof(m) <:Real
        m_over_q = ones(npart)*m/q
    elseif typeof(q) == typeof(m) <: Vector{<:Real}
        m_over_q = m ./ q
    else
        error("Typef og q and m ($(typeof(q)), $(typeof(m))) not supported")
    end
    for i = 1:npart
        B = norm(B_vec[:,i])
        b_vec = B_vec[:,i]/B
        ExBdrift =  (E_vec[:,i] × b_vec)/B
        vel_in_E_frame = vel[:,i] - ExBdrift
        # Calculate the guiding centre posistion 
        R[:,i] = pos[:,i] - m_over_q[i]/B * (vel_in_E_frame × b_vec)
        # Calculate the velocity parallell to the magnetic field -- vparal
        vparal[i] = vel[:,i] ⋅ b_vec
        # Calculate mangetic moment -- mu
        vperp = vel_in_E_frame - vparal[i]*b_vec
        mu[i] = m*norm(vperp)^2/(2B)
    end 
    return R, vparal, mu
end

