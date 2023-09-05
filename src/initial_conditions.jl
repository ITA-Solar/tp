#-------------------------------------------------------------------------------
# Created 23.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                initial_conditions.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

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

