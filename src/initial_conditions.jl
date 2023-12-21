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
function (ic::FloatIC)(_::Integer=1)
    return ic.float
end

struct VectorIC{RealT<:Real} <: AbstractVariableIC
    vec::Vector{RealT}
end
function (ic::VectorIC)(i::Integer=1)
    return ic.vec[i]
end


struct PhaseSpace1DIC <: AbstractInitialConditions
    x::AbstractVariableIC
    v::AbstractVariableIC
end
function (ic::PhaseSpace1DIC)(i::Int64=1)
    return [ic.x(i); ic.v(i)]
end

struct GCAPitchAngleScatteringIC <: AbstractInitialConditions
    rx    ::AbstractVariableIC
    ry    ::AbstractVariableIC
    rz    ::AbstractVariableIC
    vparal::AbstractVariableIC
    beta  ::AbstractVariableIC
end
function (ic::GCAPitchAngleScatteringIC)(i::Integer=1)
    return [ic.rx(i), ic.ry(i), ic.rz(i), ic.vparal(i), ic.beta(i)]
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


#-------------------------------------------------------------------------------
# Functions for the generation of initial values based on Experiment Parameters
# giving distributions and bounds
#
function get_initial_positions(
    npart      ::Integer,
    pos_distr  ::String,
    pos_xbounds::Vector{<:Real},
    pos_ybounds::Vector{<:Real},
    pos_zbounds::Vector{<:Real},
    RealT      ::DataType,
    ;
    seed::Integer=0
    )
    ndims = 3
    if pos_distr == "uniform"
        # Set default bounds if not defined
        # Default bounds are the mesh domain boundaries.
        println("           Drawing positions")
        @time pos = inituniform(
            npart,
            pos_xbounds,
            pos_ybounds,
            pos_zbounds,
            RealT
            ;
            seed
        )
    elseif pos_distr == "point"
        pos = ones(RealT, ndims, npart)
        # If particles are given by coordinates
        if (length(pos_xbounds) ==
            length(pos_ybounds) ==
            length(pos_zbounds) == npart)
            pos[1,:] = pos_xbounds
            pos[2,:] = pos_ybounds
            pos[3,:] = pos_zbounds
        else
            pos[1,:] *= pos_xbounds[1]
            pos[2,:] *= pos_ybounds[1]
            pos[3,:] *= pos_zbounds[1]
        end
    else
        error("Invalid parameter value: pos_distr")
    end
    # IMPLEMENT e.g
    #
    # position distributed according to density
    #
    return pos
end

function get_initial_velocities(
    npart      ::Integer,
    vel_distr  ::String,
    vel_xbounds::Vector{<:Real},
    vel_ybounds::Vector{<:Real},
    vel_zbounds::Vector{<:Real},
    RealT      ::DataType,
    ;
    seed::Integer=0,
    stds=nothing
    )
    ndims = 3
    vel = Array{RealT}(undef, ndims, npart)
    if vel_distr == "mb" || vel_distr == "mb-onetemp"
        println("           Drawing velocities")
        @time for i = 1:npart
            vel[:,i] = randn(0.0, stds[i], (ndims); seed=seed+i-1)
        end
    elseif vel_distr == "uniform"
        println("           Drawing velocities")
        @time vel = inituniform(
            npart,
            vel_xbounds,
            vel_ybounds,
            vel_zbounds,
            RealT
            ;
            seed=seed
            )
    elseif vel_distr == "point"
        if (length(vel_xbounds) ==
            length(vel_ybounds) ==
            length(vel_zbounds) == npart)
            vel[1,:] = vel_xbounds
            vel[2,:] = vel_ybounds
            vel[3,:] = vel_zbounds
        else
            vel = ones(ndims, npart)
            vel[2,:] *= vel_ybounds[1]
            vel[3,:] *= vel_zbounds[1]
        end
    else
        error("Invalid parameter value: vel_distr")
    end
    return vel
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

