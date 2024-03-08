#-------------------------------------------------------------------------------
# Created 29.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dataprocessing.jl
#
#-------------------------------------------------------------------------------
# Contains functions for processing results.
#-------------------------------------------------------------------------------
"""
    getinitialstate(
        u::Vector{<:ODESolution},
        )
Return the initial state (or initial condition) of a vector of `ODEsolution`'s.
"""
function getinitialstate(
    u::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    ndof = length(u[1].u[1])
    RealT = typeof(u[1].u[1][1])
    u0 = Array{RealT}(undef, ndof, npart)
    for i = 1:npart
        u0[:,i] .= u[i].u[1]
    end
    return u0
end


"""
    getinitialstate(
        u::Array{<:Real, 3},:
        )
Return the initial state (or initial condition) of a solution stored as a
3D-array.
"""
function getinitialstate(
    u::Array{<:Real, 3},
    )
    return u[1,:,:]
end


"""
    getfinalstate(
        u::Vector{<:ODESolution},
        )
Return the final state of a vector of `ODEsolution`'s.
"""
function getfinalstate(
    u::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    ndof = length(u[1].u[1])
    RealT = typeof(u[1].u[1][1])
    uf = Array{RealT}(undef, ndof, npart)
    for i = 1:npart
        uf[:,i] .= u[i].u[end]
    end
    return uf
end


"""
    getfinalstate(
        u::Array{<:Real, 3},:
        )
Return the final state of a solution stored as a 3D-array.
"""
function getfinalstate(
    u::Array{<:Real, 3},
    )
    return u[end,:,:]
end


"""
    getfinaltime(
        u::Vector{<:ODESolution},
        )
Return the final times of a vector of `ODEsolution`'s.
"""
function getfinaltime(
    u::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    RealT = typeof(u[1].u[1][1])
    tf = Vector{RealT}(undef, npart)
    for i = 1:npart
        tf[i] = u.t[end]
    end
    return tf
end


function calc_energy(args...; kwargs...)
    ek0comp, ekfcomp, v0, vf = calc_energycomponents(args...; kwargs...)
    return sum.(ek0comp), sum.(ekfcomp), v0, vf
end


function calc_energycomponents(
    solution::Vector{Tuple{Any, Any, Any, Any, Any}},
    expname ::String,
    snap    ::Integer,
    expdir  ::String,
    ;
    mass  ::Real=tp.m_e,
    charge::Real=tp.e,
    units ::String="SI",
    destagger::Bool=true,
    itp_bc=(Flat(), Flat())
    )
    itpvec = get_br_emfield_vecof_interpolators(
        expname, snap, expdir;
        units=units, destagger=destagger, itp_bc=itp_bc
        )
    calc_energycomponents(solution, itpvec; mass=mass, charge=charge)
end


function calculateenergy(
    solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}}
    ;
    mass  ::Real=tp.m_e,
    charge::Real=tp.e
    )
    nofparticles = length(solution)
    ek0 = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    ekf = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    for (part, i) in zip(solution, 1:nofparticles)
        ek0[i] = kineticenergy(
            part[4],
            part[5],
            part[1][4],
            part[3],
            mass,
            )
        ekf[i] = kineticenergy(
            part[6],
            part[7],
            part[2][4],
            part[3],
            mass,
            )
    end
    return ek0, ekf
end


function calc_energycomponents(
    solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}},
    args...
    ;
    mass  ::Real=tp.m_e,
    charge::Real=tp.e
    )
    v0, vf = calc_drifts_and_vperp(solution, args..., mass=mass, charge=charge)
    return calc_energycomponents(solution, v0, vf)..., v0, vf
end

function calc_energycomponents(
    solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}},
    v0::Vector{<:Tuple{Real, Vector{<:Vector}}},
    vf::Vector{<:Tuple{Real, Vector{<:Vector}}},
    ;
    mass  ::Real=tp.m_e,
    )
    nofparticles = length(solution)
    ek0 = Vector{Vector{typeof(solution[1][1][1])}}(undef, nofparticles)
    ekf = Vector{Vector{typeof(solution[1][1][1])}}(undef, nofparticles)
    for (part, i) in zip(solution, 1:nofparticles)
        ek0[i] = [
            part[1][4] ^ 2,
            v0[i][1] ^2,
            norm(sum(v0[i][2])) ^2
            ]
        ekf[i] = [
            part[2][4] ^ 2,
            vf[i][1] ^2,
            norm(sum(vf[i][2])) ^2
            ]
    end
    return 0.5mass*ek0, 0.5mass*ekf
end


function calc_drifts_and_vperp(
    solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}}
    ;
    mass  ::Real=tp.m_e,
    kwargs...
    )
    nofparticles = length(solution)
    RealT = typeof(solution[1][1][1])
    v0 = Vector{Tuple{RealT, Vector{Vector{RealT}}}}(undef, nofparticles)
    vf = Vector{Tuple{RealT, Vector{Vector{RealT}}}}(undef, nofparticles)
    for (part, i) in zip(solution, 1:nofparticles)
        v0[i] = ( 
            perpendicular_velocity(part[3], mass, norm(part[4])),
            [exbdrift(part[4], part[5])[1]]
            )
        vf[i] = ( 
            perpendicular_velocity(part[3], mass, norm(part[6])),
            [exbdrift(part[6], part[7])[1]]
            )
    end
    return v0, vf
end


function calc_drifts_and_vperp(
    solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}},
    emfield_itpvec::Vector{<:AbstractInterpolation}
    ;
    charge::Real=-tp.e,
    mass  ::Real=tp.m_e,
    )
    nofparticles = length(solution)
    RealT = typeof(solution[1][1][1])
    v0 = Vector{Tuple{RealT, Vector{Vector{RealT}}}}(undef, nofparticles)
    vf = Vector{Tuple{RealT, Vector{Vector{RealT}}}}(undef, nofparticles)
    for (part, i) in zip(solution, 1:nofparticles)
        v0[i] = drifts_and_vperp_2Dxz(
            part[1][1:3],
            part[1][4],
            charge,
            mass,
            part[3],
            emfield_itpvec
            )
        vf[i] = drifts_and_vperp_2Dxz(
            part[2][1:3],
            part[2][4],
            charge,
            mass,
            part[3],
            emfield_itpvec
            )
    end
    return v0, vf
end


function avgnumsteps(
    solution::Vector{Tuple{
        Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any
        }}
    )
    sum = 0
    for particle in solution
        sum += particle[8]
    end
    return sum/length(solution)
end


function findinittemp(
    solution::Vector{Tuple{
        Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any
        }},
    tempitp::AbstractInterpolation
    )
    nofparticles = length(solution)
    T = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    for (particle, i) in zip(solution, 1:nofparticles)
        T[i] = tempitp(particle[1][1], particle[1][3])
    end
    return T
end

