#-------------------------------------------------------------------------------
# Created 05.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Particles.jl
#
#-------------------------------------------------------------------------------
# Module containing the Particle structs and methods
#-------------------------------------------------------------------------------



#------------------#
# Global variables #
#------------------#------------------------------------------------------------
#              mass charge
specieTable = [m_e   -e     # Electron
               m_p    e     # Proton
               1.0   3.0    # Unit mass and charge = 3
               1.0   1.0    # Unit mass and charge
               1.0  -1.0    # Unit mass and negative unit charge
               ]

#-------------#   
# Structs     # 
#-------------#-----------------------------------------------------------------

struct ODEParticle{Equations} <: AbstractParticle
    eom  ::Equations
    ic   ::AbstractInitialConditions
    p    ::AbstractProblemParameters
    npart::Int64
end

function get_problem(part::ODEParticle, tspan)
    return ODEProblem(part.eom, part.ic(), tspan, part.p())
end

function get_prob_func(part::ODEParticle)
    return (prob, i, repeat) -> (
                #println("tp.jl: Solving particle $i");
                remake(prob, u0=part.ic(i), p=part.p(i))
                )
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For stochastic differential equations (SDEs)
#
struct SDEParticle{DriftEq, DiffusionEq} <: AbstractParticle
    eom_drift    ::DriftEq
    eom_diffusion::DiffusionEq
    ic           ::AbstractInitialConditions
    p_drift      ::AbstractProblemParameters
    p_diffusion  ::AbstractProblemParameters
    npart        ::Int64
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#------------------#
# Legacy-particles #
#------------------#------------------------------------------------------------
    
mutable struct ParticleSoA <: TraceParticle
    pos    ::Array{T, 3} where {T<:Real}
    vel    ::Array{T, 3} where {T<:Real}
    species::Vector{T} where {T<:Integer}   # Particle specie identifier
                                            # (e.g. electron, proton)
    alive  ::Vector{Bool}
    weight ::Vector{T} where {T<:Real}
    
    
    # Constructors
    #--------------------------------------------------------------------------
    """
        ParticleSoA(pos::Matrix, vel::Matrix, specie, numSteps)

    One would normally only pass initial conditions. This constructor handles
    the creation of the type accordingly, by adding the initial conditions to
    an higher order array.
    """
    # "Default" constructor
    function ParticleSoA( 
        pos    ::Array{T, 3} where {T<:Real},
        vel    ::Array{T, 3} where {T<:Real},
        species::Vector{T} where {T<:Integer},   # Particle specie ID
                                                 # (e.g. electron, proton)
        alive  ::Vector{Bool},
        weight ::Vector{T} where {T<:Real}
        )
        return new(pos, vel, species, alive, weight)
    end # constructor 
    #|
    function ParticleSoA(
        pos    ::Array{T, 3} where {T<:Real},
        vel    ::Array{T, 3} where {T<:Real},
        species::Vector{T} where {T<:Integer},   # Particle specie ID
                                                 # (e.g. electron, proton)
        ;
        alive  ::Vector{Bool}=ones(Bool, size(pos)[2]),
        weight ::Vector{T} where {T<:Real}=ones(typeof(pos[1]),
                                                size(pos)[2])
        )
        return new(pos, vel, species, alive, weight)
    end # constructor 
    #|
    function ParticleSoA(
        pos0     ::Matrix{T} where {T<:Real},
        vel0     ::Matrix{T} where {T<:Real},
        species ::Vector{T} where {T<:Integer},
        numSteps::Integer
        ;
        wfp::DataType=typeof(pos0[1])
        )
        numDims, numParticles = size(pos0)
        positions  = zeros(wfp, numDims, numParticles, numSteps + 1)
        velocities = zeros(wfp, numDims, numParticles, numSteps + 1)
        positions[:, :, 1] .= pos0
        velocities[:, :, 1] .= vel0
        alive = ones(Bool, numParticles)
        weight = ones(wfp, numParticles)
        return new(positions, velocities, species, alive, weight)
    end # constructor 
    #|
    function ParticleSoA(
        pos0    ::Vector{T} where {T<:Real},
        vel0    ::Vector{T} where {T<:Real},
        species ::Vector{T} where {T<:Integer},
        numSteps::Integer
        ;
        wfp::DataType=typeof(pos0[1])
        )
        numParticles = 1
        numDims = length(pos0)
        positions  = zeros(wfp, numDims, numParticles, numSteps + 1)
        velocities = zeros(wfp, numDims, numParticles, numSteps + 1)
        positions[:, :, 1] .= pos0
        velocities[:, :, 1] .= vel0
        alive = ones(Bool, numParticles)
        weight = ones(wfp, numParticles)
        return new(positions, velocities, species, alive, weight)
    end # constructor 
end # mutable struct ParticleSoA


mutable struct GCAParticleSoA <: TraceParticle
    R      ::Array{T, 3} where {T<:Real} # The position of the guiding centre
    vparal ::Matrix{T} where {T<:Real}   # Velocity parallel to the magnetic field
    μ      ::Vector{T} where {T<:Real}   # Magnetic moment μ of particle
    species::Vector{T} where {T<:Integer}# Particle specie ID (e.g. electron, proton)
    alive  ::Vector{Bool}
    weight ::Vector{T} where {T<:Real}
    
    
    # Constructors
    #--------------------------------------------------------------------------
    """
        GCAParticleSoA(R, vparal, μ, species, alive, weight)

    One would normally only pass initial conditions. This constructor handles
    the creation of the type accordingly, by adding the initial conditions to
    an higher order array.
    """
    function GCAParticleSoA( 
        R      ::Array{T, 3} where {T<:Real},# The position of the guiding centre
        vparal ::Matrix{T} where {T<:Real},   # The velocity parallel to the
                                              # magnetic field
        μ      ::Vector{T} where {T<:Real},   # Magnetic moment μ of particle
        species::Vector{T} where {T<:Integer},# Particle specie identifier
                                              #   (e.g. electron, proton)
        ;
        alive  ::Vector{Bool}=ones(Bool, size(pos)[2]),
        weight ::Vector{T} where {T<:Real}=ones(typeof(R[1]),
                                                size(R)[2])
        )
        return new(R, vparal, μ, species, alive, weight)
    end # constructor 
    #|
    # Given only initial conditions
    function GCAParticleSoA( 
        R0      ::Matrix{T} where {T<:Real},
        vparal0 ::Vector{T} where {T<:Real},
        μ       ::Vector{T} where {T<:Real},
        species ::Vector{T} where {T<:Integer},
        numSteps::Integer
        ;
        wfp::DataType=typeof(R0[1])
        )
        numDims, numParticles = size(R0)
        R      = zeros(wfp, numDims, numParticles, numSteps + 1)
        vparal = zeros(wfp, numParticles, numSteps + 1)
        R[:, :, 1] .= R0
        vparal[:, :, 1] .= vparal0
        alive = ones(Bool, numParticles)
        weight = ones(wfp, numParticles)
        return new(R, vparal, μ, species, alive, weight)
    end # constructor 
    #|
    function GCAParticleSoA( # Given only initial conditions
        R0      ::Vector{T} where {T<:Real},
        vparal0 ::Vector{T} where {T<:Real},
        μ       ::Vector{T} where {T<:Real},
        species ::Vector{T} where {T<:Integer},
        numSteps::Integer
        ;
        wfp::DataType=typeof(R0[1])
        )
        numParticles = 1
        numDims = length(R0)
        R      = zeros(wfp, numDims, numParticles, numSteps + 1)
        vparal = zeros(wfp, numParticles, numSteps + 1)
        R[:, :, 1] .= R0
        vparal[:, :, 1] .= vparal0
        alive = ones(Bool, numParticles)
        weight = ones(wfp, numParticles)
        return new(R, vparal, μ, species, alive, weight)
    end # constructor
end # mutable struct ParticleSoA


#----------------#
# Base functions #
#-------------------------------------------------------------------------------
"""
    Base.copy(tp::TraceParticle)
Make a deep copy of a TraceParticle-type.
"""
function Base.copy(tp::ParticleSoA)
    ParticleSoA(copy(tp.pos),
                copy(tp.vel),
                copy(tp.species),
                copy(tp.alive),
                copy(tp.weight)
                )
end # function Base.copy
#|
function Base.copy(tp::GCAParticleSoA)
    GCAParticleSoA(copy(tp.R),
                   copy(tp.vparal),
                   copy(tp.μ),
                   copy(tp.species),
                   copy(tp.alive),
                   copy(tp.weight)
                   )
end # function Base.copy
                     

function Base.Multimedia.display(tp::ParticleSoA)
    println("""Instance of mutable struct:
    Particles.ParticleSoA <: Particles.TraceParticle
        pos    ::Array{T, 3} where {T<:Real}
        vel    ::Array{T, 3} where {T<:Real}
        species::Vector{T} where {T<:Integer}   
        alive  ::Vector{Bool}
        weight ::Vector{T} where {T<:Real}
    """)
    println("................................................................")
    println("GCAPaticleSoA.pos:")
    display(tp.pos)
    println("................................................................")
    println("GCAPaticleSoA.vel:")
    display(tp.vel)
    println("................................................................")
    println("GCAPaticleSoA.species:")
    display(tp.species)
    println("................................................................")
    println("GCAPaticleSoA.alive:")
    display(tp.alive)
    println("................................................................")
    println("GCAPaticleSoA.weight:")
    display(tp.weight)
end # function Base.Multimedia.display
#|
function Base.Multimedia.display(tp::GCAParticleSoA)
    println("""Instance of mutable struct:
    Particles.ParticleSoA <: Particles.TraceParticle
        R      ::Array{T, 3} where {T<:Real} 
        vparal ::Matrix{T} where {T<:Real}   
        μ      ::Vector{T} where {T<:Real}  
        species::Vector{T} where {T<:Integer}   
        alive  ::Vector{Bool}
        weight ::Vector{T} where {T<:Real}
    """)
    println("................................................................")
    println("GCAPaticleSoA.R:")
    display(tp.R)
    println("................................................................")
    println("GCAPaticleSoA.vparal:")
    display(tp.vparal)
    println("................................................................")
    println("GCAPaticleSoA.μ:")
    display(tp.μ)
    println("................................................................")
    println("GCAPaticleSoA.species:")
    display(tp.species)
    println("................................................................")
    println("GCAPaticleSoA.alive:")
    display(tp.alive)
    println("................................................................")
    println("GCAPaticleSoA.weight:")
    display(tp.weight)
end # function Base.Multimedia.display


#------------------------#
# Particle get-functions #
#-------------------------------------------------------------------------------
function getpos(particles::ParticleSoA)
    return particles.pos
end
#|
function getpos(particles::GCAParticleSoA)
    return particles.R
end


function getvel(particles::ParticleSoA)
    return particles.vel
end
#------------------------#
# Particle set-functions #
#-------------------------------------------------------------------------------
function reset!(particles::ParticleSoA)
    n = length(particles.pos[1,1,:])
    particles.pos[:, :, 2:n] .= 0.0
    particles.vel[:, :, 2:n] .= 0.0
    particles.alive .= true
end #function reset!


function revive!(tp::TraceParticle)
    tp.alive .= true
end # function revive!


function setinitpos!(particles::ParticleSoA,
                     pos      ::Matrix{T} where {T<:Real})
    particles.pos[:, :, 1] .= pos
end # function setinitpos
#|
function setinitpos!(particles::ParticleSoA,
                     pos      ::Vector{T} where {T<:Real},
                     partIdx  ::Integer)
    particles.pos[:, partIdx, 1] .= pos
end # function setinitpos


function setinitvel!(particles::ParticleSoA,
                     vel      ::Matrix{T} where {T<:Real})
    particles.vel[:, :, 1] .= vel
end # function setinitvel
#|
function setinitvel!(particles::ParticleSoA,
                     vel      ::Vector{T} where {T<:Real},
                     partIdx  ::Integer)
    particles.vel[:, partIdx, 1] .= vel
end # function setinitvel


#---------------------------#
# Push particles a timestep #
#-------------------------------------------------------------------------------
function push!(
    tp    ::GCAParticleSoA,
    mesh  ::Mesh,
    time  ::Integer,
    dt    ::Real,
    solver::Function,
    interp::Function,
    scheme::Function,
    periodicBC::Tuple{Bool, Bool, Bool},
    )
    for j in eachindex(tp.alive)
        if tp.alive[j]
            pos = tp.R[:,j,time]
            vel = tp.vparal[j,time]
            nextpos, nextvel = solver(
                pos,
                vel,
                tp.μ[j],
                tp.species[j],
                mesh,
                dt,
                interp,
                scheme,
            )
            checkboundary!(pos, tp.alive, j, periodicBC, mesh.domain, time)
            if tp.alive[j]
                tp.R[:,j,time+1] .= nextpos
                tp.vparal[j,time+1] = nextvel
            else
                tp.R[:,j,time+1:end] .= pos
                tp.vparal[j,time+1:end] .= vel
            end
        end
    end # loop over particles
end # function push!
#|
function push!(
    tp    ::ParticleSoA,
    mesh  ::Mesh,
    time  ::Integer,
    dt    ::Real,
    solver::Function,
    interp::Function,
    scheme::Function,
    periodicBC::Tuple{Bool, Bool, Bool},
    )
    for j in eachindex(tp.alive)
        if tp.alive[j]
            pos = tp.pos[:,j,time]
            vel = tp.vel[:,j,time]
            nextpos, nextvel = solver(
                pos,
                vel,
                tp.species[j],
                mesh,
                dt,
                interp,
                scheme,
            )
            checkboundary!(nextpos, tp.alive, j, periodicBC, mesh.domain, time)
            if tp.alive[j]
                tp.pos[:,j,time+1] .= nextpos
                tp.vel[:,j,time+1] .= nextvel
            else
                tp.pos[:,j,time+1:end] .= pos
                tp.vel[:,j,time+1:end] .= vel
            end
        end
    end # loop over particles
end # function push!

function checkboundary!(
    pos       ::Vector{T} where {T<:Real},
    alive     ::Vector{Bool},
    j         ::Integer,
    periodicBC::Tuple{Bool, Bool, Bool},
    domain    ::Matrix{T} where {T<:Real},
    time      ::Integer
    )
    for k in 1:length(domain[:,1])
        if pos[k] < domain[k,1]
            if periodicBC[k] == true
                pos[k] = domain[k, 2] + (pos[k] - domain[k, 1])
            else
                alive[j] = false # kill particle
                #println("Stp $time: Particle $j killed at lower $k boundary")
                break
            end
        elseif pos[k] > domain[k, 2]
            if periodicBC[k] == true
                pos[k] = domain[k, 1] + (pos[k] - domain[k, 2])
            else
                alive[j] = false # kill particle
                #println("Step $time: Particle $j killed at upper $k boundary")
                break
            end # if particle alive
        end # if particle outside domain
    end
end # function checkboundary!

#----------------------#
# Auxiliary quantities #
#-------------------------------------------------------------------------------
function kineticenergy(particles::ParticleSoA)
    _, npart, N = size(particles.pos) 
    v = norm3(particles.vel) 
    Ek = zeros(npart, N)
    for j = 1:npart
        mass = specieTable[particles.species[j], 1]
        @. Ek[j,:] = 0.5*mass*v[j,:]^2
    end # loop j
    return Ek
end # function kineticenergy
#|
function kineticenergy(
    particles::ParticleSoA,
    mesh     ::Mesh,
    interp   ::Function
    )
    return kineticenergy(particles)
end 
#|
function kineticenergy(
    particles::GCAParticleSoA,
    mesh     ::Mesh,
    interp   ::Function
    )
    npart, N = size(particles.vparal) 
    mass = specieTable[particles.species[:], 1]
    vperp = getvperp(particles.R, particles.μ, mass, mesh, interp)
    Ek = zeros(npart, N)
    Ekparal = zeros(npart, N)
    Ekperp = zeros(npart, N)
    for j = 1:npart
        v2 = particles.vparal[j,:].^2 + vperp[j,:].^2
        @. Ek[j,:] = 0.5*mass[j]*v2
        @. Ekparal[j,:] = 0.5*mass[j]*particles.vparal[j,:]^2
        @. Ekperp[j,:] = 0.5*mass[j]*vperp[j,:]^2
    end # loop j
    return Ek, Ekparal, Ekperp
end # function kineticenergy


function getvperp(
    R     ::Array{T, 3} where {T<:Real},
    μ     ::Vector{T} where {T<:Real},
    mass  ::Vector{T} where {T<:Real},
    mesh  ::Mesh,
    interp::Function,
    )
    _, numparticles, numsteps = size(R)
    vperp = zeros(numparticles, numsteps)
    for i = 1:numsteps
        for j = 1:numparticles
            bfield, _ = gridinterp(mesh.bField, interp, R[:,j,i],
                                   mesh.xCoords, mesh.yCoords, mesh.zCoords
                                   )
            B = norm(bfield)
            vperp[j,i] = √(2μ[j]*B/mass[j])
        end # loop over j: particles
    end # loop over i: timesteps
    return vperp
end # function getvperp

"""
    computeμ(
        tp::ParticleSoA,
        mesh::Mesh,
        intep::Function,
        )
Calculates and returns μ for all particles at all time steps by
interpolating to the particle positions. This is a post-processing procudure and
is only accurate if the same interpolation method is used as in the simulation
itself.
"""
function computeμ(
    tp::ParticleSoA,
    mesh::Mesh,
    interp::Function,
    )
    m = specieTable[tp.species, 1]
    _, numparticles, numsteps = size(tp.pos)
    μ = zeros(numparticles, numsteps)
    for i = 1:numsteps
        for j = 1:numparticles
            v⃗ = tp.vel[:,j,i]
            r⃗ = tp.pos[:,j,i]
            bfield, _ = gridinterp(mesh.bField, interp, r⃗,
                                   mesh.xCoords, mesh.yCoords, mesh.zCoords
                                   )
            B = norm(bfield)
            b̂ = bfield/B
            vparal = v⃗⋅b̂
            vperp =  norm(v⃗ - vparal*b̂)
            μ[j,i] = 0.5*m[j]*vperp^2/B
        end # loop over j: particles
    end # loop over i: timesteps
    return μ
end # function computeμ

