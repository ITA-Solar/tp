#-------------------------------------------------------------------------------
# Created 17.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Utilities.jl
#
#-------------------------------------------------------------------------------
# Module containing the utility functions.
#-------------------------------------------------------------------------------


module Utilities

# Standard libraries
using LinearAlgebra:    norm
using Random:           MersenneTwister
# Internal libraries
using Constants:        k_B

export randn
export rand
export initparticlesuniform
export initparticlesmaxwellian
export norm4


#----------------#
# Linear albebra #
#----------------#--------------------------------------------------------------
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



#----------------#
# Distributions  #
#----------------#--------------------------------------------------------------
"""
    normaldistr(
        x::Vector{T} where {T<:Real},
        μ::Real,
        σ::Real
        )

Function returning the values of `x` on a 1D normalised normal distribution with
expectation value `μ` and standard deviation `σ`.

See also [`Utilities.uniformdistr`](@ref).
"""
function normaldistr(
    x::Array{T} where {T<:Real},
    μ::Real,
    σ::Real
    )
    return @.  1/(σ*√(2π))*exp(-0.5((x - μ)/σ)^2)
end # normaldistr


"""
    uniformdistr(
        x::Vector{T} where {T<:Real},
        a::Real,
        b::Real
        )
Function returning the values of `x` on a 1D normalised unifrom distribution on
the interval [`a, `b`].

See also [`Utilities.normaldistr`](@ref).
"""
function uniformdistr(
    x::Array{T} where {T<:Real},
    a::Real,
    b::Real
    )
    mask = a .<= x .<= b
    stepheight = 1.0/(b - a)
    prob = zeros(typeof(x[1]), size(x))
    prob[mask] .= stepheight
    return prob
end # function uniformdistr


function maxwellBoltzmanndistr(
    v          ::AbstractArray{T} where {T<:Real},
    temperature::Real,
    mass       ::Real
    )
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    return @.  (1/(2π))^(3/2) *σ^(-3) * exp(-0.5(v/σ)^2) * 4π*v^2
end # normaldistr
#------------------#
# Random variables #
#-------------------------------------------------------------------------------
"""
    Base.randn(μ, σ, dims)
Method which return random variables from a normal distribution with 
expectation value `μ` and standard deviation `σ`.
"""
function Base.randn(
    μ   ::Real,
    σ   ::Real,
    dims::Tuple{Vararg{Int}}
    ;
    seed::Integer=0
    )
    return μ .+ σ.*randn(MersenneTwister(seed), dims)
end # function randn
#|
function Base.randn(
    μ   ::Real,
    σ   ::Real,
    dims...
    ;
    seed::Integer=0
    )
    return μ .+ σ.*randn(MersenneTwister(seed), dims)
end # function randn


"""
    Base.rand(a, b, dims)
Method which return random variables from a uniform distribution on the interval
[a, b].
"""
function Base.rand(
    a   ::Real,
    b   ::Real,
    dims::Tuple{Vararg{Int}}
    )
    return a .+ rand(dims...) .* (b - a)
end # function rand
#|
function Base.rand(
    a   ::Real,
    b   ::Real,
    dims...
    )
    return a .+ rand(dims...) .* (b - a)
end # function rand
#|
function Base.rand(
    domain::Matrix{T} where {T<:Real}
    )
    numaxes = size(domain)[1]
    r = zeros(numaxes)
    for i = 1:numaxes
        a = domain[i,1]
        b = domain[i,2]
        r[i] = a .+ rand() .* (b - a)
    end
    return  r
end # function rand

"""
    importancesampling(
        target  ::Function, 
        proposal::Function, 
        randgen ::Function, 
        N       ::Integer,    
        )
Sample `N` points from the `proposal`-distribution and compute the importance
weights with respect to the `target`-distribution.
"""
function importancesampling(
    target  ::Function, # Target distribution
    proposal::Function, # Proposal distruv
    randgen ::Function, # Random variable generator. Following proposal pdf.
    dims    ::Tuple{Vararg{Integer}}, # Number of samples
    )
    samples = randgen(dims)
    weights = target(samples) ./ proposal(samples)
    return samples, weights
end # function importancesampling


function rejectionsampling(
    target    ::Function,
    maxvalue  ::Real,
    numsamples::Integer,
    domain    ::Matrix{T} where {T<:Real}
    )
    numdims = size(domain)[1]
    positions = zeros((numdims, numsamples))
    accepted = 0
    rejected = 0
    while accepted < numsamples
        pos     = rand(domain)
        yguess  = rand(0.0, maxvalue, 1)
        ytarget = target(pos)[1] # Target should return a single float, but in
        # case this float is in a 1-element Vector we specify the first index.
        if yguess < ytarget
            accepted += 1
            positions[:,accepted] .= pos
        else
            rejected += 1
        end
    end
    acceptencerate = accepted/rejected
    if numdims == 1
        return positions[1,:], acceptencerate
    else
        return positions, acceptencerate
    end
end # function rejection sampling


#----------------------------------------#
# Vector potential generation            #
# Mesh generation from analytical fields #
#-------------------------------------------------------------------------------
function createaxes(
    (x0, y0, z0)::Tuple{Real, Real, Real},
    (xf, yf, zf)::Tuple{Real, Real, Real},
    (nx, ny, nz)::Tuple{Integer, Integer, Integer}
    )
    xx = collect(LinRange(x0, xf, nx))
    dx = xx[2] - xx[1]
    yy = collect(LinRange(y0, yf, ny))
    dy = yy[2] - yy[1]
    # Account for single point in the z-axis
    if nz == 1
        zz = [z0]
        dz = 0.
            else
        zz = collect(LinRange(z0, zf, nz))
        dz = zz[2] - zz[1]
    end
    return xx, yy, zz, dx, dy, dz
end # function createaxes


function discretise!(
    field::Array{T, 4} where {T<:Real},
    xx   ::Vector{T} where {T<:Real},
    yy   ::Vector{T} where {T<:Real},
    zz   ::Vector{T} where {T<:Real},
    func::Function,
    args...
    )
    for i in eachindex(xx)
        for j in eachindex(yy)
            for k in eachindex(zz)
                f = func(xx[i], yy[j], zz[k], args...)
                field[:,i,j,k] .= f
            end
        end
    end
end # function discretise
#|
function discretise!(
    field::Array{T, 3} where {T<:Real},
    xx  ::Vector{T} where {T<:Real},
    yy  ::Vector{T} where {T<:Real},
    zz  ::Vector{T} where {T<:Real},
    func::Function,
    args...
    )
    for i in eachindex(xx)
        for j in eachindex(yy)
            for k in eachindex(zz)
                f = func(xx[i], yy[j], zz[k], args...)
                field[i,j,k] = f
            end
        end
    end
end # function discretise


function mirroringfield(
    x::Real,
    y::Real,
    z::Real,
    B0::Real,
    L ::Real,
    )
    a = B0*z/L^2
    return [-x*a, -y*a, B0 + z*a]
end # mirroringfield


function dipolefield(
    x::Real,
    y::Real,
    z::Real,
    M ::Real,
    )
    a = M/(x^2 + y^2 + z^2)^(5/2)
    return [3a*z*x, 3a*z*y, a*(2z^2 - x^2 - y^2)]
end # function dipolefield

"""
    normal3Donlyz((x0, y0, z0), 
                  (xf, yf, zf), 
                  (nx, ny, nz), 
                  (μx, μy), 
                  (σx, σy),
                  amplitude
                  )
According to given spatial domain, resolution, expectationvalue, standard
deviation and amplitude, creates a vector field who's z-component is normally
distributed in x and y according to the formula fz(i, j) = amplitude *
fx(i)fy(j), where fx and fy are normal distributions in x, and y with
expectation value and std equal to μx, μy, σx, σy, respetively. 
"""
function normal3Donlyz(
    (x0, y0, z0)::Tuple{Real, Real, Real},
    (xf, yf, zf)::Tuple{Real, Real, Real},
    (nx, ny, nz)::Tuple{Integer, Integer, Integer},
    (μx, μy)    ::Tuple{Real, Real},
    (σx, σy)    ::Tuple{Real, Real},
    amplitude   ::Real=1.0
    )
    # Create spatial axes and find the grid sizes
    xx, yy, zz, dx, dy, dz = createaxes((x0, y0, z0),
                                        (xf, yf, zf),
                                        (nx, ny, nz)
                                        )
    # Initialise the vector field
    ndims = 3
    A = zeros(ndims, nx, ny, nz)
    # Evaluate the z-component of the vecor field to be normally distributed in
    # the x and y dimensions.
    gaussx = normaldistr(xx, μx, σx)
    gaussy = normaldistr(yy, μy, σy)
    for i = 1:nx
        for j = 1:ny
            A[3,i,j,:] .= amplitude * gaussx[i] * gaussy[j]
        end
    end
    return (xx, yy, zz), (dx, dy, dz), A
end # function normal3donlyz


function fadeevEquilibrium(
    (x0, y0, z0)::Tuple{Real, Real, Real},
    (xf, yf, zf)::Tuple{Real, Real, Real},
    (nx, ny, nz)::Tuple{Integer, Integer, Integer},
    λ           ::Real,
    ϵ           ::Real,
    B0          ::Real
    )
    # Create spatial axes and find the grid sizes
    xx, yy, zz, dx, dy, dz = createaxes((x0, y0, z0),
                                        (xf, yf, zf),
                                        (nx, ny, nz)
                                        )
    # Initialise the vector field
    ndims = 3
    A = zeros(ndims, nx, ny, nz)
    for i = 1:nx
        for j = 1:ny
            A[3,i,j,:] .= B0 * λ * log2( ϵ * cos(xx[i]/λ) + cosh(yy[j]/λ) )
        end
    end
    return (xx, yy, zz), (dx, dy, dz), A
    
end # function FadeevEquilibrium
#-------------------------#
# Particle initialisation #
#-------------------------------------------------------------------------------
function initparticlesuniform(
    numparticles::Integer,
    pos0        ::Vector{T} where {T<:Real}, 
    posf        ::Vector{T} where {T<:Real}, 
    vel0        ::Vector{T} where {T<:Real}, 
    velf        ::Vector{T} where {T<:Real}, 
    seed        ::Integer=0  # random-seed
    )
    numdims = 3
    spatialextent = posf .- pos0
    velocityrange = velf .- vel0
    r = rand(MersenneTwister(seed),
             typeof(pos0[1]), (Int64(2numdims), Int64(numparticles)))
             # 2 for both position and
    # velocity 
    positions  = pos0 .+ (spatialextent .* r[1:numdims, :])
    velocities = vel0 .+ (velocityrange .* r[4:2numdims, :])
    return positions, velocities
end # function initparticlesuniform

function initparticlesmaxwellian(
    numparticles::Integer,
    pos0        ::Vector{T} where {T<:Real}, 
    posf        ::Vector{T} where {T<:Real}, 
    temperature ::Real, # temperature of the Maxwellian distribution
    mass        ::Real, # mass of particles
    seed        ::Integer=0  # random-seed
    )
    numdims = 3
    #
    # Velocities
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    # components 
    μ = 0.0 # Expectation-value of velocity distributions. 
    # Generate velocities from a normal distribution
    velocities = μ .+ σ .* randn(MersenneTwister(seed),
                                 typeof(pos0[1]), (numdims, numparticles))
    #
    # Posistions: Generate from a uniform distribution
    spatialextent = posf .- pos0
    positions = pos0 .+ 
        (spatialextent .* rand(MersenneTwister(seed),
                               typeof(pos[0]), (numdims, numparticles)))
    #
    return positions, velocities
end # function initparticlesmaxwellian

function inituniform(
    numparticles::Integer,
    xbounds::Vector{T} where {T<:Real},
    ybounds::Vector{T} where {T<:Real},
    zbounds::Vector{T} where {T<:Real},
    wfp::DataType,
    ;
    seed   ::Integer=0
    )
    lowerbounds = [xbounds[1], ybounds[1], zbounds[1]]
    upperbounds = [xbounds[end], ybounds[end], zbounds[end]]
    spatialextent = upperbounds .- lowerbounds
    numdims = 3
    r = rand(MersenneTwister(seed),
             wfp, (numdims, numparticles)) # 2 for both position and
    positions  = lowerbounds .+ (spatialextent .* r[1:numdims, :])
    return convert(Matrix{wfp}, positions)
end

function initparticlesmaxwellianx(
    numparticles::Integer,
    pos0        ::Vector{T} where {T<:Real}, 
    posf        ::Vector{T} where {T<:Real}, 
    temperature ::Real, # temperature of the Maxwellian distribution
    mass        ::Real, # mass of particles
    seed        ::Integer=0  # random-seed
    ;
    wfp::DataType=typeof(pos0[1])
    )
    #
    numdims = 3
    # Velocities
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    # components 
    μ = 0.0 # Expectation-value of velocity distributions. 
    # Generate velocities from a normal distribution
    velocitiesx = μ .+ σ .* randn(MersenneTwister(seed),
                                  wfp, (numparticles))
    velocities = zeros(wfp, numdims, numparticles)
    velocities[1, :] = velocitiesx
    #
    # Posistions: Generate from a uniform distribution
    spatialextent = posf .- pos0
    positions = pos0 .+ 
        (spatialextent .* rand(MersenneTwister(seed),
                               wfp, (numdims, numparticles)))
    #
    return positions, velocities
end # function initparticlesmaxwellian


""" 
    initparticlesimsam(
        proposal    ::Function,
        randgen     ::Function,
        numparticles::Integer,
        pos0        ::Vector{T} where {T<:Real}, 
        posf        ::Vector{T} where {T<:Real}, 
        temperature ::Real, # temperature of the Maxwellian distribution
        mass        ::Real  # mass of particles
        )
Initialise particle position and velocity using importance sampling of the
Maxwellian velocity.

The position of the `numparticles` particles is uniformly distributed in the
domain defined by `pos0` and `posf`.

The the three velocity components are sampled from the `proposal` distribution
using `randgen` and given an importance weight according to the corresponding
probability in appropriate normal-distribution (defined by `temperature` and
particle `mass).
"""
function initparticlesimsam(
    proposal    ::Function,
    randgen     ::Function,
    numparticles::Integer,
    pos0        ::Vector{T} where {T<:Real}, 
    posf        ::Vector{T} where {T<:Real}, 
    temperature ::Real, # temperature of the Maxwellian distribution
    mass        ::Real  # mass of particles
    ;
    wfp::DataType=typeof(pos0[1])
    )
    #
    numdims = 3
    # Velocities
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    # components 
    μ = 0.0 # Expectation-value of velocity distributions. 
    #  Define target distribution
    targetdistr(v) = normaldistr(v, μ, σ)
    N = (numdims, numparticles)
    velocities, weights = importancesampling(targetdistr,
                                             proposal,
                                             randgen,
                                             N)
    totweight = @. weights[1,:] * weights[2,:] * weights[3,:]
    #
    # Posistions: Generate from a uniform distribution
    spatialextent = posf .- pos0
    positions = pos0 .+ 
        (spatialextent .* rand(wfp, (numdims, numparticles)))
    #
    return positions, velocities, weights, totweight
end # function initparticlesmaxwellian

end # module utilities
