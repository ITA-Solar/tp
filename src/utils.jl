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

"""
    Base.dropdims(arr::AbstractArray)
Extend `dropdims` to find and drop all axes with one point when `dims` keyword 
is excluded.
"""
function Base.dropdims(arr::AbstractArray)
    meshsize = size(arr)
    idxs = findall(x -> x == 1, meshsize)
    for dim in idxs
        arr = dropdims(arr, dims=dim)
    end
    return arr
end


"""
    Base.dropdims(axes::Tuple{Vararg{Vector})
Extend `dropdims` to find and drop all axes with one point when `dims` keyword 
is excluded.
"""
function Base.dropdims(axes::Tuple{Vararg{Vector}})
    mask = findall(x -> length(x) != 1, axes)
    return axes[mask]
end


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
    ;
    wfp::DataType=typeof(pos0[1])
    )
    numdims = 3
    #
    # Velocities
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    # components 
    μ = 0.0 # Expectation-value of velocity distributions. 
    # Generate velocities from a normal distribution
    velocities = μ .+ σ .* randn(MersenneTwister(seed),
                                 wfp, (numdims, numparticles))
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

