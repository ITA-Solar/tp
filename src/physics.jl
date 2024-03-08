# Created 15.01.24 by eilfso
#
#
#-------------------------------------------------------------------------------
# Created 15.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                physics.jl
#
#-------------------------------------------------------------------------------
#  Contains physics of charge particles in plasma.
#-------------------------------------------------------------------------------


"""
    kineticenergy(velocity, mass)
Return the kinetic energy of a particle with `velocity` and `mass`.
"""
function kineticenergy(velocity, mass)
    0.5*mass*velocity^2
end


"""
    kineticenergy(velocity::Vector, mass)
Return the kinetic energy of a particle given the velocity *vector* and mass.
"""
function kineticenergy(velocity::Vector, mass)
   kineticenergy(norm(velocity), mass)
end


"""
    kineticenergy(
        parallel_velocity::Real,
        vperp            ::Real,
        exbdrift         ::Vector{<:Real},
        ∇Bdrift          ::Vector{<:Real},
        curvaturedrift   ::Vector{<:Real},
        polarisationdrift::Vector{<:Real},
        mass             ::Real,
        )
Return the kinetic energy of charged particle given the particle mass,
its drifts, and `parallel_velocity` and perpendicular (`vperp`)
velocity of its *guiding centre*. Other drifts are neglected
"""
function kineticenergy(
    parallel_velocity::Real,
    vperp            ::Real,
    exbdrift         ::Vector{<:Real},
    ∇Bdrift          ::Vector{<:Real},
    curvaturedrift   ::Vector{<:Real},
    polarisationdrift::Vector{<:Real},
    mass             ::Real,
    )
    vdrift = norm(exbdrift + ∇Bdrift + curvaturedrift + polarisationdrift)
    velsquared = parallel_velocity^2 + vperp^2 + norm(vdrift)^2
    return 0.5*mass*velsquared
end


"""
     kineticenergy(
         parallel_velocity::Real,
         vperp            ::Real,
         exbdrift         ::Vector{<:Real},
         mass             ::Real,
         )
Return the kinetic energy of charged particle given the particle mass,
its E cross B drift, and `parallel_velocity` and perpendicular (`vperp`)
velocity of its *guiding centre*. Drifts other than the E cross B drift
are neglected.
"""
function kineticenergy(
    parallel_velocity::Real,
    vperp            ::Real,
    exbdrift         ::Vector{<:Real},
    mass             ::Real,
    )
    velsquared = parallel_velocity^2 + vperp^2 + norm(exbdrift)^2
    return 0.5*mass*velsquared
end


"""
    kineticenergy(
        bfield           ::Vector{<:Real},
        efield           ::Vector{<:Real},
        parallel_velocity::Real,
        magneticmoment   ::Real,
        mass             ::Real,
        )
Return the kinetic energy of charged particle given the `parallel_velocity`
of its *guiding centre*, its `magnetic_moment`, its `mass`, and the external
electromagnetic field. Drifts other than the E cross B drift are neglected.
"""
function kineticenergy(
    bfield           ::Vector{<:Real},
    efield           ::Vector{<:Real},
    parallel_velocity::Real,
    magneticmoment   ::Real,
    mass             ::Real,
    )
    velsquared = parallel_velocity^2 + perpendicular_velocity(
        magneticmoment,
        mass,
        norm(bfield)
        )^2 + norm(exbdrift(bfield, efield)[1])^2
    return 0.5*mass*velsquared
end


"""
Kinetic energy of a charged particle based on its guiding centre, parallel
velocity, magnetic moment, and the electromagnetic field in which it is
embedded. Accounts for E cross B drift, magnetic gradient drift, curvature
drift, and polarisationdrift.
"""
function kineticenergy(
    R     ::Vector{<:Real},
    vparal::Real,
    q     ::Real, # Particle charge
    m     ::Real, # Particle mass
    μ     ::Real, # Particle magnetic moment
    itpvec::Vector{<:AbstractInterpolation}, # electromagnetic field interpolators
    )
    vperp, drifts = drifts(R, vparal, q, m, μ, itpvec)
    kineticenergy(vparal, vperp, drifts...)
end

"""
    perpendicular_velocity(magnetic_moment, mass, magneticfieldstrength)
Return the perpendicular velocity of a charged particle in a magnetic field.
"""
function perpendicular_velocity(magnetic_moment, mass, magneticfieldstrength)
    sqrt(2magnetic_moment*magneticfieldstrength/mass)
end


"""
    magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
Return the magnetic moment of a chargred particle in a magnetic field.
"""
function magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
    0.5mass*perpendicular_velocity^2/magneticfieldstrength
end


"""
exbdrift(
    b̂::Vector{<:Real},
    E_vec::Vector{<:Real},
    B_inv::Real,
    )
"""
function exbdrift(
    b̂::Vector{<:Real},
    E_vec::Vector{<:Real},
    B_inv::Real,
    )
    return (E_vec × b̂)*B_inv
end


"""
    exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
Calculate the E cross B drift in a given electromagnetic field. Also return
the magnetic field strength and mangnetic field direction.
"""
function exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
    B = norm(magneticfield) # Magnetic field strength
    B_inv = 1/B
    b̂ = magneticfield*B_inv     # Magnetic field direction (unit vector)
    return exbdrift(b̂, electricfield, B_inv), B, B_inv, b̂
end



"""
    gradbdrift(
        b̂    ::Vector{<:Real},
        ∇B   ::Vector{<:Real},
        B_inv::Real,
        μ    ::Real,
        q_inv::Real,
        )
Calculate the ∇B-drift given the inverse of the magnetic field strength `B_inv`,
the magnetic field direction `b̂` (a unit vector), the gradient of the magnetic
field strength `∇B`, inverse of particle charge `q_inv` and magnetic moment `μ`.
"""
function gradbdrift(
    b̂    ::Vector{<:Real}, # The direction of the magnetic fielda
    ∇B   ::Vector{<:Real}, # The gradient of the magnetic field strength
    B_inv::Real, # The inverse of the magnetic field strength
    μ    ::Real, # the magnetic moment of the particle
    q_inv::Real  # The inverse of the charge of the particle
    )
    return q_inv*B_inv*μ*(b̂ × ∇B)
end


"""
    curvaturedrift(
        b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
        db̂dt  ::Matrix{<:Real},# Material derivative of the magn. field direction
        vparal::Real, # Particle velocity component parallel to the magnetic field
        B_inv ::Real, # The inverse of the magnetic field strength
        q_inv ::Real, # The inverse of the charge of the particle
        mass  ::Real  # The particle mass
        )
Calculate the curvature drift of a charged particle in an electromagnetic field.
"""
function curvaturedrift(
    b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
    db̂dt  ::Vector{<:Real},# Material derivative of the magn. field direction
    vparal::Real, # Particle velocity component parallel to the magnetic field
    B_inv ::Real, # The inverse of the magnetic field strength
    q_inv ::Real, # The inverse of the charge of the particle
    mass  ::Real  # The particle mass
    )
    return q_inv*B_inv*mass*b̂ × (vparal*db̂dt)
end


"""
    polarisationdrift(
        b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
        dExBdt::Matrix{<:Real},# Material derivative of the E cross B drift
        B_inv ::Real, # The inverse of the magnetic field strength
        q_inv ::Real, # The inverse of the charge of the particle
        mass  ::Real  # The particle mass
        )
Calculate the polarisation drift of a charged particle in an electromagnetic
field.
"""
function polarisationdrift(
    b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
    dExBdt::Vector{<:Real},# Material derivative of the E cross B drift
    B_inv ::Real, # The inverse of the magnetic field strength
    q_inv ::Real, # The inverse of the charge of the particle
    mass  ::Real  # The particle mass
    )
    return q_inv*B_inv*mass*b̂ × dExBdt
end


"""
    drifts_and_vperp(
        R     ::Vector{<:Real},
        vparal::Real,
        q     ::Real, # Particle charge
        m     ::Real, # Particle mass
        μ     ::Real, # Particle magnetic moment
        itpvec::Vector{<:AbstractInterpolation}, # electromagnetic field interpolators
        )
Evaluate the drifts of a particle at location `R` in an electromagnetic field
given by a vector of interpolation objects.
"""
function drifts_and_vperp(
    R     ::Vector{<:Real},
    vparal::Real,
    q     ::Real, # Particle charge
    m     ::Real, # Particle mass
    μ     ::Real, # Particle magnetic moment
    itpvec::Vector{<:AbstractInterpolation}, # electromagnetic field interpolators
    )
    # NOTE: Comments are copied from 3D implementation. Running times are
    # not correct in this case since this is 2D version.
    #
    q_inv = 1/q # Inverse of q - to replace division with multiplication
    # Use the gyrocentre position interpolate the vectors
    # scalars from the interpolation objects.
    B_vec = [itpvec[i](R...) for i in 1:3]
    E_vec = [itpvec[i](R...) for i in 4:6]
    ExBdrift, B, B_inv, b = exbdrift(B_vec, E_vec)

    # Calculate the gradient of the magnetic field strength
    ∇B = ForwardDiff.gradient(R) do x
        sqrt(itpvec[1](x...)^2 + itpvec[2](x...)^2 + itpvec[3](x...)^2)
    end
    ∇B = [∇B[1], 0f0, ∇B[2]]
    #
    # Calculate the gradient of the magnetic field direction
    jacobian_matrix = stack(
        [Interpolations.gradient(itp, R...) for itp in itpvec],
        dims=1
        )
    # Add zeros-column representing derivatives along the y-axis
    jacobian_matrix = [
        jacobian_matrix[:,1];;
        jacobian_matrix[:,2];;
        jacobian_matrix[:,3]
        ]
    ∇B_vec = jacobian_matrix[1:3,:]
    ∇b = (∇B_vec - b * ∇B')*B_inv
    #
    # Calculate the Jacobian matrix of the ExB-drift
    ∇E_vec = jacobian_matrix[4:6,:]
    skewE = skewsymmetric_matrix(E_vec)
    skewb = skewsymmetric_matrix(b)
    ∇ExB = (-skewb*∇E_vec + skewE*∇b - ExBdrift * ∇B')*B_inv

    # Total time derivatives. Assumes ∂/∂t = 0,
    dbdt = vparal * (∇b * b) + ∇b*ExBdrift
    dExBdt = vparal * (∇ExB * b) + ∇ExB*ExBdrift

    # Calculate other drifts
    ∇Bdrift = gradbdrift(b, ∇B, B_inv, μ, q_inv)
    Rdrift = curvaturedrift(b, dbdt, vparal, B_inv, q_inv, m)
    Pdrift = polarisationdrift(b, dExBdt, B_inv, q_inv, m)
    # Compute the perpendicular velocity
    vperp = √(2B*μ/m)

    return vperp, [ExBdrift, ∇Bdrift, Rdrift, Pdrift]
end


function drifts_and_vperp_2Dxz(
    R     ::Vector{<:Real},
    vparal::Real,
    q     ::Real, # Particle charge
    m     ::Real, # Particle mass
    μ     ::Real, # Particle magnetic moment
    itpvec::Vector{<:AbstractInterpolation}, # electromagnetic field interpolators
    )
    # NOTE: Comments are copied from 3D implementation. Running times are
    # not correct in this case since this is 2D version.
    #
    Rx, Rz = R[1], R[3]
    # Extract parameters
    q_inv = 1/q # Inverse of q - to replace division with multiplication
    # Use the gyrocentre position interpolate the vectors
    # scalars from the interpolation objects.
    B_vec = [itpvec[i](Rx, Rz) for i in 1:3]
    E_vec = [itpvec[i](Rx, Rz) for i in 4:6]
    B = norm(B_vec)   # The magnetic field strength
    B_inv = 1/B       # Inverse of B - to replace divition with multiplication
    b = B_vec*B_inv   # An unit vector pointing in the direction of the
                      #  magnetic field
    ExBdrift = (E_vec × b)/B # The E cross B-drift

    # Calculate the gradient of the magnetic field strength
    ∇B = ForwardDiff.gradient([Rx, Rz]) do x
        sqrt(itpvec[1](x...)^2 + itpvec[2](x...)^2 + itpvec[3](x...)^2)
    end
    ∇B = [∇B[1], 0f0, ∇B[2]]
    #
    # Calculate the gradient of the magnetic field direction
    jacobian_matrix = stack(
        [Interpolations.gradient(itp, Rx, Rz) for itp in itpvec],
        dims=1
        )
    # Add zeros-column representing derivatives along the y-axis
    jacobian_matrix = [
        jacobian_matrix[:,1];;
        zeros(typeof(Rx), 6);;
        jacobian_matrix[:,2]
        ]
    ∇B_vec = jacobian_matrix[1:3,:]
    ∇b = (∇B_vec - b * ∇B')*B_inv
    #
    # Calculate the Jacobian matrix of the ExB-drift
    ∇E_vec = jacobian_matrix[4:6,:]
    skewE = skewsymmetric_matrix(E_vec)
    skewb = skewsymmetric_matrix(b)
    ∇ExB = (-skewb*∇E_vec + skewE*∇b - ExBdrift * ∇B')*B_inv

    # Total time derivatives. Assumes ∂/∂t = 0,
    dbdt = vparal * (∇b * b) + ∇b*ExBdrift
    dExBdt = vparal * (∇ExB * b) + ∇ExB*ExBdrift

    # Calculate other drifts
    ∇Bdrift = gradbdrift(b, ∇B, B_inv, μ, q_inv)
    Rdrift = curvaturedrift(b, dbdt, vparal, B_inv, q_inv, m)
    Pdrift = polarisationdrift(b, dExBdt, B_inv, q_inv, m)
    # Compute the perpendicular velocity
    vperp = √(2B*μ/m)

    return vperp, [ExBdrift, ∇Bdrift, Rdrift, Pdrift]
end


"""
    get_guidingcentre(
        pos          ::Vector{<:Real},
        vel          ::Vector{<:Real},
        magneticfield::Vector{<:Real},
        electricfield::Vector{<:Real},
        charge       ::Real,
        mass         ::Real
        )
Calculates and returns the guiding-centre position, parallel velocity, and
 magnetic moment of a charged particle in a electromagnetic field. 
"""
function get_guidingcentre(
    pos          ::Vector{<:Real},
    vel          ::Vector{<:Real},
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    charge       ::Real,
    mass         ::Real
    )
    ExBdrift, B, b_vec = exbdrift(magneticfield, electricfield)
    vel_in_E_frame = vel - ExBdrift
    # Calculate the guiding centre posistion 
    R = pos - mass/(charge*B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal*b_vec
    mu = magneticmoment(norm(vperp), mass, B)
    return R, vparal, mu
end


"""
    get_fullorbit(
        magneticfield::Vector{<:Real},
        electricfield::Vector{<:Real},
        R            ::Vector{<:Real},
        vparal       ::Real,
        μ            ::Real,
        mass         ::Real,
        charge       ::Real,
        phaseangle   ::Real,
        )
Get the position and velocity from the guiding centre position `R`, parallel
velocity `vparal`, and magnetic moment `μ` of a charged particle with `mass`
and `charge` in an electromagnetic field.

Requires an arbitrary `phaseangle`
of the gyromotion because we are going from 5 parameters to 6. The phase angle
is with respect to vector perpendicular to the magnetic field and guiding
centre position vector.
"""
function get_fullorbit(
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    R            ::Vector{<:Real}, # Guiding centre position
    vparal       ::Real,           # Velocity parallel to the magnetic field
    μ            ::Real,           # Magnetic moment of particle
    mass         ::Real,
    charge       ::Real,
    phaseangle   ::Real,           # Arbitrary phase angle of gyration.
    )
    v_exb, B, b = exbdrift(magneticfield, electricfield) # get E cross B drift,
    # magnetic field strength and magnetic field direction
    vperp = perpendicular_velocity(μ, mass, B)
    larmorradius = mass*vperp/(charge*B)
    e₁ = (R × b)/norm(R)
    e₂ = e₁ × b
    cθ = cos(phaseangle)
    sθ = sin(phaseangle)
    position = R + larmorradius*(cθ*e₁ + sθ*e₂)
    velocity = vparal*b + v_exb + vperp*(cθ*e₂ - sθ*e₁)
    return [position; velocity]
end


"""
    cosineof_pitchangle(
        magneticfield    ::Vector{<:Real},
        parallel_velocity::Real,
        mass             ::Real,
        magneticmoment   ::Real,
        )
Calculates the cosine of the pitch angle of a charge particle in a magnetic
field, given the particle's parallel velocity, mass and magnetic moment.
"""
function cosineof_pitchangle(
    magneticfield    ::Vector{<:Real},
    parallel_velocity::Real,
    mass             ::Real,
    magneticmoment   ::Real,
    )
    B = norm(magneticfield)
    return sqrt(1/(2B*magneticmoment/(mass*parallel_velocity^2) + 1))
end


"""
    maxwellboltzmanndistr(
    v          ::AbstractArray{<:Real},
    temperature::Real,
    mass       ::Real
    )
Return the probability of the velocity `v` of a Maxwell-Boltzmann 
distribution with given particle `temperature` and `mass`.
"""
function maxwellboltzmanndistr(
    v          ::AbstractArray{<:Real},
    temperature::Real,
    mass       ::Real
    )
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    return @.  (1/(2π))^(3/2) *σ^(-3) * exp(-0.5(v/σ)^2) * 4π*v^2
end


"""
    magneticmirrorfield(
        x::Real,
        y::Real,
        z::Real,
        B0::Real,
        L ::Real
        )
Return the magnetic field at (`x`, `y`, `z`) in a magnetic mirror/bottle with 
length `L` and strength `B0`.

B⃗(x, y, z) = -xzB₀/L² x̂  - yzB₀/L² ŷ  + (B₀ + zB₀/L²)ẑ
"""
function magneticmirrorfield(
    x::Real,
    y::Real,
    z::Real,
    B0::Real,
    L ::Real,
    )
    a = B0*z/L^2
    return [-x*a, -y*a, B0 + z*a]
end # mirroringfield


"""
    magneticdipolefield(
        x::Real,
        y::Real,
        z::Real,
        M ::Real
        )
Return the magnetic field at (`x`, `y`, `z`) in a magnetic dipole field with 
mangetic moment `M`.

B⃗(x, y, z) = 3αzx x̂  + 3αzy ŷ  + α(2z² - x² - y²)ẑ, 
where 
α = M/(x² + y² + z²)^(2.5).
"""
function magneticdipolefield(
    x::Real,
    y::Real,
    z::Real,
    M ::Real,
    )
    a = M/(x^2 + y^2 + z^2)^(5/2)
    return [3a*z*x, 3a*z*y, a*(2z^2 - x^2 - y^2)]
end


"""
    fadeevEquilibrium(
        (x0, y0, z0)::Tuple{Real, Real, Real},
        (xf, yf, zf)::Tuple{Real, Real, Real},
        (nx, ny, nz)::Tuple{Integer, Integer, Integer},
        λ           ::Real,
        ϵ           ::Real,
        B0          ::Real
        )
"""
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


"""
    gcagradients_from_emfield(
        bField ::Array{T, 4} where {T<:Real},
        eField ::Array{T, 4} where {T<:Real},
        xCoords::Vector{T} where {T<:Real},
        yCoords::Vector{T} where {T<:Real},
        zCoords::Vector{T} where {T<:Real},
        scheme ::Function
        )
Compute the gradients used in the guiding centre approximation from an
electromagnetic field. The gradients are:
- The gradient of the magnetic field strength
- the gradient of the magnetic field direction
- the gradient of the E cross B drift.
"""
function gcagradients_from_emfield(
    bField ::Array{T, 4} where {T<:Real},
    eField ::Array{T, 4} where {T<:Real},
    xCoords::Vector{T} where {T<:Real},
    yCoords::Vector{T} where {T<:Real},
    zCoords::Vector{T} where {T<:Real},
    scheme ::Function
    )
    wfp = typeof(bField[1])
    nx, ny, nz, ncomp = size(bField)
    BB = norm4(bField, axis=4)
    b̂ = zeros(wfp, nx, ny, nz, ncomp)
    ExBdrift = zeros(wfp, nx, ny, nz, ncomp)
    # I guess this could be done much more efficiently
    for i = 1:nx
        for j= 1:ny
            for k = 1:nz
                B⃗ = bField[i,j,k,:]
                E⃗ = eField[i,j,k,:]
                B = BB[i,j,k]
                b̂[i,j,k,:]  .= B⃗ ./ B
                ExBdrift[i,j,k,:] .= (E⃗ × B⃗) ./ B^2
            end
        end
    end
    ∇B = ∇(BB, xCoords, yCoords, zCoords, scheme)
    ∇b̂ = ∇(b̂,  xCoords, yCoords, zCoords, scheme)
    ∇ExBdrift= ∇(ExBdrift, xCoords, yCoords, zCoords, scheme)
    return ∇B, ∇b̂, ∇ExBdrift
end




