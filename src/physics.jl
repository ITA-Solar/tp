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
    kineticenergy(velocity::Vector{<:Vector}, mass)
Return the kinetic energy of a particle given the velocity *vector* and mass.
"""
function kineticenergy(velocity::Vector{<:Vector}, mass)
   kineticenergy.(norm.(velocity), mass)
end


"""
    kineticenergy(
        parallel_velocity::Vector{<:Real},
        magnetic_moment  ::Vector{<:Real},
        mass             ::Real,
        magneticfield    ::Vector{<:Real},
        electricfield    ::Vector{<:Real}
        )
Return the kinetic energy of charged particle given the `parallel_velocity` of
its *guiding centre*, its `magnetic_moment`, `mass`, and the external
electromagnetic field.
"""
function kineticenergy(
    parallel_velocity::Vector{<:Real},
    magnetic_moment  ::Vector{<:Real},
    mass             ::Real,
    magneticfield    ::Vector{<:Real},
    electricfield    ::Vector{<:Real}
    )
    v_exb, B, _ = exbdrift(magneticfield, electricfield)
    # Add more drifts if relevant
    vperp = perpendicular_velocity(magnetic_moment, mass, B)
    return 0.5mass*(vperp + parallel_velocity + norm(v_exb))^2
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
    exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
Calculate the E cross B drift in a given electromagnetic field. Also return
the magnetic field strength and mangnetic field direction.
"""
function exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
    B = norm(magneticfield) # Magnetic field strength
    b = magneticfield/B     # Magnetic field direction (unit vector)
    return (electricfield × b)/B, B, b
end


"""
    gradbdrift(
        b ::Vector{<:Real},
        ∇B::Vector{<:Real},
        B ::Real,
        μ ::Real,
        q ::Real,
        )
Calculate the ∇B-drift given the magnetic field strength `B`, the magnetic
field direction `b` (a unit vector), the gradient of the magnetic field
strength `∇B`, particle charge `q` and magnetic moment `μ`.
"""
function gradbdrift(
    b ::Vector{<:Real},
    ∇B::Vector{<:Real},
    B ::Real,
    μ ::Real,
    q ::Real,
    )
    μ/(q*B)*(b × ∇B)
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
