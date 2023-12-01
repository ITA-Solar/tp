#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Solvers.jl
#
#-------------------------------------------------------------------------------
# Module containing the various solvers.
#-------------------------------------------------------------------------------

"""
    fullOrbit(pos, vel, specie, bField, eField, dt, scheme)

Solves the Lorentz equation of motion using an arbitrary numerical scheme
(defined by the argument `scheme`).
"""
function fullOrbit_interstaticfield(
    pos         ::Vector{T} where {T<:Real},
    v           ::Vector{T} where {T<:Real}, # velocity
    specie      ::Integer,
    mesh        ::Mesh,
    dt          ::Real,
    interpolator::Function,
    scheme      ::Function
    )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    # Interpolate fields
    fields, _ = gridinterp(mesh,
                           interpolator,
                           pos)
    B = fields[1]
    E = fields[2]
    #
    statevector = [pos; v]
    statevectorNext = scheme(statevector,
                             dt, 
                             eomLorentzforce,
                             B, E, q, m)
    return statevectorNext[1:3], statevectorNext[4:6]
end # funcion fullOrbit_interstaticfield

function fullOrbit(pos        ::Vector{T} where {T<:Real},
                   v           ::Vector{T} where {T<:Real}, # velocity
                   specie      ::Integer,
                   mesh        ::Mesh,
                   dt          ::Real,
                   interpolator::Function,
                   scheme      ::Function
                   )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    #
    statevector = [pos; v]
    statevectorNext = scheme(statevector,
                             dt, 
                             eomLorentzforce,
                             q, m, mesh, interpolator,
                             )
    return statevectorNext[1:3], statevectorNext[4:6]
end # function fullOrbit


function eomLorentzforce(
    s::Vector{T} where {T<:Real}, # The state vector
    B::Vector{T} where {T<:Real}, # The magnetic field
    E::Vector{T} where {T<:Real}, # The electric field
    q::Real,         # Charge
    m::Real          # Mass
    )
    x = s[1:3] # The position vector
    v = s[4:6] # The velocity vector
    dvdt = q/m * (E + v × B) 
    dxdt = v
    dsdt = [dxdt; dvdt]
    return dsdt
end # function eomLorentzforce
#|
function eomLorentzforce(
    statevector ::Vector{T} where {T<:Real}, # The state vector
    q           ::Real,         # Charge
    m           ::Real,         # Mass
    mesh        ::Mesh,            # The mesh containing the magnetic field
    interpolator::Function         # Interpolation function used for evaluating
                                   #   the field at the stavector-location
    )
    x = statevector[1:mesh.numdims] # The position vector
    v = statevector[mesh.numdims+1:2mesh.numdims] # The velocity vector
    # Interpolate fields
    fields, _ = gridinterp(mesh,
                           interpolator,
                           x)
    B = fields[1]
    E = fields[2]
    dvdt = q/m * (E + v × B)
    dxdt = v
    dsdt = [dxdt; dvdt]
    return dsdt
end # function eomLorentzforce


function fieldtracingforward(
    statevector ::Vector{T} where {T<:Real}, # Should be just position
    vectorfield ::Array{T, 4} where {T<:Real},
    interpolator::Function,
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real}
    )
    interpfield, _ = gridinterp(vectorfield, interpolator, statevector,
                                xx, yy, zz)
    fieldstrength = norm(interpfield)
    fielddirection = interpfield ./ fieldstrength
    dsdt = fielddirection
    return dsdt
end # function fieldtracing

function fieldtracingbackward(
    statevector ::Vector{T} where {T<:Real}, # Should be just position
    vectorfield ::Array{T, 4} where {T<:Real},
    interpolator::Function,
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real}
    )
    interpfield, _ = gridinterp(vectorfield, interpolator, statevector,
                                xx, yy, zz)
    fieldstrength = norm(interpfield)
    fielddirection = interpfield ./ fieldstrength
    dsdt = -fielddirection
    return dsdt
end # function fieldtracing
    

function relFullOrbitExplLeapFrog(pos         ::Vector{T} where {T<:Real},
                                  vel         ::Vector{T} where {T<:Real}, 
                                  specie      ::Integer,
                                  mesh        ::Mesh,
                                  dt          ::Real,
                                  interpolator::Function,
                                  scheme      ::Function
                                  )
    # Extract particle mass and charge
    mass  = specieTable[specie, 1]
    charge = specieTable[specie, 2]

    #
    # Step 1: Evaluate half-step in time for position
    #
    posHalf = positionHalfStep(pos, vel, dt)
    # Interpolate fields to this location
    fields, _ = gridinterp(mesh,
                  interpolator,
                  posHalf)
    bField = fields[1]
    eField = fields[2]

    # 
    # Step 2: Evaluate full time step in velocity, which is shceme-dependent.
    #
    velNext = scheme(vel, bField, eField, mass, charge, dt)
    
    #
    # Step 3: Evaluate second half of time step in position
    # 
    posNext = positionHalfStep(posHalf, velNext, dt)

    return posNext, velNext
end # function relFullOrbitExpLeapFrog


function gca(pos         ::Vector{T} where {T<:Real},
             vel         ::Vector{T} where {T<:Real},
             specie      ::Integer,
             mesh        ::Mesh,
             dt          ::Real,
             interpolator::Function,
             scheme      ::Function
             )
    # Interpolate fields to this location
    fields, _ = gridinterp(mesh, interpolator, pos)
    # i, j k are cell corner indexes in mesh. Corresponding to the 
    #  position of the particle.
    bField = fields[1] # The magnetic field
    eField = fields[2] # The electric field
    ∇B = fields[3] # The gradient of the magnetic field.
    

    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    
    vperp  = vel[4]    # Particle velocity perpendicular to the magne. field
    vparal = vel[5]    # Particle velocity parallel to the magnetic field
    μ  = vel[6]      # The magnetic moment
    B = norm(bField) # The magnetic field strength
    b̂ = bField/B     # An unit vector pointing in the direction of the
                       #  magnet field
    
    # Electric field component parallel to the magnetic field
    Eparal = eField ⋅ b̂ 
    
    # Compute the acceleration 
    accparal = (q*Eparal - μ*b̂⋅∇B)/m # along the magnetic field lines
    acc = accparal*b̂              # The vector
    # With spatially changing fields, the velocity at this point will not be the
    # same as the last, independent of time.
    velHere = vparal*b̂ + b̂/B × (-c*eField + μ*c/q*∇B)
    
    # Use integration scheme to find velocities at the next time step
    vparalnext    = scheme(vparal, accparal, dt)
    posNext, v = scheme(pos, velHere, acc, dt) # Use v∥next? Will
    # essentially be used if the scheme is euler cromer since a∥ is in acc.
    # norm(velNext) - norm(velhere) should equal v∥next
    #or just? posNext, v = scheme(pos, vel[1:3], acc, dt)
    
    # Compute some auxiliary quantities
    vperpnext = √(norm(v)^2 - vparalnext^2) 
    μnext = m*vperpnext^2/(2B) #  (should be constant for all times)
    # Maybe μ should be forced constant and kept as a parameter to this solver.
    #   This would require a change in the implementation of solvers and
    #   Patch.run!, where e.g. the particle type is passed to solver. Or that
    #   run! is passed with the particle type, not the patch, such that one may
    #   define different run methods depending on the particle type. 

    velNext = [v[1], v[2], v[3], vperpnext, vparalnext, μnext]
    return posNext, velNext
end # function GCA

function gca(
    pos         ::Vector{T} where {T<:Real},
    vel         ::Real,
    μ           ::Real,
    specie      ::Integer,
    mesh        ::Mesh,
    dt          ::Real,
    interpolator::Function,
    scheme      ::Function
    )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    #
    statevector = [pos; vel]
    statevectorNext = scheme(statevector,
                             dt, 
                             eom_gca,
                             q, m, μ, mesh, interpolator
                             )
    return statevectorNext[1:3], statevectorNext[4]
    
end # function GCA

function eom_gca(
    statevector ::Vector{T} where {T<:Real},
    q           ::Real,
    m           ::Real,
    μ           ::Real,
    mesh        ::Mesh,
    interpolator::Function
    )
    R      = statevector[1:3]
    vparal = statevector[4] # Particle velocity parallel to the magnetic field

    # Interpolate fields to this location
    fields, _ = gridinterp(mesh, interpolator, R)
    # i, j k are cell corner indexes in mesh. Corresponding to the 
    #  position of the particle.
    B⃗ = fields[1] # The magnetic field
    E⃗ = fields[2] # The electric field
    ∇B = fields[3]     # The gradient of the magnetic field.
    ∇b̂ = fields[4]
    ∇ExB = fields[5]
    local B = norm(B⃗)   # The magnetic field strength
    b̂ = B⃗/B       # An unit vector pointing in the direction of the
                       #  magnetic field
    # Electric field component parallel to the magnetic field
    Eparal = E⃗⋅b̂ 
    # Calculate drifts
    ExBdrift = (E⃗ × b̂)/B
    ∇Bdrift = μ/(q*B)*(b̂ × ∇B)
    # Total time derivatives. Assumes ∂/∂t = 0,
    db̂dt = vparal * (∇b̂ * b̂) + ∇b̂*ExBdrift
    dExBdt = vparal * (∇ExB * b̂) + ∇ExB*ExBdrift
    
    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m*b̂/(q*B) × (vparal*db̂dt + dExBdt)
    #dRperpdt = b̂/B × (-E⃗ + μ/q*∇B + m/q*(vparal*db̂dt + dExBdt))  # old

    # Compute the acceleration 
    dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m # along the magnetic field lines
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the velocity
    dRdt = vparal*b̂ + dRperpdt
    
    # How to store the perpendicular velocity? Would need to know b̂ at
    #   each R calculate vperp at each R. Could be interesting to store this
    #   as an auxiliary variable somehow, to se how the drift evolves.
    dsdt = [dRdt; dvparaldt]
    return dsdt
end # function GCA


struct LorentzForce
end
function (lf::LorentzForce)(du, u, p, t)
    q, m, bx, by, bz, ex, ey, ez = p
    x = u[1:3] # The position vector
    v = u[4:6] # The velocity vector
    B = [bx(x...), by(x...), bz(x...)]
    E = [ex(x...), ey(x...), ez(x...)]
    dvdt = q/m * (E + v × B) 
    dxdt = v
    du[:] = [dxdt; dvdt]
end


struct GCA
end
function (gca::GCA)(du, u, p, _)
    R = u[1:3]
    vparal = u[4] # Particle velocity parallel to the magnetic field
    # Extract parameters
    q, m, μ, B_itp, E_itp, gradB_itp, gradb_itp, gradExB_itp = p
    B_vec = [B_itp[1](R...), B_itp[2](R...), B_itp[3](R...)]
    E_vec = [E_itp[1](R...), E_itp[2](R...), E_itp[3](R...)]
    gradB_vec = [ gradB_itp[1](R...), gradB_itp[2](R...), gradB_itp[3](R...) ] 
    gradb = [ 
        gradb_itp[1,1](R...) gradb_itp[1,2](R...) gradb_itp[1,3](R...)
        gradb_itp[2,1](R...) gradb_itp[2,2](R...) gradb_itp[2,3](R...)
        gradb_itp[3,1](R...) gradb_itp[3,2](R...) gradb_itp[3,3](R...)
        ]
    gradExB = [
        gradExB_itp[1,1](R...) gradExB_itp[1,2](R...) gradExB_itp[1,3](R...)
        gradExB_itp[2,1](R...) gradExB_itp[2,2](R...) gradExB_itp[2,3](R...)
        gradExB_itp[3,1](R...) gradExB_itp[3,2](R...) gradExB_itp[3,3](R...) 
        ]
    B = norm(B_vec)   # The magnetic field strength
    b_vec = B_vec/B       # An unit vector pointing in the direction of the
                       #  magnetic field
    # Electric field component parallel to the magnetic field
    Eparal = E_vec⋅b_vec
    # Calculate drifts
    ExBdrift = (E_vec × b_vec)/B
    ∇Bdrift = μ/(q*B)*(b_vec × gradB_vec)
    # Total time derivatives. Assumes ∂/∂t = 0,
    dbdt = vparal * (gradb * b_vec) + gradb*ExBdrift
    dExBdt = vparal * (gradExB * b_vec) + gradExB*ExBdrift
    
    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m*b_vec/(q*B) × (vparal*dbdt + dExBdt)
    #dRperpdt = b̂/B × (-E⃗ + μ/q*∇B + m/q*(vparal*db̂dt + dExBdt))  # old

    # Compute the acceleration 
    dvparaldt = (q*Eparal - μ*b_vec⋅gradB_vec)/m # along the magnetic field lines
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the velocity
    dRdt = vparal*b_vec+ dRperpdt
    
    # How to store the perpendicular velocity? Would need to know b̂ at
    #   each R calculate vperp at each R. Could be interesting to store this
    #   as an auxiliary variable somehow, to se how the drift evolves.
    du[:] = [dRdt; dvparaldt]
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic differential equations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct ConstantAdvection
end
function (eom::ConstantAdvection)(du, _, p, _)
    du .= p
end
function (eom::ConstantAdvection)(_, p, _)
    return p
end

struct ConstantDiffusion
end
function (eom::ConstantDiffusion)(du, _, p, _)
    du .= p
end
function (eom::ConstantDiffusion)(_, p, _)
    return p
end

struct ProportionalAdvection
end
function (eom::ProportionalAdvection)(du, u, p, _)
    du .= p*u
end
function (eom::ProportionalAdvection)(u, p, _)
    return p*u
end

struct ProportionalDiffusion
end
function (eom::ProportionalDiffusion)(du, u, p, _)
    du .= p*u
end
function (eom::ProportionalDiffusion)(u, p, _)
    return p*u
end

struct DuffingVanDerPolOscillatorDrift
end
function (eom::DuffingVanDerPolOscillatorDrift)(du, u, p, _)
    X, V = u
    dXdt = V
    dVdt = (X*(p - X^2) - V)
    du .= [dXdt; dVdt]
end
function (eom::DuffingVanDerPolOscillatorDrift)(u, p, _)
    X, V = u
    dXdt = V
    dVdt = (X*(p - X^2) - V)
    return [dXdt; dVdt]
end

struct DuffingVanDerPolOscillatorDiffusion
end
function (eom::DuffingVanDerPolOscillatorDiffusion)(du, u, p, _)
    X, _ = u
    dXdW = 0.0
    dVdW = p*X
    du .= [dXdW; dVdW]
end
function (eom::DuffingVanDerPolOscillatorDiffusion)(u, p, _)
    X, _ = u
    dXdW = 0.0
    dVdW = p*X
    return [dXdW; dVdW]
end

struct OrnsteinUhlenbeckDrift
end
function (eom::OrnsteinUhlenbeckDrift)(du, u, p, _)
    _, V = u
    dXdt = V
    dVdt = -V/p
    du .= [dXdt, dVdt]
end
function (eom::OrnsteinUhlenbeckDrift)(u, p, _)
    _, V = u
    dXdt = V
    dVdt = -V/p
    return [dXdt, dVdt]
end

struct OrnsteinUhlenbeckDiffusion
end
function (eom::OrnsteinUhlenbeckDiffusion)(du, _, p, _)
    dXdW = 0.0
    dVdW = p
    du .= [dXdW, dVdW]
end
function (eom::OrnsteinUhlenbeckDiffusion)(_, p, _)
    dXdW = 0.0
    dVdW = p
    return [dXdW, dVdW]
end

struct GCA_pitchAngleFriction
end
function (eom::GCA_pitchAngleFriction)(du, u, p, _)
    R = u[1:3]    # Position of the gyrocentre
    vparal = u[4] # Particle velocity parallel to the magnetic field
    beta = u[5]   # Cosine of pitch angle

    # Extract parameters:
    #
    # Particle charge and particle mass and the Coulomb logarithm
    q, m, coulomb_logarithm = p[1:3]
    # interpolation objects of
    # the magnetic vector field, the electric vector field, the magnetic
    # gradient vector field, the gradient of the magnetic direction
    # (a 3rd order tensor field), the gradient of the ExB-drift (a 3rd
    # order tensor field), the number density of electrons (a scalar
    # field), and the temperature (a scalar field)
    B_itp, E_itp, gradB_itp, gradb_itp, gradExB_itp, n_itp, T_itp = p[4:10]

    # Use the gyrocentre position interpolate the vectors, matrices and
    # scalars from the interpolation objects.
    B_vec = [B_itp[1](R...), B_itp[2](R...), B_itp[3](R...)]
    E_vec = [E_itp[1](R...), E_itp[2](R...), E_itp[3](R...)]
    gradB_vec = [ gradB_itp[1](R...), gradB_itp[2](R...), gradB_itp[3](R...) ]
    gradb = [
        gradb_itp[1,1](R...) gradb_itp[1,2](R...) gradb_itp[1,3](R...)
        gradb_itp[2,1](R...) gradb_itp[2,2](R...) gradb_itp[2,3](R...)
        gradb_itp[3,1](R...) gradb_itp[3,2](R...) gradb_itp[3,3](R...)
        ]
    gradExB = [
        gradExB_itp[1,1](R...) gradExB_itp[1,2](R...) gradExB_itp[1,3](R...)
        gradExB_itp[2,1](R...) gradExB_itp[2,2](R...) gradExB_itp[2,3](R...)
        gradExB_itp[3,1](R...) gradExB_itp[3,2](R...) gradExB_itp[3,3](R...)
        ]
    n = n_itp(R...)
    T = T_itp(R...)

    # Calculate some quantities
    B = norm(B_vec)   # The magnetic field strength
    b_vec = B_vec/B   # An unit vector pointing in the direction of the
                      #  magnetic field
    mu = m*vparal^2/2B*(1/beta^2 - 1) # The magnetic moment
    Eparal = E_vec⋅b_vec # Electric field component parallel to the
                         # magnetic field

    # Calculate drifts
    ExBdrift = (E_vec × b_vec)/B
    ∇Bdrift = mu/(q*B)*(b_vec × gradB_vec)

    # Material derivatives of the magnetic field strength, the magnetic field
    # direction, and the ExB-drift. Assumes ∂/∂t = 0
    # See Ripperda et al. (2018) and notes.
    dBdt = vparal * gradb ⋅ gradB_vec + ExBdrift ⋅ gradB_vec
    dbdt = vparal * (gradb * b_vec) + gradb*ExBdrift
    dExBdt = vparal * (gradExB * b_vec) + gradExB*ExBdrift

    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m*b_vec/(q*B) × (vparal*dbdt + dExBdt)

    # Compute the acceleration along the magnetic field
    dvparaldt = (q*Eparal - mu*b_vec⋅gradB_vec)/m
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the total velocity of the guiding centre
    dRdt = vparal*b_vec + dRperpdt

    # Compute the pitch angle rate of change
    dbetadt = (1/vparal*dvparaldt - 1/2B*dBdt)*beta*(1-beta^2)

    # Compute the electron collision frequency
    eta = 2.91e-6 * n * coulomb_logarithm * T^(-1.5)
    beta_friction = -eta*beta


    # Update the statevector
    du .= [dRdt; dvparaldt; dbetadt + beta_friction]

end

struct GCA_pitchAngleDiffusion
end
function (eom::GCA_pitchAngleDiffusion)(du, u, params, _)
    R = u[1:3]    # Position of the gyrocentre
    beta = u[5]   # Cosine of pitch angle

    # Extract parameters
    # coulomb logarithm
    coulomb_logarithm = params[3]
    # number density and temperature interpolation objects
    n_itp, T_itp = params[9:10]
    # Interpolate density to guiding centre location
    n = n_itp(R...)
    T = T_itp(R...)

    # Compute the electron collision frequency
    eta = 2.91e-6 * n * coulomb_logarithm * T^(-1.5)

    # Compute the diffusion coefficients
    dRdW = [0,0,0]
    dvparaldW = 0
    dbetadW = √(eta*(1-beta^2))

    du .= [dRdW; dvparaldW; dbetadW]
end
