```@meta
CurrentModule = TraceParticles
```

# Gallery

## [Magnetic mirroring](@id mirror-verifaction)
When a charge particle experiences increased magnetic field strength, the gyration radius decreases.
Field strength increases in the plane of the gyration give rise to a drift along $\nabla B/\times\mathbf{B}$.
If the field strength increases perpendicular to the gyration, parallel to the magnetic field direction, the paricle might be mirrored if the magnetic moment of the particle is conserved.
As the the Larmor radius decreases, conservation of the moment requires a larger perpendicular velocity, and if the total energy is conserved (for example if there are no electric field), the parallel velocity must decrease.
If the magnetic field becomes strong enough, all the energy will be transfered from the parallel component into the perpendicular, and the particle will be mirrored.
In other words, in a field with a magnet gradient that i along its direction, a charged particle will experience a force in the opposite direction, $F_\text{mirror} \propto \nabla B \cdot \hat{\mathbf{b}}$.
```@example mirror-verification
using OrdinaryDiffEq
using LinearAlgebra
using TraceParticles
using CairoMakie

#...............................................................................
# Electromagnetic field
function magneticbottle(x, y, z; B0, L)
    a = B0 * z / L^2
    return [-x * a, -y * a, B0 + z * a]
end

B0 = 10.0
L = 0.1
emfield(x, y, z, t) = zeros(3), magneticbottle(x, y, z; B0=B0, L=L)

#...............................................................................
# Particle parameters
tf = 1.5
tspan = (0, tf)
mass = 1
charge = 10

#...............................................................................
# Initial conditions
vel0 = [0.0, 0.4, 0.4]
rL = mass * √(vel0[1]^2 + vel0[2]^2) / (charge * B0)
pos0 = [-1rL, 0.0, 0.0]
E, B = emfield(pos0..., 0)
R0, vparal, μ = get_guidingcentre(pos0, vel0, B, E, charge, mass)

#...............................................................................
# Create problem
prob_FO = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, electromagneticfield=emfield),
)
prob_GCA = ODEProblem(
    guidingcentreapproximation!,
    [R0; vparal],
    tspan,
    (
        charge=charge,
        mass=mass,
        magneticmoment=μ,
        electromagneticfield=emfield
    )
)

#...............................................................................
# Run simulation
sol_FO = solve(prob_FO)
sol_GCA = solve(prob_GCA)

#-------------------------------------------------------------------------------
# Analytical solution

function z_analytical(t; μ, m, B0, L, A, ϕ)
    ω = √(2μ * B0 / (m* L^2))
    return @. A * sin(ω * t + ϕ)
end
Bmax = B0 * norm(vel0)^2 / (vel0[1]^2 + vel0[2]^2)
zmax = L * √(Bmax / B0 - 1)
ϕ = 0
times = range(0.0, tf, length=1000)
analytical = z_analytical(times; μ=μ, m=mass, B0=B0, L=L, A=zmax, ϕ=ϕ)

#-------------------------------------------------------------------------------
# Plot
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, [sol_FO(t)[3] for t in times], [sol_FO(t)[1] for t in times],
    label="Full orbit", linewidth=0.5
)
lines!(ax, sol_GCA[3,:], sol_GCA[1,:];
    label="GCA", linewidth=2.0
)
lines!(ax, analytical, [0.0 for _ in analytical];
    linestyle=:dash, label="Analytical solution", linewidth=1.0
)
ax.aspect = DataAspect()
ax.limits = ((nothing, nothing), (-0.01,0.03))
ax.yticks = ([-0.01,0.00,0.01])
axislegend(ax, position=:ct, orientation=:horizontal, framevisible=false)
```

## [Van Allen belts (magnetic dipole)](@id dipole-verification)
