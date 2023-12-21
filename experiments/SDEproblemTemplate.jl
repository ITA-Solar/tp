using tp

npart = 100
tspan = (0,5)
dt = 1e-3

#X0 = VectorIC(collect(LinRange(-4,-2,npart)))
#V0 = FloatIC(0.0)
#ic = PhaseSpace1DIC(X0, V0)
ic = FloatIC(1.0)

p_drift = SingleCoefficientParams(1.0)
p_diffusion = SingleCoefficientParams(1.0)

particles = SDEParticle(
    ProportionalAdvection(),
    ProportionalDiffusion(),
    ic,
    p_drift,
    p_diffusion,
    npart
    )

# SOLVE
u, t = solve(particles, dt, tspan, euler_maruyama)
