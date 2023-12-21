#-------------------------------------------------------------------------------
# Created 28.09.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 solve.jl
#
#-------------------------------------------------------------------------------
# Module containing the the solve-function
#-------------------------------------------------------------------------------
function solve(
    particles::SDEParticle,
    dt       ::Real,
    tspan    ::Tuple{Number, Number},
    solver   ::Function
    ;
    dW=nothing,
    callback=nothing,
    seed=0,
    )
    t = collect(tspan[1]:dt:tspan[2])
    nsteps = convert(Int64, (tspan[2] - tspan[1])/dt)
    nvar = length(particles.ic())
    u = Array{Float64}(undef, nsteps+1, nvar, particles.npart)
    if dW === nothing
        dW = randn(MersenneTwister(seed), nsteps+1, particles.npart)*√dt
    end
    if callback === nothing
    # Add solve without callbacks if not needed
        callback = DiscreteTPCallback()
    #     solve_wo_cb(particles, dt, tspan, solver; dW=dW)
    end
    for n = 1:particles.npart
        u[1,:,n] .= particles.ic(n)
        for i = 1:nsteps
            u[i+1,:,n] = u[i,:,n] .+ solver(
                u[i,:,n],
                t[i],
                dt,
                dW[i,n],
                particles.eom_drift,
                particles.eom_diffusion,
                particles.p_drift(n),
                particles.p_diffusion(n)
                )
            if callback.condition(u[i+1,:,n])
                callback.affect!(u, i, n)
            end
        end
    end
    return u, t
end


# Currently only available for scalar SDEs
function solve_stoppingtime(
    particles::SDEParticle,
    dt       ::Real,
    u_max    ::Real,
    solver   ::Function
    )
    t = zeros(particles.npart)
    nsteps = zeros(Int64, particles.npart)
    for n = 1:particles.npart
        u = particles.ic(n)
        while u < u_max
            dW = randn()*√dt
            u += solver(
                u,
                t[n],
                dt,
                dW,
                particles.eom_drift,
                particles.eom_diffusion,
                particles.p_drift(n),
                particles.p_diffusion(n)
                )
            t[n] += dt
            nsteps[n] += 1
        end
    end
    return t, nsteps
end
