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
    eom_drift    ::Function,
    eom_diffusion::Function,
    ic           ::Vector{<:Vector},
    dt           ::Real,
    tspan        ::Tuple{Number, Number},
    p_drift      ::Tuple{Vararg{Any}},
    p_diffusion  ::Tuple{Vararg{Any}},
    scheme       ::Function,
    npart        ::Integer
    ;
    dW=nothing,
    callback=nothing,
    seed=0,
    )
    t = collect(tspan[1]:dt:tspan[2])
    nsteps = convert(Int64, (tspan[2] - tspan[1])/dt)
    nvar = length(ic[1])
    u = Array{Float64}(undef, nsteps+1, nvar, npart)
    if dW === nothing
        dW = randn(MersenneTwister(seed), nsteps+1, npart)*√dt
    end
    if callback === nothing
    # Add solve without callbacks if not needed
        callback = DiscreteTPCallback()
    #     solve_wo_cb(particles, dt, tspan, solver; dW=dW)
    end

    verbose = 1
    if verbose == 1
        statement1 = string("tp.jl: Running simulation:\n",
            "\tnpart:  $(npart)\n",
            "\tnsteps: $(@sprintf("%.2e", nsteps))\n",
            "\tdt:     $(dt)\n",
            "\tNumber of iterations: ",
            "$(@sprintf("%.2e",npart*nsteps))"
            )
        println(statement1)
    end

    # Loop over particles
    for n = 1:npart
        u[1,:,n] .= ic[n]
        # Loop over timesteps
        for i = 1:nsteps
            u[i+1,:,n] = u[i,:,n] .+ scheme(
                u[i,:,n],
                t[i],
                dt,
                dW[i,n],
                eom_drift,
                eom_diffusion,
                p_drift,
                p_diffusion
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
    eom_drift    ::Function,
    eom_diffusion::Function,
    ic           ::Vector{<:Vector},
    dt           ::Real,
    u_max        ::Real,
    p_drift      ::Tuple{Vararg{Any}},
    p_diffusion  ::Tuple{Vararg{Any}},
    scheme       ::Function,
    npart        ::Integer
    )
    t = zeros(npart)
    nsteps = zeros(Int64, npart)
    for n = 1:npart
        u = ic[n]
        while u < u_max
            dW = randn()*√dt
            u += scheme(
                u,
                t[n],
                dt,
                dW,
                eom_drift,
                eom_diffusion,
                p_drift,
                p_diffusion
                )
            t[n] += dt
            nsteps[n] += 1
        end
    end
    return t, nsteps
end
