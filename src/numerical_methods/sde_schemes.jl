#-------------------------------------------------------------------------------
# Created 02.12.22. Code from back in 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 sde_schemes
#
#-------------------------------------------------------------------------------
# Contains various schemes for solving ordinary stochastic differential 
# equations
#-------------------------------------------------------------------------------

function euler_maruyama(
    X,
    t,
    dt,
    dW,
    drift,
    diffusion,
    drift_params,
    diffusion_params
    )
    # The drift coefficient
    a = drift(X, drift_params, t)
    b = diffusion(X, diffusion_params, t)
    return a*dt + b*dW
end

function milstein_central(
    X,
    t,
    dt,
    dW,
    drift,
    diffusion,
    drift_params,
    diffusion_params
    )
    # Extract the parameters:
    #   The constant drift coefficient
    #   the constan diffusion coefficient
    a = drift(X, drift_params, t)
    b = diffusion(X, diffusion_params, t)
    delta = 0.5e-6
    dbdX = (diffusion(X .+ delta, diffusion_params, t)
           .- diffusion(X .- delta, diffusion_params, t)
           )./2delta
    return a*dt + + b*dW + 0.5*b .* dbdX * (dW^2 - dt)
end
