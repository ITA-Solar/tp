
abstract type AbstractTPCallback end

struct DiscreteTPCallback{F1, F2} <: AbstractTPCallback
    condition::F1
    affect!  ::F2
    function DiscreteTPCallback(condition::F1, affect!::F2) where {F1, F2}
        new{F1, F2}(condition, affect!)
    end
end
function DiscreteTPCallback()
    condition(x) = false
    affect! = nothing
    return DiscreteTPCallback(condition, affect!)
end


function killparticle!(integrator)
    terminate!(integrator)
end

"""
    outside2dnullpointzoom(u, t, integrator)
Callback condition for particles exiting the zoomed-in domain in the 2D
nullpoint Bifrost simulation used for the results in the Japan Hinode/IRIS
poster of septermber 2023.

The effect of this callback should be termination, which means that the
callback should be created like this:
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
"""
function outside2dnullpointzoom(u,t,integrator)
    return u[1] <= 14.52e6 || u[1] >= 18.83e6 ||
            u[3] <= -7.79e6 || u[3] >= -3.88e6
end


