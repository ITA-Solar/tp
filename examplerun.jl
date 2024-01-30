using tp
using DifferentialEquations
using Distributed

eom = lorentzforce!
ic = [0.0]
p = [1.0,1.0]
tspan = (0,1)
npart = 1

prob = ODEProblem(eom, ic, p, tspan)

# For distributed parallelisation
@everywhere function prob_func_distributed(prob, i, repeat)
    remake(prob, u0=ic[i], p=p)
end
function prob_func(prob, i, repeat)
    remake(prob, u0=ic[i], p=p)
end

ensamble_prob = EnsembleProblem(
    prob, 
    prob_func=prob_func
    ;
    safetycopy=false,
)
@time sim = DifferentialEquations.solve(
    ensamble_prob, 
    EnsembleSerial()
    #EnsembleThreads()
    #EnsembleSplitThreads()
    #EnsembleDistributed()
    ;
    trajectories=npart, 
    progress=true,
    #callback=cb,
    maxiters=1000  
)


