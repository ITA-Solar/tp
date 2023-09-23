#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Patches.jl
#
#-------------------------------------------------------------------------------
# Contains the Patch structs and methods
#-------------------------------------------------------------------------------



#-------------#   
# structs     # 
#-------------#-----------------------------------------------------------------

struct DEPatch <: AbstractPatch
    tspan::Tuple{Real, Real}
    tp   ::AbstractParticle
    mesh ::AbstractMesh
end 

function run!(patch::DEPatch)
    update(patch)
end 


function update(patch::DEPatch)

    function condition(u,t,integrator) 
        return u[1] <= 14.52e6 || u[1] >= 18.83e6 ||
                u[3] <= -7.79e6 || u[3] >= -3.88e6
    end
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)

    prob = get_problem(patch.tp, patch.tspan)

    ## For distributed parallelisation
    #println("Defining prob-function")
    #@everywhere function prob_func(prob, i, repeat)
    #    remake(prob, u0=patch.tp.ic(i), p=patch.tp.p(i))
    #end
    #function prob_func(prob, i, repeat)
    #    remake(prob, u0=patch.tp.ic(i), p=patch.tp.p(i))
    #end

    ensamble_prob = EnsembleProblem(
        prob, 
        prob_func=get_prob_func(patch.tp)
        ;
        safetycopy=false,
    )
    @time sim = solve(
        ensamble_prob, 
        EnsembleSerial()
        #EnsembleThreads()
        #EnsembleSplitThreads()
        #EnsembleDistributed()
        ;
        trajectories=patch.tp.npart, # or just patch.solver?
        progress=true,
        callback=cb,
        maxiters=1000  
    )
end

#------------------#
# Legacy-patch #
#------------------#------------------------------------------------------------

mutable struct Patch <: AbstractPatch
    mesh        ::Mesh
    tp          ::TraceParticle # The trace particles
    solver      ::Function
    scheme      ::Function
    interpolator::Function
    dt          ::Real
    numSteps    ::Integer
    numParticles::Integer
    periodicBC  ::Tuple{Bool, Bool, Bool}

    # Constructors
    #--------------------------------------------------------------------------
    function Patch(mesh        ::Mesh,
                   tp          ::TraceParticle, # The trace particles
                   solver      ::Function,
                   scheme      ::Function,
                   interpolator::Function,
                   dt          ::Real,
                   numSteps    ::Integer,
                   numParticles::Integer
                   )
        new(mesh, 
            tp, 
            solver, 
            scheme, 
            interpolator, 
            dt, 
            numSteps, 
            numParticles,
            (false, false, false)
            )
    end # constructor

    function Patch(mesh        ::Mesh,
                   tp          ::TraceParticle, # The trace particles
                   solver      ::Function,
                   scheme      ::Function,
                   interpolator::Function,
                   dt          ::Real,
                   numSteps    ::Integer,
                   numParticles::Integer,
                   periodicBC  ::Tuple{Bool, Bool, Bool}
                   )
        new(mesh, 
            tp, 
            solver, 
            scheme, 
            interpolator, 
            dt, 
            numSteps, 
            numParticles,
            periodicBC
            )
    end # constructor

    function Patch(
        expname     ::String,
        snap        ::Int,
        expdir      ::String,
        tp          ::TraceParticle, # The trace particles
        solver      ::Function,
        scheme      ::Function,
        interpolator::Function,
        dt          ::Real,
        numSteps    ::Integer,
        numParticles::Integer,
        periodicBC  ::Tuple{Bool, Bool, Bool}
        )
        mesh = Mesh(expname, snap, expdir)
        new(mesh, 
            tp, 
            solver, 
            scheme, 
            interpolator, 
            dt, 
            numSteps, 
            numParticles,
            periodicBC
            )
    end # constructor
end # mutable struct Patch

"""
    PatchDiffEq
AbstractPatch subtype for using the DifferentialEquations.jl julia-package.
"""
#struct PatchDiffEq <: AbstractPatch
#    
#end


#---------#
# Methods #
#---------#---------------------------------------------------------------------
function run!(patch::Patch)
    onepercentofsnap = ceil(Int64, patch.numSteps/100)
    for i = 1:patch.numSteps # Over timesteps
        if i%onepercentofsnap == 0
            print("Stepping progress: $(i/onepercentofsnap)% \r")
            flush(stdout)
        end
        push!(patch.tp,
                        patch.mesh,
                        i,
                        patch.dt,
                        patch.solver,
                        patch.interpolator,
                        patch.scheme,
                        patch.periodicBC,
                        )
    end # loop over timesteps (i)
end # function: run


#----------------#
# Base functions #
#-------------------------------------------------------------------------------
"""
    Base.copy(p::Patch)
Make a deep copy of a Patch-type.
"""
function Base.copy(p::Patch)
    Patch(p.mesh,
          copy(p.tp),
          p.solver,
          p.scheme,
          p.interpolator,
          p.dt,
          p.numSteps,
          p.numParticles,
          p.periodicBC
          )
end # function Base.copy


function Base.Multimedia.display(p::Patch)
    println("""Instance of mutable struct:
    Patches.Patch
        mesh        ::Meshes.Mesh
        tp          ::Particles.TraceParticle 
        solver      ::Function
        scheme      ::Function
        interpolator::Function
        dt          ::Real
        numSteps    ::Integer
        numParticles::Integer""")
end # function Base.Multimedia.display
#-------------------------------------------------------------------------------

