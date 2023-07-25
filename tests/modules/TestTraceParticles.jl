#-------------------------------------------------------------------------------
# Created 01.07.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestTraceParticles.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the TraceParticles-module
#-------------------------------------------------------------------------------

module TestTraceParticles

using Test

using TraceParticles

export test_tp_save_and_tp_load

function test_tp_save_and_tp_load(verbose::Bool)
    @testset verbose=verbose "tp_save & tp_load" begin
	# Set simulation parameters
	include("../experiments/PlasmoidParameters.jl")
	time = 1e-1
	params = Base.invokelatest(PlasmoidParameters.plasmoidParameters,
				    (100,100,2), 1e-3
				   )
	params.tp_expdir = "tests/modules"

    	# Initialise experiment
    	exp = tp_init!(params);
    	
    	# Run experiment
    	tp_run!(exp)
    	
    	# Save simulation
    	tp_save(exp)
    	
    	# Load simulation parameters by using params-file
	include("modules/$(params.tp_expname)_params.jl")
    	
    	# Load particles, mesh and background
    	
    	exp2 = tp_load(params)
    	
    	# Check that variables of exp and exp2 is equal
    	@testset verbose = true "tp" begin
    	    @test exp.patch.tp.pos == exp2.patch.tp.pos
    	    @test !(exp.patch.tp.pos === exp2.patch.tp.pos)
    	    @test exp.patch.tp.vel == exp2.patch.tp.vel
    	    @test !(exp.patch.tp.vel === exp2.patch.tp.vel)
    	end
    	@testset verbose = true "mesh" begin
    	    @test exp.patch.mesh.xCoords == exp2.patch.mesh.xCoords
    	    @test !(exp.patch.mesh.xCoords === exp2.patch.mesh.xCoords)
    	    @test exp.patch.mesh.yCoords == exp2.patch.mesh.yCoords
    	    @test !(exp.patch.mesh.yCoords === exp2.patch.mesh.yCoords)
    	    @test exp.patch.mesh.zCoords == exp2.patch.mesh.zCoords
    	    @test !(exp.patch.mesh.zCoords === exp2.patch.mesh.zCoords)
    	end
    	@testset verbose = true "bg" begin
    	    @test exp.patch.mesh.bField == exp2.patch.mesh.bField
    	    @test !(exp.patch.mesh.bField === exp2.patch.mesh.bField)
    	    @test exp.patch.mesh.eField == exp2.patch.mesh.eField
    	    @test !(exp.patch.mesh.eField === exp2.patch.mesh.eField)
    	    @test exp.patch.mesh.∇B == exp2.patch.mesh.∇B
    	    @test !(exp.patch.mesh.∇B === exp2.patch.mesh.∇B)
    	    @test exp.patch.mesh.∇b̂ == exp2.patch.mesh.∇b̂
    	    @test !(exp.patch.mesh.∇b̂ === exp2.patch.mesh.∇b̂)
    	    @test exp.patch.mesh.∇ExB == exp2.patch.mesh.∇ExB
    	    @test !(exp.patch.mesh.∇ExB === exp2.patch.mesh.∇ExB)
    	end
    	end
    end # testset
end  # function 
