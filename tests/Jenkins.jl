#-------------------------------------------------------------------------------
# Created 10.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running tests checked by Jenkins after it registrers activity on 
# repository.
#-------------------------------------------------------------------------------

using Test
using TestInterpolations
using TestSchemes
using TestSolvers

verbose = 3

@testset verbose = verbose ≥ 1 "All tests" begin
    @testset verbose = verbose ≥ 2 "Unit tests" begin
        @testset verbose = verbose ≥ 3 "Interpolations" begin
            testtrilinear(verbose ≥ 4)
            testlocateCell(verbose ≥ 4)
        end
        @testset verbose = verbose ≥ 3 "Schemes" begin
            testeuler(verbose ≥ 4)
            testeulerCromer(verbose ≥ 4)
        end
        @testset verbose = verbose ≥ 3 "Solvers" begin
            testfullOrbit(verbose ≥ 4)
            testrelfullOrbitExplLeapFrog(verbose ≥ 4)
        end
    end # testset Unit tests
    
    @testset verbose = verbose ≥ 3 "Experiments" begin
        include("experiments/testExBdrift.jl")
    end # testset Experiments

end # testset All test
