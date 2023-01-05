#-------------------------------------------------------------------------------
# Created 13.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running tests.
#-------------------------------------------------------------------------------

using Test
using TestInterpolations
using TestSchemes
using TestSolvers

verbose = 3
    
@testset verbose = verbose ≥ 1 "Modules" begin
    @testset verbose = verbose ≥ 2 "Interpolations" begin
        testtrilinear(verbose ≥ 1)
        testlocateCell(verbose ≥ 1)
    end
    @testset verbose = verbose ≥ 2 "Schemes" begin
        testeuler(verbose ≥ 3)
        testeulerCromer(verbose ≥ 3)
    end
    @testset verbose = verbose ≥ 2 "Solvers" begin
        testfullOrbit(verbose ≥ 3)
        testrelfullOrbitExplLeapFrog(verbose ≥ 3)
    end
end # testset Modules
