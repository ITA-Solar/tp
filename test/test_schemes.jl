#-------------------------------------------------------------------------------
# Created 02.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestSchemes.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Schemes-module.
#-------------------------------------------------------------------------------



"""
"""


function testeuler(verbose::Bool)
    # Test parameters
    pos = [1, 2, 2, 2.5]
    vel = [1, 2, 3, 1.5]
    acc = [1, 3, 4, 10.5]
    f(args...) = [vel;  acc]
    dt = 0.1
    @testset verbose=verbose "euler" begin
        posAnswer = [1.1, 2.2, 2.3, 2.65]
        velAnswer = [1.1, 2.3, 3.4, 2.55]
        svNext = euler([pos; vel], dt, f)
        @test posAnswer == svNext[1:4]
        @test velAnswer == svNext[5:8]
    end # testset euler
end # function testeuler


function testeulerCromer(verbose::Bool)
    # Test parameters
    pos = [1, 2, 2, 2.5]
    vel = [1, 2, 3, 1.5]
    acc = [1, 3, 4, 10.5]
    f(args...) = [vel;  acc]
    dt = 0.1
    @testset verbose=verbose "eulerCromer" begin
        posAnswer = [1.11, 2.23, 2.34, 2.755]
        velAnswer = [1.1,   2.3,  3.4,  2.55]
        svNext = eulerCromer([pos; vel], dt, f)
        @test posAnswer == svNext[1:4]
        @test velAnswer == svNext[5:8]
    end # testset euler
    end # function testeuler
