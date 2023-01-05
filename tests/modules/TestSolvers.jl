#-------------------------------------------------------------------------------
# Created 02.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestSolvers.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Particles-module.
#-------------------------------------------------------------------------------

module TestSolvers

using Test

using WorkingPrecision: wpInt, wpFloat
using Solvers
using Schemes
using Interpolations
using Meshes

"""
"""

export testfullOrbit
export testrelfullOrbitExplLeapFrog

# Test parameters
pos = [1.0, 2.0, 2.5] # Arbitrary position
vel = [1.0, 2.0, 1.5] # Arbitrary velocity
relVel = 6e7 .* vel   # m/s -> γ ≈ 1.19
specie = 3 # Should yield mass = 1 and charge = 3
dt = 0.1
interpolator = Interpolations.trilinear
scheme = Schemes.euler

#------------------------------#
# Answer for the various tests #
#-------------------------------------------------------------------------------
# No fields
posAnswer1 = [1.1, 2.2, 2.65]
# Fields, but no velocity
B = [1.0, 2.0, 3.0]
E = [1.0, 2.0, 3.0]
posAnswer2, velAnswer2 = scheme(pos, [0.0, 0.0, 0.0], 3E, dt)
# Fields, non-relativistic velocity
accAnswer3 = 3E .+ 3*[3.0, -1.5, 0.0]
posAnswer3, velAnswer3 = scheme(pos, vel, accAnswer3, dt)

function testfullOrbit(verbose::Bool)
    @testset verbose=verbose "fullOrbit" begin
        #
        # No fields
        #
        B = zeros(3, 3, 3, 3)
        E = zeros(3, 3, 3, 3)
        mesh = Mesh(B, E)
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             vel,
                                             specie,
                                             mesh,
                                             dt,
                                             interpolator,
                                             scheme)
        @test nextPos == posAnswer1 
        @test nextVel == vel
        
        #
        # Fields, but no velocity
        #
        B[1, :, :, :] .= 1.0
        B[2, :, :, :] .= 2.0
        B[3, :, :, :] .= 3.0
        E[1, :, :, :] .= 1.0
        E[2, :, :, :] .= 2.0
        E[3, :, :, :] .= 3.0
        mesh = Mesh(B, E)
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             [0.0, 0.0, 0.0],
                                             specie,
                                             mesh,
                                             dt,
                                             interpolator,
                                             scheme)
        @test nextPos == posAnswer2 
        @test nextVel == velAnswer2

        #
        # B-field, but no E-field
        #
        B = zeros(3, 3, 3, 3)
        B[3, :, :, :] .= 1.0
        E = zeros(3, 3, 3, 3)
        mesh = Mesh(B, E)
        nextPos, nextVel = Solvers.fullOrbit([0.0, 0.0, 0.0],
                                             [1.0, 0.0, 1.0],
                                             specie,
                                             mesh,
                                             dt,
                                             interpolator,
                                             scheme)
        @test nextPos ≈ [0.1,   0.0, 0.1]
        @test nextVel ≈ [1.0,  -0.3, 1.0]

        #
        # Fields, non-relativistic velocity
        #
        B[1, :, :, :] .= 1.0
        B[2, :, :, :] .= 2.0
        B[3, :, :, :] .= 3.0
        E[1, :, :, :] .= 1.0
        E[2, :, :, :] .= 2.0
        E[3, :, :, :] .= 3.0
        mesh = Mesh(B, E)
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             vel,
                                             specie,
                                             mesh,
                                             dt,
                                             interpolator,
                                             scheme)
        @test nextPos == posAnswer3 
        @test nextVel == velAnswer3
        
        #
        # Fields, relativistic velocity (not really necessary)
        #
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             relVel,
                                             specie,
                                             mesh,
                                             dt,
                                             interpolator,
                                             scheme)
        accAnswer = 3E[:, 1, 1, 1] .+ 3*[1.8e8, -9.0e7, 0.0]
        posAnswer, velAnswer = scheme(pos, relVel, accAnswer, dt)
        @test nextPos == posAnswer 
        @test nextVel == velAnswer

    end # testset fullOrbit
end # function testfullOrbit


function testrelfullOrbitExplLeapFrog(verbose::Bool)
    @testset verbose=verbose "relFullOrbitExplLeapFrog" begin
        #
        # No fields
        #
        B = zeros(3, 3, 3, 3)
        E = zeros(3, 3, 3, 3)
        mesh = Mesh(B, E)
        posAnswer = [1.1, 2.2, 2.65]

        nextPosVay1, nextVelVay1 = Solvers.relFullOrbitExplLeapFrog(pos,
                                                            vel,
                                                            specie,
                                                            mesh,
                                                            dt,
                                                            interpolator,
                                                            Schemes.vay)
        nextPosBoris1, nextVelBoris1 = Solvers.relFullOrbitExplLeapFrog(pos,
                                                            vel,
                                                            specie,
                                                            mesh,
                                                            dt,
                                                            interpolator,
                                                            Schemes.boris)

        #
        # Fields, but no velocity
        #
        B[1, :, :, :] .= 1.0
        B[2, :, :, :] .= 2.0
        B[3, :, :, :] .= 3.0
        E[1, :, :, :] .= 1.0
        E[2, :, :, :] .= 2.0
        E[3, :, :, :] .= 3.0
        mesh = Mesh(B, E)

        nextPosVay2, nextVelVay2 = Solvers.relFullOrbitExplLeapFrog(pos,
                                                            [0.0, 0.0, 0.0],
                                                            specie,
                                                            mesh,
                                                            dt,
                                                            interpolator,
                                                            Schemes.vay)
        nextPosBoris2, nextVelBoris2 = Solvers.relFullOrbitExplLeapFrog(pos,
                                                            [0.0, 0.0, 0.0],
                                                            specie,
                                                            mesh,
                                                            dt,
                                                            interpolator,
                                                            Schemes.boris)
        @testset verbose=verbose "Vay" begin
            @test nextPosVay1 ≈ posAnswer1
            @test nextVelVay1 ≈ vel
        # Vay (and the other explicit leap-frog methods that integrates the
        # relativistic Lorentz force) is semi-implicit in the way that the 
        # position is advanced by an velocity equal to the average of the
        # current velocity and the velocity in the next time step.
            @test nextPosVay2 ≈ posAnswer2+0.5dt*([0.0,0.0,0.0] .+ nextVelVay2)
            @test nextVelVay2 ≈ velAnswer2
        end # testset Vay

        # Test if the Boris scheme gave the same result as Vay.
        @testset verbose=verbose "Boris" begin
            @test nextPosVay1 ≈ nextPosBoris1
            @test nextVelVay1 ≈ nextVelBoris1
            @test nextPosVay2 ≈ nextPosBoris2
            @test nextVelVay2 ≈ nextVelBoris2
        end # testset Boris
        
        # The following tests currently fail. I'm not yet sure whether they
        # should pass actually or not.
        """
        #
        # B-field, but no E-field
        #
        nextPos, nextVel = Solvers.vay([0.0, 0.0, 0.0],
                                       [1.0, 0.0, 1.0],
                                       specie,
                                       [0.0, 0.0, 1.0],
                                       [0.0, 0.0, 0.0],
                                       dt,
                                       scheme)
        @test nextPos ≈ [0.1, -0.03, 0.1]
        @test nextVel ≈ [1.0,  -0.30, 1.0]
        
        #
        # Fields, non-relativistic velocity
        #
        nextPos, nextVel = Solvers.vay(pos,
                                       vel,
                                       specie,
                                       B,
                                       E,
                                       dt,
                                       scheme)
        @test nextPos ≈ posAnswer3 
        @test nextVel ≈ velAnswer3
        """
    end # testset vay
end # function testvay

end # module TestSchemes
        
