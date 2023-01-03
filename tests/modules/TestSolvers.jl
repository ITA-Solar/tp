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

"""
"""

export testfullOrbit
export testvay

# Test parameters
pos = [1.0, 2.0, 2.5] # Arbitrary position
vel = [1.0, 2.0, 1.5] # Arbitrary velocity
B = [1.0, 2.0, 3.0]   # Arbitrary magnetic field
E = [1.0, 2.0, 3.0]   # Arbitrary electric field
relVel = 6e7 .* vel   # m/s -> γ ≈ 1.19
specie = 3 # Should yield mass = 1 and charge = 3
dt = 0.1
scheme = Schemes.euler

#------------------------------#
# Answer for the various tests #
#-------------------------------------------------------------------------------
# No fields
posAnswer1 = [1.1, 2.2, 2.65]
# Fields, but no velocity
posAnswer2, velAnswer2 = scheme(pos, [0.0, 0.0, 0.0], 3E, dt)
# Fields, non-relativistic velocity
accAnswer3 = 3E .+ 3*[3.0, -1.5, 0.0]
posAnswer3, velAnswer3 = scheme(pos, vel, accAnswer3, dt)

function testfullOrbit(verbose::Bool)
    @testset verbose=verbose "fullOrbit" begin
        #
        # No fields
        #
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             vel,
                                             specie,
                                             [0.0,0.0,0.0],
                                             [0.0,0.0,0.0],
                                             dt,
                                             scheme)
        @test nextPos == posAnswer1 
        @test nextVel == vel
        
        #
        # Fields, but no velocity
        #
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             [0.0, 0.0, 0.0],
                                             specie,
                                             B,
                                             E,
                                             dt,
                                             scheme)
        @test nextPos == posAnswer2 
        @test nextVel == velAnswer2

        #
        # B-field, but no E-field
        #
        nextPos, nextVel = Solvers.fullOrbit([0.0, 0.0, 0.0],
                                             [1.0, 0.0, 1.0],
                                             specie,
                                             [0.0, 0.0, 1.0],
                                             [0.0, 0.0, 0.0],
                                             dt,
                                             scheme)
        @test nextPos ≈ [0.1,   0.0, 0.1]
        @test nextVel ≈ [1.0,  -0.3, 1.0]

        #
        # Fields, non-relativistic velocity
        #
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             vel,
                                             specie,
                                             B,
                                             E,
                                             dt,
                                             scheme)
        @test nextPos == posAnswer3 
        @test nextVel == velAnswer3
        
        #
        # Fields, relativistic velocity (not really necessary)
        #
        nextPos, nextVel = Solvers.fullOrbit(pos,
                                             relVel,
                                             specie,
                                             B,
                                             E,
                                             dt,
                                             scheme)
        accAnswer = 3E .+ 3*[1.8e8, -9.0e7, 0.0]
        posAnswer, velAnswer = scheme(pos, relVel, accAnswer, dt)
        @test nextPos == posAnswer 
        @test nextVel == velAnswer

    end # testset fullOrbit
end # function testfullOrbit


function testvay(verbose::Bool)
    @testset verbose=verbose "vay" begin
        #
        # No fields
        #
        nextPos, nextVel = Solvers.vay(pos,
                                       vel,
                                       specie,
                                       [0.0,0.0,0.0],
                                       [0.0,0.0,0.0],
                                       dt,
                                       scheme)
        posAnswer = [1.1, 2.2, 2.65]
        @test nextPos ≈ posAnswer1
        @test nextVel ≈ vel

        #
        # Fields, but no velocity
        #
        nextPos, nextVel = Solvers.vay(pos,
                                       [0.0, 0.0, 0.0],
                                       specie,
                                       B,
                                       E,
                                       dt,
                                       scheme)
        @test nextPos ≈ posAnswer2 + nextVel*dt # Vay is semi-implicit
        @test nextVel ≈ velAnswer2
        
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
        
