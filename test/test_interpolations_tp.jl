#-------------------------------------------------------------------------------
# Created 13.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestInterpolations_tp.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Interpolations_tp-module.
#-------------------------------------------------------------------------------


"""
    testtrilinear(verbose::Bool)

Testing trilinear interpolation function `trilinear` in a static electromagnetic
mesh and in a linear electromagnetic mesh.
In practice, this test also tests the Mesh-constructors: 
    Mesh(bField, eField)
    Mesh(bField, eField, xCoords, yCoords, zCoords)
"""
function testtrilinear(verbose::Bool)
    @testset verbose=verbose "trilinear" begin
        #
        # Test interpolation of uniform field
        #
        B0 = ones(3,3,3,3) # uniform magnetic field
        E0 = ones(3,3,3,3) # uniform electric field
        mesh0 = Mesh(B0, E0) # Create mesh
        cellIndices00 = (1, 1, 1)
        cellIndices01 = (2, 2, 2)
        point00 = (0.0, 0.0, 0.0)
        point01 = (2/3, 2/3, 2/3)
        answer0 = ([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])

        result00 = trilinear(mesh0, cellIndices00, point00)
        result01 = trilinear(mesh0, cellIndices01, point01)

        @test result00[1] ≈ answer0[1]
        @test result00[2] ≈ answer0[2]
        @test result01[1] ≈ answer0[1]
        @test result01[2] ≈ answer0[2]

        #
        # Test interpolation of linear field gradient
        #
        N = 3 # Number of vertices (grid points) on each axis of the mesh
        dx = 1/(N-1) # grid spacing

        # Field function
        linearField(x,y,z,a,b,c,d) = a*x + b*y + c*z .+ d
        a = b = c = d = 1 # Constants
        # Generate cartesian axes
        xx = collect(LinRange(0,1,N))
        yy = collect(LinRange(0,1,N))
        zz = collect(LinRange(0,1,N))
        # Evaluate fields on grid points
        B1 = zeros(3, N, N, N)
        E1 = zeros(3, N, N, N)
        # Broadcast into meshgrid?
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    for dim = 1:3
                        B1[dim, i,j,k] = linearField(xx[i],
                                                     yy[j],
                                                     zz[k],
                                                     a, b, c, d)
                        E1[dim, i,j,k] = linearField(xx[i],
                                                     yy[j],
                                                     zz[k],
                                                     2a, 2b, 2c, 2d)
                    end # loop over dim
                end # loop k
            end # loop j
        end # loop i

        # Initialise mesh
        mesh1 = Mesh(B1, E1, xx, yy, zz)
        # Define points to interpolate.
        #   Here the middle point of all cells.
        pointx = collect(0:N-1) .- dx/2
        pointy = collect(0:N-1) .- dx/2
        pointz = collect(0:N-1) .- dx/2
        # Evaluate the correct field values for these points
        answer = zeros(3, N-1, N-1, N-1) 
        # Loop through all cells of the mesh.
        #    Each index triple represent the corner of a celle
        for i = 1:N-1
            for j = 1:N-1
                for k = 1:N-1
                    answerB = linearField(pointx[i],
                                          pointy[j],
                                          pointz[k],
                                          a, b, c, d)
                    answerE = linearField(pointx[i],
                                          pointy[j],
                                          pointz[k],
                                          2a, 2b, 2c, 2d)
                    # Interpolate
                    result = trilinear(mesh1, (i,j,k), (pointx[i],
                                                        pointy[j],
                                                        pointz[k]))
                    # Check result
                    @test result[1] ≈ [answerB, answerB, answerB]
                    @test result[2] ≈ [answerE, answerE, answerE]
                end # over k
            end # over j
        end # over i
    end # testset
end # function testtrilinear


"""
    testlocateCell(verbose::bool)
Tests the grid search function `locateCell` which uses a binary search
algorithm. 
"""
function testlocateCell(verbose)
    @testset verbose=verbose "locateCell" begin
        N = 5
        coords = collect(LinRange(0, 1, N))
        dx = 1/(N-1)
        point = 2/3
        answer1 = ceil(Int64, point/dx)
        answer2 = 1
	answer3 = 2
	answer4 = 3
	answer5 = 4
        @test locateCell(coords, point) == answer1
        @test locateCell(coords, dx   ) == answer2
        @test locateCell(coords, 0.50 ) == answer3
        @test locateCell(coords, 0.75 ) == answer4
        @test locateCell(coords, 1.00 ) == answer5
    end # testset
end # function testLocateCell
