# Created 21.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dipoleloop.jl
#
#-------------------------------------------------------------------------------
# Runs a series of magnetic dipole runs using 'dipole.jl' to reproduce fig. 19
# in Ripperda et al., 2018. The tests uses a previous result to compare
# with. The previous result was empirically validated by making a plot of
# 'phisNew' and comparing it with fig 19. 
#-------------------------------------------------------------------------------

qMmvec = collect(1:10)*10.0
phis = zeros(3, length(qMmvec))

for i = 1:length(qMmvec)
    global qMm = qMmvec[i]
    include("dipole.jl")
    phis[1,i] = ϕ[end]
    phis[2,i] = ϕ_FO[end]
    phis[3,i] = ϕ_GCA[end]
end

# The calculation of the angle phi is by evaluating atan(y/x).
# Rotations more than π/2 degrees thus needs to be adjusted.
# The adjustments above are the ones that make the result most
# equal to Ripperda et al., 2018. All of them are pretty certain
# except maybe the first angle (qMm = 10) where the analytical
# approximation fail.
phisNew = copy(phis)
phisNew[2:3,4:end] .= phis[2:3,4:end] .+ π
phisNew[2:3,3] .= phis[2:3,3] .+ 2π
phisNew[2:3,2] .= phis[2:3,2] .+ 3π
phisNew[2:3,1] .= phis[2:3,1] .+ 4π

 phisResult = [18.1437  9.07184  6.04789  4.53592  3.62873  3.02395  2.59195  2.26796  2.01596  1.81437
               11.619   8.8902   6.00361  4.50694  3.61259  2.99498  2.57028  2.24829  1.9973   1.79954
               11.8231  9.05323  6.03546  4.5266   3.62127  3.01772  2.58665  2.26328  2.01182  1.81067]

@testset verbose = true "GCA: RK4" begin
    @test all(@. isapprox(phisNew[3,:], phisResult[3,:], rtol=5e-6))
end
@testset verbose = true "Full orbit: RK4" begin
    @test all(@. isapprox(phisNew[2,:], phisResult[2,:], rtol=3e-6))
end
@testset verbose = true "An approx.: RK4" begin
    @test all(@. isapprox(phisNew[1,:], phisResult[1,:], rtol=2e-6))
end


             
