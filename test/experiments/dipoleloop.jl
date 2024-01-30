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
    phis[1,i] = ϕ
    phis[2,i] = ϕ_FO
    phis[3,i] = ϕ_GCA
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

# Old (without DifferentialEquations)
# phisResult = [18.1437  9.07184  6.04789  4.53592  3.62873  3.02395  2.59195  2.26796  2.01596  1.81437
#               11.619   8.8902   6.00361  4.50694  3.61259  2.99498  2.57028  2.24829  1.9973   1.79954
#               11.8231  9.05323  6.03546  4.5266   3.62127  3.01772  2.58665  2.26328  2.01182  1.81067]
phisResult = [18.143673169089446 9.071836584544723 6.047891056363149 4.5359182922723615 3.6287346338178885 3.0239455281815744 2.5919533098699206 2.2679591461361808 2.0159636854543828 1.8143673169089443; 13.477450313435256 10.122700032698448 6.02968383872634 4.0884423647442665 3.082113161073851 2.416777925984258 2.0094530756929636 1.7722890869213554 1.701666083090213 1.59628860675075; 11.730322020262847 9.041947510562848 6.054674192829948 4.469158187455404 3.5419696390159485 2.9781025502692096 2.5541303382812446 2.253375420931734 1.982363732483801 1.7898986454849026]

@testset verbose = true "GCA: RK4" begin
    @test all(@. isapprox(phisNew[3,:], phisResult[3,:], rtol=5e-6))
end
@testset verbose = true "Full orbit: RK4" begin
    @test all(@. isapprox(phisNew[2,:], phisResult[2,:], rtol=3e-6))
end
@testset verbose = true "An approx.: RK4" begin
    @test all(@. isapprox(phisNew[1,:], phisResult[1,:], rtol=2e-6))
end


             
