using TraceParticles
using Test

# Set simulation parameters
params = Parameters(
    dt=1e-8,
    nsteps=Int64(1e3),
    npart=5,
    tp_expname="nullpoint2k",
    tp_expdir="/Users/eilifo/repos/tp/tests/experiments",
    solver="full-orbit",
    scheme="RK4",
    interp="bilinear_xz",
    br_expname="nullpoint2k",
    br_expdir="/Users/eilifo/data/ohf/nullpoint2k",
    br_isnap=700,
    wp_part=Float64,
    wp_snap=Float32,
    pos_distr="point",
    posxbounds=[15.547e6, 15.739e6, 15.846e6, 15.947e6, 16.226e6],
    posybounds=[1e6, 1e6, 1e6, 1e6, 1e6],
    poszbounds=[-6.709e6, -6.519e6, -6.410e6, -6.314e6, -6.048e6],
    vel_distr="mb",
    pbc = (false, true, false),
    specie=[1,1,1,1,1],
)

# Initialise experiment
exp = tp_init!(params);

# Run experiment
tp_run!(exp)

# Save simulation
tp_save(exp)

# Load simulation parameters by using params-file
include("nullpoint2k_params.jl")

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
