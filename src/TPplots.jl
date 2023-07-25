#-------------------------------------------------------------------------------
# Created 19.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                TPplots.jl
#
#-------------------------------------------------------------------------------
# Module containing plotting functions
#-------------------------------------------------------------------------------

module TPplots

# External libraries
import Plots
import PyPlot
const plt = PyPlot
using LinearAlgebra:    norm
using LaTeXStrings
# Internal libraries
using Meshes
using Patches:          Patch
using Particles
using Utilities:        rejectionsampling, norm4, norm2
using Interpolations_tp
using Schemes
using Solvers

export plt
export plot
export plotKE
export plot3D
export plottraj
export streamplotslice!
export pcolormeshslice!
export pcolormesh_xz!
export pcolormesh_xy!
export pcolormesh_yz!
export trajectoryslice!
export trajectory3D!
export quiverslice!

#------#
# Mesh #
#-------------------------------------------------------------------------------
function quiverslice!(
    ax  ::plt.PyCall.PyObject,
    mesh::Mesh, 
    normal::String, 
    point::Integer,
    field::String="B"
    ;
    step=1,
    yflip=false,
    )
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        ax.quiver(mesh.yCoords[begin:step:end],
               mesh.zCoords[begin:step:end],
               transpose(f[2, point, begin:step:end, begin:step:end]),
               transpose(f[3, point, begin:step:end, begin:step:end])
               )
    elseif normal == "y"
        ax.quiver(mesh.xCoords[begin:step:end],
               mesh.zCoords[begin:step:end],
               transpose(f[1, begin:step:end, point, begin:step:end]),
               transpose(f[3, begin:step:end, point, begin:step:end])  
               )
    elseif normal == "z"
        ax.quiver(mesh.xCoords[begin:step:end],
               mesh.yCoords[begin:step:end],
               transpose(f[1, begin:step:end, begin:step:end, point]),
               transpose(f[2, begin:step:end, begin:step:end, point])  
               )
    else
        println("Error: Plane not valid.")
    end
    if yflip
        ax.invert_yaxis()
    end
end # function quiverslice


function streamplotslice!(
    ax    ::plt.PyCall.PyObject,
    mesh  ::Mesh, 
    normal::String, 
    point ::Integer,
    field ::String="B")
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        xx = mesh.yCoords
        yy = mesh.zCoords
        uu = transpose(f[2, point, :, :])
        vv = transpose(f[3, point, :, :])
    elseif normal == "y"
        xx = mesh.xCoords
        yy = mesh.zCoords
        uu = transpose(f[1, :, point, :])
        vv = transpose(f[3, :, point, :])
    elseif normal == "z"
        xx = mesh.xCoords
        yy = mesh.yCoords
        uu = transpose(f[1, :, :, point])
        vv = transpose(f[2, :, :, point])
    else
        println("Error: Plane not valid.")
    end
    ax.streamplot(
        xx,
        yy,
        uu,
        vv,
        linewidth=0.3,
        arrowsize=0.6,
        color="black"
    )
end # function streamplotslice


function pcolormesh_xz!(
    ax     ::plt.PyCall.PyObject,
    A      ::Matrix{T} where {T<:Real},
    mesh   ::Mesh
    ;
    yflip  ::Bool=false,
    cmlabel::String="",
    title  ::String=""
    )
    pfabs = copy(transpose(A))
    uuu = mesh.xCoords' .* ones(length(mesh.zCoords))
    vvv = ones(length(mesh.xCoords))' .* mesh.zCoords
    pcm = ax.pcolormesh(uuu, vvv, pfabs,
                        alpha=1.0,
                        cmap=plt.get_cmap("plasma"))
    cb = plt.colorbar(pcm,
                      label=cmlabel,
                      ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(title)
    if yflip
        ax.invert_yaxis()
    end
end # function pcolormesh_xz!
#|
function pcolormesh_xz!(
    ax     ::plt.PyCall.PyObject,
    A      ::Matrix{T} where {T<:Real},
    x      ::Vector{T} where {T<:Real},
    y      ::Vector{T} where {T<:Real},
    z      ::Vector{T} where {T<:Real},
    ;
    yflip  ::Bool=false,
    cmlabel::String="",
    title  ::String=""
    )
    pfabs = copy(transpose(A))
    uuu = x' .* ones(length(z))
    vvv = ones(length(x))' .* z
    pcm = ax.pcolormesh(uuu, vvv, pfabs,
                        alpha=1.0,
                        cmap=plt.get_cmap("plasma"))
    cb = plt.colorbar(pcm,
                      label=cmlabel,
                      ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(title)
    if yflip
        ax.invert_yaxis()
    end
end # pcolormesh_xz!

function pcolormesh_xy!(
    ax     ::plt.PyCall.PyObject,
    A      ::Matrix{T} where {T<:Real},
    mesh   ::Mesh
    ;
    yflip  ::Bool=false,
    cmlabel::String="",
    title  ::String=""
    )
    pfabs = copy(transpose(A))
    uuu = mesh.xCoords' .* ones(length(mesh.yCoords))
    vvv = ones(length(mesh.xCoords))' .* mesh.yCoords
    pcm = ax.pcolormesh(uuu, vvv, pfabs,
                        alpha=1.0,
                        cmap=plt.get_cmap("plasma"))
    cb = plt.colorbar(pcm,
                      label=cmlabel,
                      ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    if yflip
        ax.invert_yaxis()
    end
end # function pcolormesh_xy!


function pcolormesh_yz!(
    ax     ::plt.PyCall.PyObject,
    A      ::Matrix{T} where {T<:Real},
    mesh   ::Mesh
    ;
    yflip  ::Bool=false,
    cmlabel::String="",
    title  ::String=""
    )
    pfabs = copy(transpose(A))
    uuu = mesh.yCoords' .* ones(length(mesh.zCoords))
    vvv = ones(length(mesh.yCoords))' .* mesh.zCoords
    pcm = ax.pcolormesh(uuu, vvv, pfabs,
                        alpha=1.0,
                        cmap=plt.get_cmap("plasma"))
    cb = plt.colorbar(pcm,
                      label=cmlabel,
                      ax=ax)
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_title(title)
    if yflip
        ax.invert_yaxis()
    end
end # function pcolormesh_yz!


function pcolormeshslice!(
    ax    ::plt.PyCall.PyObject,
    mesh  ::Mesh, 
    normal::String, 
    point ::Integer,
    field ::String="B")
    if field == "B"
        f = mesh.bField
        label = latexstring("\$|\\mathbf{B}|\$")
    elseif field == "E"
        f = mesh.eField
        label = latexstring("\$|\\mathbf{E}|\$")
    else
        println("Error: Invalid field.")
    end
    fabs = norm4(f)
    if normal == "x"
        pfabs = copy(transpose(fabs[point,:,:]))
        uuu = mesh.yCoords' .* ones(length(mesh.zCoords))
        vvv = ones(length(mesh.yCoords))' .* mesh.zCoords
    elseif normal == "y"
        pfabs = copy(transpose(fabs[:,point,:]))
        uuu = mesh.xCoords' .* ones(length(mesh.zCoords))
        vvv = ones(length(mesh.xCoords))' .* mesh.zCoords
    elseif normal == "z"
        pfabs = copy(transpose(fabs[:,:,point]))
        uuu = mesh.xCoords' .* ones(length(mesh.yCoords))
        vvv = ones(length(mesh.xCoords))' .* mesh.yCoords
    else
        println("Error: Plane not valid.")
    end
    pcm = ax.pcolormesh(uuu, vvv, pfabs,
                        alpha=1.0,
                        cmap=plt.get_cmap("Greys"))
    cb = plt.colorbar(pcm,
                      label=label,
                      ax=ax)
end # function colormeshslice


function trajectoryslice!(
    ax    ::plt.PyCall.PyObject,
    patch ::Patch,
    normal::String, 
    ;
    linestyle::String="",
    labelling::Bool=false,
    periodic::Bool=false
    )
    # Set variables depending on slice orientation
    if normal == "x"
        k = 1
    elseif normal == "y"
        k = 2
    elseif normal == "z"
        k = 3
    else
        println("Error: Plane not valid.")
    end
    pos = getpos(patch.tp)[1:end .!= k, :, :]
    domain = patch.mesh.domain[1:end .!= k, :]
    # Set colour-scheme
    cm = plt.get_cmap(:tab20)
    colorrange = (0:patch.numParticles) ./ patch.numParticles
    for j = 1:patch.numParticles
        if periodic
            plotperiodictrajectory(ax, pos[:,j,:],
                                   j, domain,
                                   cm, colorrange,
                                   linestyle)
        else
            ax.plot(pos[1,j,:], pos[2,j,:])
        end
    end

    # Mark initial positions of particles
    for j = 1:patch.numParticles
        if j == 1
            ax.plot(pos[1,1,1], pos[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
        else
            ax.plot(pos[1,j,1], pos[2,j,1], marker=".", color="Black")
        end
    end
end # function trajectoryslice!


function trajectory3D!(
    ax       ::plt.PyCall.PyObject,
    patch    ::Patch,
    linestyle::String="-",
    labelling::Bool=false,
    )
    # It seems that a 3D-projected PyPlot does not support more than ~16300
    # points per graph. Hence we have to skip some data
    skip = trunc(Integer, ceil(patch.numSteps/16_000))
    if skip == 0
        skip = 1
    end
    pos = getpos(patch.tp)
    for j = 1:patch.numParticles
        ax.plot(pos[1, j, 1:skip:end],
                pos[2,j,1:skip:end],
                pos[3,j,1:skip:end],
                linestyle=linestyle)
    end # loop over particles
    # Mark initial positions of particles
    for j = 1:patch.numParticles
        if j == 1
            ax.plot(pos[1,1,1], pos[2,1,1], pos[3,1,1], marker=".",
                    color="Black", label=latexstring("\$t_0\$"),
                    linestyle="None")
        else
            ax.plot(pos[1,j,1], pos[2,j,1], pos[3,j,1],
                    marker=".", color="Black")
        end
    end
end #function trajectory3D!

#-----------#
# Particles #
#-------------------------------------------------------------------------------
function plotenergydistr(particles::Particles.ParticleSoA, 
                         snap     ::Integer,
                         numbins  ::Integer,
                         title
                         )
    absvel = norm2(particles.vel[:, :, snap])
    binrange = range(minimum(absvel), maximum(absvel), length=numbins)
    Plots.histogram(absvel, bins=binrange)
    Plots.xlabel!("Absolute velocity, m/s")
    Plots.ylabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr
#
function plotenergydistr(particles::Particles.GCAParticleSoA, 
                         snap     ::Integer,
                         numbins  ::Integer,
                         title
                         )
    vparal = particles.vparal[:, snap]
    binrange = range(minimum(vparal), maximum(vparal), length=numbins)
    Plots.histogram(vparal, bins=binrange)
    Plots.xlabel!("Parallel velcity, m/s")
    Plots.ylabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr
function plotenergydistr(absvel ::Vector{T} where {T<:Real},
                         numbins::Integer,
                         title
                         )
    binrange = range(minimum(absvel), maximum(absvel), length=numbins)
    Plots.histogram(absvel, bins=binrange)
    Plots.xlabel!("Absolute velocity, m/s")
    Plots.ylabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr



function plotperiodictrajectory(
    ax::plt.PyCall.PyObject,
    pos::Matrix{T} where {T<:Real},
    partidx::Integer,
    domain ::Matrix{T} where {T<:Real},
    cm,
    colorrange,
    linestyle::String,
    )
    extent = domain[:,2] .- domain[:,1]
    posjumps = diff(pos, dims=2)
    posjumps = @. abs(posjumps)
    mask = @. isapprox(posjumps, extent, rtol=0.01)
    indices = findall(mask[1:2, :])
    numindices = length(indices)
    if numindices > 0
        j = 1
        for s = 1:length(indices)
            i = indices[s][2] # We only care about what time step the particle
            # hit the boundary, not which boundary was crossed.
            ax.plot(pos[1,j:i], pos[2,j:i], color=cm(colorrange[partidx]),
                    label="Particle $partidx", linestyle=linestyle)
            j = i+1
            break
        end
        ax.plot(pos[1,j:end], pos[2,j:end], color=cm(colorrange[partidx]))
    else
        ax.plot(pos[1,:], pos[2,:], color=cm(colorrange[partidx]),
                    label="Particle $partidx", linestyle=linestyle)
    end
end 

#---------#
# Patches #
#-------------------------------------------------------------------------------
function plot(
    patch    ::Patch,
    normal   ::String="z",
    ;
    point    ::Integer=1,
    labelling::Bool=false,
    periodic  ::Bool=false,
    streamplot::Bool=false,
    step      ::Integer=1,
    yflip     ::Bool=false,
    efield    ::Bool=false,
    )
    #
    if efield
        fig, axes = plt.subplots(1,2)
        ax = axes[1]
        eax = axes[2]
    else
        fig, ax = plt.subplots(1,1)
    end
   
    plot(ax, patch, normal, point, labelling, "B";
         periodic=periodic,
         streamplot=streamplot,
         step=step,
         yflip=yflip
         )
    if efield
        plot(eax, patch, normal, point, labelling, "E";
             periodic=periodic,
             streamplot=streamplot,
             step=step,
             yflip=yflip
             )
    end
    # Make energy distribution-plots
    #numbins = 20
    #p1 = plotenergydistr(patch.tp, 1, numbins, "Initial energy")
    #p2 = plotenergydistr(patch.tp, patch.numSteps+1, numbins, "final energy")
    #Plots.plot(p1,p2)#p3,p4)
    return fig, ax
end # function plot


function plot(
    ax       ::plt.PyCall.PyObject,
    patch    ::Patch,
    normal   ::String="z",
    point    ::Integer=1,
    labelling::Bool=false,
    field    ::String="B",
    ;
    periodic  ::Bool=false,
    streamplot::Bool=false,
    step      ::Integer=1,
    yflip     ::Bool=false,
    )   
    # Make pcolormesh of magnetic field strength
    pcolormeshslice!(ax, patch.mesh, normal, point, field)
    # Make streamplot or quiverplot of magnetic field
    if streamplot
        streamplotslice!(ax, patch.mesh, normal, point, field)
    else
        quiverslice!(ax, patch.mesh, normal, point, field,
                     step=step,
                     yflip=yflip
                     )
    end
    # Plot particle trajectories
    trajectoryslice!(ax, patch, normal, linestyle="-", periodic=periodic)

    # Set title and labelling
    if typeof(patch.tp) == ParticleSoA
        ax.set_title("Full-orbit")
    else
        ax.set_title("GCA")
    end
    if labelling
        ax.legend()
    end
    setcartesianaxes!(ax, normal)
    setaxislimits!(ax, patch.mesh.domain, normal)
end


function plot(
    tp::ParticleSoA,
    normal::String="z"
    ;
    step  ::Integer=1
    )
    _, npart, _ = size(tp.pos)
    fig, ax = plt.subplots(1, npart)
    for i = 1:npart
        if normal == "x"
            ax[i].plot(tp.pos[2,i,begin:step:end], tp.pos[3,i,begin:step:end])
        elseif normal == "y"
            ax[i].plot(tp.pos[1,i,begin:step:end], tp.pos[3,i,begin:step:end])
        elseif normal == "z" 
            ax[i].plot(tp.pos[1,i,begin:step:end], tp.pos[2,i,begin:step:end])
        end
        ax[i].set_title("Trajectory of particle $i")
        setcartesianaxes!(ax[i], normal)
    end
    return fig, ax
end
    

function plotKE(
    patch    ::Patch,
    labelling::Bool=false,
    )
    # Calculate kinetic energy at each timestep
    times = collect(range(0.0, step=patch.dt, length=patch.numSteps+1))
    Ek = kineticenergy(patch.tp, patch.mesh, trilinear)
    # Create figure and plot energies
    fig, ax = plt.subplots(1,1)
    plotKE(ax, Ek, times)
    ax.set_xlabel("Time, s")
    ax.set_ylabel("Energy, J")
    ax.set_title("Kinetic energy")
    if labelling
        ax.legend()
    end
end # function plotKE
#|
function plotKE(
    ax   ::plt.PyCall.PyObject,
    Ek   ::Matrix{T} where {T<:Real},
    times::Vector{T} where {T<:Real},
    )
    nj, ni = size(Ek)
    for j = 1:nj
        ax.plot(times, Ek[j,:])
    end 
end 
#|
function plotKE(
    ax   ::plt.PyCall.PyObject,
    Ek   ::Tuple{Matrix{T} where {T<:Real}, Matrix{T} where {T<:Real}, Matrix{T} where {T<:Real}},
    times::Vector{T} where {T<:Real},
    )
    Ekk, Ekparal, Ekperp = Ek
    nj, ni = size(Ekk)
    cm = plt.get_cmap(:tab20)
    colorrange = (0:nj) ./ nj
    for j = 1:nj
        ax.plot(times, Ekk[j,:], color=cm(colorrange[j]))
        ax.plot(times, Ekparal[j,:], linestyle="--", color=cm(colorrange[j]))
        ax.plot(times, Ekperp[j,:], linestyle=":", color=cm(colorrange[j]))
    end 
end 


function plottraj(
    patch::Patch,
    normal   ::String="z",
    linestyle::String="solid",
    labelling::Bool=false
    )
    fig, axes = plt.subplots(1,1)
    trajectoryslice!(axes, patch, normal, linestyle)
    # Set title and labelling
    if typeof(patch.tp) == ParticleSoA
        axes.set_title("Full-orbit")
    else
        axes.set_title("GCA")
    end
    if labelling
        axes.legend()
    end
    setcartesianaxes!(axes, normal)
end 


function plot3D(
    patch::Patch,
    linestyle::String="solid",
    labelling::Bool=false
    )
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel(latexstring("\$x\$"))
    ax.set_ylabel(latexstring("\$y\$"))
    ax.set_zlabel(latexstring("\$z\$"))
    ax.set_title("Particle trajectories")
    trajectory3D!(ax, patch, linestyle, labelling)
    return fig, ax
end # function plot3D


#-----------#
# Particles #
#-------------------------------------------------------------------------------
function plotKE(
    tp       ::ParticleSoA,
    times    ::Vector{T} where {T<:Real},
    mesh     ::Mesh,
    labelling::Bool=false,
    )
    Ek = kineticenergy(tp)
    
end # function plotKE

#---------------------#
# Field line tracing  #
#-------------------------------------------------------------------------------
function generatefieldline(
    mesh,
    initpos,
    stepsize,
    numsteps,
    interpolator,
    scheme
    )

    # Find out how many steps are needed in both directions
#    pos = zeros(mesh.numdims)
#    # Forward
#    pos .= initpos
#    stepsforward = -1
#    # While the position is till in the domain
#    while all(pos .>= mesh.domain[:, 1]) & all(pos .<= mesh.domain[:, 2])
#        # Get magnetic field at position by interpolation
#        fields, _ = grid(mesh, interpolator, pos)
#        bfield, _ = fields
#        # The movement should be dependent on the field direction, not the field
#        # strength 
#        B = norm(bfield) # field strength
#        b = bfield/B     # unit vector pointing i the magnetic field direction
#        @. pos = pos + b * stepsize
#        stepsforward += 1
#    end
#    # Backward
#    pos .= initpos
#    stepsbackward = -1 # We don't want to step if the the fist step is out of
#    # domain 
#    while all(pos .>= mesh.domain[:, 1]) & all(pos .<= mesh.domain[:, 2])
#        fields, _ = grid(mesh, interpolator, pos)
#        bfield, _ = fields
#        B = norm(bfield)
#        b = bfield/B
#        @. pos = pos - b * stepsize
#        stepsbackward += 1
#    end

    # Store line
    stepsforward = numsteps
    stepsbackward = numsteps
    lineforward  = zeros((mesh.numdims, stepsforward + 1))
    linebackward = zeros((mesh.numdims, stepsbackward + 1))
    lineforward[:, 1] .= initpos
    linebackward[:, 1] .= initpos

    # Follow field line forward
    for i = 1:stepsforward
        pos = lineforward[:, i]
        fields, _ = gridinterp(mesh, interpolator, pos)
        bfield, _ = fields
        B = norm(bfield)
        b = bfield/B
        #@. lineforward[:, i + 1] = lineforward[:, i] + b * stepsize
        lineforward[:, i + 1] .= scheme(pos,
                                        stepsize,
                                        Solvers.fieldtracingforward,
                                        mesh.bField,
                                        interpolator,
                                        mesh.xCoords,
                                        mesh.yCoords,
                                        mesh.zCoords
                                        )
    end
    # Follow field line backward
    for i = 1:stepsbackward
        pos = linebackward[:, i]
        fields, _ = gridinterp(mesh, interpolator, pos)
        bfield, _ = fields
        B = norm(bfield)
        b = bfield/B
        #@. linebackward[:, i + 1] = linebackward[:, i] - b * stepsize
        linebackward[:, i + 1] .= scheme(pos,
                                        stepsize,
                                        Solvers.fieldtracingbackward,
                                        mesh.bField,
                                        interpolator,
                                        mesh.xCoords,
                                        mesh.yCoords,
                                        mesh.zCoords
                                        )
    end

    # Fix points outside the domain
    for i = 1:stepsforward
        if any(lineforward[:,i] .< mesh.domain[:,1]) |
            any(lineforward[:,i] .> mesh.domain[:,2])
            if i == 1
                lineforward[:, i:end] .= lineforward[:,i]
            else
                lineforward[:, i:end] .= lineforward[:,i-1]
            end
            break
        end
    end
    for i = 1:stepsbackward
        if any(linebackward[:,i] .< mesh.domain[:,1]) |
            any(linebackward[:,i] .> mesh.domain[:,2])
            if i == 1
                linebackward[:, i:end] .= linebackward[:,i]
            else
                linebackward[:, i:end] .= linebackward[:,i-1]
            end
            break
        end
    end
    
    # Reverse linebackward and contatenate the lines
    fieldline = [linebackward[:, end:-1:1];; lineforward]
    return fieldline, lineforward, linebackward
end 


function plotfieldlines(mesh,
                        numlines::Integer,
                        stepsize,
                        interpolator,
                        scheme;
                        uniformsampling=false)
    fieldstrength = norm4(mesh.bField)
    maxfieldstrength = maximum(fieldstrength)
    # Create a function which returns the field strength at a given position.
    B(pos) = gridinterp(fieldstrength,
                        interpolator,
                        pos,
                        mesh.xCoords,
                        mesh.yCoords,
                        mesh.zCoords)
    # Sample initial position for magnetic field lines based on the field
    # strength 
    if uniformsampling
        # Hardcoded for nullpoint2k-experiment
        positions = zeros(3, numlines)
        positions[1,:] .= rand(mesh.xCoords[1], mesh.xCoords[end], numlines)
        positions[2,:] .= ones(numlines) .* 0.5
        positions[3,:] .= ones(numlines) .* (
            (mesh.zCoords[end] - mesh.zCoords[1])/2)
        display(positions)
    else
        positions, ratio = rejectionsampling(B,
                                             maxfieldstrength,
                                             numlines,
                                             mesh.domain)
    end
    # Assumes the lines are no greater than twice the maximum domain extent.
    maxextent = maximum(mesh.domain[:, 2] .- mesh.domain[:, 1]) 
    numsteps = Integer(2maxextent/stepsize)
    lines = zeros((mesh.numdims, numlines, 2*numsteps+2))
    for i = 1:numlines
        l, _, _ = generatefieldline(mesh,
                                    positions[:,i],
                                    stepsize,
                                    numsteps,
                                    interpolator,
                                    scheme
                                    )
        lines[:,i,:] .= l
    end
    plt.figure()
    plt.plot(lines[1,1,:], lines[2,1,:], color="black", linewidth=0.5)
    for i = 2:numlines
        plt.plot(lines[1,i,:], lines[2,i,:], color="black", linewidth=0.5)
    end
    return lines, positions
end
#|
function plotfieldlines(
    mesh,
    initpos::Matrix{T} where {T<:Real},
    stepsize,
    interpolator,
    scheme
    )
    numlines = size(initpos)[2]
    positions = initpos
    maxextent = maximum(mesh.domain[:, 2] .- mesh.domain[:, 1]) 
    numsteps = Integer(2maxextent/stepsize)
    lines = zeros((mesh.numdims, numlines, 2*numsteps+2))
    for i = 1:numlines
        l, _, _ = generatefieldline(mesh,
                                    positions[:,i],
                                    stepsize,
                                    numsteps,
                                    interpolator,
                                    scheme
                                    )
        lines[:,i,:] .= l
    end
    plt.figure()
    plt.plot(lines[1,1,:], lines[2,1,:], color="black", linewidth=0.5)
    for i = 2:numlines
        plt.plot(lines[1,i,:], lines[2,i,:], color="black", linewidth=0.5)
    end
    return lines, positions
end




#-----------#
# Utilities #
#-------------------------------------------------------------------------------
function setcartesianaxes!(
    ax    ::plt.PyCall.PyObject,
    normal::String="z",
    )
    if normal == "x"
        ax.set_xlabel(latexstring("\$y\$"))
        ax.set_ylabel(latexstring("\$z\$"))
    elseif normal == "y"
        ax.set_xlabel(latexstring("\$x\$"))
        ax.set_ylabel(latexstring("\$z\$"))
    elseif normal == "z"
        ax.set_xlabel(latexstring("\$x\$"))
        ax.set_ylabel(latexstring("\$y\$"))
    end
end # function setcartesianaxes


function setaxislimits!(
    ax    ::plt.PyCall.PyObject,
    domain::Matrix{T} where {T<:Real},
    normal::String="z",
    )
    if normal == "x"
        ax.set_xlim([domain[2,1], domain[2,2]])
        ax.set_ylim([domain[3,1], domain[3,2]])
    elseif normal == "y"
        ax.set_xlim([domain[1,1], domain[1,2]])
        ax.set_ylim([domain[3,1], domain[3,2]])
    elseif normal == "z"
        ax.set_xlim([domain[1,1], domain[1,2]])
        ax.set_ylim([domain[2,1], domain[2,2]])
    end
end # function setaxlimits!

#---------------------------------------------------#
# Messy plotting function hardcoded for RoCMI 2023  #
#-------------------------------------------------------------------------------
"""
    plotRoCMI(patch, labelling)
Plotting function used to make RoCMI-poster figures.
"""
function plotRoCMI(
    patch    ::Patch,
    labelling::Bool=false
    )
    times = collect(range(0.0, step=patch.dt, length=patch.numSteps+1))
    Ek = kineticenergy(patch.tp)
    plotplasmoid(
        patch.tp,
        times,
        Ek,
        patch.mesh.xCoords, patch.mesh.yCoords, patch.mesh.zCoords,
        patch.mesh.bField,
        patch.numParticles,
        labelling
    )
    numbins = 20
    p1 = plotenergydistr(patch.tp, 1, numbins, "Initial energy")
    p2 = plotenergydistr(patch.tp, patch.numSteps+1, numbins, "final energy")
    #p3 = plotenergydistr(patch.tp.vel[1, :, 1], numbins, "Vx-initial")
    #p4 = plotenergydistr(patch.tp.vel[1, :, end], numbins, "Vx-final")
    Plots.plot(p1,p2)#p3,p4)
end


"""
    plotRoCMI(patch, labelling)
Function for attempt to make nice subplots for the RoCMI-poster (Did not finish).
"""
function plotRoCMIsubplot(
    fbpatch  ::Patch,
    GCApatch ::Patch,
    labelling::Bool=false
    )
    times = collect(range(0.0, step=fbpatch.dt, length=fbpatch.numSteps+1))
    Ek = kineticenergy(fbpatch.tp)
    plotplasmoid(
        fbpatch.tp,
        GCApatch.tp,
        times,
        Ek,
        fbpatch.mesh.xCoords, fbpatch.mesh.yCoords, fbpatch.mesh.zCoords,
        fbpatch.mesh.bField,
        fbpatch.numParticles,
        labelling
    )
end


"""
    plotplasmoid(tp::TraceParticle, ...)
Function for plotting GCA or full orbit trajectories of the plasmoid-experiment
used to make figures for the RoCMI 2023 poster.
"""
function plotplasmoid(
    tp::TraceParticle,
    times,
    Ek,
    xx,
    yy, 
    zz,
    B,
    numparticles,
    labelling
    )

    pos = getpos(tp)
    # Energy plot
    plt.figure()
    for i = 1:numparticles
        if i == 2
            plt.plot(times[1:1904], Ek[i, 1:1904], label="Particle $i")
        elseif i == 3
            plt.plot(times[1:1744], Ek[i, 1:1744], label="Particle $i")
        elseif i == 4
            plt.plot(times[1:793], Ek[i, 1:793], label="Particle $i")
        else
            plt.plot(times, Ek[i, :], label="Particle $i")
        end
    end
    if typeof(tp) == ParticleSoA
        plt.title("Full-orbit: Kinetic energy")
    else
        plt.title("GCA: kinetic energy")
    end
    plt.xlabel("Time, s")
    plt.ylabel("Energy, J")
    if labelling
        plt.legend()
    end
    
    # Streamplot of magnetic field with initial positions and velocity
    #plt.figure()
    #plt.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
    #for i = 1:numparticles
    #    plt.plot(pos[1,i,1], pos[2,i,1], marker="o")
    #end
    #if typeof(tp) == ParticleSoA
    #    # Plot starting velocity as an arrow
    #    vel = getvel(tp)
    #    plt.quiver(pos[1, :, 1], pos[2, :, 1],
    #                  vel[1, :, 1], vel[2, :, 1],
    #                  width=0.003)
    #    plt.title("Full-orbit trajectories")
    #else
    #    plt.title("GCA trajectories")
    #end
    
    # Streamplot of magnetic field with particle trajectories
    babs = norm4(B) # Magnetic field strength
    plt.figure()
    #plt.axis("equal")
    #plt.contourf(xx, yy, transpose(babs[:,:,1]), levels=70,
    #                alpha=1.0,
    #                cmap=plt.get_cmap("Greys"))
    pbabs = copy(transpose(babs[:,:,1]))
    uuu = xx' .* ones(length(yy))
    vvv = ones(length(xx))' .* yy
    plt.pcolormesh(uuu, vvv, pbabs,
                    alpha=1.0,
                    cmap=plt.get_cmap("Greys"))
    cb = plt.colorbar(label=latexstring("\$|\\mathbf{B}|\$"),
                         ticks=[minimum(pbabs), maximum(pbabs)]
                         )
    println(minimum(pbabs))
    cb.ax.set_yticklabels(
        [latexstring("\$B_{min}\$"),latexstring("\$B_{max}\$")])
    plt.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]),
                      linewidth=0.3,
                      arrowsize=0.6,
                      color="black")
    plt.xlabel(latexstring("\$x\$"))
    plt.ylabel(latexstring("\$y\$"))
    cm = plt.get_cmap(:tab20)
    colorrange = (0:numparticles) ./ numparticles

    for i = 1:numparticles
        plotperiodictrajectory(pos[:,i,:], i, cm, colorrange)
    end
    for i = 1:numparticles
        # Mark position after t=1.5
        if i == 1
            #plt.plot(pos[1,1,1501], pos[2,1,1501], marker="^",
            #            color="black",label="t = 1.5s", linestyle="None") 
            plt.plot(pos[1,1,1], pos[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
        else
            #plt.plot(pos[1,i,1501], pos[2,i,1501], marker="^",
            #            color="black") 
            plt.plot(pos[1,i,1], pos[2,i,1], marker=".", color="Black")
        end
    end

    if typeof(tp) == ParticleSoA
        # Plot starting velocity as an arrow
        #vel = getvel(tp)
        #plt.quiver(pos[1, :, 1], pos[2, :, 1],
        #              vel[1, :, 1], vel[2, :, 1],
        #              width=0.003)
        plt.title("Full-orbit")
    else
        plt.title("GCA")
    end

    if labelling
        plt.legend()
    end

end # function plotplasmoid


"""
    plotplasmoid(tp::ParticleSoA, tp::GCAParticleSoA, ...)
Subplot both full orbit and GCA for comparison in RoCMI-poster.
"""
function plotplasmoid(
    fbtp::ParticleSoA,
    gcatp::GCAParticleSoA,
    times,
    Ek,
    xx,
    yy, 
    zz,
    B,
    numparticles,
    labelling
    )
    #fig = plt.figure(figsize=(10,6))
    fig = plt.figure()
    energyplot = fig.add_subplot(3,4,(9,12))
    # Energy plot
    for i = 1:numparticles
        if i == 2
            energyplot.plot(times[1:1904], Ek[i, 1:1904], label="Particle $i")
        elseif i == 3
            energyplot.plot(times[1:1744], Ek[i, 1:1744], label="Particle $i")
        elseif i == 4
            energyplot.plot(times[1:793], Ek[i, 1:793], label="Particle $i")
        else
            energyplot.plot(times, Ek[i, :], label="Particle $i")
        end
    end
    energyplot.set_title("Full-orbit: Kinetic energy")
    energyplot.set_xlabel("Time, s")
    energyplot.set_ylabel("Energy, J")
    if labelling
        energyplot.legend()
    end
    
    # Streamplot of magnetic field with particle trajectories
    fbplot = fig.add_subplot(3,4,(3))
    pos = fbplot.get_position()
    fbplot.set_position([pos.x0*2, pos.y0*2, pos.width, pos.height/2])
    gcaplot = fig.add_subplot(3,4,(1))
    pos = gcaplot.get_position()
    fbplot.set_title("Full-orbit")
    gcaplot.set_title("GCA")
    fbplot.set_xlabel(latexstring("\$x\$"))
    fbplot.set_ylabel(latexstring("\$y\$"))
    gcaplot.set_xlabel(latexstring("\$x\$"))
    # Full orbit
    babs = norm4(B) # Magnetic field strength
    pcm1 = fbplot.contourf(xx, yy, transpose(babs[:,:,1]),
                    alpha=1.0,
                    cmap=plt.get_cmap("Greys"))
    fbcb = plt.colorbar(pcm1,
#                         label=latexstring("\$|\\mathbf{B}|\$"),
                         ticks=[0, maximum(babs)],
                         ax=fbplot
                         )
    fbcb.ax.set_yticklabels(["0",latexstring("\$B_{max}\$")])
    fbplot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]),
                      linewidth=0.3,
                      arrowsize=0.6,
                      color="black")
    # GCA
    pcm2 = gcaplot.contourf(xx, yy, transpose(babs[:,:,1]),
                           alpha=1.0,
                           cmap=plt.get_cmap("Greys"))
    
    gcacb = plt.colorbar(pcm2,
#                         label=latexstring("\$|\\mathbf{B}|\$"),
                         ticks=[0, maximum(babs)],
                         ax=gcaplot
                         )
    gcacb.ax.set_yticklabels(["0",latexstring("\$B_{max}\$")])
    gcaplot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]),
                      linewidth=0.3,
                      arrowsize=0.6,
                      color="black")

    # Plot particle trajectories
    cm = plt.get_cmap(:tab20)
    colorrange = (0:numparticles) ./ numparticles
    for i = 1:numparticles
        plotperiodictrajectory(fbplot, fbtp.pos[:,i,:], i, cm, colorrange)
        plotperiodictrajectory(gcaplot, gcatp.R[:,i,:], i, cm, colorrange)
    end

    for i = 1:numparticles
        # Mark position after t=1.5
        if i == 1
            fbplot.plot(fbtp.pos[1,1,1], fbtp.pos[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
            gcaplot.plot(gcatp.R[1,1,1], gcatp.R[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
        else
            fbplot.plot(fbtp.pos[1,i,1], fbtp.pos[2,i,1], marker=".", color="Black")
            gcaplot.plot(gcatp.R[1,i,1], gcatp.R[2,i,1], marker=".", color="Black")
        end
    end

    if labelling
        fbplot.legend()
        gcaplot.legend()
    end

end # function plotplasmoid


function plotperiodictrajectoryRoCMI(pos, partidx, cm, colorrange)
    extent = [1.0, 1.0, 1.0] # Hardcoded for now
    posjumps = diff(pos, dims=2)
    posjumps = @. abs(posjumps)
    mask = @. isapprox(posjumps, extent, rtol=0.01)
    indices = findall(mask[1:2, :])
    numindices = length(indices)
    if numindices > 0
        j = 1
        for s = 1:length(indices)
            i = indices[s][2] # We only care about what time step the particle
            # hit the boundary, not which boundary was crossed.
            plt.plot(pos[1,j:i], pos[2,j:i], #color=cm(colorrange[partidx]),
                        label="Particle $partidx")                       
            j = i+1
            break
        end
    else
        plt.plot(pos[1,:], pos[2,:], #color=cm(colorrange[partidx]),
                    label="Particle $partidx")
    end
end 
#|
function plotperiodictrajectoryRoCMI(pyaxes, pos, partidx, cm, colorrange)
    extent = [1.0, 1.0, 1.0] # Hardcoded for now
    posjumps = diff(pos, dims=2)
    posjumps = @. abs(posjumps)
    mask = @. isapprox(posjumps, extent, rtol=0.01)
    indices = findall(mask[1:2, :])
    numindices = length(indices)
    if numindices > 0
        j = 1
        for s = 1:length(indices)
            i = indices[s][2] # We only care about what time step the particle
            # hit the boundary, not which boundary was crossed.
            pyaxes.plot(pos[1,j:i], pos[2,j:i], #color=cm(colorrange[partidx]),
                        label="Particle $partidx")                       
            j = i+1
            break
        end
    else
        pyaxes.plot(pos[1,:], pos[2,:], #color=cm(colorrange[partidx]),
                    label="Particle $partidx")
    end
end 



end # module TPplots
