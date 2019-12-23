import Makie: Contour, Scatter, Lines
using Colors

function visualize_RampE(xrange, parameters)#::TWOTEMPParameters)
    sub = Scene(resolution = (1000, 1000))
    x   = collect(xrange)
    #sub  = Scene(scene, IRect(50, 50, 900, 900))
    _visualizeRameE(x1) = RampE(x1, parameters)
    lines!(sub, xrange, _visualizeRameE.(x), color = :red, linewidth = 2)
    xmin = minimum(x)
    xmax = maximum(x)
    points = [Point2f0(xmin, 0) => Point2f0(xmax, 0)]
    linesegments!(sub, points, color = :black, linewidth = 2, limits = FRect(xmin, 0, xmax-xmin, 3))

    axis = sub[Axis]
    axis[:names, :axisnames] = ("x", "E (K)")
    axis[:names, :font] = ("Latin Modern Math", "Latin Modern Math")
    axis[:ticks, :font] = ("Latin Modern Math", "Latin Modern Math")
    ΔT = round(parameters.T₁ - parameters.T₂, digits = 3)
    text!(sub,
    "ΔT = $(ΔT), λ1 = $(parameters.λ₁), λ = $(parameters.λ), n = $(parameters.n) ",
    position = Point2f0(xmax, 3),
    align = (:right,  :bottom),
    textsize = 1,
    font = "Latin Modern Math", limits = FRect(xmin, 0, xmax-xmin, 3))
    sub
end

function visualize(xrange, yrange, parameters::TWOTEMPParameters) #This is the function that will visualize the background potential
    @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = parameters
    _visualizeHarmonicsE(x1, x2) = HarmonicsE(x1, x2, parameters)
    _visualizeFENEE(x1, x2)      = FENEE(x1, x2, parameters)
    _visualizeRameE(x1)          = RampE(x1, parameters)
    _visualizeRampE(x1, x2)      = RampE(x1, parameters) + RampE(x2, parameters)

    scene    = Scene(resolution = (1000, 1000))
    subRE    = Scene(scene, IRect(10, 100, 400, 400))
    subRE1   = Scene(scene, IRect(95,  50, 280, 50))
    subRE2   = Scene(scene, IRect(430, 175, 50, 280))

    #visualize the RampE
    if α  != 0
        contour!(subRE, xrange, yrange, _visualizeRampE.(xrange, yrange'), levels = 0, color = :viridis, linewidth = 0, fillrange = true)
        lines!(subRE1, xrange, _visualizeRameE.(xrange), color = :red, show_axis = false)
        lines!(subRE2, _visualizeRameE.(xrange), xrange, color = :red, show_axis = false)
        axis = subRE[Axis]
        axis[:names, :axisnames] = ("x₁", "x₂")
    end

    #visualize the HarmonicsE
    #=if r₀ == 0
        contour!(sub2, xrange, yrange, visualizeHarmonicsE.(xrange, yrange'), levels = 0, color = :viridis, linewidth = 0, fillrange = true)
        axis = sub2[Axis]
        axis[:names, :axisnames] = ("x₁", "x₂")
    end=#
    display(scene)
end

function visualize_paths(paths; N = 1) #visualize named tupule of the form (t, anything₁, anything₂)
    colors = :tomato
    scene    = Scene(resolution = (500, 500))
    for i = 1:N
        lines!(scene, paths.x₁[:,i], paths.x₂[:,i], color = colors, linewidth = 2)
    end
    axis = scene[Axis]
    axis[:names, :axisnames] = ("x₁ (μm)", "x₂ (μm)")
    display(scene)
end

function visualize_xv(v, x; N = 1)
    colors = :tomato
    scene  = Scene(resolution = (1000, 1000))
    subv₁  = Scene(scene, IRect( 50,   50, 400, 400))
    subv₂  = Scene(scene, IRect( 50,  550, 400, 400))
    subx₁  = Scene(scene, IRect( 550,  50, 400, 400))
    subx₂  = Scene(scene, IRect( 550, 550, 400, 400))

    for i = 1:N
        lines!(subv₁, v.t[:,i], v.v₁[:,i], color = colors, linewidth = 2)
        lines!(subv₂, v.t[:,i], v.v₂[:,i], color = colors, linewidth = 2)
        lines!(subx₁, x.t[:,i], x.x₁[:,i], color = colors, linewidth = 2)
        lines!(subx₂, x.t[:,i], x.x₂[:,i], color = colors, linewidth = 2)
    end

    axis = subv₁[Axis]
    axis[:names, :axisnames] = ("t (s)", "v₁ (μm/s)")
    axis = subv₂[Axis]
    axis[:names, :axisnames] = ("t (s)", "v₂ (μm/s)")
    axis = subx₁[Axis]
    axis[:names, :axisnames] = ("t (s)", "x₁ (μm)")
    axis = subx₂[Axis]
    axis[:names, :axisnames] = ("t (s)", "x₂ (μm)")

    display(scene)

end

function visualize_vparametric(v, solution) #for now only work on 1 trajectory per prblem instance

    scene     = Scene(resolution = (650, 550))
    main      = Scene(scene, IRect(   25,  0, 500, 500))
    scalebar  = Scene(scene, IRect( 550,  125, 50, 400))
    a = reverse(to_colormap(:RdBu , size(solution, 1)))

    for i = 1:size(solution, 1)
        lines!(main, v[i][1], v[i][2][:,1], color = a[i], linewidth = 2)[end]
    end

    ΔT  = [solution[i].ΔT for i = 1:size(solution, 1)]
    hbar(scalebar, ΔT; colormaps = reverse(to_colormap(:RdBu, size(solution, 1))));
    scene.children[2].plots[1][:ticks, :textsize] = (5, 10)
    scene
end

function visualize_drift(drift)
    scene     = Scene(resolution = (650, 1100))

    main1     = Scene(scene, IRect(  25,   25, 500, 500))
    scalebar1 = Scene(scene, IRect( 550,  200,  50, 400))

    main2     = Scene(scene, IRect(   25, 575, 500, 500))
    scalebar2 = Scene(scene, IRect(  550, 725,  50, 400))

    #a = reverse(to_colormap(:RdBu , size(, 1)))

    ΔT = [drift[i].ΔT for i = 1:size(drift, 1)]
    v₁ = [drift[i].vstat₁[1] for i = 1:size(drift, 1)]
    v₂ = [drift[i].vstat₂[1] for i = 1:size(drift, 1)]
    tmin = minimum(Δ)
    tmax = maximum(Δ)
    scatter!(main1, Δ, (v₁ .+ v₂)/2, markersize = 10, color = reverse(to_colormap(:RdBu, size(drift, 1))), limits =  FRect(tmin, -0.25, tmax-tmin, 0.5))
    axis = main1[Axis]
    axis[:names, :axisnames] = ("ΔT", "⟨dR/dt⟩")
    axis[:names, :font] = ("Latin Modern Math", "Latin Modern Math")
    axis[:ticks, :font] = ("Latin Modern Math", "Latin Modern Math")
    #hbar(scalebar1, ΔT; colormaps = reverse(to_colormap(:RdBu, size(drift, 1))));

    scatter!(main2, Δ, (v₁ .- v₂)/2, markersize = 10, color = reverse(to_colormap(:RdBu, size(drift, 1))), limits =  FRect(tmin, -0.5, tmax-tmin, 1))
    axis = main2[Axis]
    axis[:names, :axisnames] = ("ΔT", "⟨dr/dt⟩")
    axis[:names, :font] = ("Latin Modern Math", "Latin Modern Math")
    axis[:ticks, :font] = ("Latin Modern Math", "Latin Modern Math")
    hbar(scalebar2, ΔT; colormaps = reverse(to_colormap(:RdBu, size(drift, 1))));




    scene.children[2].plots[1][:ticks, :textsize] = (5, 10)
    scene.children[2].plots[1][:ticks, :font] = ("Latin Modern Math", "Latin Modern Math")
    scene
end

function hbar(scene, value; colormaps = :viridis)
    x = LinRange(0, 0.1, 10)
    y = LinRange(minimum(value), maximum(value), length(value))
    z = ones(size(x)) .* y'
    v = heatmap!(scene, x, y, z, interpolate = true, fxaa = false, colormap = colormaps, axis = (names = (axisnames = ("", ""),),))
    v[1][:ticks, :textcolor] = (:white, :black)
    return v
end

#animatation_realspace(solutionIndependentParticle[], "IndependentParticle.gif", solutionHarmonics[1].parameter )
function animatation_realspace(solution, filename, parameter; skip = 0)
    mytime = Node(1)
    if solution.T₁ >= solution.T₂
        x₁ = solution.x₁
        x₂ = solution.x₂
    else #This is just for colorschemewise
        x₁ = solution.x₂
        x₂ = solution.x₁
    end
    maxlim = 2*floor(maximum([abs(maximum(x₁)), abs(maximum(x₂)), abs(minimum(x₁)), abs(minimum(x₂))]))
    x = collect(-maxlim:0.01:maxlim)

    scene  = visualize_RampE(x, parameter)
    linesegments!(scene, lift(t -> [Point2(x₁[t],0) => Point2(x₂[t],0)], mytime), color = :tan1, linewidth = 3, raw = true, limits = scene[:limits].val)
    scatter!(scene, lift(t -> [x₁[t]], mytime), [0], color = :tomato, markersize = 0.5, limits = scene[:limits].val)
    scatter!(scene, lift(t -> [x₂[t]], mytime), [0], color = :royalblue4, markersize = 0.5, limits = scene[:limits].val)
    #Project Upward
    scatter!(scene, lift(t -> [x₁[t]], mytime), lift(t -> [RampE(x₁[t], parameter)], mytime), color = :tomato, markersize = 0.5 , limits = scene[:limits].val)
    scatter!(scene, lift(t -> [x₂[t]], mytime), lift(t -> [RampE(x₂[t], parameter)], mytime), color = :royalblue4, markersize = 0.5, limits = scene[:limits].val)

    if skip != 0
        record(scene, filename, collect(1:skip:length(x₁))) do i
            mytime[] = i
        end
    else
        record(scene, filename, collect(1:length(x₁))) do i
            mytime[] = i
        end
    end

end
