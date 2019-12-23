using Makie

λ = slider(LinRange(  0.0,  5.0, 50), raw = true, camera = campixel!, start = 0.01)
δ = slider(LinRange(-0.99, 0.99, 50), raw = true, camera = campixel!)
n = slider(LinRange(    1,  100, 50), raw = true, camera = campixel!)

data = lift(n[end][:value]) do v
     map(LinRange(0, 2pi, 100)) do x
         4f0 .* Point2f0(sin(x) + (sin(x * v) .* 0.1), cos(x) + (cos(x * v) .* 0.1))
     end
end
p = scatter(data, markersize = s1[end][:value])

RecordEvents(
     hbox(p, vbox(s1, s2), parent = Scene(resolution = (500, 500))),
     "./docs/media/sliders"
)=#


s3[end][:value]
