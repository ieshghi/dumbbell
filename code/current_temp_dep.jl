module current_temp_dep
using HDF5
using Statistics
include("dual_current.jl")
function get_temps()
	l = [[] for i = 1:10]
	push!(l[1],h5read("../data/hists/collected_01.h5","l"))
	push!(l[2],h5read("../data/hists/collected_02.h5","l"))
	push!(l[3],h5read("../data/hists/collected_03.h5","l"))
	push!(l[4],h5read("../data/hists/collected_04.h5","l"))
	push!(l[5],h5read("../data/hists/collected_05.h5","l"))
	push!(l[6],h5read("../data/hists/collected_06.h5","l"))
	push!(l[7],h5read("../data/hists/collected_07.h5","l"))
	push!(l[8],h5read("../data/hists/collected_08.h5","l"))
	push!(l[9],h5read("../data/hists/collected_09.h5","l"))
	push!(l[10],h5read("../data/hists/collected_10.h5","l"))
	
	x = h5read("../data/hists/collected_10.h5","x")
	y = h5read("../data/hists/collected_10.h5","y")

	ts = 0.1*(1:10)
	jvals = zeros(10)
	for i = 1:10
		pars = [0.1,ts[i],1,0.7,0]
		jx,jy,ln,xn,y,fx,fy,lx,ly = dual_current.calc_curr(l[i][1],pars,x,y,1)
		jxo,xo,yo = dual_current.coarse_grain_2(jx,x,y,4)
		sumj = zeros(size(jxo)[1])
		dy = yo[2]-yo[1]
		for j = 1:size(jxo)[1]
			sumj[j] = dy*traprule(jxo[j,:])
		end
		jvals[i] = median(sumj)
	end

	return jvals
end

function traprule(a::Array{Float64,1})
	s = (a[1]+a[end])/2
	n = length(a)
	for i = 2:(n-1)
		s += a[i]
	end

	return s
end

end
