idmin = parse(Int,ARGS[1])
idmax = parse(Int,ARGS[2])

using HDF5

for i = idmin:idmax
	id = string(i)
	lt = h5read(string("../../data/hists/run",i,".h5"),"l")
	if @isdefined l
		l += lt
	else
		global l = deepcopy(lt)
	end
	
end
l /= length(idmin:idmax)
x = h5read(string("../../data/hists/run",idmin,".h5"),"x")
y = h5read(string("../../data/hists/run",idmin,".h5"),"y")

h5open("../../data/hists/collected.h5","w") do file
	write(file,"l",l)
	write(file,"x",x)
	write(file,"y",y)
end
