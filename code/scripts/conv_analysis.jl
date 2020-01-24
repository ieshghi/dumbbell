using HDF5

ifext = 0
convvals = [1000,500,200,100,50];
for nx = convvals
for i = idmin:idmax
	id = string(i)
	lt = h5read(string("../../data/hists/",nx,"run",i,"ifext",ifext,".h5"),"l")
	if @isdefined l
		l += lt
	else
		global l = deepcopy(lt)
	end
	
end
l /= length(idmin:idmax)
x = h5read(string("/home/ie355/code/dumbbell/data/hists/run",nx,"run",idmin,"ifext",ifext,".h5"),"x")
y = h5read(string("/home/ie355/code/dumbbell/data/hists/run",nx,"run",idmin,"ifext",ifext,".h5"),"y")

h5open("/home/ie355/code/dumbbell/data/hists/collected",nx,"ifext",ifext,".h5","w") do file
	write(file,"l",l)
	write(file,"x",x)
	write(file,"y",y)
end
end
