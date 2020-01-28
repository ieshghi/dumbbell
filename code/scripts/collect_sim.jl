idmin = parse(Int,ARGS[1])
idmax = parse(Int,ARGS[2])

using HDF5
err = 0
for i = idmin:idmax
	id = string(i)
	try
		lt = h5read(string("../../data/hists/run",i,".h5"),"l")
		if @isdefined l
			l += lt
		else
			global l = deepcopy(lt)
		end
	catch
		print(string("run ",id," failed.\n"))
		global err+=1
	end
end
l ./= (length(idmin:idmax)-err)
x = h5read(string("/home/ie355/code/dumbbell/data/hists/run",idmin,".h5"),"x")
y = h5read(string("/home/ie355/code/dumbbell/data/hists/run",idmin,".h5"),"y")

h5open("/home/ie355/code/dumbbell/data/hists/collected.h5","w") do file
	write(file,"l",l)
	write(file,"x",x)
	write(file,"y",y)
end
