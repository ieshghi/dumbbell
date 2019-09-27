module histutils

function genmesh(x,y)
	n = size(x)[1];
	m = size(y)[1];

	pairs = zeros(n*m,2);
	k = 1
	for i = 1:n
		for j = 1:m
			pairs[k,:] = [x[i],y[j]];
			k += 1;
		end
	end
	hist = zeros(n*m);
	return pairs,hist
end

function placeinmesh(p,pairs,hist,l)
	n = size(pairs)[1];
	dx = pairs[2,1]-pairs[1,1];
	dy = pairs[2,2]-pairs[1,2];
	for i = 1:n
		a = mod(p[1]+1-l,1);
		b = p[2]-p[1];
		if ((a-pairs[i,1])<=dx)&&((b-pairs[i,2])<=dy)
			hist[i] += 1;
			break
		end
	end
	return hist
end

function shape_prob(xi,yi,z)
	x = unique(xi);
	y = unique(yi);
	xmesh = [i for j=y,i=x];
	ymesh = [j for j=y,i=x];
	zmesh = zeros(size(xmesh));
	for i = 1:(size(x)[1])
		for j = 1:(size(y)[1])
			zmesh[j,i] = z[(xi.==x[i]).*(yi.==y[j])][1];
		end
	end
	return xmesh,ymesh,zmesh
end

end
