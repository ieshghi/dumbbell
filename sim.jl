module sim
using DistributedArrays
using Distributed
using Statistics

#Parameters are to be fed into the pars array in the following order (already unitless)
# - T_1 (T_1 / U_0)
# - T_2 (T_2 / U_0)
# - k (k*l^2/U_0)
# - lambda (lambda/l)
# - gravity (g/(u0/gamma)^2/(l))
# - potential choice (1 = sawtooth, 2 = parabolas, ...)

function make_unitless(t1,t2,k,gamma,u0,lam,l,grav,pot,dt) #If you don't like the unitless guys this converts them for you
  t1bar = t1/u0;
  t2bar = t2/u0;
  kbar = k*l^2/u0;
  lambar = lam/l;
  gravbar = grav/(u0/(gamma*lam)); #gravity is basically just imposing a speed

  pars = [t1bar,t2bar,kbar,lambar,gravbar,pot];
  dtbar = dt*u0/(gamma*l^2);

  return pars,dtbar
end

function parallelcall(N,x0,dt,nsteps,pars,undersamp = 1) #simulates N particles in parallel and returns their average path
  x = @DArray[evolve(x0,dt,nsteps,pars,undersamp)[1] for j = 1:N];
  x_av = sum(x)./N;
  x_out = mean(x_av,dims = 2);
  t = LinRange(0,dt*nsteps,length(x_out));
  return t[1:end],x_out[1:end] #use undersamp if you want to plot subset of points
end

function evolve(x0,dt,nsteps,pars,undersamp) #Feed in unitless parameters here!
  nkeep = nsteps/undersamp;
  x = zeros(nkeep,2);
  t = zeros(nkeep);
  x[1,:] = x0;
  t[1] = 0;
  b = [sqrt.(2*dt*pars[1]),sqrt.(2*dt*pars[2])];
  dw = randn(2,nsteps);
  x_curr = x[1,:];
  j = 1;  
  for i = 2:nsteps
   a = rhs(x_curr,pars);
   x_curr += a*dt + b.*dw[:,i];
   if mod(i-1,undersamp)==0
       j +=1;
       t[j] = (i-1)*dt;
       x[j,:] = x_curr;
   end
  end
  return x,t
end

function rhs(x,pars)
  lam = pars[4];
  pot = pars[6];
  extforce = gpot.(x,lam,pot);
  k = pars[3];
  g = pars[5];
  sprforce = [k*(x[1]-x[2]),-k*(x[1]-x[2])];
  gforce = g*ones(2);

  return -(extforce + sprforce + gforce)
end

function gpot(x,l,pot)
  z = x-floor(x);
  z = ((x%1+1)%1);
  if pot == 1
    if z<l
      return 1/l;
    else
      return -1/(1-l);
    end
  elseif pot == 2
    if z < l
      return 2*z/l^2
    else
      return 2*(z-1)/(1-l)^2
    end
  else
    return nan
    print("Wrong choice of potential! Choose 1 or 2")
  end
end

end
