module sim_dumbbell
using DistributedArrays
using Distributed
using Statistics
using NPZ
    function gpot(x1,x2,a,l)
        z1 = x1-floor(x1);
        z2 = x2-floor(x2);
        
        if z1<=l
            p1 = -a/l;
        else
            p1 = a/(1-l);
        end
        if z2<=l
            p2 = -a/l;
        else
            p2 = a/(1-l);
        end

        return [p1,p2]
    end

    function rhs(x,pars)
        x1 = x[1];
        x2 = x[2];
        kbar = pars[1];
        abar = pars[2];
        lbar = pars[3];
        delta = pars[4];
        grav = pars[5];
        
        p = gpot(x1,x2,abar,lbar);
        return [-kbar*(x1-x2)-p[1],kbar*(x1-x2)-p[2]]-grav.*ones(2);
    end

    function evolve(x1i,x2i,dt,nsteps,pars)
        delta = pars[4];
 
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];

        b = [sqrt(1+delta),sqrt(1-delta)]; 
        for i = 2:nsteps
            a = rhs(x[i-1,:],pars);
            dw = randn(2);
            x[i,:] = x[i-1,:] + a*dt + b.*dw;            
        end
        return x,t
    end

    function manypart(N,x1i,x2i,dt,nsteps,pars)
        xs = @DArray[evolve(x1i,x2i,dt,nsteps,pars)[1] for j=1:N];
        q = sum(xs)./N;
        return q,xs
    end
    
    function linfit(xdat,ydat)
        X = zeros(size(xdat,1),2);
        X[:,1] = xdat;  X[:,2] = 1.0*ones(size(xdat,1),1);
        coeff_pred = X\ydat;

        slope = coeff_pred[1];
        intercept = coeff_pred[2];

        return slope, intercept
    end
    
    function heatmap(a,k)
        n = 19;
        darr = LinRange(.05,.95,n);
        larr = LinRange(.05,.95,n);
        nsteps = 1e6;
        dt = .001;

        npart = 20;
        xinit = 0;
        dx_init = 1;

        out = zeros(n,n);

        t = LinRange(0,nsteps*dt,nsteps);

        for i = 1:n
            for j = i:n
                d = darr[i];
                l = darr[j];
                pars = [k,a,l,d,0];
                q,xs = manypart_euler(npart,xinit,dx_init,dt,nsteps,pars);
                out[i,j],inter = linfit(t,(q[:,1]+q[:,2])/2);
                print("Heatmap: delta = $(d), l = $(l)");
            end
        end
        return out
    end
end
