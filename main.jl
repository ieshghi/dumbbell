module sim_dumbbell
using DistributedArrays
using Distributed
using Statistics
using NPZ
    export rk_evolve,manypart

    function extpot(x,u,pottype = 0)
        #Potential types are as follows:
        # 0: sawtooth
        # 1: four sines that look like a sawtooth
        # 2: quadratics that look like a sawtooth
        # 3: flat potential

        if pottype == 0
            z = ((x%10+10)%10);
            if z<9
                out = u*z/9;
            else
                out = -u*(z-10);
            end
        elseif pottype == 1
            out = (sin(pi*(x)/5)+sin(2*pi*(x)/5)/4+sin(3*pi*(x)/5)/9+sin(4*pi*(x)/5)/16+sin(5*pi*x/5)/25+1)/2;
        elseif pottype == 2
            z = ((x%10+10)%10);
            if z<9
                out = u/81*(z-9)^2;
            else
                out = u*(z-9)^2;
            end
        elseif pottype == 3
            out = 0;
        end

        return out
    end
    
    function extgpot(x,u,pottype=0) #returns GRADIENTS of potentials
        #Potential types are as follows:
        # 0: sawtooth
        # 1: four sines that look like a sawtooth
        # 2: quadratics that look like a sawtooth

        if pottype == 0
            z = ((x%10+10)%10);
            if z<9
                out = u/9;
            else
                out = -u;
            end
        elseif pottype == 1
            out = (-1)*u*(cos(pi*(x)/5)+cos(2*pi*(x)/5)/2+cos(3*pi*(x)/5)/3+cos(4*pi*(x)/5)/4+cos(5*pi*x/5)/5)/2;
        elseif pottype == 2
            z = ((x%10+10)%10);
            if z<9
                out = -2*u*(z-9)/81;
            else
                out = -2*u*(z-9);
            end
        elseif pottype == 3
            out = 0;
        end
        return out
    end

    function sprgpot(x1,x2,k,l=0) #potential gradient ON PARTICLE 1. For particle 2 take the negative of this
        if l == 0 #l sets the maximum size of the spring, if it's 0 then the spring is infinitely flexible
            return k*(x1-x2);
        else        
            return k*(x1-x2)/abs((1-((x1-x2)/l)^2));
        end
    end

    function xdot(x1,x2,g,k,u,grav,pottype,l)
        force = [-sprgpot(x1,x2,k,l) - extgpot(x1,u,pottype) ,sprgpot(x1,x2,k,l)-extgpot(x2,u,pottype) ]-grav.*ones(2);
        return force./g;
    end

    function rhs(x,t,grav,pottype,g,t1,t2,l,k)
        x1 = x[1];
        x2 = x[2];
        tr = (t1+t2)/2;
        k = k*tr/50; #make the spring such that if l (l=10) is the period of the potential, then .5kl^2 = Tr
        
        u = 1.;

        #                                    #
        #   Change physical parameters here  #
        #                                    #
 
        return xdot(x1,x2,g,k,u,grav,pottype,l);
    end

    function rkint(x,t,dt,grav,pottype,g,l,k) 
        h = dt;
        tm = t;
        k1 = h.*rhs(x,tm,grav,pottype,g,t1,t2,l,k);
        k2 = h.*rhs((x+k1*1/2),(tm+h*1/2),grav,pottype,g,t1,t2,l,k);
        k3 = h.*rhs((x+k2*1/2),(tm+h*1/2),grav,pottype,g,t1,t2,l,k);
        k4 = h.*rhs((x+k3),(tm+h),grav,pottype,g,t1,t2,l,k);
        x += 1/6*(k1 + 2*k2 + 2*k3 + k4);
        tm = tm+h;
        ynew = x;
        tnew = tm;
        return ynew,tnew
    end


    function rk_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];

        #                                   #
        #             And Here              #
        #                                   #

        g = .1;
        b = [sqrt(2*t1*dt/g),sqrt(2*t2*dt/g)];
        for i = 2:nsteps
            dw = randn(2);
            x[i,:],t[i] = rkint(x[i-1,:],t[i-1],dt,grav,pottype,g,t1,t2,l,k) ;
            x[i,:] += b.*dw;
        end
        return x,t
    end
    function euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];

        #                                   #
        #             And Here              #
        #                                   #

        g = .1;
        b = [sqrt(2*t1*dt/g),sqrt(2*t2*dt/g)]; 
        for i = 2:nsteps
            a = rhs(x[i-1,:],t[i-1],grav,pottype,g,t1,t2,l,k);
            dw = randn(2);
            x[i,:] = x[i-1,:] + a*dt + b.*dw;            
        end
        return x,t
    end

    function manypart(N,x1i,x2i,dt,nsteps,t1,t2,grav,pottype = 0,l=0,k=1)
        xs = @DArray[rk_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k)[1] for j=1:N];
        q = sum(xs)./N;
        return q,xs
    end

    function manypart_euler(N,x1i,x2i,dt,nsteps,t1,t2,grav,pottype = 0,l=0,k=1)
        xs = @DArray[euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k)[1] for j=1:N];
        q = sum(xs)./N;
        return q,xs
    end

    function manypart_euler_dnv(N,x1i,x2i,dt,nsteps,t1,t2,grav,pottype = 0,l=0,k=1) #decorrelated noise variables
        xs = @DArray[euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k)[1] for j=1:N];
        q = sum(xs)./N;
        alpha = t2/(t1+t2);
        beta = t1/(t1+t2);
        
        xa = q[:,1];
        xb = q[:,2];

        r = xa-xb;
        R = alpha*xa+beta*xb;

        q_out = hcat(r,R);

        return q_out,xs
    end

    function linfit(xdat,ydat)
        X = zeros(size(xdat,1),2);
        X[:,1] = xdat;  X[:,2] = 1.0*ones(size(xdat,1),1);
        coeff_pred = X\ydat;

        slope = coeff_pred[1];
        intercept = coeff_pred[2];

        return slope, intercept
    end

    function gettscale(tmin)
#        dts = LinRange(0,1,100);
        dts = 1.;
        sl = zeros(100,1);
        t1 = tmin;
        t2 = tmin + dts;
        tim = LinRange(0,10000000*0.001,10000000);
        for i = 1:size(t2,1);
            print("T1 = $(t1), T2 = $(t2[i]), finding drift.\n");
            q,xs = manypart_euler(20,0,1,.001,10000000,t1[i],t2[i],0,0);
            sl[i],inter = linfit(tim[100:end],q[100:end,2]);
        end
        return sl
    end

    function getk(kvals)
        sl = zeros(size(kvals));
        t1 = .05;
        t2 = .9;
        tim = LinRange(0,10000000*0.001,10000000);
        for i = 1:size(kvals,1);
            print("T1 = $(t1), T2 = $(t2), k = $(kvals[i]).\n");
            q,xs = manypart_euler(50,0,1,.001,10000000,t1,t2,0,0,0,kvals[i]);
            sl[i],inter = linfit(tim[100:end],q[100:end,2]);
        end
        return sl

    end

    function getgrav(pot,tt,tmin,grange,whatpart=0)
        gi = grange[1];
        gf = grange[2];
        gs = LinRange(gi,gf,60);
        sl = zeros(60,1);
        
        maxtime = 10000;
        timestep = 0.001;
        nsteps = Int(maxtime/timestep);
        npart = 20;
        xinit = 0;
        dx_init = 1;
        t1 = tmin;
        t2 = tmin + tt;

        t = LinRange(0,maxtime,nsteps);
        for i = 1:size(gs,1)
            print("Potential = $(pot), dT = $(tt), Tm = $(tmin), Gravity = $(gs[i])\n");
            if whatpart == 0
                g = gs[i]
            elseif whatpart == 1
                g = gs[i].*[1,0]
            elseif whatpart == 2
                g = gs[i].*[0,1];
            end
            q,xs = manypart_euler(npart,xinit,dx_init,timestep,nsteps,t1,t2,g,pot);
            sl[i],inter = linfit(t[100:end],q[100:end,2]);
            tmin_wr = Int(floor(tmin*100));
            #npzwrite("./tempdata/pot$(pot)tm$(tmin_wr)dt$(tt)grav$(i).npz",q);
        end
        return sl
    end

    function heatmap(tmin,tmax,n,pot)
        t1arr = LinRange(tmin,tmax,n);
        t2arr = LinRange(tmin,tmax,n);
        out1 = zeros(n,n);        
        out2 = zeros(n,n);        
        maxtime = 10000;
        timestep = 0.001;
        nsteps = Int(maxtime/timestep);
        npart = 20;
        xinit = 0;
        dx_init = 1;
        t = LinRange(0,maxtime,nsteps);

        for i = 1:n
            for j = i:n
                t1 = t1arr[i];
                t2 = t2arr[j];
                q,xs = manypart_euler(npart,xinit,dx_init,timestep,nsteps,t1,t2,0,pot);
                out1[i,j],inter = linfit(t,q[:,1]);
                out2[i,j],inter = linfit(t,q[:,2]);
                print("Heatmap: Potential = $(pot), T1 = $(t1), T2 = $(t2)\n");
            end
        end
        return out1,out2
    end

    function runsim()
        sl_stiffspr = getk(LinRange(1,40,200));
        sl_softspr = getk(LinRange(.01,1,200));
        
        npzwrite("./data/sl_stiffspr.npz",sl_stiffspr);
        npzwrite("./data/sl_softspr.npz",sl_softspr);
    end

end
