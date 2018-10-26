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
                out = -u*z/9+1;
            else
                out = u*(z-10)+1;
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

    function sprgpot(x1,x2,k) #potential gradient ON PARTICLE 1. For particle 2 take the negative of this
        return k*(x1-x2);
    end

    function xdot(x1,x2,g,k,u,grav,pottype)
        force = [-sprgpot(x1,x2,k) - extgpot(x1,u,pottype) ,sprgpot(x1,x2,k)-extgpot(x2,u,pottype) ]-grav.*ones(2);
        return force./g;
    end

    function rhs(x,t,grav,pottype,g)
        x1 = x[1];
        x2 = x[2];
        k = .01;
        u = 1.;

        #                                    #
        #   Change physical parameters here  #
        #                                    #
 
        return xdot(x1,x2,g,k,u,grav,pottype);
    end

    function rkint(x,t,dt,grav,pottype,g) #Naive Runge-Kutta. Seems to be unreliable.
        h = dt;
        tm = t;
        k1 = h.*rhs(x,tm,grav,pottype,g);
        k2 = h.*rhs((x+k1*1/2),(tm+h*1/2),grav,pottype,g);
        k3 = h.*rhs((x+k2*1/2),(tm+h*1/2),grav,pottype,g);
        k4 = h.*rhs((x+k3),(tm+h),grav,pottype,g);
        x += 1/6*(k1 + 2*k2 + 2*k3 + k4);
        tm = tm+h;
        ynew = x;
        tnew = tm;
        return ynew,tnew
    end


    function rk_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];

        #                                   #
        #             And Here              #
        #                                   #

        g = .01;
        #b = [sqrt(2*t1*dt/g),sqrt(2*t2*dt/g)]; This is what it should be 
        b = [sqrt(t1*dt/g),sqrt(t2*dt/g)]; #This gives the right answer...
        for i = 2:nsteps
            #a = rhs(x[i-1,:],t[i-1],grav,pottype,g);
            dw = randn(2);
            #x[i,:] = x[i-1,:] + a*dt + b.*dw;            
            x[i,:],t[i] = rkint(x[i-1,:],t[i-1],dt,grav,pottype,g) ;
            x[i,:] += b.*dw;
        end
        return x,t
    end
    function euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];

        #                                   #
        #             And Here              #
        #                                   #

        g = .01;
        b = [sqrt(t1*dt/g),sqrt(t2*dt/g)]; #This gives the right answer...
        for i = 2:nsteps
            a = rhs(x[i-1,:],t[i-1],grav,pottype,g);
            dw = randn(2);
            x[i,:] = x[i-1,:] + a*dt + b.*dw;            
        end
        return x,t
    end

    function manypart(N,x1i,x2i,dt,nsteps,t1,t2,grav,pottype = 0)
        xs = @DArray[rk_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype)[1] for j=1:N];
        q = sum(xs)./N;
        return q,xs
    end

    function manypart_euler(N,x1i,x2i,dt,nsteps,t1,t2,grav,pottype = 0)
        xs = @DArray[euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype)[1] for j=1:N];
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

    function getgrav(pot,tt,tmin,grange)
        gi = grange[1];
        gf = grange[2];
        gs = LinRange(gi,gf,40);
        sl = zeros(40,1);
        
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
            q,xs = manypart_euler(npart,xinit,dx_init,timestep,nsteps,t1,t2,gs[i],pot);
            sl[i],inter = linfit(t[100:end],q[100:end,2]);
            tmin_wr = Int(tmin*10);
            npzwrite("./tempdata/pot$(pot)tm$(tmin_wr)dt$(tt)grav$(i).npz",q);
        end
        return sl
    end

    function heatmap(tmin,tmax,n,pot)
        t1arr = LinRange(tmin,tmax,n);
        t2arr = LinRange(tmin,tmax,n);
        out = zeros(n,n);        
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
                out(i,j) = linfit(t,mean(q,2));
                print("Heatmap: Potential = $(pot), T1 = $(t1), T2 = $(t2)\n");
            end
        end
    end

    function runsim()

    sl1_hot_smooth = zeros(40,4);
    sl1_cold_smooth = zeros(40,4);
    tts = [0,2,4,6];

    for i = 1:4
        sl1_hot_smooth[:,i] = getgrav(0,tts[i],.5,[0,.007]);
        sl1_cold_smooth[:,i] = getgrav(0,tts[i],0.1,[0,.02]);
    end
    npzwrite("./data/sl1_hot_smooth.npz",sl1_hot_smooth);
    npzwrite("./data/sl1_cold_smooth.npz",sl1_cold_smooth);
    heatmap_1 = heatmap(0,1.2,10,0);
    npzwrite("./data/heatmap_1.npz",heatmap_1);
    end

    function plotres(pot,hot,tt)
     if pot == 0
         a = npzread("./data/sl1_hot.npz");
         b = npzread("./data/sl1_cold.npz");
     elseif pot == 1
         a = npzread("./data/sl2_hot.npz");
         b = npzread("./data/sl2_cold.npz");
     elseif pot == 2
         a = npzread("./data/sl3_hot.npz");
         b = npzread("./data/sl3_cold.npz");
     end
     if hot == 1
        data = a;
     elseif hot == 0
        data = b;
     end

     gs = LinRange(0,0.3,40);

     if tt == 0
         ind = 1;
     elseif tt == 1
         ind = 2;
     elseif tt == 2
         ind = 3;
     elseif tt == 4
         ind = 4;
     elseif tt == 5
         ind = 5;
     elseif tt == 6
         ind = 6;
     end
     plot(gs,data[:,ind],label="dT = %.2f"%(tt))
     ylabel("Slope")
     xlabel("Gravity")
     
    end

end
