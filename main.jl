module sim_dumbbell
using DistributedArrays
using Distributed
using Statistics
using DelimitedFiles
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

    function sprgpot(x1,x2,k,l) #potential gradient ON PARTICLE 1. For particle 2 take the negative of this
            return k*(x1-x2-l)#l sets the rest size of the spring, if it's 0 then the spring is infinitely flexible
    end

    function xdot(x1,x2,g,k,u,grav,pottype,l)
        force = [-sprgpot(x1,x2,k,l) - extgpot(x1,u,pottype) ,sprgpot(x1,x2,k,l)-extgpot(x2,u,pottype) ]-grav.*ones(2);
        return force./g;
    end

    function rhs(x,t,grav,pottype,g,t1,t2,l,k,u)
        x1 = x[1];
        x2 = x[2];
        tr = (t1+t2)/2;
        k = k*tr/50; #make the spring such that if l (l=10) is the period of the potential, then .5kl^2 = Tr

        #                                    #
        #   Change physical parameters here  #
        #                                    #
 
        return xdot(x1,x2,g,k,u,grav,pottype,l*10);#multiply l by 10 because we've set the potential period to 10
    end

    function euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k,u)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];

        #                                   #
        #             And Here              #
        #                                   #

        g = .1;
        b = [sqrt(2*t1*dt/g),sqrt(2*t2*dt/g)]; 
        for i = 2:nsteps
            a = rhs(x[i-1,:],t[i-1],grav,pottype,g,t1,t2,l,k,u);
            dw = randn(2);
            x[i,:] = x[i-1,:] + a*dt + b.*dw;            
        end
        return x,t
    end

    function manypart_euler(N,x1i,x2i,dt,nsteps,t1,t2,grav = 0,pottype = 0,l=0,k=1,u=1)
        xs = @DArray[euler_evolve_sde(x1i,x2i,dt,nsteps,t1,t2,grav,pottype,l,k,u)[1] for j=1:N];
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

    function getuscale(umax,kmult = 1)
        u = LinRange(0,umax,100);
#        dts = 1.;
        sl = zeros(100,1);
        t1 = 0.1;
        t2 = 1.;
        tim = LinRange(0,1000000*0.01,1000000);
        for i = 1:size(u,1);
            print("U = $(u[i]), finding drift.\n");
            q,xs = manypart_euler(20,0,1,.01,1000000,t1,t2,0,2,0,kmult,u[i]);
            sl[i],inter = linfit(tim[100:end],q[100:end,2]);
        end
        return sl
    end

    function gettscale(tmin)
        dts = LinRange(0,1,100);
#        dts = 1.;
        sl = zeros(100,1);
        t1 = tmin;
        t2 = tmin + dts;
        tim = LinRange(0,10000000*0.001,10000000);
        for i = 1:size(t2,1);
            print("T1 = $(t1), T2 = $(t2[i]), finding drift.\n");
            q,xs = manypart_euler(20,0,1,.001,10000000,t1,t2[i],0,2);
            sl[i],inter = linfit(tim[100:end],q[100:end,2]);
        end
        return sl
    end

    function gettscale_grav(tmin,g)
        dts = LinRange(0,1,100);
#        dts = 1.;
        sl = zeros(100,1);
        t1 = tmin;
        t2 = tmin + dts;
        tim = LinRange(0,10000000*0.001,10000000);
        for i = 1:size(t2,1);
            print("T1 = $(t1), T2 = $(t2[i]), finding drift.\n");
            q,xs = manypart_euler(20,0,1,.001,10000000,t1,t2[i],g,2);
            sl[i],inter = linfit(tim[100:end],q[100:end,2]);
        end
        return sl
    end

    function gettscale_spring(tmin,kmult)
        dts = LinRange(0,1,100);
#        dts = 1.;
        sl = zeros(100,1);
        t1 = tmin;
        t2 = tmin + dts;
        tim = LinRange(0,10000000*0.001,10000000);
        for i = 1:size(t2,1);
            print("T1 = $(t1), T2 = $(t2[i]), finding drift.\n");
            q,xs = manypart_euler(50,0,1,.001,10000000,t1,t2[i],0,2,0,kmult);
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
            q,xs = manypart_euler(50,0,1,.001,10000000,t1,t2,0,2,0,kvals[i]);
            sl[i],inter = linfit(tim[100:end],q[100:end,2]);
        end
        return sl

    end

    function getl(lvals)
        sl = zeros(size(lvals));
        t1 = .1;
        t2 = 1;
        tim = LinRange(0,10000000*0.001,10000000);
        for i = 1:size(lvals,1);
            print("T1 = $(t1), T2 = $(t2), k = $(lvals[i]).\n");
            q,xs = manypart_euler(50,0,1,.001,10000000,t1,t2,0,2,lvals[i],1);
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
        maxtime = 10000;
        timestep = 0.001;
        
        #want to test dependencies at high kmult, so that the patching solution is correct:
#        sl_dt4 = gettscale_spring(0.1,4);
#        writedlm("./data/parab_strongk/sl_dt4",sl_dt4);
#        sl_dt3 = gettscale_spring(0.1,3);
#        writedlm("./data/parab_strongk/sl_dt3",sl_dt3);
#        sl_dt2 = gettscale_spring(0.1,2);
#        writedlm("./data/parab_strongk/sl_dt2",sl_dt2);
#        sl_dt1 = gettscale_spring(0.1,1);
#        writedlm("./data/parab_strongk/sl_dt1",sl_dt1);

        sl_u5 = getuscale(5,5);
        writedlm("./data/parab_strongk/sl_u5",sl_u5);
        sl_u4 = getuscale(5,4);
        writedlm("./data/parab_strongk/sl_u4",sl_u4);
        sl_u3 = getuscale(5,3);
        writedlm("./data/parab_strongk/sl_u3",sl_u3);
        sl_u2 = getuscale(5,2);
        writedlm("./data/parab_strongk/sl_u2",sl_u2);
        sl_u1 = getuscale(5,1);
        writedlm("./data/parab_strongk/sl_u1",sl_u1);
        
    end

end
