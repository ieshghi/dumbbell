module sim_dumbbell
using DistributedArrays
using Distributed
using Statistics
using NPZ
using PyPlot
    export rk_evolve,manypart

    function extpot(x,u,pottype = 0)
        #Potential types are as follows:
        # 0: sawtooth
        # 1: four sines that look like a sawtooth
        # 2: quadratics that look like a sawtooth

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
            out = (-1)*(cos(pi*(x)/5)+cos(2*pi*(x)/5)/2+cos(3*pi*(x)/5)/3+cos(4*pi*(x)/5)/4+cos(5*pi*x/5)/5)/2;
        elseif pottype == 2
            z = ((x%10+10)%10);
            if z<9
                out = -2*u*(z-9)/81;
            else
                out = -2*u*(z-9);
            end
        end
        return out
    end

    function sprgpot(x1,x2,k) #potential gradient ON PARTICLE 1. For particle 2 take the negative of this
        return k*(x1-x2);
    end

    function xdot_atherm(x1,x2,g,k,u)
        force = [-sprgpot(x1,x2,k) - extgpot(x1,u),sprgpot(x1,x2,k)-extgpot(x2,u)];
        return force./g;
    end
    function xdot(x1,x2,g,k,u,t1,t2,grav,dt,pottype)
        thermvec = sqrt(2g/dt).*[sqrt(t1).*randn(1)[1],sqrt(t2).*randn(1)[1]]
        force = [-sprgpot(x1,x2,k) - extgpot(x1,u,pottype) + thermvec[1],sprgpot(x1,x2,k)-extgpot(x2,u,pottype) + thermvec[2]]-grav.*ones(2);
        return force./g;
    end

    function rhs(x,t,dt,tt,grav,pottype,tmin)
        x1 = x[1];
        x2 = x[2];
        g = 0.01;
        k = .01;
        u = 1.;
        t2 = tmin;
        t1 = tmin + tt;

        #                                    #
        #   Change physical parameters here  #
        #                                    #

        return xdot(x1,x2,g,k,u,t1,t2,grav,dt,pottype)
        #return xdot_atherm(x1,x2,g,k,u,x0)
    end

    function rk_evolve(x1i,x2i,dt,nsteps,tt,grav,pottype,tmin)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];
        for i = 2:nsteps
            x[i,:],t[i] = rkint(x[i-1,:],t[i-1],dt,tt,grav,pottype,tmin);
        end
        return x,t
    end

    function rkint(x,t,dt,tt,grav,pottype,tmin)
        h = dt;
        tm = t;
        k1 = h.*rhs(x,tm,h,tt,grav,pottype,tmin);
        k2 = h.*rhs((x+k1*1/2),(tm+h*1/2),h,tt,grav,pottype,tmin);
        k3 = h.*rhs((x+k2*1/2),(tm+h*1/2),h,tt,grav,pottype,tmin);
        k4 = h.*rhs((x+k3),(tm+h),h,tt,grav,pottype,tmin);
        x += 1/6*(k1 + 2*k2 + 2*k3 + k4);
        tm = tm+h;
        ynew = x;
        tnew = tm;
        return ynew,tnew
    end
    
    function manypart(N,x1i,x2i,dt,nsteps,tt,grav,pottype = 0,tmin = 1)
        xs = @DArray[rk_evolve(x1i,x2i,dt,nsteps,tt,grav,pottype,tmin)[1] for j=1:N];
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

    function getgrav(pot,tt,tmin)
        gs = LinRange(0,.05,40);
        sl = zeros(40,1);
        t = LinRange(0,10000,10000000);
        for i = 1:size(gs,1)
            print("Potential = $(pot), dT = $(tt), Tm = $(tmin), Gravity = $(gs[i])\n");
            q,xs = manypart(20,0,10,0.001,10000000,tt,gs[i],pot,tmin);
            sl[i],inter = linfit(t[1000:end],q[1000:end,2]);
        end
        return sl
    end


    function runsim()

    sl1_cold_smooth = zeros(40,4);

    tts = [0,2,4,6];

    for i = 1:4
        sl1_cold_smooth[:,i] = getgrav(0,tts[i],.1);
    end
    npzwrite("./data/sl1_cold_smooth.npz",sl1_cold_smooth);
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
