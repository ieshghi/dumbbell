module sim_dumbbell
    export rk_evolve
    
    function extgpot(x,u)
        if x%10<9
            out = u/9;
        else
            out = -u;
        end
        return out
    end

    function sprgpot(x1,x2,k,x0) #potential gradient ON PARTICLE 1. For particle 2 take the negative of this
        return k*(x1-x2-x0);
    end


    function xdot_atherm(x1,x2,g,k,u,x0)
        force = [-sprgpot(x1,x2,k,x0) - extgpot(x1,u),sprgpot(x1,x2,k,x0)-extgpot(x2,u)];
        return force./g;
    end
    function xdot(x1,x2,g,k,u,t1,t2,x0,grav,dt)
        thermvec = sqrt(2g/dt).*[sqrt(t1).*randn(1)[1],sqrt(t2).*randn(1)[1]]
        force = [-sprgpot(x1,x2,k,x0) - extgpot(x1,u) + thermvec[1],sprgpot(x1,x2,k,x0)-extgpot(x2,u) + thermvec[2]]-grav.*ones(2);
        return force./g;
    end

    function rhs(x,t,dt,tt,grav)
        x1 = x[1];
        x2 = x[2];
        x0 = 0.1;
        g = 0.1;
        k = 0.1;
        u = 1;
        t2 = 1;
        t1 = 1 + tt*t2;

        #                                    #
        #   Change physical parameters here  #
        #                                    #

        return xdot(x1,x2,g,k,u,t1,t2,x0,grav,dt)
        #return xdot_atherm(x1,x2,g,k,u,x0)
    end


    function rk_evolve(x1i,x2i,dt,nsteps,tt,grav,substeps = 1)
        x = zeros(nsteps,2);
        t = zeros(nsteps);
        x[1,:] = [x1i,x2i];
        for i = 2:nsteps
            x[i,:],t[i] = rkint(x[i-1,:],t[i-1],dt,substeps,tt,grav);
        end
        return x,t
    end

    function rkint(x,t,dt,n,tt,grav)
        h = dt/n;
        tm = t;
        for i = 1:n
            k1 = h.*rhs(x,tm,h,tt,grav);
            k2 = h.*rhs((x+k1*1/2),(tm+h*1/2),h,tt,grav);
            k3 = h.*rhs((x+k2*1/2),(tm+h*1/2),h,tt,grav);
            k4 = h.*rhs((x+k3),(tm+h),h,tt,grav);
            x += 1/6*(k1 + 2*k2 + 2*k3 + k4);
            tm = tm+h;
        end
        ynew = x;
        tnew = tm;
        return ynew,tnew
    end

end

