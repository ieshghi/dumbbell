module plotting
using Plots
pyplot()

function maketuples(a1,a2)
    s = size(a1);
    a = Array{Tuple{Float64,Float64},2}(undef,s[1],s[2]);
    for i = 1:s[1]
        for j = 1:s[2]
            a[i,j] = (a1[i,j],a2[i,j]);
        end
    end
    return a
end

function quiv_under(x,y,jx,jy,u,m,rat = 1)
    x = x[1:u:end,1:u:end];
    y = y[1:u:end,1:u:end];
    jx = m*jx[1:u:end,1:u:end];
    jy = m*jy[1:u:end,1:u:end];
    return quiver(x,y,quiver=maketuples(jx,jy),aspect_ratio = rat,arrow=arrow(.3,.1),linecolor=:steelblue);
end







end
