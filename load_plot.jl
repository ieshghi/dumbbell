module load_plot
using JLD
using PyPlot
using SpecialFunctions

#Define font sizes to be uniform
labelfont = 15;
titlefont = 20;
legendfont = 13;
fontnm = "Liberation Sans"
rc("font",family = fontnm)

function make_plots()

#First, make the individual paths plot

dat = load("/home/ieshghi/Documents/code/dumbbell/data/individual_paths/basic.jld");
temp_vals = dat["parvals"].-0.1;
time = dat["t"];
trajs = -dat["xvals"];
figure(1);

for i = 1:length(temp_vals)
	lab = string("\$\\delta T = ",string(temp_vals[i]),"\$");
	plot(time,trajs[:,i],"-",label = lab);
end
xlabel("t",fontsize = labelfont,fontname = fontnm)
ylabel("x",fontsize = labelfont,fontname = fontnm)
text(250,840,L"T_{min} = 0.1",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
xticks(0:3000:maximum(time)+1)
yticks(0:400:maximum(trajs)+1)
tight_layout();

#Next, the large kappa plot
dat = load("/home/ieshghi/Documents/code/dumbbell/data/k_runs/stiff.jld");
kappa = dat["parvals"];
trajs = -dat["vvals"];

figure(2);
showvals = kappa.>10;
plot(kappa,trajs,".",label = "data");
plot(kappa[showvals],kappa[showvals].^(-1/2)/2.5,label = L"\kappa^{-1/2}");

xlabel(L"\kappa",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
text(0.95,0.105,"\$ T_{min}  = 0.1\$ \n\$T_{max} = 1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
xscale("log")
yscale("log")
xticks([1,10])
yticks([0.1])
minorticks_off()
tight_layout();

#Next, the small kappa plot
dat = load("/home/ieshghi/Documents/code/dumbbell/data/k_runs/soft.jld");
kappa = dat["parvals"];
trajs = -dat["vvals"];

figure(3);
showvals = kappa.>0.5;
plot(kappa,trajs,".",label = "data");
plot(kappa[showvals],kappa[showvals].^(1)/30,label = L"\kappa^{1}");

xlabel(L"\kappa",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
text(0.047,0.015,"\$ T_{min}  = 0.1\$ \n\$T_{max} = 1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
xscale("log")
yscale("log")
#xticks([1,10])
#yticks([0.1])
minorticks_off()
tight_layout();

#Next, the temperature difference plot
dat = load("/home/ieshghi/Documents/code/dumbbell/data/temp_runs/stiffer.jld");
temp = dat["parvals"].-0.1;
trajs = -dat["vvals"];

figure(4);
plot(temp,trajs,".",label = "data");
plot(temp,patching_soln.(0.1,(temp.+0.1),10,0.1));

xlabel(L"\delta T",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
#text(0,.0265,"\$ T_{min}  = 0.1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
#xticks([1,10])
#yticks([0,.01,.02])
tight_layout();

end

function patching_soln(t1,t2,k,lam)
	k = k*(t1+t2)/2
	j = 1/(pi*erf(2/sqrt(t1+t2)))*2*(exp(-4/(t1+t2))*k*(t1-t2))*((lam-1)*sqrt((1+k*(lam-1)^2)/(16*t1*t2+16*k*t1*t2*(1-lam)^2+k^2*(t1+t2)^2*(1-lam)^4))+lam*sqrt((1+k*lam^2)/(16*t1*t2+16*k*t1*t2*lam^2+k^2*(t1+t2)^2*lam^4)))
	return j
end

end
