module load_plot
using JLD
using DelimitedFiles
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
#savefig("figures/figure_1.eps", bbox_inches="tight");
close(1)

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
#savefig("figures/figure_2.eps", bbox_inches="tight")
close(2)

#Next, the small kappa plot
dat = load("/home/ieshghi/Documents/code/dumbbell/data/k_runs/soft.jld");
kappa = dat["parvals"];
trajs = -dat["vvals"];

figure(3);
showvals = (kappa.>0.1).*(kappa.<0.8);
plot(kappa,trajs,".",label = "data");
#plot(kappa[showvals],kappa[showvals].^(1)/30,label = L"\kappa^{1}");
plot(kappa[showvals],kappa[showvals].^(3/2)/20,label = L"\kappa^{3/2}");

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
#savefig("figures/figure_3.eps", bbox_inches="tight")
close(3)

#Next, the temperature difference plot
#dat = readdlm("/home/ieshghi/Documents/code/dumbbell/data/archive/parab_strongk/sl_dt1");
#temp = LinRange(0,1,100);

dat1 = load("/home/ieshghi/Documents/code/dumbbell/data/temp_runs/colder.jld");
dat2 = load("/home/ieshghi/Documents/code/dumbbell/data/temp_runs/stiffer.jld");
dat3 = load("/home/ieshghi/Documents/code/dumbbell/data/temp_runs/stifferer.jld");
temp1 = dat1["parvals"].-0.1;
temp2 = dat2["parvals"].-0.1;
temp3 = dat3["parvals"].-0.1;
trajs1 = -dat1["vvals"];
trajs2 = -dat2["vvals"];
trajs3 = -dat3["vvals"];

figure(4);
plot(temp1,trajs1,".",label = L"\kappa = 1");
#plot(temp,dat,".",label = L"\kappa = 1");
plot(temp2,trajs2,".",label = L"\kappa = 10");
plot(temp3,trajs3,".",label = L"\kappa = 20");
plot(temp1,patching_soln.(0.1,(temp1.+0.1),1,0.1));
plot(temp2,patching_soln.(0.1,(temp2.+0.1),50,0.1)*3);
plot(temp3,patching_soln.(0.1,(temp3.+0.1),20,0.1));

xlabel(L"\delta T",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
text(0.01,.09,"\$ T_{min}  = 0.1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
#xticks([1,10])
#yticks([0,.01,.02])
#xscale("log")
#yscale("log")
tight_layout();
#savefig("figures/figure_4.eps", bbox_inches="tight")
#close(4)

#next, cold gravity plot with small gravity values
dat1 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_cold_dt0.jld");
dat2 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_cold_dt1.jld");
dat3 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_cold_dt2.jld");
dat4 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_cold_dt3.jld");
grav1 = dat1["parvals"];
grav2 = dat2["parvals"];
grav3 = dat3["parvals"];
grav4 = dat4["parvals"];
trajs1 = -dat1["vvals"];
trajs2 = -dat2["vvals"];
trajs3 = -dat3["vvals"];
trajs4 = -dat4["vvals"];

figure(5);
plot(grav1,trajs1,".",label = L"\delta T = 0");
plot(grav2,trajs2,".",label = L"\delta T = 1");
plot(grav3,trajs3,".",label = L"\delta T = 2");
plot(grav4,trajs4,".",label = L"\delta T = 3");

xlabel(L"g",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
text(0.01,.37,"\$ T_{min}  = 0.1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
#xticks([1,10])
#yticks([0,.01,.02])
tight_layout();
#savefig("figures/figure_5.eps", bbox_inches="tight")
close(5)

#Next, cold gravity plot with large gravity values

dat1 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/large_cold_dt0.jld");
dat2 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/large_cold_dt1.jld");
dat3 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/large_cold_dt2.jld");
dat4 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/large_cold_dt3.jld");
grav1 = dat1["parvals"];
grav2 = dat2["parvals"];
grav3 = dat3["parvals"];
grav4 = dat4["parvals"];
trajs1 = -dat1["vvals"];
trajs2 = -dat2["vvals"];
trajs3 = -dat3["vvals"];
trajs4 = -dat4["vvals"];

figure(6);
plot(grav1,trajs1,"-",label = L"\delta T = 0");
plot(grav2,trajs2,"-",label = L"\delta T = 1");
plot(grav3,trajs3,"-",label = L"\delta T = 2");
plot(grav4,trajs4,"-",label = L"\delta T = 3");
plot(grav4,grav4,"k--")
xlabel(L"g",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
text(-0.2,5.8,"\$ T_{min}  = 0.1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
#xticks([1,10])
#yticks([0,.01,.02])
#xlim([-0.5,6]);
tight_layout();
#savefig("figures/figure_6.eps", bbox_inches="tight")
close(6)

#Finally, warm gravity plot with small gravity values

dat1 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_hot_dt0.jld");
dat2 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_hot_dt1.jld");
dat3 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_hot_dt2.jld");
dat4 = load("/home/ieshghi/Documents/code/dumbbell/data/grav_runs/small_hot_dt3.jld");
grav1 = dat1["parvals"];
grav2 = dat2["parvals"];
grav3 = dat3["parvals"];
grav4 = dat4["parvals"];
trajs1 = -dat1["vvals"];
trajs2 = -dat2["vvals"];
trajs3 = -dat3["vvals"];
trajs4 = -dat4["vvals"];

figure(7);
plot(grav1,trajs1,".",label = L"\delta T = 0");
plot(grav2,trajs2,".",label = L"\delta T = 1");
plot(grav3,trajs3,".",label = L"\delta T = 2");
plot(grav4,trajs4,".",label = L"\delta T = 3");
xlabel(L"g",fontsize = labelfont,fontname = fontnm)
ylabel("v",fontsize = labelfont,fontname = fontnm)
text(0.00025,.014,"\$ T_{min}  = 0.1\$",verticalalignment = "top",fontsize = legendfont)
legend(fontsize = legendfont,frameon = 0)
tick_params(labelsize = labelfont)
#xticks([1,10])
#yticks([0,.01,.02])
ylim([-.023,.03])
tight_layout();
#savefig("figures/figure_7.eps", bbox_inches="tight")
close(7)

end

function patching_soln(t1,t2,k,lam)
	t2 = t2
	j = 1/(pi*erf(2/sqrt(t1+t2)))*2*(exp(-4/(t1+t2))*k*(t1-t2))*((lam-1)*sqrt((1+k*(lam-1)^2)/(16*t1*t2+16*k*t1*t2*(1-lam)^2+k^2*(t1+t2)^2*(1-lam)^4))+lam*sqrt((1+k*lam^2)/(16*t1*t2+16*k*t1*t2*lam^2+k^2*(t1+t2)^2*lam^4)))
	return j
end

#function patching_soln(t1,t2,k,lam)
#	j = -exp(-4/(t1+t2))*(t1-t2)*(t1^1-14*t1*t2+t2^2)*(2*lam-1)/(pi*sqrt(k)*erf(2/sqrt(t1+t2))*(t1+t2)^3*(lam-1)^2*lam^2);
#	return j
#end

end
