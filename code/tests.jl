module tests
include("base_sim.jl")
include("dual_current.jl")

function equil_sim_test(nvals) #currents should be 0 in equilibrium
	pars = [0.1,0.1,1.,0.7,0]
	err_sim = zeros(size(nvals))
	err_calc = zeros(size(nvals))
	err_all = zeros(size(nvals))
	@inbounds for i = 1:length(nvals)
		print(i)
		c,x,y = simbase.run_and_bin([0.,0.],.01,pars,10^7,nvals[i],nvals[i],5.0)
		jx,jy,cn,a,b,c,d = dual_current.calc_curr(c,pars,x,y)
		cb_tru = dual_current.boltzmann(x,y,pars)
		cb_tru = cb_tru.*maximum(cn)
		jx_calc,jy_calc,cb_calc,a,b,c,d = dual_current.calc_curr(cb_tru,pars,x,y)
		jx_tru = zeros(size(cb_calc))
		jy_tru = zeros(size(cb_calc))
	
		err_sim[i] = sqrt(sum((cn-cb_tru).^2)/sum(ones(size(cn))))
		err_calc[i] = sqrt(sum((jx_calc-jx_tru).^2+(jy_calc-jy_tru).^2)/sum(ones(size(cn))))
		err_all[i] = sqrt(sum((jx-jx_tru).^2+(jy-jy_tru).^2)/sum(ones(size(cn))))
	end
	return err_sim,err_calc,err_all
end



end
