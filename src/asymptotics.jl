#####
# Some simple functions to assess the asymptotic 
# performance
#####


#Convenience wrapper
#probably exists somewhere in base.
const sqrt_2 = sqrt(2)
invnorm(x) = sqrt_2 * erfinv(2x-1)

#asymptotic var given by gaussian approximation
function asymptotic_var(vs, alphas, eps_)
	const phat, alpha0 = calc_phat(alphas)
	dot(vs, phat) + 
		invnorm(1-eps_)*sqrt(vs' * sigma(alphas) * vs)
end

#Ratios to optimality for near-optimal sets
kl_const(eps_) = sqrt(2*log(1/eps_)) / invnorm(1-eps_)
chisq_const(eps_) = sqrt(1/eps_ - 1) / invnorm(1-eps_)

#Ratios for the credible regions
#VG change these to all reference same functin...
kl_cov_const(eps_, d) = 
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)

chisq_cov_const(eps_, d) = 
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)

#the general covering ratio
cov_const(eps_, d) = 
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)
