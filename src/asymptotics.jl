#VG This should likely exist somewhere in Base
invnorm(x) = sqrt(2) * erfinv(2x-1)

function asymptotic_var(vs, alphas, eps_)
	const phat, alpha0 = calc_phat(alphas)
	dot(vs, phat) + 
		invnorm(1-eps_)*sqrt(vs' * sigma(alphas) * vs)
end

function kl_const(eps_)
	log_eps = log(1/eps_)
	denom = invnorm(1-eps_)
	sqrt(2log_eps) / denom
end

function chisq_const(eps_)
	kappa = sqrt(1/eps_ - 1)
	denom = invnorm(1-eps_)
	kappa / denom
end

#this is the ratio for the BenTal SEt
function kl_cov_const(eps_, d)
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)
end

function chisq_cov_const(eps_, d)
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)
end


#the general covering ratio
function cov_const(eps_, d)
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)
end
