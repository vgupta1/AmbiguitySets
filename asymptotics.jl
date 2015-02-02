
function calc_phat(alphas) 
	const alpha0 = sum(alphas)
	alphas/alpha0, alpha0
end

#This is sigma from document
function unscaled_sigma_(phat)
	const d = length(phat)
	sigma = zeros(Float64, d, d)
	for i = 1:d
		sigma[i, i] = phat[i]*(1-phat[i])
		for j = i+1:d
			sigma[i, j] = -phat[i] * phat[j]
			sigma[j, i] = sigma[i, j]
		end
	end
	sigma
end

#this is the finite sample scaling, equiv to assymptotic
#preferred usage
function sigma(alphas)
	const phat, alpha0 = calc_phat(alphas)
	unscaled_sigma_(phat) / (alpha0 + 1)
end

#VG This should likely exist somewhere
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

function mom_const(eps_)
	kappa = sqrt(1/eps_ - 1)
	denom = invnorm(1-eps_)
	kappa / denom
end

#this is the ratio for the BenTal SEt
function kl_chisq_const(eps_, d)
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)
end

function chisq_const(eps_, d)
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)
end
