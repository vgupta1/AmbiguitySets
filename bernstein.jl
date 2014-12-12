#########
# Using the benrstein type stuff
##########
using JuMP

mgf_p(vbars, alphas, lam, eps_) = lam * log(1/eps_) - lam * dot(alphas,  log(1-vbars ./ lam))
deriv(vbars, alphas, lam, eps_) = log(1/eps_) + dot(alphas, vbars./(vbars - lam)) - dot(alphas, log(1- vbars / lam))
lammin(vbars; TOL=1e-8) = max(maximum(vbars), 0) + TOL

function mgf_p(vbars, alphas, eps_)
	lamstar = fzero(l->deriv(vbars, alphas, l, eps_), lammin(vbars),  1e2)
	lammin(vbars), lamstar, mgf_p(vbars, alphas, lamstar, eps_)
end

function bernVar(vs, alphas, eps_)
	fzero(t-> mgf_p(vs-t, alphas, eps_)[3], -norm(vs, Inf), norm(vs, Inf))
end

#chernoff  bound on P(vbars^T p >=0)
function chernoff(vbars, alphas)
	#minimize the log for kicks...
	f(lam) = dot(-alphas, log(1 - lam * vbars))
	if maximum(vbars) < 0
		l_ = 1e2
	else
		l_ = minimum([1/vi for vi in vbars[vbars .> 0]])
	end
	res = optimize(f, 1e-10, l_)
	return exp(res.f_minimum)
end

#computes the finite sample sigma
function sig_finite(alphas)
	const d = length(alphas)
	sigma = zeros(d, d)
	const alpha0 = sum(alphas)
	const denom = (alpha0 + 1)alpha0^2
	for i = 1:d
		sigma[i, i] = (alpha0 - alphas[i]) * alphas[i] / denom
		for j=i+1:d
			sigma[i, j] = -alphas[i] * alphas[j]/denom
			sigma[j, i] = sigma[i, j] 
		end
	end
	sigma
end

#Cheap 2 moment VaR bound
function secondMomentVar(vs, alphas, eps_; gauss_approx=false)
	const kappa = gauss_approx ? quantile(Normal(), 1-eps_) : sqrt(1/eps_ - 1)
	const mu = dot(vs, alphas)/sum(alphas)
	const sig = sqrt(vs' * sig_finite(alphas) * vs)[1]
	vmin = minimum(vs)
	vmax = maximum(vs)
	min(max(mu + kappa * sig, vmin), vmax)
end
