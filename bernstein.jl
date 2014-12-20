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