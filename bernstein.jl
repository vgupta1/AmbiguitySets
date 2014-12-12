#########
# Using the benrstein type stuff
##########
using JuMP

mgf_p(vbars, alphas, lam, eps_) = lam * log(1/eps_) - lam * dot(alphas,  log(1-vbars ./ lam))
deriv(vbars, alphas, lam, eps_) = log(1/eps_) + dot(alphas, vbars./(vbars - lam)) - dot(alphas, log(1- vbars / lam))
lammin(vbars; TOL=1e-8) = max(maximum(vbars), 0) + 1e-8

function mgf_p(vbars, alphas, eps_)
	lamstar = fzero(l->deriv(vbars, alphas, l, eps_), lammin(vbars),  1e2)
	lammin(vbars), lamstar, mgf_p(vbars, alphas, lamstar, eps_)
end

function bernVar(vs, alphas, eps_)
	fzero(t-> mgf_p(vs-t, alphas, eps_)[3], -norm(vs, Inf), norm(vs, Inf))
end

#computes the best chernoff style bound on P(vbars^T p >=0)
function chernoff(vbars, alphas)
	#minimize the log for kicks...
	f(lam) = dot(-alphas, log(1 - lam * vbars) )
	if maximum(vbars) < 0
		l_ = 1e2
	else
		l_ = minimum([ 1/ abs(vi) for vi in vbars[vbars .> 0]])
	end
	res = optimize(f, 1e-10, l_)
	return exp(res.f_minimum)
end

function mgf_d(vbars, alphas, lam, eps_)
	m = Model()
	@defVar(m, ys[1:d] >= 0)
	@defVar(m, s)

	@addNLConstraint(m, 
		sum{ (ys[i] - alphas[i])+ alphas[i]*log(alphas[i]/ys[i]), i=1:d}
			 <= log(1/eps_) + s )
	@setObjective(m, Max, sum{vbars[i]*ys[i], i=1:d} - lam*s)

	solve(m)
	getObjectiveValue(m), getValue(ys), getValue(s)
end

# The super cheap 2 moment VaR bound
function secondMomentVar(vs, alphas, eps_)
	const kappa = sqrt(1/eps_ - 1)
	const phat = alphas/sum(alphas)
	const d = length(alphas);  @assert d == length(vs)
	const N = sum(alphas)
	sigma = sigStar(phat) / N
	dot(vs, phat) + kappa * sqrt( vs' * sigma * vs )[1]
end
