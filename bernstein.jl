#########
# Using the benrstein type stuff
##########
using JuMP
function supp_rep(vs, phat, Gamma)
	m = Model()
	@defVar(m, theta)
	@defVar(m, lam >=0)
	@addConstraint(m, theta >= maximum(vs))

	const d = length(vs)
	#define the inner bits
	@defVar(m, inner[1:d] >= 0)	
	for i = 1:d
		@addNLConstraint(m, theta - vs[i] == lam * inner[i])
	end

	@defVar(m, term)
	@addNLConstraint(m, term <= lam * sum{phat[i] * log(inner[i]), i=1:d})

	@setNLObjective(m, Min, (Gamma-1)*lam + theta - term)
	solve(m)
	getObjectiveValue(m), getValue(theta), getValue(lam)
end

function bernVar(vs, phat, Gamma; TOL = 1e-8, factor=10)
	if maximum(vs) - minimum(vs) < TOL
		println("Degenerate")
		return maximum(vs)
	end

	#fcn arises from KKT conditions
	function f(theta)
		lam = 1/sum(phat ./ (theta-vs))
		p = lam*phat ./ (theta-vs)
		dot(phat, log(phat./p)) - Gamma, p
	end
	low, high = mult_bracket(theta->f(theta)[1] ,maximum(vs) + TOL)
	theta_opt = fzero(theta->f(theta)[1], low, high)
	p_opt = f(theta_opt)[2]
	dot(vs, p_opt), p_opt
end

#computes a valid bracket by multiplicative search
function mult_bracket(fun, val_; factor=2, max_iter=50)
	val    = val_
	f_init = fun(val)
	iter   = 0
	while (f_init*fun(val) > 0) && (iter < max_iter)
		val = val * factor
		iter = iter + 1
	end
	iter == max_iter && error("Bracketing Failed")
	(val/factor, val)
end

# mgf_p(vs, phat, lam, Gamma) = lam*Gamma - lam*dot(phat, log(1 - vs./lam))
# deriv(vs, phat, lam, Gamma) = Gamma + dot(phat, vs./(vs-lam)) - dot(phat, log(1 - vs/lam))
# lammin(vbars; TOL=1e-8) = max(maximum(vbars), 0) + TOL

# function mgf_p(vs, phat, Gamma)
# 	lamstar = fzero(l->deriv(vs, phat, l, Gamma), lammin(vs),  1e2)
# 	lammin(vs), lamstar, mgf_p(vs, phat, lamstar, Gamma)
# end



