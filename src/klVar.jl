#########
# The KL Set (Bernstein type inequality)
##########
using JuMP
function KLVarSol(vs, phat, Gamma; TOL = 1e-8, factor=10)
	#degeneracy
	if maximum(vs) - minimum(vs) < TOL		
		return maximum(vs), phat
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

#For the near optimal variant
function KLVar(vs, alphas, eps_)
	phat = alphas/sum(alphas)
	Gamma = log(1/eps_)/sum(alphas)  #approximates N = sum(alphas)
	KLVarSol(vs, phat, Gamma)[1]
end

#For the coverage variant
function KLVarCov(vs, alphas, eps_)
	phat = alphas/sum(alphas)
	Gamma = Distributions.Chisq(length(alphas)-1, 1-eps_)/2/sum(alphas)
	KLVarSol(vs, phat, Gamma)[1]
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