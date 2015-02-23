module Dir

using Distributions, Roots

include("asymptotics.jl")
include("chiSqVar.jl")
include("klVar.jl")
include("oracles.jl")

#the relevant integrand
#vbars = vs - te
fint(z, vbars, alphas) = prod((1 + z*vbars).^-alphas)/z
fint_deriv(z, vbars, alphas) = prod((1 + z*vbars).^-alphas)

#prob_fun(vbars, alphas) => Prob( X(vbars, alphas) >= 0 )
function VaR(vs, alphas, eps_, prob_fun = prob_direct; 
				tol=1e-6)
	#rescale vs to make life easier
	scale = mean(abs(vs))
	vs = vs/scale
	tmin = minimum(vs) + tol
	tmax = maximum(vs) - tol
	
	if tmax-tmin <= 2tol #degenerate case where all vs equal
		return scale
	else
		scale * fzero(t->prob_fun(vs-t, alphas)-1+eps_, tmin, tmax)
	end
end

#finds the optimal positive contour for fundamental integral
#assumes vbars has at least one pos component
function calc_shift(vbars, alphas_)
	ps = -1./vbars
	alphas = alphas_[ps .> 0]
	ps = ps[ps .> 0]
	ix = indmin(ps)
	ps[ix]/(1 + alphas[ix])
end

#finds optimal pos contour for the deriv integral
#assumes vbars has both pos and neg elements
function calc_shift_deriv(vbars, alphas)
	imin = indmin(vbars)	
	imax = indmax(vbars)
	pneg = -1/vbars[imax]
	ppos = -1/vbars[imin]
	@assert pneg*ppos < 0
	(pneg*alphas[imin] + ppos*alphas[imax])/(alphas[imin] + alphas[imax])
end

#Computes Prob(X(vbars, alphas) <= 0) = Prob(v p <= t)
#includes pi scaling
function prob_direct(vbars, alphas)
	#handle some degnerate cases gracefully
	if minimum(vbars) > 0
		return 0.
	elseif maximum(vbars) <= 0
		return 1.
	end
	const a = calc_shift(vbars, alphas)
	real(quadgk(s->fint(a + s * 1im, vbars, alphas)/pi, 0, Inf)[1])
end
prob_direct(vs, alphas, t) = prob_direct(vs-t, alphas)

#the gradient of P(vbars'p <= 0) with respect to vbars
function calc_prob_grad(vbars, alphas)
	d = length(vbars)
	@assert minimum(alphas) > 0 "alphas not positive $alphas"

	#degenerate case
	if (maximum(vbars) < 0) || (minimum(vbars) > 0)
		return zeros(Float64, d)
	end

	grad = zeros(Float64, d)
	alphas = copy(alphas) #safety over efficiency
	for i = 1:d
		alphas[i] +=1
		const a = calc_shift_deriv(vbars, alphas)
		grad[i] = real(quadgk(s->fint_deriv(a + s * 1im, 
											vbars, alphas)/pi, 
								0, Inf)[1])				
		alphas[i] -=1
		grad[i] *= -alphas[i]
	end
	grad
end
calc_prob_grad(vs, alphas, t) = calc_prob_grad(vs-t, alphas)

#gradient of VaR
function grad_VaR(vs, alphas, var)
	TOL = 1e-10
	#handle the degenerate case
	if maximum(abs(vs-var)) <= TOL
		#any element should work?
		return alphas/sum(alphas)
	end

	#first compute the unscaled partial derivs of P(v'p <=t)
	const d = length(vs)
	grad = calc_prob_grad(vs, alphas, var)
	@assert abs(sum(grad)) > TOL "gradient is zero for var $var \n v $vs \n alphas $alphas \t $grad"
	grad/sum(grad)
end

function kl_ratio(vs, alphas, eps_)
	phat, alpha0 = calc_phat(alphas)
    Gamma = log(1/eps_) / alpha0
    bval = KLVarSol(vs, phat, Gamma)[1] - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    bval/val
end

function kl_cov_ratio(vs, alphas, eps_, d, N)
	phat, alpha0 = calc_phat(alphas)
    Gamma = quantile(Chisq(d-1), 1-eps_)/2N
		    bval = KLVarSol(vs, phat, Gamma)[1] - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    bval/val
end

function chisq_ratio(vs, alphas, eps_)
	phat, alpha0 = calc_phat(alphas)
    sm_val = chiSqVar(vs, alphas, eps_) - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    sm_val/val
end

function chisq_cov_ratio(vs, alphas, eps_, d, N)
	phat, alpha0 = calc_phat(alphas)
	sm_val = chiSqVarCov(vs, alphas, eps_, N) - dot(phat, vs)
	val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
	sm_val / val
end

end #ends module

