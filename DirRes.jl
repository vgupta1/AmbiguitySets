module Dir

using Distributions, Roots

include("asymptotics.jl")
include("secMom.jl")
include("bernstein.jl")

#the relevant integrand
#vbar = v - te
fint(z, vbars, alphas) = prod((1 + z*vbars).^-alphas)/z

#prob_fun(vbars, alphas) => Prob( X(vbars, alphas) > 0 )
#VG Revisit this to add a timing function
function VaR(vs, alphas, eps_, prob_fun; 
				tol=1e-6, 
				tmax=maximum(vs) - tol, 
				tmin=minimum(vs) + tol, 
)
	fzero(t->prob_fun(vs -t, alphas)-1 +eps_, tmin, tmax)
end

#finds the optimal positive contour
function calc_shift(vbars_, alphas_)
	ps = -1./vbars_
	alphas =alphas_[ps .> 0]
	ps = ps[ps .> 0]
	ix = indmin(ps)
	ps[ix]/(1 + alphas[ix])
end

function prob_gauss_approx(vbars, alphas)
	const mu = dot(vbars, alphas)
	const sig = sqrt(dot(alphas, vbars.^2))
	cdf(Normal(), -mu/sig)
end

#Computes Prob(X(vbars, alphas) <= 0)
#Equivalent ot computing the fundamental integral
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

function kl_ratio(vs, alphas, eps_, N)
	phat, alpha0 = calc_phat(alphas)  #Remember, N and alpha0 may be different
    Gamma = log(1/eps_) / N
    bval = bernVar(vs, phat, Gamma)[1] - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    bval/val
end

function kl_ratio_bt(vs, alphas, eps_, d, N)
	phat, alpha0 = calc_phat(alphas)  #Remember, N and alpha0 may be different
    Gamma = quantile(Chisq(d-1), 1-eps_)/2N
    bval = bernVar(vs, phat, Gamma)[1] - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    bval/val
end

function mom_ratio(vs, alphas, eps_)
	phat, alpha0 = calc_phat(alphas)
    sm_val = secondMomentVar(vs, alphas, eps_) - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    sm_val/val
end

function mom_ratio_bt(vs, alphas, eps_, d, N)
	phat, alpha0 = calc_phat(alphas)
	sig = sqrt(vs'*sigma(alphas)*vs * alpha0/N )[1] #corrects for def
	sm_val = sqrt(quantile(Chisq(d-1), 1-eps_)) * sig 
	val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
	sm_val / val
end

end #ends module

