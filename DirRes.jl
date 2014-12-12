module Dir

using Distributions, Roots

include("bernstein.jl")
include("asymptotics.jl")

#the relevant integrand
#vbar = v - te
fint(z, vbars, alphas) = prod((1 + z*vbars).^-alphas)/z

#prob_fun(vs, alphas) => Prob( X(vs, alphas) > 0 )
#VG Revisit this to add a timing function
function VaR(vs, alphas, eps_, prob_fun; 
				tmax=bernVar(vs, alphas, eps_), 
				tmin=dot(vs, alphas)/sum(alphas), 
				tol=1e-6)
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

function probNeg_direct(vbars, alphas)
	#handle some degnerate cases gracefully
	if minimum(vbars) > 0
		return 0.
	elseif maximum(vbars) <= 0
		return 1.
	end
	const a = calc_shift(vbars, alphas)
	real(quadgk(s->fint(a + s * 1im, vbars, alphas)/pi, 0, Inf)[1])
end

#evaluates Dir tail Prob for 3dimensions only
function Prob3d(alphas, vs, t)
	@assert length(alphas) == 3 "Only defined for 3 dimensions"
	if t >= norm(vs, 1)
		return 0
	elseif t <= -norm(vs, 1)
		return 1
	end

	##VG a Cludge
	if abs(vs[2] - vs[3]) < 1e-10
		vs[2] = vs[2] + 1e-8
	end

	beta1 = Beta(alphas[1], alphas[2] + alphas[3])
	beta2 = Beta(alphas[2], alphas[3])
	#the integrand
	function f(p1)
		val = (t - vs[1]*p1)/(1-p1)/(vs[2]-vs[3]) - vs[3]/(vs[2]-vs[3])
	 	cdf(beta2, val) * pdf(beta1, p1)
	end

	I, E = Base.quadgk(f, 1e-10, 1-1e-10)
	if vs[2] > vs[3]
	 	I = 1 - I
	end
	 I
end

end #ends module

