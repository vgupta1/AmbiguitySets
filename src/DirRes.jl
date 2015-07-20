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
	if maximum(vs)-tmin <= 2tol #degenerate case where all vs equal
		return scale
	end

	tmax = KLVar(vs, alphas, eps_) #a more expensive, but better starting pt
	scale * fzero(t->prob_fun(vs-t, alphas)-1+eps_, tmin, tmax)
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

	#check the asymptotic result
	a = dot(alphas, vbars)/ dot(alphas, vbars.^2)
	if pneg <= a <= ppos
		return a
	else
		println("old way")
		(pneg*alphas[imin] + ppos*alphas[imax])/(alphas[imin] + alphas[imax])
	end
end

#Computes Prob(X(vbars, alphas) <= 0) = Prob(v p <= t)
#includes pi scaling
function prob_direct(vbars, alphas)
	#handle some degnerate cases gracefully
	if minimum(vbars) >= 0
		return 0.
	elseif maximum(vbars) <= 0
		return 1.
	end
	const a = calc_shift(vbars, alphas)
	I, E = quadgk(s->fint(a + s * 1im, vbars, alphas)/pi, 0, Inf)
	if E > 1e-6
		println("Bad Integral:\t", real(I), "\t", E)
		println("vs:\t", vbars)
		println("alphas:\t", alphas)
		println()
	end
	real(I)
	# if abs(E/real(I)) > 1e-5
	# 	println("using Big Floats")
	# 	I, E = quadgk(s->fint(a + s * 1im, vbars, alphas)/pi, zero(BigFloat), BigFloat(Inf))
	# end
	# float(real(I))
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
		#try it once with ordinary floats
		I, E = quadgk(s->fint_deriv(a + s * 1im, vbars, alphas)/pi, 
								0., 3.)
		if E > 1e-6
			println("Bad Deriv Integral:\t", i, "\t", real(I), "\t", E, "\t", E/real(I))
			println("vs:\t", vbars)
			println("alphas:\t", alphas)
			println("hderiv(a):\t", fint_deriv(a, vbars, alphas)/pi)
			println()

		end

		#Use exact arithmetic to avoid the numerical stability
		#very slow...
		# if abs(E/real(I)) > 1e-5
		# 	println("using Big Floats")
		# 	I, E = quadgk(s->fint_deriv(a + s * 1im, vbars, alphas)/pi, 
		# 						zero(BigFloat), BigFloat(Inf))
		# end

		grad[i] = real(I)				
		alphas[i] -=1
		grad[i] *= -alphas[i]
	end
	min(0., grad)  #this is to correct for small overflow errors
end
calc_prob_grad(vs, alphas, t) = calc_prob_grad(vs-t, alphas)

#gradient of VaR
function grad_VaR(vs, alphas, var)
	TOL = 1e-10
	#rescale for nicety
	scale = mean(abs(vs))
	vs = vs/scale
	var = var/scale 

	#handle the degenerate case
	if maximum(abs(vs-var)) <= TOL
		#any element should work
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

