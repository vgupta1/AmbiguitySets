module Dir

using Distributions, Roots, Optim, JuMP

#Approximate the residue of fun at z0
function numRes(fun, z0; r = .01)
	f(t) = r * exp(2im * pi * t) * fun(z0 + r * exp(2im * pi * t) )
	quadgk(f, 0., 1.)
end

#fun takes a complex argument
function numRes(fun, z0, r, m; old=false)
	if old
		return numRes(fun, z0, r=r)
	end
	f(t) = r * exp(2im * pi * t) * fun(z0 + r * exp(2im * pi * t) )
	sum(map(f, linspace(0, 1, m+1)[1:end-1]))/m
end

#inserts zero proxies and retains ordering
function sort_poles(vbars_, alphas_)
	ps    = [-1./vbars_, 0.0]
	alphas    = [alphas_, 1]
	indxs = sortperm(ps)
	[vbars_, 1][indxs], ps[indxs], alphas[indxs]
end


#approximates a good raduii for numRes and conditioning constants
#assumes that vbars, ps, alphas are ordered now and proxied
function optRs(vbars, ps, alphas)
	rs    = zeros(length(ps))  			#radius for integration
	lams    = zeros(length(ps))			#conditioning of the integrand

	zstar_l = Inf
	for ix = 1:length(ps)
		#compute the right "midpoint"
		if ix == length(ps)
			zstar_r = Inf
		else
			zstar_r = ps[ix+1] - ps[ix]
		end

		#maximum size of the annulus
		if zstar_l < zstar_r
			rmax, alpha_max = zstar_l, alphas[ix-1]
			vmax = vbars[ix-1]
		else
			rmax, alpha_max = zstar_r, alphas[ix + 1]
			vmax = vbars[ix + 1]
		end

		#locally optimum r
		rs[ix] = alphas[ix] * rmax / (alphas[ix] + alpha_max)
		lams[ix] = abs((vbars[ix] * rs[ix])^(2*alphas[ix] - 1))

		zstar_l = zstar_r
	end
	rs, lams
end

#the relevant integrand
#vbar = v - te
fint(z, vbars, alphas) = 1/ z / prod((1 + z * vbars) .^ alphas)

#defunct....
function eq6(vbars, alphas)
	warn("Deprecated")
	f(z) = fint(z, vbars, alphas)
	I1, E1 = quadgk(f, -10im, -.01im)
	I2, E2 = quadgk(f, .01im, 10im)
	.5 + 1/2pi * (I1 + I2) / 1im
end

function probXNeg(vbars_, alphas_; method = :PosContour, m=100)
	vbars, ps, alphas = sort_poles(vbars_, alphas_)
	rs, lams          = optRs(vbars, alphas)
	f(z)              = fint(z, vbars, alphas)
	ps_rs_ls = zip(ps, rs, lams)

	if method == :PosContour || method == :NegContour
		if method == :PosContour
			ps_rs_ls = filter(p_r_l -> p_r_l[1] > 0, ps_rs_ls)
		elseif method ==:NegContour
			ps_rs_ls = filter(p_r_l -> p_r_l[1] < 0, ps_rs_ls)
		else
			error()
		end
		residues = zeros(length(ps_rs_ls))
		for (ix, (p, r, l)) in enumerate(ps_rs_ls)
			residues[ix] = real(numRes(z->l * f(z), p, r, m, old=false)[1])
		end

		inv_ls = [1/l for (p, r, l) in ps_rs_ls]

		if method == :PosContour
			return -sum_kbn(inv_ls .* residues)
		elseif method ==:NegContour
			return 1 + sum_kbn(inv_ls .*residues)
		else
			error()
		end
	elseif method == :Shift
		if maximum(ps) < 0 
			#VG Cludge
			kappa = 1
		else
			kappa = .5 * minimum(ps[ ps .> 0])
		end
		println(kappa)
		I,E = quadgk(f, kappa - 10im, kappa + 10im)  #VG more principled solution needed
		println("Error \t $E")
		return I / (2im * pi)
	else
		error("Method not recognzied")
	end
end

function Var(alphas, vs, eps_; method=:PosContour)
	exp_val = dot(vs, alphas ./sum(alphas))
	tmax = bernVar(vs, alphas, eps_)
	#VG Wrap the function to do some tracing
	ix = 0
	function f(t)
		const val = @time probXNeg(vs - t, alphas, method=method)
		println("$ix \t $t \t $val")
		ix = ix + 1
		real(val) -1. + eps_
	end
	fzero(f, exp_val, tmax)
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
	 I, E
end

function Var3d(alphas, vs, eps_)
	fzero(t->Prob3d(alphas, vs, t)[1] -eps_, -norm(vs, Inf), norm(vs, Inf))
end

function convexTest(v1, v2; eps_=.1, alphas=[2 3 10])
	for lam = linspace(0, 1, 20)
		v = lam * v1 + (1-lam)*v2
		println("$lam \t $(Var3d(alphas, v, eps_))")
	end
end

testProbs(vbars, alphas) = 
	(Prob3d(alphas, vbars, 0)[1], real(1-probXNeg(vbars, alphas)[1]), 
								real(1-probXNeg(vbars, alphas, method=:NegContour)[1]), 
								1-probXNeg(vbars, alphas, method=:Shift)[1] )

#Seems to be a numerical stability problem for some values of vbars alphas
# vbars -0.925296  -0.264872  -0.905954
# alphas 3.0  4.0  5.0

#also the negative computation seems to be off by 1 (too large) sometimes, and sometimes weird
#This example, it is 1 too small.
# vbars = -0.669911 -0.59414 0.930822
# alphas = 4 5 7


#why not integrate the function directly?  The singularity only seems to affect
#the imaginary part... And the function (seems) odd.... Maybe fold it and evaluate 
# as a real integral?

# You are Free to scale vbar.... does this help?


#########
# Using the benrstein type stuff
##########
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

#Computes the Sigma Star limit for Dirichlet
function sigStar(pstar)
	const d = length(pstar)
	sigma = zeros(d, d)
	for i = 1:d
		sigma[i, i] = pstar[i]*(1-pstar[i])
		for j = i+1:d
			sigma[i, j] = -pstar[i] * pstar[j]
			sigma[j, i] = sigma[i, j]
		end
	end
	sigma
end

end #ends module

