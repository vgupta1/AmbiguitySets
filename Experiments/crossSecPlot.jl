## Generates the cross-section plots

using JuMP, Roots
include("../DirRes.jl")

####
#Calls used to generate the plots
#crossSec_plot("Results/pstar50.csv", 50, Dir.VaR, "pstar")
#crossSec_plot("Results/secMom50.csv", 50, Dir.chiSqVar, "ChiSq")
#crossSec_plot("Results/bernVar50.csv", 50, Dir.KLVar, "ChiSq")

function crossSec_dir(dir, d, eps_, alphas, varfun)
	phat = alphas/sum(alphas)

	#set up the basic model
	m = Model()
	@defVar(m, p[1:d] >=0)
	@addConstraint(m, sum{p[i], i=1:d}==1)
	for j = 4:d
		@addConstraint(m, p[j] == phat[j])
	end
	@setObjective(m, Max, sum{dir[i]*p[i], i=1:d})
	lazy_P!(m, p, alphas, eps_, phat, varfun)
end

#varfun should take vs, alphas, eps_ and return var
function crossSec_plot(path, N, varfun, tag, useSkew=false)
	d = 5; eps_ = .1
	if useSkew
		alphas = [d, 1, 2, 2, 2]
		alphas = ifloor(N/sum(alphas)) * alphas
	else
		alphas = fill(ifloor(N/d), d)
	end

	fp = open(path, "w")
	writecsv(fp, ["NumPts" "Type" "x" "y" "p1" "p2"])
	for theta = linspace(0, 2pi, 100)
		dir = [cos(theta), sin(theta), zeros(d-2)]
		pstar = crossSec_dir(dir, d, eps_, alphas, varfun)[:]
		writecsv(fp, [N tag cos(theta) sin(theta) pstar[1:2]'])
	end
	close(fp)
end

#repeatedly solve given model adding constraints corresponding to p \in P
#phat must be a point strictly interior to P
function lazy_P!(m, p, alphas, eps_, phat, varfun)
	TOL = 1e-6; MAXITER = 100; iter = 1

	const d = length(phat)
	prod = 2.
	var = 1.
	vs = ones(Float64, d)
	while (prod > var + TOL) && (iter < MAXITER)
		@addConstraint(m, sum{vs[i]*p[i], i=1:d} <= var)

		solve(m)
		pstar = getValue(p)
		vs = pstar - phat
		var = varfun(vs, alphas, eps_)
		prod = dot(vs, pstar)
		iter +=1
	end
	iter == MAXITER && error("Max iterations reached")
	return getValue(p)
end