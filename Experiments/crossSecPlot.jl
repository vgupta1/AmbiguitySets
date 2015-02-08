## Generates the cross-section plots

using JuMP, Roots
include("../DirRes.jl")

function crossSec_dir(dir, d, eps_, alphas)
	phat = alphas/sum(alphas)

	#set up the basic model
	m = Model()
	@defVar(m, p[1:d] >=0)
	@addConstraint(m, sum{p[i], i=1:d}==1)
	for j = 4:d
		@addConstraint(m, p[j] == phat[j])
	end
	@setObjective(m, Max, sum{dir[i]*p[i], i=1:d})
	lazy_P!(m, p, alphas, eps_, phat)
end

function crossSec_plotPstar(path, N)
	d = 5; eps_ = .1
	alphas = fill(ifloor(N/d), d)

	fp = open(path, "w")
	writecsv(fp, ["NumPts" "Type" "x" "y" "p1" "p2"])
	for theta = linspace(0, 2pi, 100)
		dir = [cos(theta), sin(theta), zeros(d-2)]
		pstar = crossSec_dir(dir, d, eps_, alphas)[:]
		writecsv(fp, [N "pstar" cos(theta) sin(theta) pstar[1:2]'])
	end
	close(fp)
end


#repeatedly solve given model adding constrains corresponding to p \in P
#phat must be a point strictly interior to P
function lazy_P!(m, p, alphas, eps_, phat)
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
		var = Dir.VaR(vs, alphas, eps_)
		prod = dot(vs, pstar)
		iter +=1
	end
	iter == MAXITER && error("Max iterations reached")
	return getValue(p)
end