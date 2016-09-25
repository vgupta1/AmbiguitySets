###
# Testing the shift hypothesis
###

include("DirRes.jl")
using JuMPeR

n = 5
alphas = fill(10, n)
w = Dir.OptDirichlet(alphas, .1)

println("W Details:\t", w.n)

srand(8675309)

for numRuns = 1:100
	vs = 1000 * randn(n)
	var = Dir.VaR(vs, w.alphas, w.eps_)
	grad = Dir.grad_VaR(vs, w.alphas, var)
	shift = var - dot(vs, grad)
	abs(shift) > 1e-8 && println("Shift:\t", shift)
	# println("Grad:\t", minimum(grad), "\t", sum(grad))
	# println()
end

# #the straightup VaR
# println("VaR:\t", Dir.VaR(vs, w.alphas, w.eps_))

# #the suppfcn
# println("Supp:\t", Dir.suppFcn(vs, w))

# m = RobustModel()
# setDefaultOracle!(m, w)
# @defVar(m, xs[1:n] >=0)
# @addConstraint(m, sum{xs[i], i=1:n} == 1)

# @defVar(m, t)
# @defUnc(m, ps[1:n])
# addConstraint(m, t>= sum([ps[i]*xs[i] for i =1:n]))

# @setObjective(m, Min, t)
# solveRobust(m)

# println("xs:\t", getValue(xs))
# println("obj value:\t", getObjectiveValue(m))

