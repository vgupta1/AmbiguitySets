###
# Test the oracles
###

include("DirRes.jl")
using JuMPeR

n = 5
alphas = fill(10, n)
eps_ = .1
w = Dir.KLSet(alphas, eps_)

srand(8675309)
vs = randn(n)

#thesuppFcn for min
var, grad = Dir.suppFcn(vs, w, :Min)
println("Supp:\t", var, "\t", grad)

dir = randn(n) * .01
var2 = -Dir.VaR(-vs, alphas, eps_)
var_ = -Dir.VaR(-vs - dir, alphas, eps_)
println("Var Comparison:\t", var, "\t", var2)
println("Diff1:\t", (var_ - var2))
println("Exact:\t", dot(dir, grad))




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

