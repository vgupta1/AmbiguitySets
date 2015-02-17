###
# Brief test suite
###
#Currently values not fixed... just checking for failures.  
include(../"DirRes.jl")

using FactCheck

facts("simpleVarCalcs") do
n = 5; eps_ = .1
alphas = fill(10, n)
srand(8675309)
vs = randn(n)

@fact error("checking")
#Dir.VaR(vs, alphas, eps_)
end

# function VaR(vs, alphas, eps_, prob_fun = prob_direct; 
# 				tol=1e-6, 
# 				tmax=maximum(vs) - tol, 
# 				tmin=minimum(vs) + tol, 
# )
# 	fzero(t->prob_fun(vs -t, alphas)-1 +eps_, tmin, tmax)
# end

FactCheck.exitstatus()