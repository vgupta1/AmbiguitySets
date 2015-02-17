###
# Brief test suite to allow code refactoring
###
include("../DirRes.jl")

using FactCheck

facts("simpleVarCalcs") do
	n = 5; eps_ = .1
	alphas = fill(10, n)
	srand(8675309)
	vs = randn(n)

	var = Dir.VaR(vs, alphas, eps_)
	@fact var =>roughly(0.6911925526746056)

	dir = randn(n)
	grad = Dir.grad_VaR(vs, alphas, var)	
	@fact dot(dir, grad) => roughly(0.03525674543867546)
end

facts("cross_sec Plot") do
	n = 5; eps_ = .1
	alphas = fill(10, n)
	srand(8675309)
	vs = randn(n)

	var = Dir.VaR(vs, alphas, eps_)
	@fact var =>roughly(0.6911925526746056)

	klvar = Dir.KLVar(vs, alphas, eps_)
	@fact klvar =>roughly(0.7966523725611896)

	momvar = Dir.chiSqVar(vs, alphas, eps_)
	@fact momvar =>roughly(0.9007022382437936)
end


facts("asymptotics") do
	eps_ = .1
	d = 5
	@fact Dir.kl_const(eps_) => roughly(1.6745061876440455)
	@fact Dir.kl_cov_const(eps_, d) => roughly(2.1763968656568484)
	@fact Dir.chisq_const(eps_) => roughly(2.340912438217137)
	@fact Dir.chisq_cov_const(eps_, d) => roughly(2.1763968656568484)
end

facts("relative size") do
	n = 5; eps_ = .1
	alphas = fill(10, n)
	N = sum(alphas)
	srand(8675309)
	vs = randn(n)

	@fact Dir.kl_ratio(vs, alphas, eps_) => roughly(1.6715439839816046)
	@fact Dir.chisq_ratio(vs, alphas, eps_) => roughly(2.334109702455164)
	@fact Dir.kl_cov_ratio(vs, alphas, eps_, n, N) => roughly(2.1596923050794143)
	@fact Dir.chisq_cov_ratio(vs, alphas, eps_, n, N) => roughly(2.1916655054582344)
end

facts("support Fcns") do
	n = 5; eps_ = .1
	alphas = fill(10, n)
	srand(8675309)
	vs = randn(n)

	wstar = Dir.OptDirichlet(alphas, eps_)
	var, pstar = Dir.suppFcn(vs, wstar, :Min)
	@fact var => roughly(0.37703095047921426)
	@fact pstar[1] => roughly(0.24713880804852523)

	wKL = Dir.KLSet(alphas, eps_)
	klvar, klpstar = Dir.suppFcn(vs, wKL, :Min)
	@fact klvar =>roughly(0.27280988597450007)
	@fact klpstar[1] => roughly(0.2883035532766228)

	wChiSq = Dir.ChiSqSet(alphas, eps_) 
	chi_var, chi_pstar = Dir.suppFcn(vs, wChiSq, :Min)
	@fact chi_var =>roughly(0.34690230081245177)
	@fact chi_pstar[1] => roughly(0.2544700296340291)
end

FactCheck.exitstatus()