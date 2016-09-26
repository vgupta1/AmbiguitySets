###
# Portfolio Experiment for the paper based on kenneth french data
###
include("../src/DirRes.jl")
module port
using Distributions, Gurobi, JuMP, JuMPeR, Dir

#creates supp of distribution most recent numPts of KF dataset
function getMktSupp(numPts, path="Data//12_Industry_Portfolios_clean.csv")
	#load up the data set
	dat, header = readdlm(path, ',', '\r', header=true)
	dat[end-numPts+1:end, 2:end]
end

#Compute the CVaR of each asset in the dataset based on the previous numPts months
function mktCVars(numPts, eps_=.1, path="Data//12_Industry_Portfolios_clean.csv")
	dat = getMktSupp(numPts, path)
	d = size(dat, 1)
	phat = ones(d)/d
	cvars = zeros(size(dat, 2))
	for i = 1:size(dat, 2)
		cvars[i] = cvar_sort(dat[:, i], phat, eps_)
	end
	cvars
end

#create N synthetic observations, drawn equally likely from support
function buildSynthMkt(N, supp)
	#sample counts
	#now sample counts from this distribution
	const d = size(supp, 1)
	dist = Categorical(ones(d)/d)
	sample = rand(dist, N)
	cnts = zeros(d)
	for j = 1:N
		cnts[sample[j]] += 1
	end
	cnts
end

#assess out-of-sample performance
#return the mean and cvar at eps_
#VG uses assumption that pstar = 1/N
function out_perf(x, eps_, supp)
	d = size(supp, 1)
	pstar = ones(d)/d
	rets = supp * x
	(pstar' * supp * x)[1], Dir.cvar_sort(rets, pstar, eps_)
end

# # VG drop this
# # this is defunct...  a truly synthetic version
# function buildMkt2(d, N; numAssets=5)
# 	#sample supp from an interesting continuous distribution
# 	supp = randn(d, numAssets)
# 	supp = supp * diagm(linspace(.1, .3, numAssets))
# 	supp = broadcast(+, linspace(0, .1, numAssets)', supp)

# 	#now sample counts from this distribution	
# 	dist = Categorical(ones(d)/d)	
# 	sample = rand(dist, N)
# 	cnts = zeros(d)
# 	for j = 1:N
# 		cnts[sample[j]] += 1
# 	end
# 	supp, cnts
# end

###
# CVAR-constrained portfolio allocation problems for various methods
###

#1/N Diversification
#Rescaled so that CVaR_eps is close to budget
function naive_port(eps_, budget, supp, counts)
	d, numAssets = size(supp)
	phat, alpha0 = Dir.calc_phat(counts)
	xs = ones(numAssets)
	cvar = Dir.cvar_sort(supp * xs, phat, eps_)
	xs /= cvar/budget  #ensures new cvar close to budget
	phat' * supp * xs, xs
end

#The minimum variance portfolio
function min_var_port(eps_, budget, supp, counts)
	m = Model(solver=GurobiSolver(OutputFlag=false))
	d, numAssets = size(supp)
	phat, alpha0 = Dir.calc_phat(counts)
	sigma = Dir.unscaled_sigma_(phat)

	#M^T Sigma Mkt is the relevant quadratic form
	Q = supp' * sigma * supp

	@defVar(m, xs[1:numAssets] >= 0)
	@addConstraint(m, sum(xs) == 1)
	@setObjective(m, Min, sum{xs[ix]*xs[jx]*Q[ix,jx], ix=1:numAssets, jx=1:numAssets})
	solve(m)
	xs = getValue(xs[:])

	#rescale to get CVaR close to budget
	cvar = Dir.cvar_sort(supp * xs, phat, eps_)
	xs /= cvar/budget  #ensures new cvar close to budget
	phat' * supp * xs, xs
end

function saa_port(eps_, budget, supp, counts)
	# solve for the SAA portfolio
	m = Model(solver=GurobiSolver(OutputFlag=false))
	d, numAssets = size(supp)
	N = sum(counts); phat = counts / N

	@defVar(m, xs[1:numAssets]>=0)
	@addConstraint(m, sum(xs) <= 1)
	@setObjective(m, Max, sum{phat[j] * dot(vec(supp[j, :]), xs), j=1:d})

	@defVar(m, Beta)
	@defVar(m, zs[1:d] >= 0)
	#define CVAR
	for j=1:d
		@addConstraint(m, -dot(vec(supp[j, :]), xs[:])  - Beta <= zs[j])
	end
	@addConstraint(m, sum{phat[j]*zs[j]/eps_, j=1:d} + Beta <= budget)
	solve(m)
	getObjectiveValue(m), getValue(xs[:])
end

function chisq_port(eps_, budget, supp, cnts, useCover)
	m = RobustModel(solver=GurobiSolver(OutputFlag=false))
	oracle = Dir.ChiSqSet(cnts + 1, eps_, useCover; debug_printcut=false) 
	setDefaultOracle!(m, oracle)

	d, numAssets = size(supp)
	@defVar(m, xs[1:numAssets]>=0)
	@addConstraint(m, sum(xs) <= 1)

	@defUnc(m, ps[1:d])

	@defVar(m, zs[1:d] >= 0) #
	@defVar(m, Beta)
	@defVar(m, obj)
	@addConstraint(m, sum([ps[j]*dot(vec(supp[j, :]), xs) for j=1:d]) >= obj)
	@setObjective(m, Max, obj)

	for j=1:d
		@addConstraint(m, -dot(vec(supp[j, :]), xs[:])  - Beta <= zs[j])
	end
	@addConstraint(m, sum([ps[j]*zs[j]/eps_ for j=1:d]) + Beta <= budget)	
	solve(m, prefer_cuts=true)
	getObjectiveValue(m), getValue(xs[:])
end

function kl_port(eps_, budget, supp, cnts, useCover)
	m = RobustModel(solver=GurobiSolver(OutputFlag=false))
	oracle = Dir.KLSet(cnts + 1, eps_, useCover; debug_printcut=false) 
	setDefaultOracle!(m, oracle)

	d, numAssets = size(supp)
	@defVar(m, xs[1:numAssets]>=0)
	@addConstraint(m, sum(xs) <= 1)
	@defUnc(m, ps[1:d])

	@defVar(m, zs[1:d] >= 0)
	@defVar(m, Beta)
	@defVar(m, obj)
	@addConstraint(m, sum([ps[j]*dot(vec(supp[j, :]), xs) for j=1:d]) >= obj)
	@setObjective(m, Max, obj)

	for j=1:d
		@addConstraint(m, -dot(vec(supp[j, :]), xs[:])  - Beta <= zs[j])
	end
	@addConstraint(m, sum([ps[j]*zs[j]/eps_ for j=1:d]) + Beta <= budget)	
	solve(m, prefer_cuts=true)
	getObjectiveValue(m), getValue(xs[:])
end

####
# The actual tests
####

#Simple test for debugging.
function test(d, N; budget=.2, seed=nothing)
	seed != nothing && srand(seed)
	supp = getMktSupp(d)
	cnts = buildSynthMkt(N, supp)

	zsaa, xsaa = saa_port(.1, budget, supp, cnts)
	zchisq, xchisq = chisq_port(.1, budget, supp, cnts, false)
	zkl, xkl = kl_port(.1, budget, supp, cnts, false)
	zchisq_cov, xchisq_cov = chisq_port(.1, budget, supp, cnts, true)
	zkl_cov, xkl_cov = kl_port(.1, budget, supp, cnts, true)

	println("Saa:\t $zsaa")
	show(xsaa'); println()
	println("Out:\t", out_perf(xsaa, .1, supp))

	println()
	println("ChiSq:\t $zchisq")
	show(xchisq'); println()

	println("Out:\t", out_perf(xchisq, .1, supp))

	println()
	println("KL:\t $zkl")
	show(xkl'); println()
	println("Out:\t", out_perf(xkl, .1, supp))

	println()
	println("Chisq_cov:\t $zchisq_cov")
	show(xchisq_cov'); println()
	println("Out:\t", out_perf(xchisq_cov, .1, supp))

	println()
	println("KL_cov:\t $zkl_cov")
	show(xkl_cov'); println()
	println("Out:\t", out_perf(xkl_cov, .1, supp))

end

## dump some stuff to files to look at plots
function test2(d, numRuns=100; budget=3, seed=nothing, eps_=.1)
#	f = open("Results/portExp2a_$(d)_$(budget).csv", "w")
	f = open("Results/portExp2b_$(d)_$(budget).csv", "w")
	writecsv(f, ["Run" "N" "Method" "inReturn" "outReturn" "CVaR" "X_Norm"])
	seed != nothing && srand(seed)
	supp = getMktSupp(d)

	for N in 100:100:1000
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp)

			zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
			outsaa, cvarsaa = out_perf(xsaa, eps_, supp)
			writecsv(f, [iSim N "SAA" zsaa outsaa cvarsaa norm(xsaa)])

			zchisq, xchisq = chisq_port(eps_, budget, supp, cnts, false)
			outchisq, cvarchisq = out_perf(xchisq, eps_, supp)
			writecsv(f, [iSim N "Chisq" zchisq outchisq cvarchisq norm(xchisq)])

			zchisq_cov, xchisq_cov = chisq_port(eps_, budget, supp, cnts, true)
			out_chi_cov, cvar_chi_cov = out_perf(xchisq_cov, eps_, supp)
			writecsv(f, [iSim N "ChisqCov" zchisq_cov out_chi_cov cvar_chi_cov norm(xchisq_cov)])

			zkl, xkl = kl_port(eps_, budget, supp, cnts, false)
			outkl, cvarkl = out_perf(xkl, eps_, supp)
			writecsv(f, [iSim N "KL" zkl outkl cvarkl norm(xkl)])

			zkl_cov, xkl_cov = kl_port(eps_, budget, supp, cnts, true)
			outkl_cov, cvarkl_cov = out_perf(xkl_cov, eps_, supp)
			writecsv(f, [iSim N "KLCov" zkl_cov outkl_cov cvarkl_cov norm(xkl_cov)])

			z_naive, x_naive = naive_port(eps_, budget, supp, cnts)
			out_naive, cvar_naive = out_perf(x_naive, eps_, supp)
			writecsv(f, [iSim N "Naive" z_naive out_naive cvar_naive norm(x_naive)])

			z_minvar, x_minvar = min_var_port(eps_, budget, supp, cnts)
			out_minvar, cvar_minvar = out_perf(x_minvar, eps_, supp)
			writecsv(f, [iSim N "MinVar" z_minvar out_minvar cvar_minvar norm(x_minvar)])
		end
	end
	close(f)
end

function test_incr_d(N, numRuns; budget=3, eps_=.1, seed=nothing)
	f = open("Results/incrDExp_a_$(N)_$(budget).csv", "w")
	writecsv(f, ["Run" "d" "Method" "inReturn" "outReturn" "CVaR"])
	seed != nothing && srand(seed)

	for numPts = 24:12:120
		supp = getMktSupp(numPts)
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp)
			zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
			zchisq_cov, xchisq_cov = chisq_port(eps_, budget, supp, cnts, true)
			zchisq, xchisq = chisq_port(eps_, budget, supp, cnts, false)
			zkl_cov, xkl_cov = kl_port(eps_, budget, supp, cnts, true)
			zkl, xkl = kl_port(eps_, budget, supp, cnts, false)

			outsaa, cvarsaa = out_perf(xsaa, eps_, supp)
			writecsv(f, [iSim numPts "SAA" zsaa outsaa cvarsaa])

			outchisq, cvarchisq = out_perf(xchisq, eps_, supp)
			writecsv(f, [iSim numPts "Chisq" zchisq outchisq cvarchisq])

			out_chi_cov, cvar_chi_cov = out_perf(xchisq_cov, eps_, supp)
			writecsv(f, [iSim numPts "ChisqCov" zchisq_cov out_chi_cov cvar_chi_cov])

			outkl, cvarkl = out_perf(xkl, eps_, supp)
			writecsv(f, [iSim numPts "KL" zkl outkl cvarkl])

			outkl_cov, cvarkl_cov = out_perf(xkl_cov, eps_, supp)
			writecsv(f, [iSim numPts "KLCov" zkl_cov outkl_cov cvarkl_cov])
		end
	end
	close(f)
end

#Generate data and a wrong prior and then assess
function test_wrongPrior(d, N, numRuns; budget=3, seed=8675309)
	f = open("Results/randWrongPrior_$(N)_$(d)_$(seed).csv", "w")
	writecsv(f, ["Run" "Dist" "Scale" "Method" "inReturn" "outReturn" "CVaR"])
	seed != nothing && srand(seed)
	eps_ = .1
	scale_grid = collect(.1:.1:3.)

	supp = getMktSupp(d)
	cnts = buildSynthMkt(N, supp)
	pstar = ones(d)/d

	for scale in scale_grid
		for iSim = 1:numRuns
			pseudo_cnts = rand(d)
			pseudo_cnts /= sum(pseudo_cnts)
			#calculate the distance
			dist = norm(pseudo_cnts - pstar)

			#scale up
			pseudo_cnts *= N * scale

			zchisq, xchisq = chisq_port(eps_, budget, supp, cnts+ pseudo_cnts, false)
			zkl, xkl = kl_port(eps_, budget, supp, cnts + pseudo_cnts, false)

			outchisq, cvarchisq = out_perf(xchisq, eps_, supp)
			writecsv(f, [iSim dist scale "Chisq" zchisq outchisq cvarchisq])

			outkl, cvarkl = out_perf(xkl, eps_, supp)
			writecsv(f, [iSim dist scale "KL" zkl outkl cvarkl])
		end
	end
	close(f)
end


#Generate data and a wrong prior and then assess
function test_wrongPriorScale(numRuns; budget=3, d=72, N=300, seed=nothing)
	f = open("Results/randWrongPriorScale_$(N)_$(d).csv", "w")
	writecsv(f, ["Run" "Scale" "Method" "inReturn" "outReturn" "CVaR"])
	seed != nothing && srand(seed)
	eps_ = .1

	supp = getMktSupp(d)
	pstar = ones(d)/d

	for tau_0 = 1:25:301
		pseudo_cnts = ones(d)
		pseudo_cnts[1] = tau_0
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp)
			zkl, xkl = kl_port(eps_, budget, supp, cnts + pseudo_cnts, false)
			outkl, cvarkl = out_perf(xkl, eps_, supp)
			writecsv(f, [iSim tau_0 "KL" zkl outkl cvarkl])
		end
	end
	close(f)
end

function test_wrongPriorConv(d, numRuns; budget=3, seed=nothing, scale=20)
	f = open("Results/randWrongPrior2_$(scale).csv", "w")
	writecsv(f, ["Run" "N" "Method" "inReturn" "outReturn" "CVaR"])
	seed != nothing && srand(seed)
	eps_ = .1

	supp = getMktSupp(d)

	prior1 = collect(1:d)/d 
	prior1 /= sum(prior1)
	prior1 *= scale

	prior2 = exp((1:d)/d)
	prior2 /= sum(prior2) 
	prior2 *= scale

	prior3 = exp(5(1:d)/d)
	prior3 /= sum(prior3) 
	prior3 *= scale

	for N = 10:50:500
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp)

			zkl1, xkl1 = kl_port(eps_, budget, supp, cnts + prior1, false)
			zkl2, xkl2 = kl_port(eps_, budget, supp, cnts + prior2, false)
			zkl3, xkl3 = kl_port(eps_, budget, supp, cnts + prior3, false)

			outkl1, cvarkl1 = out_perf(xkl1, eps_, supp)
			writecsv(f, [iSim N "prior1" zkl1 outkl1 cvarkl1])

			outkl2, cvarkl2 = out_perf(xkl2, eps_, supp)
			writecsv(f, [iSim N "prior2" zkl2 outkl2 cvarkl2])

			outkl3, cvarkl3 = out_perf(xkl3, eps_, supp)
			writecsv(f, [iSim N "prior3" zkl3 outkl3 cvarkl3])
		end
	end
	close(f)
end

function fullSim(d; budget=3, seed=nothing)
	f = open("Results/real_sim_$(d)_$(budget).csv", "w")
	writecsv(f, ["ix" "Method" "inReturn" "outReturn"])
	seed != nothing && srand(seed)
	T = 200 + d + 1
	eps_ = .1
	mkt = getMktSupp(T)
	cnts = ones(d)

	turnover = zeros(Float64, 5)
	prev_saa = Float64[]
	prev_chisq = Float64[]
	prev_kl = Float64[]
	prev_chisq_cov = Float64[]
	prev_kl_cov = Float64[]

	for ix = d:T-1
		supp = mkt[(ix-d+1):ix, :]

		zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
		zchisq_cov, xchisq_cov = chisq_port(eps_, budget, supp, cnts, true)
		zchisq, xchisq = chisq_port(eps_, budget, supp, cnts, false)
		zkl_cov, xkl_cov = kl_port(eps_, budget, supp, cnts, true)
		zkl, xkl = kl_port(eps_, budget, supp, cnts, false)

		outsaa =dot(vec(mkt[ix+1, :]), xsaa)
		outchisq = dot(vec(mkt[ix+1, :]), xchisq)
		outchisq_cov = dot(vec(mkt[ix+1, :]), xchisq_cov)
		outkl = dot(vec(mkt[ix+1, :]), xkl)
		outkl_cov = dot(vec(mkt[ix+1, :]), xkl_cov)

		writecsv(f, [ix "SAA" zsaa outsaa])
		writecsv(f, [ix "Chisq" zchisq outchisq])
		writecsv(f, [ix "ChisqCov" zchisq_cov outchisq_cov])
		writecsv(f, [ix "KL" zkl outkl])
		writecsv(f, [ix "KLCov" zkl_cov outkl_cov])

		if ix > d
			turnover[1] += norm(xsaa - prev_saa, 1)
			turnover[2] += norm(xchisq - prev_chisq, 1)
			turnover[3] += norm(xkl - prev_kl, 1)
			turnover[4] += norm(xchisq_cov - prev_chisq_cov, 1)
			turnover[5] += norm(xkl_cov - prev_kl_cov, 1)
		end
		prev_saa = xsaa
		prev_chisq = xchisq
		prev_kl = xkl
		prev_chisq_cov = xchisq_cov
		prev_kl_cov = xkl_cov
	end

	#rescale the turnovers
	turnover /= (T-d-1)
	println("Turnover:")
	show( [["SAA", "ChiSq", "KL", "Chisq_cov", "KL_cov"] turnover] )
	println()

	close(f)
end

end #module