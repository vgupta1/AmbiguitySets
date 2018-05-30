###
# Portfolio Experiment for the paper based on kenneth french data
###
include("../src/DirRes.jl")
module port
using Distributions, Gurobi, JuMP, JuMPeR, Dir, Iterators

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
function buildSynthMkt(N, supp, probs=Float64[])
	#sample counts
	#now sample counts from this distribution
	const d = size(supp, 1)
	if isempty(probs)
		probs = ones(d)/d
	end
	dist = Categorical(probs)
	sample = rand(dist, N)
	cnts = zeros(d)
	for j = 1:N
		cnts[sample[j]] += 1
	end
	cnts
end

#get cluster Support and probs
function clusterMkt(path="Data//ClusterScenarios.csv")
	dat, header = readcsv(path, header=true)
	#drop first column from R's weird row naming,1 for cluster name
	supp = dat[:, 3:end-1]
	probs= dat[:, end] 
	supp, probs
end

#assess out-of-sample performance
#return the mean and cvar at eps_
function out_perf(x, eps_, supp, pstar=Float64[])
	d = size(supp, 1)
	if isempty(pstar)		
		pstar = ones(d)/d
	end
	rets = supp * x
	(pstar' * supp * x)[1], Dir.cvar_sort(rets, pstar, eps_)
end

###
# CVAR-constrained portfolio allocation problems for various methods
###

#1/N Diversification
#Rescaled so that CVaR_eps is close to budget
#if budget < 0, no scaling
function naive_port(eps_, budget, supp, counts)
	d, numAssets = size(supp)
	phat, alpha0 = Dir.calc_phat(counts)
	xs = ones(numAssets)/numAssets
	if budget > 0 
		cvar = Dir.cvar_sort(supp * xs, phat, eps_)
		if cvar >= budget 
			xs /= cvar/budget 
		end
	end
	phat' * supp * xs, xs
end

#The minimum variance portfolio
#If budget < 0, no scaling of CVAR
function min_var_port(eps_, budget, supp, counts)
	m = Model(solver=GurobiSolver(OutputFlag=false))
	d, numAssets = size(supp)
	phat, alpha0 = Dir.calc_phat(counts)

	@defVar(m, xs[1:numAssets] >= 0)
	@addConstraint(m, sum(xs) == 1)
	@defVar(m, mu)
	@setObjective(m, Min, sum{phat[ix]*(dot(vec(supp[ix, :]), xs) - mu)^2, ix = 1:d})
	solve(m)
	xs = getValue(xs[:])

	#rescale to get CVaR close to budget
	if budget > 0
		cvar = Dir.cvar_sort(supp * xs, phat, eps_)
		if cvar >= budget 
			xs /= cvar/budget  #ensures new cvar close to budget
		end
	end
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
	#supp = getMktSupp(d)
	#cnts = buildSynthMkt(N, supp)

	supp, probs = clusterMkt()
	cnts = buildSynthMkt(N, supp, probs)

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

# Generates data for plots for dependence in N plots
function test2(d=72, numRuns=1000; budget=3, seed=nothing, eps_=.1)
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

## Generates the dependence in N, synthetic cluster-based market
function test3(numRuns=100; budget=3, seed=nothing, eps_=.1)
	f = open("Results/portExp3_$(budget).csv", "w")
	writecsv(f, ["Run" "N" "Method" "inReturn" "outReturn" "CVaR" "X_Norm"])
	seed != nothing && srand(seed)
	supp, probs = clusterMkt()

	for N in 100:100:1000
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp, probs)

			zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
			outsaa, cvarsaa = out_perf(xsaa, eps_, supp, probs)
			writecsv(f, [iSim N "SAA" zsaa outsaa cvarsaa norm(xsaa)])

			zchisq, xchisq = chisq_port(eps_, budget, supp, cnts, false)
			outchisq, cvarchisq = out_perf(xchisq, eps_, supp, probs)
			writecsv(f, [iSim N "Chisq" zchisq outchisq cvarchisq norm(xchisq)])

			zchisq_cov, xchisq_cov = chisq_port(eps_, budget, supp, cnts, true)
			out_chi_cov, cvar_chi_cov = out_perf(xchisq_cov, eps_, supp, probs)
			writecsv(f, [iSim N "ChisqCov" zchisq_cov out_chi_cov cvar_chi_cov norm(xchisq_cov)])

			zkl, xkl = kl_port(eps_, budget, supp, cnts, false)
			outkl, cvarkl = out_perf(xkl, eps_, supp, probs)
			writecsv(f, [iSim N "KL" zkl outkl cvarkl norm(xkl)])

			zkl_cov, xkl_cov = kl_port(eps_, budget, supp, cnts, true)
			outkl_cov, cvarkl_cov = out_perf(xkl_cov, eps_, supp, probs)
			writecsv(f, [iSim N "KLCov" zkl_cov outkl_cov cvarkl_cov norm(xkl_cov)])

			z_naive, x_naive = naive_port(eps_, budget, supp, cnts)
			out_naive, cvar_naive = out_perf(x_naive, eps_, supp, probs)
			writecsv(f, [iSim N "Naive" z_naive out_naive cvar_naive norm(x_naive)])

			z_minvar, x_minvar = min_var_port(eps_, budget, supp, cnts)
			out_minvar, cvar_minvar = out_perf(x_minvar, eps_, supp, probs)
			writecsv(f, [iSim N "MinVar" z_minvar out_minvar cvar_minvar norm(x_minvar)])
		end
	end
	close(f)
end

#Creates the Dependence in D plots for N = 300 and N = 700  
function test_incr_d(; N=300, numRuns=1000, budget=3, eps_=.1, seed=nothing)
	f = open("Results/incrDExp_a_$(N)_$(budget).csv", "w")
	writecsv(f, ["Run" "d" "Method" "inReturn" "outReturn" "CVaR"])
	seed != nothing && srand(seed)

	for numPts = 24:12:120
		supp = getMktSupp(numPts)
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp)
			zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
			outsaa, cvarsaa = out_perf(xsaa, eps_, supp)
			writecsv(f, [iSim numPts "SAA" zsaa outsaa cvarsaa])

			zchisq, xchisq = chisq_port(eps_, budget, supp, cnts, false)
			outchisq, cvarchisq = out_perf(xchisq, eps_, supp)
			writecsv(f, [iSim numPts "Chisq" zchisq outchisq cvarchisq])

			zchisq_cov, xchisq_cov = chisq_port(eps_, budget, supp, cnts, true)
			out_chi_cov, cvar_chi_cov = out_perf(xchisq_cov, eps_, supp)
			writecsv(f, [iSim numPts "ChisqCov" zchisq_cov out_chi_cov cvar_chi_cov])

			zkl, xkl = kl_port(eps_, budget, supp, cnts, false)
			outkl, cvarkl = out_perf(xkl, eps_, supp)
			writecsv(f, [iSim numPts "KL" zkl outkl cvarkl])

			zkl_cov, xkl_cov = kl_port(eps_, budget, supp, cnts, true)
			outkl_cov, cvarkl_cov = out_perf(xkl_cov, eps_, supp)
			writecsv(f, [iSim numPts "KLCov" zkl_cov outkl_cov cvarkl_cov])

			z_minvar, x_minvar = min_var_port(eps_, budget, supp, cnts)
			out_minvar, cvar_minvar = out_perf(x_minvar, eps_, supp)
			writecsv(f, [iSim numPts "MinVar" z_minvar out_minvar cvar_minvar])

			z_naive, x_naive = naive_port(eps_, budget, supp, cnts)
			out_naive, cvar_naive = out_perf(x_naive, eps_, supp)
			writecsv(f, [iSim numPts "Naive" z_naive out_naive cvar_naive])
		end
	end
	close(f)
end

##VG Seems like we can drop
#Generate data and a wrong prior and then assess
# function test_wrongPrior(d, N, numRuns; budget=3, seed=8675309)
# 	f = open("Results/randWrongPrior_$(N)_$(d)_$(seed).csv", "w")
# 	writecsv(f, ["Run" "Dist" "Scale" "Method" "inReturn" "outReturn" "CVaR"])
# 	seed != nothing && srand(seed)
# 	eps_ = .1
# 	scale_grid = collect(.1:.1:3.)

# 	supp = getMktSupp(d)
# 	cnts = buildSynthMkt(N, supp)
# 	pstar = ones(d)/d

# 	for scale in scale_grid
# 		for iSim = 1:numRuns
# 			pseudo_cnts = rand(d)
# 			pseudo_cnts /= sum(pseudo_cnts)
# 			#calculate the distance
# 			dist = norm(pseudo_cnts - pstar)

# 			#scale up
# 			pseudo_cnts *= N * scale

# 			zchisq, xchisq = chisq_port(eps_, budget, supp, cnts+ pseudo_cnts, false)
# 			zkl, xkl = kl_port(eps_, budget, supp, cnts + pseudo_cnts, false)

# 			outchisq, cvarchisq = out_perf(xchisq, eps_, supp)
# 			writecsv(f, [iSim dist scale "Chisq" zchisq outchisq cvarchisq])

# 			outkl, cvarkl = out_perf(xkl, eps_, supp)
# 			writecsv(f, [iSim dist scale "KL" zkl outkl cvarkl])
# 		end
# 	end
# 	close(f)
# end


#An increasingly incorrect prior
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

function test_wrongPrior_IncrN(;d=72, numRuns=1000, budget=3, seed=8675309, 
								eps_=.1, tau0=175, path="Results/wrongPrior_IncrN")
	f = open("$(path)_$(d)_$(budget).csv", "w")
	writecsv(f, ["Run" "N" "Method" "inReturn" "outReturn" "CVaR" "X_norm"])
	seed != nothing && srand(seed)
	supp = getMktSupp(d)
	prior = ones(d)
	prior[1] = tau0

	for N in 100:100:1000
		for iSim = 1:numRuns
			cnts = buildSynthMkt(N, supp)

			zkl, xkl = kl_port(eps_, budget, supp, cnts + prior, false)
			outkl, cvarkl = out_perf(xkl, eps_, supp)
			writecsv(f, [iSim N "KL" zkl outkl cvarkl norm(xkl)])
		end
	end
	close(f)
end

##VG Seems like we can drop
# function test_wrongPriorConv(d, numRuns; budget=3, seed=nothing, scale=20)
# 	f = open("Results/randWrongPrior2_$(scale).csv", "w")
# 	writecsv(f, ["Run" "N" "Method" "inReturn" "outReturn" "CVaR"])
# 	seed != nothing && srand(seed)
# 	eps_ = .1

# 	supp = getMktSupp(d)

# 	prior1 = collect(1:d)/d 
# 	prior1 /= sum(prior1)
# 	prior1 *= scale

# 	prior2 = exp((1:d)/d)
# 	prior2 /= sum(prior2) 
# 	prior2 *= scale

# 	prior3 = exp(5(1:d)/d)
# 	prior3 /= sum(prior3) 
# 	prior3 *= scale

# 	for N = 10:50:500
# 		for iSim = 1:numRuns
# 			cnts = buildSynthMkt(N, supp)

# 			zkl1, xkl1 = kl_port(eps_, budget, supp, cnts + prior1, false)
# 			zkl2, xkl2 = kl_port(eps_, budget, supp, cnts + prior2, false)
# 			zkl3, xkl3 = kl_port(eps_, budget, supp, cnts + prior3, false)

# 			outkl1, cvarkl1 = out_perf(xkl1, eps_, supp)
# 			writecsv(f, [iSim N "prior1" zkl1 outkl1 cvarkl1])

# 			outkl2, cvarkl2 = out_perf(xkl2, eps_, supp)
# 			writecsv(f, [iSim N "prior2" zkl2 outkl2 cvarkl2])

# 			outkl3, cvarkl3 = out_perf(xkl3, eps_, supp)
# 			writecsv(f, [iSim N "prior3" zkl3 outkl3 cvarkl3])
# 		end
# 	end
# 	close(f)
# end

#Distance metrics to consider for assessing prior misspec
# Exp ChiSq distance and Exp Modified Chi Sq distance 
# MSE of the performance of full-info portfolio
# Suboptimality of the full-info portfolio relative to optimal portfolio for that prior

#expected mse between prior realization and pstar
mse(d_prior, pstar) = norm(mean(d_prior) - pstar)^2 + sum(var(d_prior))

#expected mse of performance of xstar under prior and pstar
function port_mse(xstar, supp, taus, pstar)
	mu = Dir.calc_phat(taus)[1]
	Sigma = Dir.sigma(taus)
	xstar' * supp' * (Sigma + (mu-pstar)*(mu-pstar)') * supp * xstar	
end

#upperbound on the exp 2 norm of diff between mean perf of assests
function asset_mse(supp, taus, pstar)
	mu = Dir.calc_phat(taus)[1]
	Sigma = Dir.sigma(taus)
	sqrt(dot(vec((Sigma + (mu-pstar)*(mu-pstar)')), vec(supp*supp')))
end

#randomly generate priors and then assess their performance
function randomWrongPriors(; d=72, numSims=300, numPriors=100, budget=3, eps_ = .1, 
							seed=8675309, path = "Results/random_wrong_priors",
							strengths = [.1, .25, .5, .75, 1., 1.25, 1.5], 
							N_grid = [300])
	srand(seed)
	f = open("$(path)_$(d)_$(budget).csv", "w")
	writecsv(f, ["iPrior" "iSim" "N" "strength" "MSE" "portMSE" "AssetMSE" "Method" "outReturn" "outCVaR"])
	const pstar = ones(d)/d
	supp = getMktSupp(d)

	#compute the full-information portfolio for comparisons
	zstar, xstar = saa_port(eps_, budget, supp, pstar)

	#We will generate priors uniformly on the simplex with strength given by a fraction of N
	prior_generator = Dirichlet(ones(d))

	for (N, strength) in product(N_grid, strengths)
		tau0 = int(N * strength) + d
		for iPrior = 1:numPriors
			taus = rand(prior_generator) * tau0
			prior = Dirichlet(taus)
			mode = (taus - 1)/(sum(taus) - d)

			#calc the distances
			MSE = mse(prior, pstar)
			PMSE = port_mse(xstar, supp, taus, pstar)
			AMSE = asset_mse(supp, taus, pstar)

			t = tic()
			for iSim = 1:numSims
				#generate some data
				cnts = buildSynthMkt(N, supp)

				#Form portfolio with prior and record
				zchisq, xchisq = chisq_port(eps_, budget, supp, cnts + taus, false)			
				ret_chisq, cvar_chisq = out_perf(xchisq, eps_, supp)
				writecsv(f, [iPrior iSim N strength MSE PMSE AMSE "ChiSq" ret_chisq cvar_chisq])
	
				zkl, xkl = kl_port(eps_, budget, supp, cnts + taus, false)			
				ret_kl, cvar_kl = out_perf(xkl, eps_, supp)
				writecsv(f, [iPrior iSim N strength  MSE PMSE AMSE  "KL" ret_kl cvar_kl])
			end
			toc()
		end #end priors loop
		flush(f)
	end
	close(f)
end

function fullSim(d; budget=4)
	f = open("Results/real_sim_$(d)_$(budget).csv", "w")
	writecsv(f, ["ix" "Method" "inReturn" "outReturn"])
	T = 200 + d + 1
	eps_ = .1
	mkt = getMktSupp(T)
	cnts = ones(d)

	turnover = zeros(Float64, 9)
	prev_saa = Float64[]
	prev_chisq = Float64[]
	prev_kl = Float64[]
	prev_chisq_cov = Float64[]
	prev_kl_cov = Float64[]
	prev_naive = Float64[]
	prev_minvar = Float64[]
	prev_naive_unsc = Float64[]
	prev_minvar_unsc = Float64[]

	for ix = d:T-1
		supp = mkt[(ix-d+1):ix, :]

		zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
		outsaa =dot(vec(mkt[ix+1, :]), xsaa)

		tic()
		zchisq, xchisq = chisq_port(eps_, budget, supp, cnts, false)
		outchisq = dot(vec(mkt[ix+1, :]), xchisq)
		toc()

		zkl, xkl = kl_port(eps_, budget, supp, cnts, false)
		outkl = dot(vec(mkt[ix+1, :]), xkl)

		zchisq_cov, xchisq_cov = chisq_port(eps_, budget, supp, cnts, true)
		outchisq_cov = dot(vec(mkt[ix+1, :]), xchisq_cov)

		zkl_cov, xkl_cov = kl_port(eps_, budget, supp, cnts, true)
		outkl_cov = dot(vec(mkt[ix+1, :]), xkl_cov)

		z_naive, x_naive = naive_port(eps_, budget, supp, cnts)
		out_naive = dot(vec(mkt[ix+1, :]), x_naive)

		z_naive_unsc, x_naive_unsc = naive_port(eps_, -1, supp, cnts)
		out_naive_unsc = dot(vec(mkt[ix+1, :]), x_naive_unsc)

		z_minvar, x_minvar = min_var_port(eps_, budget, supp, cnts)
		out_minvar = dot(vec(mkt[ix+1, :]), x_minvar)

		z_minvar_usc, x_minvar_unsc = min_var_port(eps_, -1, supp, cnts)
		out_minvar_unsc = dot(vec(mkt[ix+1, :]), x_minvar_unsc)

		writecsv(f, [ix "SAA" zsaa outsaa])
		writecsv(f, [ix "Chisq" zchisq outchisq])
		writecsv(f, [ix "ChisqCov" zchisq_cov outchisq_cov])
		writecsv(f, [ix "KL" zkl outkl])
		writecsv(f, [ix "KLCov" zkl_cov outkl_cov])
		writecsv(f, [ix "Naive" z_naive out_naive])
		writecsv(f, [ix "NaiveUnsc" z_naive_unsc out_naive_unsc])
		writecsv(f, [ix "MinVar" z_minvar out_minvar])
		writecsv(f, [ix "MinVarUnsc" z_minvar_usc out_minvar_unsc])

		if ix > d
			turnover[1] += norm(xsaa - prev_saa, 1)
			turnover[2] += norm(xchisq - prev_chisq, 1)
			turnover[3] += norm(xkl - prev_kl, 1)
			turnover[4] += norm(xchisq_cov - prev_chisq_cov, 1)
			turnover[5] += norm(xkl_cov - prev_kl_cov, 1)
			turnover[6] += norm(x_naive - prev_naive, 1)
			turnover[7] += norm(x_minvar - prev_minvar, 1)			
			turnover[8] += norm(x_naive_unsc - prev_naive_unsc, 1)
			turnover[9] += norm(x_minvar_unsc - prev_minvar_unsc, 1)
		end
		prev_saa = xsaa
		prev_chisq = xchisq
		prev_kl = xkl
		prev_chisq_cov = xchisq_cov
		prev_kl_cov = xkl_cov
		prev_naive = x_naive
		prev_naive_unsc = x_naive_unsc
		prev_minvar = x_minvar
		prev_minvar_unsc = x_minvar_unsc
	end

	#rescale the turnovers
	turnover /= (T-d-1)
	println("Turnover:")
	show( [["SAA", "ChiSq", "KL", "Chisq_cov", "KL_cov", "Naive", "MinVar", "NaiveUnsc", "MinVarUnsc"] turnover] )
	println()

	close(f)
end

end #module