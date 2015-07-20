###
# Portfolio Experiment for the presentation
###
include("../src/DirRes.jl")
module port
using Distributions, Gurobi, JuMP, JuMPeR, Dir

## build a dataset
## outputs a matrix of the various support pts, and a vector of counts
function buildMkt(N; numAssets=5)
    #the marginal distn of each asset
    ps = [.2, .2, .2, .2, .2]
    #srand(8675309)
    dist = Categorical(ps)

    #generate the raw returns
    vals = Dict()
    for i = 1:N
        rets = rand(dist, numAssets)
        if haskey(vals, rets)
            vals[rets] += 1/N
        else
            vals[rets] = 1/N
        end
    end

    #map the raw returns to supp and counts
    marg_supp = linspace(-.1, .1, 5)
    supp = zeros(length(vals), numAssets)
    cnts = zeros(length(vals))
    for (i, k) in enumerate(keys(vals))
        for j = 1:numAssets
            supp[i, j] = marg_supp[k[j]]
        end
        cnts[i] = vals[k]
    end
    supp, cnts
end


function buildMkt2(d, N; numAssets=5)
	#sample supp from an interesting continuous distribution
	supp = randn(d, numAssets)
	supp = supp * diagm(linspace(.1, .3, numAssets))
	supp = broadcast(+, linspace(0, .1, numAssets)', supp)

	#now sample counts from this distribution	
	dist = Categorical(ones(d)/d)	
	sample = rand(dist, N)
	cnts = zeros(d)
	for j = 1:N
		cnts[sample[j]] += 1
	end
	supp, cnts
end


# function saa_port(eps_, budget, supp, counts)
# 	# solve for the SAA portfolio
# 	m = Model(solver=GurobiSolver())
# 	d, numAssets = size(supp)
# 	N = sum(counts); phat = counts / N
# 	rmax = maximum(supp)
# 	@defVar(m, xs[1:numAssets]>=0)
# 	@addConstraint(m, sum(xs) <= 1)
# 	@defVar(m, zs[1:d], Bin) #whether realization exceeds budget
# 	@setObjective(m, Max, sum{phat[j] * dot(vec(supp[j, :]), xs), j=1:d})

# 	for j=1:d
# 		#if ret < budget, then z = 1
# 		@addConstraint(m, dot(vec(supp[j, :]), xs) >= budget - zs[j]*rmax)
# 	end
# 	@addConstraint(m, sum{phat[j]*zs[j], j=1:d} <= eps_)
# 	solve(m)
# 	getObjectiveValue(m), getValue(xs[:])
# end

#CVaR constraint
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
	#define CvAr
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
	addConstraint(m, sum([ps[j]*dot(vec(supp[j, :]), xs) for j=1:d]) >= obj)
	@setObjective(m, Max, obj)

	for j=1:d
		@addConstraint(m, -dot(vec(supp[j, :]), xs[:])  - Beta <= zs[j])
	end
	addConstraint(m, sum([ps[j]*zs[j]/eps_ for j=1:d]) + Beta <= budget)	
	solveRobust(m, prefer_cuts=true)
	getObjectiveValue(m), getValue(xs[:])
end

#proxy CVaR by sorting
function cvar_sort(rets, phat, eps_)
	inds = sortperm(rets)
	prob = 0.; 
	cvar = 0
	for ix = 1:length(rets)
		if prob >= eps_
			break
		else
			prob += phat[inds[ix]]
			cvar += phat[inds[ix]] * rets[inds[ix]]
		end
	end
	-cvar / prob
end

#assess out-of-sample performance
#return the mean and cvar at eps_
#VG uses assumption that pstar = 1/N
function out_perf(x, eps_, N_out, supp)
	d = size(supp, 1)
	pstar = ones(d)/d
	rets = supp * x
	(pstar' * supp * x)[1], cvar_sort(rets, pstar, eps_)
end


function test(d, N, N_out; budget=.2, seed=nothing)
	seed != nothing && srand(seed)
	
	supp, cnts = buildMkt2(d, N)
	zsaa, xsaa = saa_port(.1, budget, supp, cnts)
	zrob, xrob = chisq_port(.1, budget, supp, cnts, false)
	zbad, xbad = chisq_port(.1, budget, supp, cnts, true)

	println("Saa:\t $zsaa")
	show(xsaa'); println()
	println("Out:\t", out_perf(xsaa, .1, N_out, supp))

	println()
	println("Rob:\t $zrob")
	show(xrob'); println()
	println("Out:\t", out_perf(xrob, .1, N_out, supp))

	println()
	println("Bad:\t $zbad")
	show(xbad'); println()
	println("Out:\t", out_perf(xbad, .1, N_out, supp))
end


## dump some stuff to files to look at plots
function test2(d, N_out, numRuns=100; budget=.2, seed=nothing, eps_=.1)
	f = open("portExp_$(d)_$(budget).csv", "w")
	writecsv(f, ["Run" "N" "Method" "inReturn" "outReturn" "CVaR"])

	srand(8675309)
	for N in 100:100:1000

		for iSim = 1:numRuns
			supp, cnts = buildMkt2(d, N)
			zsaa, xsaa = saa_port(eps_, budget, supp, cnts)
			zrob, xrob = chisq_port(eps_, budget, supp, cnts, false)
			zbad, xbad = chisq_port(eps_, budget, supp, cnts, true)

			outsaa, cvarsaa = out_perf(xsaa, eps_, N_out, supp)
			writecsv(f, [iSim N "SAA" zsaa outsaa cvarsaa])

			outrob, cvarrob = out_perf(xrob, eps_, N_out, supp)
			writecsv(f, [iSim N "Rob" zrob outrob cvarrob])

			outbad, cvarbad = out_perf(xbad, eps_, N_out, supp)
			writecsv(f, [iSim N "Bad" zbad outbad cvarbad])
		end
	end
	close(f)
end

end #moduleBadbad