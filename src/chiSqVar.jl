#####
# Sets based on the ChiSq (2nd Moment)
#####

#2 moment VaR bound
#thresh is either (1/eps -1)/(alpha0 +1) or else Chisq(1-eps)/N
function chiSqVarSol(vs, alphas, thresh)
	const sqrt_thresh = sqrt(thresh)
	const phat, alpha0 = calc_phat(alphas)
	const TOL = 1e-8 

	#degenerate case since phat sums to one
	if maximum(vs) - minimum(vs)<= TOL
		return maximum(vs), phat
	end

	#first form a solution ignoring the support constranits
	sigmavs = unscaled_sigma_(phat) * vs
	const sig = sqrt(dot(vs, sigmavs))
	pstar = phat + sqrt_thresh * sigmavs / sig
	if minimum(pstar) >= 0
		#for debugging only
		mu = dot(vs, phat)
		n = length(phat)
		var = mu + sqrt_thresh * sig
		@assert abs(dot(pstar, vs) - var) <= TOL
		return dot(vs, pstar), pstar
	else
		#need to solve the actual optimization
		m = Model(solver=GurobiSolver(OutputFlag=0))
		n = length(phat)
		@defVar(m, ps[1:n] >= 0)
		@addConstraint(m, sum{ps[i], i=1:n} == 1)
		@addConstraint(m, sum([(ps[i]-phat[i])^2/phat[i] for i=1:n]) <= thresh)
		@setObjective(m, Max, dot(vs, ps))
		solve(m)
		return getObjectiveValue(m), getValue(ps)
	end
end

chiSqVar(vs, alphas, eps_) = chiSqVarSol(vs, alphas, (1/eps_-1)/(sum(alphas)+1))[1]

#VG Change this so it doesn't require N
function chiSqVarCov(vs, alphas, eps_, N=sum(alphas))
	thresh = quantile(Distributions.Chisq(length(alphas)-1), 1-eps_)/N
	chiSqVarSol(vs, alphas, thresh)[1]
end