#####
# The second moment VaR
#####
using JuMP

function calc_phat(alphas) 
	const alpha0 = sum(alphas)
	alphas/alpha0, alpha0
end

#Sigma in document
function unscaled_sigma_(phat)
	const d = length(phat)
	sigma = zeros(Float64, d, d)
	for i = 1:d
		sigma[i, i] = phat[i]*(1-phat[i])
		for j = i+1:d
			sigma[i, j] = -phat[i] * phat[j]
			sigma[j, i] = sigma[i, j]
		end
	end
	sigma
end

#Finite sample scaling, equiv to assymptotic
function sigma(alphas)
	const phat, alpha0 = calc_phat(alphas)
	unscaled_sigma_(phat)/(alpha0+1)
end

#Inverse of the unscaled SUBmatrix
function invsigma_(alphas)
	phat, alpha0 = calc_phat(alphas)
	n = length(phat)
	diagm(1./phat[1:end-1]) + fill(1./phat[end], n-1, n-1)
end


#2 moment VaR bound
#thresh is either (1/eps -1)/(alpha0 +1) or else Chisq(1-eps)/N
function chiSqVarSol(vs, alphas, thresh)
	const sqrt_thresh = sqrt(thresh)
	const phat, alpha0 = calc_phat(alphas)

	#first form a solution ignoring the support constranits
	sigmavs = unscaled_sigma_(phat) * vs
	const sig = sqrt(dot(vs, sigmavs))
	pstar = phat + sqrt_thresh * sigmavs / sig
	if minimum(pstar) >= 0
		#for debugging only
		mu = dot(vs, phat)
		n = length(phat)
		var = mu + sqrt_thresh * sig
		@assert abs(dot(pstar, vs) - var) <= 1e-8
		return dot(vs, pstar), pstar
	else
		#need to solve the actual optimization
		m = Model()
		n = length(phat)
		@defVar(m, ps[1:n] >= 0)
		@addConstraint(m, sum{ps[i], i=1:n} == 1)
		addConstraint(m, sum([(ps[i]-phat[i])^2/phat[i] for i=1:n]) <= thresh)
		@setObjective(m, Max, dot(vs, ps))
		solve(m)
		return getObjectiveValue(m), getValue(ps)
	end
end

chiSqVar(vs, alphas, eps_) = chiSqVarSol(vs, alphas, (1/eps_-1)/(sum(alphas)+1))[1]
function chiSqVarCov(vs, alphas, eps_, N)
	thresh = quantile(Distributions.Chisq(length(alphas)-1), 1-eps_)/N
	chiSqVarSol(vs, alphas, thresh)[1]
end