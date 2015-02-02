#####
# The second moment VaR
#####

#Cheap 2 moment VaR bound
function secondMomentVar(vs, alphas, eps_)
	const kappa = sqrt(1/eps_-1)
	const phat, alpha0 = calc_phat(alphas)
	const mu = dot(vs, phat)
	const sig = sqrt(vs'*sigma(alphas)*vs)[1]
	vmin = minimum(vs)
	vmax = maximum(vs)
	min(max(mu + kappa * sig, vmin), vmax)
end
