#####
# The second moment VaR
#####

#computes the finite sample sigma
function sig_finite(alphas)
	const d = length(alphas)
	sigma = zeros(d, d)
	const alpha0 = sum(alphas)
	const denom = (alpha0 + 1)alpha0^2
	for i = 1:d
		sigma[i, i] = (alpha0 - alphas[i]) * alphas[i] / denom
		for j=i+1:d
			sigma[i, j] = -alphas[i] * alphas[j]/denom
			sigma[j, i] = sigma[i, j] 
		end
	end
	sigma
end

#Cheap 2 moment VaR bound
function secondMomentVar(vs, alphas, eps_; gauss_approx=false)
	const kappa = gauss_approx ? quantile(Normal(), 1-eps_) : sqrt(1/eps_ - 1)
	const mu = dot(vs, alphas)/sum(alphas)
	const sig = sqrt(vs' * sig_finite(alphas) * vs)[1]
	vmin = minimum(vs)
	vmax = maximum(vs)
	min(max(mu + kappa * sig, vmin), vmax)
end
