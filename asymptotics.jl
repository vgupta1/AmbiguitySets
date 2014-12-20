#Asymptotic Sigma (without sqrt(alpha_0) scaling)
function asymptotic_sigma(phat)
	const d = length(phat)
	sigma = zeros(d, d)
	for i = 1:d
		sigma[i, i] = phat[i]*(1-phat[i])
		for j = i+1:d
			sigma[i, j] = -phat[i] * phat[j]
			sigma[j, i] = sigma[i, j]
		end
	end
	sigma
end

function asymptotic_var(vs, alphas, eps_)
	const alpha0 = sum(alphas)
	const phat = alphas / alpha0
	(dot(vs, phat) + quantile(Normal(), 1-eps_) / sqrt(alpha_0) * sqrt(vs'*asymptotic_sigma(phat)*vs)
end
