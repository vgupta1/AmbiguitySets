#####
# Some simple functions to assess the asymptotic 
# performance
#####

#Convenience wrapper
#probably exists somewhere in base.
const sqrt_2 = sqrt(2)
invnorm(x) = sqrt_2 * erfinv(2x-1)

#asymptotic var given by gaussian approximation
function asymptotic_var(vs, alphas, eps_)
	const phat, alpha0 = calc_phat(alphas)
	dot(vs, phat) + 
		invnorm(1-eps_)*sqrt(vs' * sigma(alphas) * vs)
end

#Ratios to optimality for near-optimal sets
kl_const(eps_) = sqrt(2*log(1/eps_)) / invnorm(1-eps_)
chisq_const(eps_) = sqrt(1/eps_ - 1) / invnorm(1-eps_)

#Ratios for the credible regions
#VG change these to all reference same functin...
kl_cov_const(eps_, d) = 
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)

chisq_cov_const(eps_, d) = 
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)

#the general covering ratio
cov_const(eps_, d) = 
	sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)) / invnorm(1-eps_)


#Approximate the CVaR of rets by sorting
#Only accurate up to discretization error
function cvar_sort(rets, probs, eps_)
    inds = sortperm(rets)
    cum_prob = 0.; 
    cvar = 0
    for ix = 1:length(rets)
        if cum_prob >= eps_
            break
        else
            cum_prob += probs[inds[ix]]
            cvar += probs[inds[ix]] * rets[inds[ix]]
        end
    end
    -cvar / cum_prob
end

#compute the 1-eps_ value at risk of the (univariate) random variable 
#taking values vals with probability probs
#only accurate up to discretization
function VaR(vals, probs, eps_)
    const d = length(probs)
    @assert 0 < eps_ < 1
    indx_sort = sortperm(vals)
    cum_prob = 0.
    ix = 0
    while cum_prob <1 - eps_
        @assert ix < d
        ix += 1
        cum_prob += probs[indx_sort[ix]]
    end
    vals[indx_sort[ix]]
end

#Compute 1-eps_ CVaR by sorting for r.v. as above
# Again only accurate up to discretization
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

function calc_phat(alphas) 
    const alpha0 = sum(alphas)
    alphas/alpha0, alpha0
end

#Unscaled Sigma is diag(p) - p*pT
#VG replace this function using kron
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

