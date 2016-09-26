###
# Top level module for all the functionality
###
module Dir
using Distributions, Roots, JuMP, Gurobi

include("asymptotics.jl")
include("chiSqVar.jl")
include("klVar.jl")
include("oracles.jl")
include("dirichletProb.jl")

function kl_ratio(vs, alphas, eps_)
    phat, alpha0 = calc_phat(alphas)
    Gamma = log(1/eps_) / alpha0
    bval = KLVarSol(vs, phat, Gamma)[1] - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    bval/val
end

function kl_cov_ratio(vs, alphas, eps_, d, N)
    phat, alpha0 = calc_phat(alphas)
    Gamma = quantile(Chisq(d-1), 1-eps_)/2N
            bval = KLVarSol(vs, phat, Gamma)[1] - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    bval/val
end

function chisq_ratio(vs, alphas, eps_)
    phat, alpha0 = calc_phat(alphas)
    sm_val = chiSqVar(vs, alphas, eps_) - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    sm_val/val
end

function chisq_cov_ratio(vs, alphas, eps_, d, N)
    phat, alpha0 = calc_phat(alphas)
    sm_val = chiSqVarCov(vs, alphas, eps_, N) - dot(phat, vs)
    val = VaR(vs, alphas, eps_, prob_direct) - dot(phat, vs)    
    sm_val / val
end

end #ends module

