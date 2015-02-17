####
# oracles for the sets
####

# The JuMPeR oracle interface
import JuMPeR
import JuMPeR: AbstractOracle, registerConstraint, setup, generateCut, generateReform

abstract AmbiguitySet <: AbstractOracle

#implement this for each ambiguity set
supp_size(w::AmbiguitySet) = error("supp_size not implemented for $(typeof(w))")
suppFcn(vs, w::AmbiguitySet, cut_sense) = error("suppFcn not implemented for $(typeof(w))")

#####
# Common boiler plate for the oracles
# If you fold this into DDUS, you can tie into DDUSoracles
#####
# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::AmbiguitySet, rm::Model, ind::Int, prefs) =
    !get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")

# JuMPeR wants us to generate a constraint for every uncertain constraint in inds
function generateCut(w::AmbiguitySet, m::Model, rm::Model, inds::Vector{Int}, active=false)
    new_cons = {}
    rd = JuMPeR.getRobust(rm)
    for ind in inds
        con = JuMPeR.get_uncertain_constraint(rm, ind)
        cut_sense, xs, lhs_const = JuMPeR.build_cut_objective(rm, con, m.colVal)
        zstar, ustar = suppFcn(xs, w, cut_sense)
        lhs_of_cut = zstar + lhs_const

        # SUBJECT TO CHANGE: active cut detection
        if active
            push!(rd.activecuts[ind], 
                JuMPeR.cut_to_scen(ustar, 
                    JuMPeR.check_cut_status(con, lhs_of_cut, w.cut_tol) == :Active))
            continue
        end

        # Check violation
        if JuMPeR.check_cut_status(con, lhs_of_cut, w.cut_tol) != :Violate
            w.debug_printcut && JuMPeR.debug_printcut(rm ,m,w,lhs_of_cut,con,nothing)
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(m, con, ustar)
        w.debug_printcut && JuMPeR.debug_printcut(rm, m, w, lhs_of_cut, con, new_con)
        push!(new_cons, new_con)
    end
    return new_cons
end

# JuMPeR asking us for any reformulations we might want to make - we make none
generateReform(w::AmbiguitySet, m::Model, rm::Model, inds::Vector{Int}) = 0

function setup(w::AmbiguitySet, rm::Model, prefs)
    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, w.cut_tol)

    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == supp_size(w)) "Num Uncertainties $(rd.numUncs) doesn't match columns in data $(supp_size(w))"

    #ignore any additional constraints on uncertainties for now
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"
end

##########
# Implements the actual support functions for each of the sets
##########
type OptDirichlet <:AmbiguitySet
	alphas
	eps_

	# Assumed to be members of all ambiguity sets
	debug_printcut #false
	cut_tol #1e-6
end
supp_size(w::OptDirichlet) = length(w.alphas)
OptDirichlet(alphas, eps_) = OptDirichlet(alphas, eps_, false, 1e-6)

function suppFcn(vs, w::OptDirichlet, cut_sense)
	toggle = 1
    if cut_sense == :Min
        vs = -vs
        toggle = -1
    end
    var = VaR(vs, w.alphas, w.eps_)
    grad = grad_VaR(vs, w.alphas, var)

    return var*toggle, grad
end

type KLSet <:AmbiguitySet
	phat
	thresh

	# Assumed to be members of all ambiguity sets
	debug_printcut #false
	cut_tol #1e-6
end
KLSet(alphas, thresh) = KLSet(alphas/sum(alphas), thresh, false, 1e-6)
function KLSet_eps(alphas, eps_, useCover=false) 
	if !useCover
		KLSet(alphas, log(1/eps_)/sum(alphas))
	else
		n = length(alphas)
		KLSet(alphas, quantile(Distributions.Chisq(n-1), 1-eps_)/sum(alphas)/2)
	end
end

supp_size(w::KLSet) = length(w.phat)

function suppFcn(vs, w::KLSet, cut_sense)
	toggle = 1
    if cut_sense == :Min
        vs = -vs
        toggle = -1
    end
    var, pstar = KLVarSol(vs, w.phat, w.thresh)
    var*toggle, pstar
end

type ChiSqSet <:AmbiguitySet
	phat
	thresh

	# Assumed to be members of all ambiguity sets
	debug_printcut #false
	cut_tol #1e-6
end
ChiSqSet(alphas, thresh) = ChiSqSet(alphas/sum(alphas), thresh, false, 1e-6)
function ChiSqSet_eps(alphas, eps_, useCover=false) 
	if !useCover
		ChiSqSet(alphas, log(1/eps_)/sum(alphas))
	else
		n = length(alphas)
		ChiSqSet(alphas, quantile(Distributions.Chisq(n-1), 1-eps_)/sum(alphas))
	end
end

function suppFcn(vs, w::ChiSqSet, cut_sense)
	toggle = 1
    if cut_sense == :Min
        vs = -vs
        toggle = -1
    end
    var, pstar = chiSqVarSol(vs, w.phat, w.thresh)
    var*toggle, pstar
end
	

