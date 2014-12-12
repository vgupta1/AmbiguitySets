module Dir

using Distributions, Roots, Optim

include("bernstein.jl")
include("asymptotics.jl")

#Problematic Test Cases
#vbars [0.24501246861638076,-0.16745207467440992,-0.11692969664839702]
#alphas = 10 * [3, 3, 3]

#Even better.. Neither pos_congtour nor neg works
#vbars [-0.3129438949072508,0.6320289286393632,0.6443383495423485,-0.8299381500899872,-0.3285773443404617]
#balphas = 20 * ones(5)

#Approximate the residue of fun at z0
#fun takes a complex argument
function numRes(fun, z0, r, m, useAdaptiveInt = false)
	f(t) = r * exp(2im * pi * t) * fun(z0 + r * exp(2im * pi * t) )
	if useAdaptiveInt
		quadgk(f, 0., 1.)[1]
	else
		t_grid = map(BigFloat, linspace(0, 1, m+1)[1:end-1] )
		t = map(f, linspace(0, 1, m+1)[1:end-1])
		sum(t)/m
	end
end

#the relevant integrand
#vbar = v - te
fint(z, vbars, alphas) = prod((1 + z*vbars).^-alphas)/z

#prob_fun(vs, alphas) => Prob( X(vs, alphas) > 0 )
#VG Revisit this to add a timing function
function Var(alphas, vs, eps_, prob_fun; 
				tmax=bernVar(vs, alphas, eps_), 
				tmin=dot(vs, alphas)/sum(alphas), 
				tol=1e-6)
	fzero(t->prob_fun(vs -t, alphas)-eps_, [tmin, tmax], tol=tol)
end

##############
function sort_poles(vbars, alphas)
	ps        = -1./vbars
	indxs     = sortperm(ps)
	vbars[indxs], ps[indxs], alphas[indxs]
end

#finds the optimal positive contour
function calc_shift(vbars_, alphas_)
	ps = -1./vbars_
	alphas =alphas_[ps .> 0]
	ps = ps[ps .> 0]
	ix = indmin(ps)
	ps[ix]/(1 + alphas[ix])
end

#finds the index of the closest pole to ps[ix].  Returns 0 if it's the origin
function closest_pole(ps, ix)
	dist = abs(ps[ix]); min_indx = 0
	#check the ordinary poles
	for jx = 1:length(ps)
		if ix == jx
			continue
		elseif abs(ps[ix] - ps[jx]) <= dist
			dist = abs(ps[ix] - ps[jx])
			min_indx = jx
		end
	end
	return dist, min_indx
end

#approximates a good raduii for numRes and conditioning constants
function optRs_NumRes(vbars, ps, alphas)
	rs      = zeros(length(ps))  			#radius for integration
	lams    = zeros(length(ps))			#conditioning of the integrand

	#This is the inefficient way....  O(n^2) operation instead of O(n)
	const n = length(ps)
	for ix = 1:n
		rmax, jx = closest_pole(ps, ix)
		if jx > 0 
			rs[ix] = alphas[ix] * rmax / (alphas[ix] + alphas[jx])
		else
			rs[ix] = alphas[ix] * rmax/ (alphas[ix] + 1)
		end
		lams[ix] = abs((vbars[ix] * rs[ix])^(2*alphas[ix] - 1))
	end	
	rs, lams
end

#evaluates the critical integral via numeric residues
function integral_residue(vbars_, alphas_; usePosContour=true, m=100, rs =nothing, useAdaptiveInt=false)
	vbars, ps, alphas = sort_poles(vbars_, alphas_)
	
	if rs == nothing
		rs, lams = optRs_NumRes(vbars, ps, alphas)
	end

	residues = BigFloat[]
	for (p, v, alpha, r) in zip(ps, vbars, alphas, rs)
		if (usePosContour && p > 0) || (!usePosContour && p < 0)
			push!(residues, real(numRes(z-> fint(z, vbars, alphas), p, r, m, useAdaptiveInt)))
		end
	end
	println("Residues:\n", residues)
	usePosContour ? -sum_kbn(residues) : 1 + sum_kbn(residues)
end


#compute the sum S_m from notes... an approximation to the residue at -1/v_k
function numResFourier(vbars, ps, alphas, k, r, m=100)
	const n = length(ps)
	indx_not_k = [1:n .!= k]
	function g(z, vbars, alphas, k) 
		terms = (1 + z*vbars).^(-alphas)
		prod(terms[indx_not_k]) / z
	end

	t_grid = map(BigFloat, linspace(0, 1, m+1)[1:end-1])
	g_out = map(t-> g(ps[k] + r * exp(2pi * 1im * t), vbars, alphas, k) * exp(-2pi * 1im * (alphas[k] - 1) * t), 
				t_grid)
	unscaled_int = sum_kbn(real(g_out))/m
	coef = r^(1-alphas[k]) * vbars[k]^(-alphas[k])
	unscaled_int * coef
end


function integral_fourier(vbars_, alphas_; m=100, usePosContour=true, rs=nothing)
	vbars, ps, alphas=sort_poles(vbars_, alphas_)
	if rs == nothing
		rs, lams = optRs_NumRes(vbars, ps, alphas)
	end
	residues = BigFloat[]
	for k = 1:length(ps)
		if (usePosContour && ps[k] > 0) || (!usePosContour && ps[k] < 0)
			push!(residues, real(numResFourier(vbars, ps, alphas, k, rs[k], m)))
		end
	end
	usePosContour ? -sum_kbn(residues) : 1 + sum_kbn(residues)
end

function prob_gauss_approx(vbars, alphas)
	const mu = dot(vbars, alphas)
	const sig = sqrt(dot(alphas, vbars.^2))
	cdf(Normal(), -mu/sig)
end


function probNeg_direct(vbars, alphas)
	#handle some degnerate cases gracefully
	if minimum(vbars) > 0
		return 0.
	elseif maximum(vbars) <= 0
		return 1.
	end

	# const a = calc_shift(vbars, alphas)
	# coef = prod((1 + a * vbars) .^(-alphas))/pi
	# vhat = vbars ./ (1 + a*vbars)
	# f(s) = coef * prod((1+ s*1im*vhat).^(-alphas))/(a + s*1im)
	# real(quadgk(f, 0, Inf)[1])

	const a = calc_shift(vbars, alphas)
	# coef = prod((1 + a * vbars) .^(-alphas))/pi
	# vhat = vbars ./ (1 + a*vbars)
	# f(s) = coef * prod((1+ s*1im*vhat).^(-alphas))/(a + s*1im)
	real(quadgk(s->fint(a + s * 1im, vbars, alphas)/pi, 0, Inf)[1])

end

function probXNeg(vbars_, alphas_; method = :PosContour, m=100)
	vbars, ps, alphas = sort_poles(vbars_, alphas_)
	rs, lams          = optRs(vbars, alphas)
	f(z)              = fint(z, vbars, alphas)
	ps_rs_ls = zip(ps, rs, lams)

	if method == :PosContour || method == :NegContour
		if method == :PosContour
			ps_rs_ls = filter(p_r_l -> p_r_l[1] > 0, ps_rs_ls)
		elseif method ==:NegContour
			ps_rs_ls = filter(p_r_l -> p_r_l[1] < 0, ps_rs_ls)
		else
			error()
		end
		residues = zeros(length(ps_rs_ls))
		for (ix, (p, r, l)) in enumerate(ps_rs_ls)
			residues[ix] = real(numRes(z->l * f(z), p, r, m, old=false)[1])
		end

		inv_ls = [1/l for (p, r, l) in ps_rs_ls]

		if method == :PosContour
			return -sum_kbn(inv_ls .* residues)
		elseif method ==:NegContour
			return 1 + sum_kbn(inv_ls .*residues)
		else
			error()
		end
	elseif method == :Shift
		if maximum(ps) < 0 
			#VG Cludge
			kappa = 1
		else
			kappa = .5 * minimum(ps[ ps .> 0])
		end
		println(kappa)
		I,E = quadgk(f, kappa - 10im, kappa + 10im)  #VG more principled solution needed
		println("Error \t $E")
		return I / (2im * pi)
	else
		error("Method not recognzied")
	end
end

#evaluates Dir tail Prob for 3dimensions only
function Prob3d(alphas, vs, t)
	@assert length(alphas) == 3 "Only defined for 3 dimensions"
	if t >= norm(vs, 1)
		return 0
	elseif t <= -norm(vs, 1)
		return 1
	end

	##VG a Cludge
	if abs(vs[2] - vs[3]) < 1e-10
		vs[2] = vs[2] + 1e-8
	end

	beta1 = Beta(alphas[1], alphas[2] + alphas[3])
	beta2 = Beta(alphas[2], alphas[3])
	#the integrand
	function f(p1)
		val = (t - vs[1]*p1)/(1-p1)/(vs[2]-vs[3]) - vs[3]/(vs[2]-vs[3])
	 	cdf(beta2, val) * pdf(beta1, p1)
	end

	I, E = Base.quadgk(f, 1e-10, 1-1e-10)
	if vs[2] > vs[3]
	 	I = 1 - I
	end
	 I
end



end #ends module

