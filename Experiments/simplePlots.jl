# Asymptotic Plots for Amb Sets Paper


using Distributions
include("../DirRes.jl")


function table_in_d()
	#creates the table in d
	eps_ = .1
	d_grid = [3, 5:5:25, 50:25:100]
	
	println("d \t KL \t KLBT \t SM \t Chisq")
	for d in d_grid
		print(d, "\t", round(Dir.kl_const(eps_), 2), "\t")
		print(round(Dir.kl_chisq_const(eps_, d), 2), "\t")
		print(round(Dir.mom_const(eps_), 2), "\t")
		print(round(Dir.chisq_const(eps_, d), 2), "\n")
	end
end


#dumps a file with several random vs
function plots_in_N(dir_path, path)
	d = 15
	eps_ = .1
	N_grid = map(ifloor, logspace(1, 6, 20)) 
	srand(8675309)	
	
	#Simulate some random vs
	vs = randn(15, 5)
	writecsv(dir_path, vs)

	fp = open(path, "w")
	writecsv(fp, ["N" "Type" "Direction" "Ratio"])

	#add the asymptotic values
	writecsv(fp, [0 "KL" 0 Dir.kl_const(eps_)])
	writecsv(fp, [0 "ChiSq" 0 Dir.mom_const(eps_)])
	writecsv(fp, [0 "KL_C" 0 Dir.kl_chisq_const(eps_, d)])
	writecsv(fp, [0 "ChiSq_C" 0 Dir.chisq_const(eps_, d)])

	#simulate a long data-stream to reuse
	dat = rand(DiscreteUniform(1, d), N_grid[end])

	for N in N_grid
	    alphas = Base.hist(dat[1:N], d)[2]
	    alphas += ones(d)

	    for ix = 1:size(vs, 2)
	    	writecsv(fp, [N "KL" ix Dir.kl_ratio(vs[:, ix], alphas, eps_, N)])
	    	writecsv(fp, [N "ChiSq" ix Dir.mom_ratio(vs[:, ix], alphas, eps_)])
	    	writecsv(fp, [N "KL_C" ix Dir.kl_ratio_bt(vs[:, ix], alphas, eps_, d, N)])
	    	writecsv(fp, [N "ChiSq_C" ix Dir.mom_ratio_bt(vs[:, ix], alphas, eps_, d, N)])
	    end
	end
	close(fp)
end

function plot_in_eps(path)
	eps_grid = linspace(.01, .4, 20)
	kl_vals = map(Dir.kl_const, eps_grid)
	mom_vals = map(Dir.mom_const, eps_grid)
	fp = open(path, "w")
	writecsv(fp, ["eps" "KL" "ChiSq"])
	writecsv(fp, [eps_grid kl_vals mom_vals])
	close(fp)
end

# creates a plot
function comp_eps_kl(path)
	eps_grid = linspace(.001, .5, 30)
	fp = open(path, "w")
	writecsv(fp, ["d" "Cov_eps" "Real_eps"])
	for d in [3, 5, 10]
		kl_cov_vals = map(eps_-> quantile(Distributions.Chisq(d-1), 1-eps_)/2, 
							eps_grid)
		real_eps_vals = 1 ./exp(kl_cov_vals)
		writecsv(fp, [d*ones(length(eps_grid)) eps_grid real_eps_vals])
	end
	close(fp)
end

