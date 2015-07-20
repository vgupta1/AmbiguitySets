# Asymptotic Plots for Amb Sets Paper


using Distributions
include("../src/DirRes.jl")

function table_in_d(eps_=.1)
	#creates the table in d
	d_grid = [3, 5:5:25, 50:25:100]
	
	println("d \t KL \t KL_Cov \t Chisq \t Chisq_Cov \t Cov")
	for d in d_grid
		print(d, "\t", round(Dir.kl_const(eps_), 2), "\t")
		print(round(Dir.kl_cov_const(eps_, d), 2), "\t")
		print(round(Dir.chisq_const(eps_), 2), "\t")
		print(round(Dir.chisq_cov_const(eps_, d), 2), "\n")
	end
end

function table_cov_d_eps()
	eps_grid = [.1, .05, .01]
	d_grid = [3, 5:5:25, 50:25:100]
	print("d/eps \t")
	map(s->print("$s \t"), eps_grid)
	print("\n")
	for d in d_grid
		print(d, "\t")
		for eps_ in eps_grid
			print(round(Dir.cov_const(eps_, d), 2), "\t")
		end
		print("\n")
	end
end


#dumps a file with several random vs
function plots_in_N(dir_path, path)
	d = 15
	eps_ = .1
	N_grid = map(ifloor, logspace(1, 4, 20)) 
	println(N_grid)

	vs = randn(15, 10)
	srand(8675309)	
	
	#Simulate some random vs
	writecsv(dir_path, vs)

	fp = open(path, "w")
	writecsv(fp, ["N" "Type" "Direction" "Ratio"])

	#add the asymptotic values
	writecsv(fp, [0 "KL" 0 Dir.kl_const(eps_)])
	writecsv(fp, [0 "ChiSq" 0 Dir.chisq_const(eps_)])
	writecsv(fp, [0 "KL_C" 0 Dir.kl_cov_const(eps_, d)])
	writecsv(fp, [0 "ChiSq_C" 0 Dir.chisq_cov_const(eps_, d)])

	#simulate a long data-stream to reuse
	dat = rand(DiscreteUniform(1, d), N_grid[end])

	for N in N_grid
	    alphas = Base.hist(dat[1:N], d)[2]
	    alphas += ones(d)

	    for ix = 1:size(vs, 2)
	    	writecsv(fp, [N "KL" ix Dir.kl_ratio(vs[:, ix], alphas, eps_)])
	    	writecsv(fp, [N "ChiSq" ix Dir.chisq_ratio(vs[:, ix], alphas, eps_)])
	    	writecsv(fp, [N "KL_C" ix Dir.kl_cov_ratio(vs[:, ix], alphas, eps_, d, N)])
	    	writecsv(fp, [N "ChiSq_C" ix Dir.chisq_cov_ratio(vs[:, ix], alphas, eps_, d, N)])
	    end
	end
	close(fp)
end

function plot_in_eps(path)
	eps_grid = linspace(.01, .4, 20)
	kl_vals = map(Dir.kl_const, eps_grid)
	mom_vals = map(Dir.chisq_const, eps_grid)

	#not a clever way to do this
	kl_cov_vals5 = map(eps_ -> Dir.kl_cov_const(eps_, 5), eps_grid)
	kl_cov_vals10 = map(eps_ -> Dir.kl_cov_const(eps_, 10), eps_grid)
	kl_cov_vals20 = map(eps_ -> Dir.kl_cov_const(eps_, 20), eps_grid)
	fp = open(path, "w")
	writecsv(fp, ["eps" "KL" "ChiSq" "KL_C5" "KL_C10" "KL_C20"])
	writecsv(fp, [eps_grid kl_vals mom_vals kl_cov_vals5 kl_cov_vals10 kl_cov_vals20])
	close(fp)
end

# creates a plot
function comp_eps_kl(path)
	eps_grid = linspace(.001, .5, 30)
	fp = open(path, "w")
	writecsv(fp, ["d" "Cov_eps" "Real_eps"])
	for d in [3, 5, 10]
		kl_cov_vals = map(eps_-> sqrt(quantile(Distributions.Chisq(d-1), 1-eps_)/2), 
							eps_grid)
		real_eps_vals = 1- cdf(Distributions.Normal(), kl_cov_vals)
		writecsv(fp, [d*ones(length(eps_grid)) eps_grid real_eps_vals])
	end

	#add one more for the chisqs
	chisq_vals = map(eps_ -> sqrt(1/eps_ - 1), eps_grid)
	real_eps_vals = 1 - cdf(Distributions.Normal(), chisq_vals)
	writecsv(fp, [fill("Chisq", length(eps_grid)) eps_grid real_eps_vals])

	close(fp)
end

