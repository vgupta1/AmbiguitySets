


#Computes the Sigma Star limit for Dirichlet
function sigStar(pstar)
	const d = length(pstar)
	sigma = zeros(d, d)
	for i = 1:d
		sigma[i, i] = pstar[i]*(1-pstar[i])
		for j = i+1:d
			sigma[i, j] = -pstar[i] * pstar[j]
			sigma[j, i] = sigma[i, j]
		end
	end
	sigma
end
