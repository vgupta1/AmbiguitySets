### Series Acceleration


####
# Adapted from Cohen, et. al, 
#"Convergence Acceleration Alternating Series"
#####
#seq = (-1)^k a_k  from k = 0 to ...
#computes sum(a_k)
function cohen_alt(seq)
	const n = length(seq)
	d = (3 + sqrt(8))^n ; d = (d + 1/d)/2
	b = -1.; c = -d; s = 0.; sgn = 1
	for k = 0:n-1
		c = b - c;
		s = s + c * sgn * seq[k + 1]  #julia indexes form 1
		b = (k+n)*(k-n)*b / (k+.5)/(k+1)
		sgn = sgn * -1
	end
	s/d
end



