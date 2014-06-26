./fisher_distance << eof
4       # non-linear P(k) in red space
3	# marginalize over 3 parameters
1 2 3	# marginalize over ln(amp), beta, and FoG
2	# use my favorite file
number_z=1.9-3.5_two_bins.txt # filename
2.0 0.3 200 # bias=2.0,kmax=0.3hMpc^-1,delta_v=200km/s for the first bin
2.5 0.3 200 # bias=2.5,kmax=0.3hMpc^-1,delta_v=200km/s for the second bin
eof
