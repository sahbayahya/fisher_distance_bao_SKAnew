./fisher_distance_bao << eof
0          # enter redshift range and # of objects as we go
400. 1.5 3.5 1  # 400deg^2 zmin=1.5 zmax=3.5 Ngal=1 million for the first bin
2 0.43 0   # bias=2,kmax=0.43hMpc^-1,delta_v=0 km/s for the first bin
1000. 3.5 6.5 10 # 1000deg^2 zmin=3.5 zmax=6.5 Ngal=10 million for the second bin
4 1 180    # bias=4,kmax=1hMpc^-1,delta_v=180 km/s for the second bin
0 0 0 0    # finish entering data
eof
