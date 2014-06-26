*** Fisher matrix calculation using Baryon Acoustic Oscillations only ***
September 17, 2008: E.Komatsu
Ref: Seo & Eisenstein, ApJ, 665, 14 (2007)

Here we provide a program for computing the Fisher matrix
of the angular diameter distance, Da, and the Hubble expansion rate, H,
using the information of Baryon Acoustic Oscillations (BAOs) ONLY.
I.e., the information in the overall shape is not used at all.

Users may choose to use any of the following elements:

- Arbitrary number of redshift bins
- Time-dependent linear galaxy bias, the maximum wavenumber used in 
the analysis (kmax), and the redshift accuracy in units of km/s

For the computation with the linear spectrum, we provide the data for 
the linear power spectrum at z=30, 
"wmap5baosn_max_likelihood_matterpower_at_z=30.dat," which was
generated using CAMB code for the maximum likelihood parameters
given in Table I of Komatsu et al.(2008) [WMAP 5-year interpretation paper] 
with "WMAP5+BAO+SN". The input file for CAMB is also provided 
(wmap5baosn_max_likelihood_params.ini). 

The program reports on the errors in ln(Da) and ln(H) in percent
as well as the correlation coefficient between them.

- To compile and use the program, edit Makefile
and simply "./make"
- It will generate executables called "fisher_distance_bao"

For convenience, we provide a sample script, "sample.cmd", which would generate
the output like this:

 For the area, the redshift range and the number of galaxies:
 1) Read in "number.txt"
 2) Read in your favorite file
 0) Enter your choice of area (in deg^2), zmin, zmax, # of galaxies in millions
    (A new file called "number.txt" will be created.)
 read in wmap5baosn_max_likelihood_matterpower_at_z=30.dat
 Enter area (in deg^2), zmin, zmax, # of galaxies (in millions)
 in the redshift bin#=( 1 )
 (Enter 0,0,0,0 when you are done.)
 Enter the values of
 bias, kmax [h Mpc^-1], and redshift error [km/s] in 
 1.500 < z < 3.500
 === redshift bin#=( 1 ) ===
 1.500 < z < 3.500
 Vsur = 3.48185 h^-3 Gpc^3
 Ngal = 1.00000 millions
 ngal = 0.28720 10^-3 h^3 Mpc^-3
 bias = 2.00000
 beta = 0.48413
 f    = 0.96826
 sigz = 0.00000 h^-1 Mpc
 spara= 6.26203 h^-1 Mpc
 sperp= 3.18151 h^-1 Mpc
 kmax = 0.43000 h Mpc^-1
 Err[lnDa](%) =  1.62498
 Err[lnH](%)  =  2.48605
 r(lnDa,lnH)  =  0.41611
 Err[lnR](%)  =  1.05243
 
 Enter area (in deg^2), zmin, zmax, # of galaxies (in millions)
 in the redshift bin#=( 2 )
 (Enter 0,0,0,0 when you are done.)
 Enter the values of
 bias, kmax [h Mpc^-1], and redshift error [km/s] in 
 3.500 < z < 6.500
 === redshift bin#=( 2 ) ===
 3.500 < z < 6.500
 Vsur =11.01239 h^-3 Gpc^3
 Ngal =10.00000 millions
 ngal = 0.90807 10^-3 h^3 Mpc^-3
 bias = 4.00000
 beta = 0.24837
 f    = 0.99347
 sigz = 1.38787 h^-1 Mpc
 spara= 3.73173 h^-1 Mpc
 sperp= 1.87198 h^-1 Mpc
 kmax = 1.00000 h Mpc^-1
 Err[lnDa](%) =  0.41523
 Err[lnH](%)  =  0.70073
 r(lnDa,lnH)  =  0.41007
 Err[lnR](%)  =  0.27941
 
 Enter area (in deg^2), zmin, zmax, # of galaxies (in millions)
 in the redshift bin#=( 3 )
 (Enter 0,0,0,0 when you are done.)
 === combined ===
 Vsur =14.49424 h^-3 Gpc^3
 Ngal =11.00000 millions
 ngal = 0.75892 10^-3 h^3 Mpc^-3
 Err[lnDa](%) =  0.40228
 Err[lnH](%)  =  0.67440
 r(lnDa,lnH)  =  0.41035
 Err[lnR](%)  =  0.27005
 
 *** Error on w, ignoring Omega_k and Omega_m ***
 Err[w]            =  0.02057
 *** Error on w, with Omega_k and Omega_m marginalized over ***
 Err[w]            =  0.06380
 Err[Omega_k]      =  0.03000
 Err[lnOmega_m](%) =  2.49328
 
 << CMB Priors Added (Edit da_accuracy & om_accuracy in "add_cmb") >>
 Err[lnDa(z=1090)](%) =  0.20000
 Err[lnOmega_m](%)    =  2.00000
 *** Error on w, ignoring Omega_k and Omega_m ***
 Err[w]            =  0.01465
 *** Error on w, with Omega_k and Omega_m marginalized over ***
 Err[w]            =  0.04251
 Err[Omega_k]      =  0.00218
 Err[lnOmega_m](%) =  1.24709
