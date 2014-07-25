from math import *from numpy import *from scipy import *import cosmolopy.constants as ccimport cosmolopy.distance as cdimport cosmolopy.perturbation as cpimport matplotlib.pyplot as pltimport matplotlib.ticker as ticimport matplotlib.axis from scipy.integrate import quadfrom scipy import specialfrom scipy import interpolatedef dp_ov_p(k,kmin, kmax,V_sur, pk,ngal):	dlnk=log(kmax) - log(kmin)	factor=(k)**3*dlnk/(2.**3 *3.1415926535**2) # Mpc^-3 h^3	wkmu=V_sur*(ngal*pk/(1.+ngal*pk))**2 # Mpc^3 h^-3	fis_7 = factor * wkmu	dp_7= (sqrt((fis_7)**-1))	return k, dp_7def dvdz(z):	# Planck best-fit parameters	cosmo = {'omega_M_0':        0.316,			 'omega_lambda_0':   0.684,    			'omega_b_0':        0.049,    			'N_eff':            3.046,   			 'h':                0.67,   			 'ns':               0.962,   			 'sigma_8':          0.834,    			'gamma':            0.55,   			 'w0':               -1.,    			'wa':               0.,   			 'sigma_nl':         7.}	cosmo = cd.set_omega_k_0(cosmo)	Vc = cd.diff_comoving_volume(z, **cosmo)	return  Vcdef da(z):	# Planck best-fit parameters	cosmo = {'omega_M_0':        0.316,			 'omega_lambda_0':   0.684,    			'omega_b_0':        0.049,    			'N_eff':            3.046,   			 'h':                0.67,   			 'ns':               0.962,   			 'sigma_8':          0.834,    			'gamma':            0.55,   			 'w0':               -1.,    			'wa':               0.,   			 'sigma_nl':         7.}	cosmo = cd.set_omega_k_0(cosmo)	d_a = cd.angular_diameter_distance(z, **cosmo)/(h)	print "Angular diameter distance = %.1f Mpc" % (d_a) 	return  d_adef V_sur(z, zmin, zmax,area):	vol = quad(dvdz, zmin, zmax)[0]		vol = area*(3.1415/180.)**2.*vol	return vol#===================Read the file==================================(k0_ov_h, k_7, pk) = loadtxt('cosmic_v_k_Euclid.txt', unpack=True)(kh, fis) = loadtxt('cosmic_v_k_Euclid_fis.txt', unpack=True)(k_Euclid,  dP_Euclid) = loadtxt('Euclid_dP_z1.dat', unpack=True)(kwmap, pkwmap)= loadtxt('pk_Euclid_z=30_wmap5.txt',unpack='True')(kh_p, pk_p)= loadtxt('pk_Euclid.txt',unpack='True')(kh_p2, pk_p2)= loadtxt('cosmic_v_k_Euclid_fis_test.txt',unpack='True')(kh_p3, pk_p3)= loadtxt('cosmic_v_k_Euclid_fis_testwmap.txt',unpack='True')#============= inputs ===========================================print 'bin size Phil', len(k_Euclid)print 'bin size me', len(kh_p)print 'width me', 2.0049310292162170E-002print 'test binning', (max(log(k_7))- max(log(k0_ov_h)))/0.02h = 0.67zmax=1.050000000000000044e+00zmin = 9.499999999999999556e-01z = (zmin+ zmax)/2.0area = 15000.0# dvdz_test = dvdz(1.)V_survey2= V_sur(z, zmin,zmax, area) *h**-3#7753233188.7160435 # output from volume.f90ngal = 1.7378236097892251E-003 # Mpc^-3 h^3V_survey3 = 25804351542.240719#29718387630.625401da_1= da(z)#==========convert Unit =======================k_Euclid = k_Euclid*h # convert units from Mpc^-1 to Mpc^-1 h V_survey = 13.38 * 10**9*(h**-3) # convert V_sur phil  to Mpc^3 h^-3#print 'V_sur ratio % = ', (V_survey2/V_survey)*100print'%.1e' % V_survey , '%.1e' % V_survey2,  '%.1e' % V_survey3#============interpolate ============================k0_new =6.8492999999999996E-005 k_new =9.5589999999999998E-005print'width Phil', log(k0_new)- log(k_new)f = interpolate.interp1d(k_7, pk)pk_new = f(k_Euclid)#========call dp_ov_p=========================k0_7= 6.8492999999999996E-005 km_7=  9.5589999999999998E-005 (k_7new ,dp_7new) = dp_ov_p(k_Euclid,k0_new, k_new,V_survey3, pk_new,ngal)(k_7, dp_7) =  dp_ov_p(k_7,k0_7,km_7 ,V_survey3, pk,ngal)#=========plotting =================================fig = plt.figure()ax = fig.add_subplot(111, yscale='log')ax.set_xscale('log')ax.plot(k_7,dp_7, color= 'green', linewidth=2.5, linestyle="--", label ='Euclid - Plank parameters')ax.plot(kh,fis, color= 'blue', linewidth=2.5, linestyle="--", label ='Euclid - WMAP5 parameters')ax.plot(k_7new,dp_7new, color='orange', linewidth=3.5, linestyle="--", label ='Euclid_binned')ax.plot(kh_p2,pk_p2, color= 'green', linewidth=2.5, linestyle="-", label ='Plank2test ')ax.plot(kh_p3,pk_p3, color= 'blue', linewidth=2.5, linestyle="-", label ='WMAP2test ')ax.plot(k_Euclid, dP_Euclid, color= 'red',  linewidth=2.5, linestyle="--", label ='Euclid_Phil')#====================== Labels ==========================ax.set_xlabel(r"$k \ [{\rm Mpc}^{-1} h ]$", fontsize=15)ax.set_ylabel(r"$\delta P/P$ ", fontsize= 17)plt.ylim(0.001, 10**4)plt.xlim(0.0001, 1.0)ax.legend(loc='upper right')plt.savefig('/Users/sahba/Dropbox/SKA Survey/Ska Survey follow up/Notes/Error on the power spectrum/deltaP_ov_p_Euclid.eps')#plt.show()fig2 = plt.figure()ax2 = fig2.add_subplot(111, yscale='log')ax2.set_xscale('log')ax2.plot(kh_p,pk_p, color= 'green', linewidth=2.5, linestyle="-", label ='Plank ')ax2.plot(kwmap,pkwmap, color= 'blue', linewidth=2.5, linestyle="-", label ='WMAP 5 ')ax2.set_ylabel(r"$P_g(k) \ [{\rm Mpc}^{3} h^{-3} ]$", fontsize=15)ax2.set_xlabel(r"$k \ [{\rm Mpc}^{-1} h ]$", fontsize=15)ax2.legend(loc='upper right')plt.savefig('/Users/sahba/Dropbox/SKA Survey/Ska Survey follow up/Notes/Error on the power spectrum/pk.eps')plt.show()print '============== The program executed successfully!=============================='