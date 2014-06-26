import math
from numpy import*
from scipy import*
from scipy.integrate import quad
import pylab as plt
import cosmolopy.distance as cd
# to find the interal of 1/H(z) with respect to z when w(z) = -1

def convert_R_to_Beta(R1,beta):

	mu=0.0000001
	dmu=1e-1
	sigma_b = zeros(len(R1))
	test_sigma =  zeros(len(R1))
	#print R1/100
	for i in range (len(beta)):
		while (mu<=1):
        		sigma_b = (R1/100) *((1.+ beta[i]* mu**2)**2)/(2*beta[i]*mu**2)
        		test_sigma += (1.+ beta[i]**2 * mu**2)/ (2*beta[i]*mu**2)
        		mu = mu + dmu
        		#print sigma_b[i]
        #sigma_b += sigma_b 
        	print sigma_b[i]/beta[i], sigma_b[i], beta[i]
        	
	return  (sigma_b)	


def Hz_func(x):
    # print 'H0', 0, h0kms * (math.sqrt(omegam*(1.0 + 0)**3 + omegal))
     return h0kms*(math.sqrt(omegam*(1.0 + x)**3 + omegal))

def integrand(x):
	return dhmpc/(math.sqrt(omegam*(1.0 + x)**3 + omegal))
def inth(x):
	return quad(integrand, 0, x )[0]

def Dv(z):
	'''__________ Dv_________________________'''
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 1.0 - 0.24, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	Dv = cd.diff_comoving_volume(z, **cosmo)
	return Dv
def err_Dv(z):
	'''___________DA__________________________'''
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 0.76, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)
	'''___________Hz___________________________'''
	H_z = cd.hubble_distance_z(z, **cosmo)
	'''________________The error on Dv___________________'''
	part1 = ( dasigma/ d_a ) **2 
	part2 = (Hsigma/ H_z)**2
	part3 = 0.0 #(cov_DaH/ (d_a* H_z))
	sigma_Dv = sqrt(Dv(Z)**2 * (part1 + part2 + part3 ))
	return  sigma_Dv
def DA_cosmo(z):
	'''___________DA__________________________'''
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 0.76, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)
	return d_a
	
def Hz_cosmo(z):
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 0.76, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	'''___________Hz___________________________'''
	H_z = cd.hubble_distance_z(z, **cosmo)
	#print z,  H_z
	#print 0,  cd.hubble_distance_z(0, **cosmo)
	return H_z
def wz(x):
	w0 =-0.82
	wa = 0.58
	w = w0 + wa*x/(1.+x)
	return (1. + w)/1. + x

if __name__ == "__main__":
    (z1, Err_lnda, Err_lnH_1,R1,B1) = loadtxt('output_1mJy_diff_14bins.txt', unpack=True)
    beta = B1
    print len(beta)
    zmin,errDa, errH,errR,B= loadtxt("output_1mJy_diff_14bins.txt",unpack='True')    
    err_beta = convert_R_to_Beta(errR,B)
    n = len(zmin)  # number of data points 
    ckms = 3.0e5
    h0kms = 73.0
    omegam = 0.24
    omegal = 1. - omegam
    #Hubble distance in Mpc
    dhmpc = ckms / h0kms
    meansigma = 0.1 # mean error of the data 
    f = open("DA_mockdatasigmainth.dat" ,"w")
    q = 0
    #for i in range(n):
    Z = zmin #random.uniform(zmin[i],zmax[i],zbin[i]) 
    q = q + len(Z)
    dasigma = zeros(len(Z))
    Hsigma= zeros(len(Z))
    beta_sigma= zeros(len(Z))
    Y = zeros(len(Z)) 
    Hz = zeros(len(Z))
    betaz = zeros(len(Z))
    f0= zeros(len(Z))
    for j in range (len(Z)): 
                dasigma[j] = errDa[j]/100.0*DA_cosmo(Z[j])
 		Hsigma[j] = errH[j]/100.0*Hz_func(Z[j])
 		#betaz[j] = beta[j]
 		#beta_sigma[j] = B1[j]*beta
         	Y[j] =(DA_cosmo(Z[j]) + dasigma[j]*random.normal())
		Hz[j]= Hz_func(Z[j]) + Hsigma[j] * random.normal()
		#f0[j] = bias[j]*beta[j] + beta_sigma[j]/100* random.normal()
         	f.write(str(Z[j]) +'\t'+ str(Y[j]) +'\t'+ str(abs(dasigma[j])) +'\n')
    from pylab import *
    p1= errorbar(Z,Y,yerr = dasigma,color='red',ecolor='black')     
    xlabel(r"redshift ($z$)",fontsize=15)
    ylabel(r"$D_A(z) $  $\rm{Mpc}$",fontsize=20) 
    legend([p1,] ,[' $1\mu$Jy with w = -1'], loc='best')
    savefig('DA_errs_0mJy_mario.eps')
    show()	
    
    
    
    p2= errorbar(Z,Hz,yerr = Hsigma,color= 'red',ecolor= 'black')     
    xlabel(r"redshift ($z$)",fontsize=15)
    ylabel(r"$H(z)$ ${\rm Mpc}^{-1}$ ${\rm  km } s^{-1}$",fontsize=20) 
    legend([p2,] ,[' $1\mu$Jy with w =-1'], loc='best')
    savefig('Hz_errs_0mJy_mario.eps') 
    show()
