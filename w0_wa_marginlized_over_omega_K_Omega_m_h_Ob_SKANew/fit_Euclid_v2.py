import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
from scipy.integrate import quad
#=============Functions ====================
def dvdz(z):
	# Planck best-fit parameters
	cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
	cosmo = cd.set_omega_k_0(cosmo)
	Vc = cd.diff_comoving_volume(z, **cosmo)
	return  Vc
def da(z):
	# Planck best-fit parameters
	cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)/(h)
	print "Angular diameter distance = %.1f Mpc)" % (d_a) 
	return  d_a
def V_sur(z, zmin, zmax,area):
	vol = quad(dvdz, zmin, zmax)[0]	
	vol = area*(3.1415/180.)**2.*vol
	return vol
def func(p,x):
   w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x) 
   print w.size
   return w
def func1(p,x):
   w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x**1.5) 
   print w.size
   return w

def residuals(p,x,y):
   w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
def D(z):
     return cp.fgrowth(z,omega_M_0, unnormed=False)
def ErrorZ(z):
	return  0.001*(1+z)
def Bias(z):
	return  sqrt(1+z)
#def func(p,x):
 #  w=p[0]*np.exp(-p[1]*x)
  # print w.size
   #return w

#def residuals(p,x,y):
 #  w=p[0]* np.exp(-p[1]*x)
  # err=w-y
  # err=err**2
  # B=sum(err)
   #return B
#================================================
xrange4 = np.linspace(0.6, 2.0, 100)
omega_M_0 = 0.27
#xrange0 = np.array([ 0.8,  0.9 ,  1.0, 1.1 , 1.2 ,1.3  , 1.4 ,1.5 , 1.6 ,1.7 , 1.8 ,1.9, 2.0])
xrange = np.array([ 0.7, 0.8,  0.9 ,  1.0, 1.1 , 1.2 ,1.3  , 1.4 ,1.5 , 1.6 ,1.7 , 1.8 ,1.9, 2.0])
#dndzrange = np.array([1.75, 2.68,  2.56 ,  2.35, 2.12 , 1.88 ,1.68  , 1.40 ,1.12 , 0.81 ,0.53 , 0.49 ,0.29, 0.16])
dndzrange_ref = np.array([ 1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11])
print 'len z', len(xrange)
print 'len dn/dz', len(dndzrange_ref)
kmax = np.linspace(0.16004, 0.2, 14)#np.array([0.16004, 0.165, 0.1700,0.1750, 0.1800, 0.1850, 0.1900, 0.1950, 0.20, 0.20, 0.20, 0.20, 0.2, 0.2])
kmin = np.linspace(0.00435,0.00334, 14)
xmin =  xrange -0.05
xmax  =xrange+0.05
#print xmin
#print xmax
#========== Fit dndz==================================
#p0=[5.52,  0.6, 4.6]
#p04=[5.74, 1.14, 3.95]
#plsqtot= opt.fmin(residuals, p0, args=(xrange0, dndzrange_ref), maxiter=10000, maxfun=10000)
#print ' |   c1  | ',       '|         c2  |',        '|         c3  |'
#print  '0 muJy ',plsqtot[0],   plsqtot[1], plsqtot[2]
#y0=func(plsqtot,xrange)
#========V survey =================================
h= 0.67
area = 15000.0
Vsurvey= []
for i in range(len(xrange)):
          Vsurvey.append(V_sur(xrange[i], xmin[i],xmax[i], area)*(h**-3))
#========================== Save================

#kmax  =  np.empty(len(xrange)); kmax.fill(0.2)
area =np.empty(len(xrange)); area.fill(15000)
ErrorZ= 0.001*(1. + xrange)
data= concatenate((reshape(area,(len(xrange),1)),reshape( xrange,(len(xrange),1)) , reshape(dndzrange_ref*10**(-3),(len(xrange),1)),reshape(Bias(xrange),(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(kmin,(len(xrange),1)), reshape(ErrorZ,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)
#datafit= concatenate((reshape(xrange4,(len(xrange4),1)),reshape(y0*10**(-3),(len(xrange4),1))),axis=1)
print 'here is clear'
#savetxt('dndz_Euclid_ref_2.txt' , datafit)
#data4= concatenate((reshape(volume,(len(xrange4),1)),reshape(xmin,(len(xrange4),1)),reshape(xmax,(len(xrange4),1)),reshape(y0*10**(-3),(len(xrange4),1)),reshape(Bias(xrange4),(len(xrange4),1)),reshape(kmax,(len(xrange4),1)),reshape(ErrorZ(xrange4),(len(xrange4),1))),axis=1)
#savetxt('number_EuclidmJy_ref.txt' , data4)
savetxt('number_EuclidmJy_ref.txt' , data)
# PLOTING
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#ax.set_yscale('log')
#p0,  = ax.plot(xrange, y0, 'ro')
#plt.xlim(0.,2.5,0.5)
#plt.ylim(1, 2)
p1, = ax.plot(xrange, dndzrange_ref , 'bo')
#p2, =  plt.plot(x, bias_rms7_3 , 'ro')
#p3, =  plt.plot(x, bias_rms23 , 'go')
#p01,  = plt.plot(xrange, y0, color='#BA5F5F')
#p11, = plt.plot(xrange, Bias(xrange), color='#0404B4')
#p21, = plt.plot(xrange, y2, color='#B40431')
#p31, = plt.plot(xrange, y3, color='#04B404')
#plt.legend([p11] ,['BIAS'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('Euclid_dndz.eps')
plt.show()

#OUTPUT FOR dNdz fitting using Obreschkow function, for 0 muJy
#6.33839678537 2.17980796636 1.3923452325

#OUTPUT FOR dNdz fitting using Obreschkow function, for  1muJy
#7.08990441653 2.84040926362 5.21751576048


#OUTPUT FOR dNdz fitting using Obreschkow function, for 7.3 muJy
#5.3074425085 0.72760414229 4.11718714615
##OUTPUT FOR dNdz fitting using Obreschkow function, for  23 muJy

#4.9727473436 0.53260579938 6.66294780323

##OUTPUT FOR dNdz fitting using Obreschkow function, for 100  muJy
#8.28343815119 4.44979424928 23.9380642576

