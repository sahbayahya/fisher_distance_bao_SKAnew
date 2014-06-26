from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
from scipy.integrate import quad
from scipy import special


#(k, p01) = loadtxt('Pk.txt', unpack=True)
(z_3, p01_3, ngal_3) = loadtxt('Pk_3.txt', unpack=True)
(z_7, p01_7, ngal_7) = loadtxt('Pk_77.txt', unpack=True)
(z_23, p01_23, ngal_23) = loadtxt('Pk_23.txt', unpack=True)
(z_100, p01_100, ngal_100) = loadtxt('Pk_100.txt', unpack=True)
#(k3, p3) = loadtxt('Pk_wmap30.txt', unpack=True)
#(k_7, p01_7) = loadtxt('Pk_linear.txt', unpack=True)
#(k_wmap, p30) = loadtxt('wmap5baosn_max_likelihood_matterpower_at_z=30.dat', unpack=True)
#(k_wmap_1, p) = loadtxt('wmap5baosn_max_likelihood_matterpower.dat', unpack=True)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#yticks = [0.1, 1, 10]
#ax.set_yscale('log')
#ax.set_title('P(0.1) $h$Mpc$^{-1}$')
#ax.plot(k, p01, color= 'red', label ='100 $\mu$Jy')
#ax.plot(k_7, p01_7, color= 'black', label =' 7 $\mu$Jy')
#ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=15)
#ax.set_ylabel(r"$P(0.1)$ ", fontsize= 15)
#ax.get_legend_handles_labels()
#ax.legend(loc='upper left')
#plt.savefig('Pk01.eps')
#plt.show()
#ay.set_yscale('log')
ax.set_title('Using WMAP 3 cosmology')
ax.set_yscale('log')
#ax.set_xscale('log')
ax.plot(z_3, p01_3, color= 'orange',  linestyle="--", linewidth=2.5,label ='3 $\mu$Jy ')
ax.plot(z_3, ngal_3, color= 'orange',  linestyle="-", linewidth=2.5)
ax.plot(z_7, p01_7, color= 'green',  linestyle="--", linewidth=2.5,label =' 7 $\mu$Jy ')
ax.plot(z_7, ngal_7, color= 'green',  linestyle="-", linewidth=2.5)
ax.plot(z_23, p01_23, color= 'black',  linestyle="--", linewidth=2.5,label ='23 $\mu$Jy ')
ax.plot(z_23, ngal_23, color= 'black',  linestyle="-", linewidth=2.5)
ax.plot(z_100, ngal_100, color= 'magenta', linestyle="-", linewidth=2.5)
ax.plot(z_100, p01_100, color= 'magenta', linestyle="--", linewidth=2.5,  label ='100 $\mu$Jy ')
#ax.plot(k_wmap_1, p, color= 'blue', label ='P$_{lin}$(k) of WMAP3')
#ax.plot(kt, p01t, color= 'cyan', label ='P(k) of 7 $\mu$Jy Using Eq.24')
#ax.plot(k3, p3,color= 'red', label='P(k) of 7 $\mu$Jy with WMAP3  at z=30')
#ax.plot(k_7, p01_7, color= 'green', label ='P$_{lin}$(k) of WMAP3 at z= 30')
#ax.plot(k_wmap, p30, color= 'blue', label ='P(k), wmap at z=30')
#ax.plot(k, p01, color= 'red', label ='1/n, 100 $\mu$Jy')
#ax.plot(k_7, p01_7, color= 'black', label =' 1/n, 7 $\mu$Jy')
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=15)
ax.set_ylabel(r"$1/n(z)$  [$h^{3} {\rm Mpc}^{-3}$] ", fontsize= 15)
ax.get_legend_handles_labels()
ax.set_ylim(100, 10**6)
ax.set_xlim(0.2, 2)
#ax.legend(loc='upper right')
#plt.legend([p1, p2] ,['$ 7 \mu$Jy',' $ 100 \mu$Jy'], loc='best')
plt.savefig('Pk_ngal.eps')
plt.show()






