from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib.ticker as tic
from scipy import special



(z, dn,dvdz)= loadtxt('no_galaxies_7point3_diff_14bins_S3.txt', unpack=True)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(z, dvdz,color='blue',linewidth=2.0, linestyle="--")
ax.set_yscale('log')
ax.set_ylabel(r'$dV/dz \  [h^{-3} \ Mpc^3 \  {\rm deg}^{-2}]$', fontsize= 20)
ax.set_xlabel('z')
plt.savefig('dVdz.eps')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(z, dn,color='red',linewidth=2., linestyle="--")
ax.set_yscale('log')
ax.set_ylabel(r'$n(z) \  [ h^{3} \ Mpc^{-3} ]$', fontsize= 20)
ax.set_xlabel('z')
plt.savefig('nz.eps')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(z, dn*dvdz/30000.0,color='green',linewidth=2., linestyle="--")
ax.set_yscale('log')
ax.set_ylabel(r'$dN/dz (z)\  [{\rm deg}^{-2}]$', fontsize= 20)
ax.set_xlabel('z')
plt.savefig('nzdvdz.eps')
plt.show()
