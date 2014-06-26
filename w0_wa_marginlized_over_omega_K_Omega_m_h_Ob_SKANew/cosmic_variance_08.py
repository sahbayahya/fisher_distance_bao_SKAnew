from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import matplotlib.axis 
from scipy.integrate import quad
from scipy import special
#===================Read the file==================================
#z, dn/dz, cosmic varaince = (1/p02)*dV/dz
(x, dNofdz_rms_0muJy,dNofdz_rms_1muJy,dNofdz_rms_3muJy,dNofdz_rms_73muJy , dNofdz_rms_23muJy, dNofdz_rms_70muJy ,dNofdz_rms_100muJy,dNofdz_rms_200muJy) = loadtxt ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)
(z_3, p01_3, ngal_3) = loadtxt('z_vs_dndz_3_08.txt', unpack=True)
(z_7, p01_7, ngal_7) = loadtxt('z_vs_dndz_7_08.txt', unpack=True)
(z_100, p01_100, ngal_100) = loadtxt('z_vs_dndz_100_08.txt', unpack=True)
#=============== plot z vs dn/dz and cosmic varaince=================
#minorLocator = MultipleLocator(10)
#yticks = [1, 10, 20,30,40,50, 60, 70, 80 ,90,100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
yticks = [1, 10,100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
xticks = arange(min(x), max(x)+1, 0.2)
#====================================================================
fig = plt.figure()
ax = fig.add_subplot(111, yscale='log')
#=========Set log scale on the yaxis========================
#========= defult plotting =================================
ax.plot(x, dNofdz_rms_3muJy, color= 'orange', linewidth=2.5,linestyle="-", label ='$3\ \mu$Jy')
ax.plot(z_3, p01_3, color= 'orange', linestyle="--", linewidth=2.5)
ax.plot(z_7, p01_7, color= 'green', linestyle="--",linewidth=2.5)
ax.plot(x, dNofdz_rms_73muJy, color= 'green', linewidth=2.5, linestyle="-", label ='  $7\ \mu$Jy')
ax.plot(x, dNofdz_rms_100muJy, color= 'magenta', linewidth=2.5,linestyle="-", label ='$100\ \mu$Jy ')
ax.plot(z_100, p01_100, color= 'magenta', linestyle="--", linewidth=2.5)
#====================== Labels ==========================
#plt.annotate('', xy=(1.39, 1992), xycoords='data', xytext=(+10, +30), textcoords='offset points', fontsize=16, arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.1"))
plt.annotate(r'${\rm shot \ noise} =   \ {\rm cosmic\ variance}$', xy=(0.8, 260), xycoords='data', xytext=(+50, +100), textcoords='offset points', fontsize=16 , arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.3"))
#plt.annotate(r'$ \frac{1}{P_{0.2}}\frac{dV}{dz} \ \frac{1}{[{\rm deg}^{2}]}$', xy=(1.39, 1992), xycoords='data', xytext=(+10, +30), textcoords='offset points', fontsize=16 , arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=15)
ax.set_ylabel(r"$\frac{d{\rm N}}{dz}(z) \ [ {\rm deg}^{-2} \ {\rm per} \ {\rm unit} \ z ]$ ", fontsize= 17)
#ax.set_ylabel(r"$\frac{d{\rm N}}{dz}(z) \ \frac{1}{[{\rm deg}^{2}]}$ ", fontsize= 17)
#==============Set Tickers =============================
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=5)
plt.tick_params(which='minor', length=2, color='g')
#======================================================
#================= yaxis
#plt.yticks(yticks,[r'$1$', r'$10$', '','','','','','','','', r'$10^2$', r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$', r'$10^7$' ])
plt.yticks(yticks,[r'$1$', r'$10$',r'$10^2$', r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$', r'$10^7$' ])
#ax.set_yscale('log')
#==================xaxis
plt.xticks(xticks,[r'$0$', r'$0.2$',r'$0.4$', r'$0.6$',r'$0.8$',r'$1.0$',r'$1.2$', r'$1.4$',r'$1.6$',r'$1.8$', r'$2.0$' ])
plt.ylim(1, 10**7)
plt.xlim(0.0, 2.0)
ax.legend(loc='upper right')
plt.savefig('cosmic_Variance_08.eps')
plt.show()

print '============== The program executed successfully!=============================='



