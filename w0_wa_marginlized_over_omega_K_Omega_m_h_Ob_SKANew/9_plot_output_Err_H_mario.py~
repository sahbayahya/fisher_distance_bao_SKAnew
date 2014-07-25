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
#******************************************************

def D(z):
     h = 0.72
     return cd.angular_diameter_distance(z , **cosmo) * h


def H(z):
	return  cd.hubble_distance_z(z, **cosmo) 
     
     

(z0, Err_lnda0, Err_lnH0,R0,B0) = loadtxt('output_0NmJy_diff_14bins_S3.txt', unpack=True); z0 = z0[6:23] ;  Err_lnH0 = Err_lnH0[6:23]
(z1, Err_lnda1, Err_lnH1,R1,B1) = loadtxt('output_Euclid_diff_14bins_S3.txt', unpack=True);z1 = z1[0:13] ;  Err_lnH1 = Err_lnda1[0:13]
(z3, Err_lnda3, Err_lnH3, R3,B3) = loadtxt('output_3mJy_diff_14bins_S3.txt', unpack=True); z3 = z3[6:23] ;  Err_lnH3 = Err_lnH3[6:23]
(z7point3, Err_lnda7point3, Err_lnH7point3,R7point3,B7point3) = loadtxt('output_7point3mJy_diff_14bins_S3.txt', unpack=True); z7point3 = z7point3[6:20] ;  Err_lnH7point3 = Err_lnH7point3[6:20]
(z23, Err_lnda23, Err_lnH23,R23,B23) = loadtxt('output_23mJy_diff_14bins_S3.txt', unpack=True);z23 = z23[6:14] ;  Err_lnH23 = Err_lnH23[6:14]
(z70, Err_lnda70, Err_lnH70,R70,B70) = loadtxt('output_70mJy_diff_14bins_S3.txt', unpack=True) ; z70 = z70[4:9] ; Err_lnH70 = Err_lnH70[4:9]
(z100, Err_lnda100, Err_lnH100,R100,B100) = loadtxt('output_100mJy_diff_14bins_S3.txt', unpack=True); z100 = z100[4:8] ; Err_lnH100 =Err_lnH100[4:8] 
(z200, Err_lnda200, Err_lnH200,R200,B200) = loadtxt('output_200mJy_diff_14bins_S3.txt', unpack=True); z200 = z200[4:7] ;  Err_lnH200 =Err_lnH200[4:7]  
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')

fontsize = 20
ax.plot(z0, Err_lnH0,color='blue',linewidth=1.5, linestyle="--", label="$0 N \mu$Jy")
ax.plot(z3, Err_lnH3 , color = 'red',linewidth=1.5, linestyle="--", label="$ 3 \mu$Jy" )
ax.plot(z7point3, Err_lnH7point3 , color = 'green',linewidth=1.5, linestyle="--", label="$ 7.3 \mu$Jy" )
ax.plot(z23, Err_lnH23, color = 'black',linewidth=1.5, linestyle="--", label="$ 23 \mu$Jy")
ax.plot(z70, Err_lnH70, color = 'cyan',linewidth=1.5, linestyle="--", label="$ 70\mu$Jy")
ax.plot(z100, Err_lnH100, color = 'magenta',linewidth=1.5, linestyle="--", label="$ 100 \mu$Jy")
ax.plot(z200, Err_lnH200, color = '#00FF00',linewidth=1.5, linestyle="--", label="$ 200 \mu$Jy")
ax.plot(z1,  Err_lnH1 , color = 'darkorange',linewidth=1.0, linestyle="--", label=r"${\rm Euclid}$")
ax.get_legend_handles_labels()
ax.legend(loc='upper right')
ax.set_title('')
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=fontsize)
ax.set_ylabel(r"$\sigma_{H}/H\%$", fontsize=fontsize)
#ax.set_fontsize(20)
#plt.rc('font', **font)
plt.plot(z0, Err_lnH0 , 'bo')
ax.errorbar(z1, Err_lnH1 , color = 'darkorange',fmt= 'o')
plt.errorbar(z3, Err_lnH3  , color = 'red',fmt= 'o')
plt.plot(z7point3, Err_lnH7point3  , 'go')
plt.plot(z23, Err_lnH23, 'ko')
plt.plot(z70, Err_lnH70, 'co')
plt.plot(z100, Err_lnH100, 'mo')
ax.errorbar(z200, Err_lnH200 , color = '#00FF00',fmt= 'o')
ax.get_yaxis().set_major_formatter(tic.ScalarFormatter())
ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%f'))
ax.yaxis.set_major_formatter(tic.FuncFormatter(lambda x, pos: str(int(round(x)))))
#======= y axis
yticks = [0, 0.2, 0.5, 1, 2,5 , 10,20]
plt.yticks(yticks,[r'$0$',r'$0.2$', r'$0.5$',r'$1$', r'$ 2$',r'$5$',r'$10$',r'$20$' ])
plt.ylim(0.1,  20)
#======= x axis 
plt.xlim(0.,3.0)
xticks = arange(0, 3.3, 0.3)
plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])
plt.savefig('output_lnH_mario_bias.eps')
plt.show()
