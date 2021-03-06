import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt

#================= Functions ======================
def func(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w

def residuals(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
   #=============Input data ============================
(z3, bias_rms0,bias_rms1,bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True) 
#(x5,rs0,rs1,rs3, rms5) =  np.loadtxt('HIBias_5muJy.txt', unpack= True)  
(z70, rms70) =  np.loadtxt('HIBias_70muJy.txt', unpack= True)  
(z100, r, rms100) =  np.loadtxt('HIBias_100muJy.txt', unpack= True)    
(x, dNofdz_rms_0muJy,dNofdz_rms_1muJy  ,dNofdz_rms_3muJy  ,  dNofdz_rms_73muJy , dNofdz_rms_23muJy,dNofdz_rms_70muJy  ,dNofdz_rms_100muJy,dNofdz_rms_200muJy ) = np.loadtxt   ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)
(xdn, dn0muJy,dn1muJy, dn3muJy, dn5muJy, dn6muJy, dn73muJy, dn10muJy, dn23muJy, dn40muJy, dn70muJy, dn100muJy, dn150muJy, dn200muJy) = loadtxt('n_table.txt', unpack=True)
#================ sqrt of the bias ====================
rms100 = sqrt(rms100)
rms70 = sqrt(rms70)
#x= x[0:10]
#dNofdz_rms_70muJy = dNofdz_rms_70muJy[0:10]
#dNofdz_rms_100muJy = dNofdz_rms_100muJy[0:10]
#dNofdz_rms_200muJy = dNofdz_rms_200muJy[0:10]
#==============Fitting ===============================
p0=[6.3,2.]
p04=[5.74, 1.14]
plsq70 =  opt.fmin(residuals, p0, args=(z70,rms70), maxiter=10000, maxfun=10000)
plsq100 = opt.fmin(residuals, p0, args=(z100,rms100), maxiter=10000, maxfun=10000)
print ' |   c1  | ',       '|         c2  |'
print x
#=============Simulate data ===========================
xrange = x
#xmin = xrange -0.1 #[0.01, 0.031 ,0.054,  0.016, 0.044,0.10749 , 0.3078, 0.52360, 0.7277, 0.888 ,  1.2734, 1.2897,  1.5303, 1.6666, 1.969999]     #
#xmax  = xrange + 0.1 # [0.03, 0.051,0.074, 0.216, 0.244, 0.3075,0.5079, 0.7236, 0.9277, 1.088, 1.0734,  1.4897, 1.6303,1.8666, 2.17]#
xmin = [ 0.01,  0.03 , 0.05,  0.07, 0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
xmax = [ 0.03,   0.05,  0.07,  0.09,  0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]
print xmin
print xmax



#==================Save Results =======================
print '========== 5, 23 , 70 , 100 and 200 muJy ==============='

x70 = x; x100 = x ; x200 = x


#============ 70 uJy ====================================== 
#x70 = x70[4:9]; dNofdz_rms_70muJy = dNofdz_rms_70muJy[4:9]
print x70
y70=func(plsq70,x70)
xmin70 = xmin; xmax70 = xmax
#xmin70 = xmin70[4:9]; xmax70= xmax70[4:9]
kmax  =  np.empty(len(x70)); kmax.fill(0.5)
err_z =np.empty(len(x70)); err_z.fill(0.0)
volume =np.empty(len(x70)); volume.fill(5000.0)
data70 = concatenate((reshape(volume,(len(x70),1)),reshape(xmin70,(len(x70),1)),reshape(xmax70,(len(x70),1)),reshape( dNofdz_rms_70muJy,(len(x70),1)),reshape(y70,(len(x70),1)),reshape(kmax,(len(x70),1)),reshape(err_z,(len(x70),1))),axis=1)
#========= 100 uJy ==========================================
#x100 = x100[4:8]; dNofdz_rms_100muJy = dNofdz_rms_100muJy[4:8]
print x100
y100=func(plsq100,x100)
xmin100 = xmin; xmax100 = xmax
#xmin100 = xmin100[4:8]; xmax100= xmax100[4:8]
kmax  =  np.empty(len(x100)); kmax.fill(0.5)
err_z =np.empty(len(x100)); err_z.fill(0.0)
volume =np.empty(len(x100)); volume.fill(5000.0)
data100 = concatenate((reshape(volume,(len(x100),1)),reshape(xmin100,(len(x100),1)),reshape(xmax100,(len(x100),1)),reshape( dNofdz_rms_100muJy,(len(x100),1)),reshape(y100,(len(x100),1)),reshape(kmax,(len(x100),1)),reshape(err_z,(len(x100),1))),axis=1)

#======== 200 uJy ============================================
#x200 = x200[4:7]; dNofdz_rms_200muJy = dNofdz_rms_200muJy[4:7]
print x200
y200=func(plsq100,x200)
xmin200 = xmin; xmax200 = xmax
#xmin200 = xmin200[4:7]; xmax200= xmax200[4:7]
kmax  =  np.empty(len(x200)); kmax.fill(0.5)
err_z =np.empty(len(x200)); err_z.fill(0.0)
volume =np.empty(len(x200)); volume.fill(5000.0)
print 'clear'
#data5 = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(xmax,(len(xrange),1)),reshape(dNofdz_rms_3muJy ,(len(xrange),1)),reshape(y3,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1))),axis=1)
data200 = concatenate((reshape(volume,(len(x200),1)),reshape(xmin200,(len(x200),1)),reshape(xmax200,(len(x200),1)),reshape( dNofdz_rms_200muJy,(len(x200),1)),reshape(y200,(len(x200),1)),reshape(kmax,(len(x200),1)),reshape(err_z,(len(x200),1))),axis=1)
#savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_5_mJy_SKANew.txt' , data5)
#===============Save the data ============================
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_70_mJy_SKANew_S3_2.txt' , data70)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_100_mJy_SKANew_S3_2.txt' , data100)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_200_mJy_SKANew_S3_2.txt' , data200)

# =================== Plotting ===========================

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#plt.xlim(0.,2.5,0.5)
#plt.ylim(1, 5)
p70, = plt.plot(x70, y70, color='m')
plt.plot(z70, rms70 , 'mo')
p100, = plt.plot(x100, y100, color='b')
plt.plot(z100, rms100 , 'bo')
p200, = plt.plot(x200, y200, color='r')
plt.plot(z100, rms100 , 'ro')
plt.legend([p70, p100,p200] ,[ ' $70\mu$Jy', ' $100\mu$Jy', ' $200\mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('fittingMario_bias_using_ObreschkowFunc_70_100_200.eps')
plt.show()
