import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
#=============Functions ====================
def func(p,x):
   w=p[0]*np.exp(p[1]*x)
   #print w.size
   return w

def residuals(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
#============== Input files ====================   
(x1, bias_rm0,bias_rm1,bias_rm3) =  np.loadtxt('HIBias_0muJ_3muJy_modified.txt', unpack= True) 
(x, bias_rms0,bias_rms1,bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True) 
#(x2, dNofdz_rms_0muJy,dNofdz_rms_1muJy  , dNofdz_rms_3muJy ,dNofdz_rms_73muJy , dNofdz_rms_23muJy, dNofdz_rms_70muJy,dNofdz_rms_100muJy, dNofdz_rms_200muJy) = np.loadtxt('data_all_NOfz_SAX3_diff_14bin_new.txt', unpack=True)
(x2, dNofdz_rms_0muJy,dNofdz_rms_1muJy 
,dNofdz_rms_3muJy,dNofdz_rms_73muJy , dNofdz_rms_23muJy, dNofdz_rms_70muJy ,dNofdz_rms_100muJy,dNofdz_rms_200muJy) = np.loadtxt ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)
(xdn, dn0muJy,dn1muJy, dn3muJy, dn5muJy, dn6muJy, dn73muJy, dn10muJy, dn23muJy, dn40muJy, dn70muJy, dn100muJy, dn150muJy, dn200muJy) = loadtxt('n_table.txt', unpack=True)

bias_rms0 = sqrt(bias_rms0)
bias_rms1 = sqrt(bias_rms1)
bias_rms3 = sqrt(bias_rms3)
#==============Fitting =======================


p0=[6.3,2.]
p04=[5.74, 1.14]
plsq0= opt.fmin(residuals, p0, args=(x,bias_rms0), maxiter=10000, maxfun=10000)
plsq1 = opt.fmin(residuals, p0, args=(x,bias_rms1), maxiter=10000, maxfun=10000)
plsq3 = opt.fmin(residuals, p0, args=(x,bias_rms3), maxiter=10000, maxfun=10000)

#==========Print results =======================
print '__________________________________________________'
print 'Rms in muJy ', '     c1   '      ,       '       |    c2  '
print '__________________________________________________'
print'0 muJy Bias = ' ,  plsq0[0], '  |  ', plsq0[1]
print'1 muJy Bias = ',  plsq1[0],  ' |  ', plsq1[1]
print'3 muJy Bias = ',  plsq3[0],  ' |  ', plsq3[1]
#=================Input data ===================
#print 'len x = ', len(x2)
xrange =  array([ 0.02, 0.04, 0.06, 0.08,  0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8,  0.9 ,  1.0,  1.1 , 1.2, 1.3  , 1.4 , 1.5 , 1.6 ,1.7 , 1.8, 1.9, 2.0])#x2
xmin = xrange -0.1 # [0.01, 0.031 ,0.054,  0.016, 0.044,0.10749 , 0.3078, 0.52360, 0.7277, 0.888 ,  1.2734, 1.2897,  1.5303, 1.6666, 1.969999]     #
xmax  =xrange + 0.1 #  [0.03, 0.051,0.074, 0.216, 0.244, 0.3075,0.5079, 0.7236, 0.9277, 1.088, 1.0734,  1.4897, 1.6303,1.8666, 2.17]# 
xmin = [ 0.01,  0.03 , 0.05,  0.07, 0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
xmax = [ 0.03,   0.05,  0.07,  0.09,  0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]
#=============costumized list  ==================
#===== comment the next few lines if you want the whole range of the input =======
#dNofdz_rms_0muJy = dNofdz_rms_0muJy[5:24]
#dNofdz_rms_1muJy = dNofdz_rms_1muJy[5:24]
#dNofdz_rms_3muJy = dNofdz_rms_3muJy[5:24]
#x2 = x2[5:24]
xmin0 = xmin ; xmax0 = xmax
xmin3 =xmin ; xmax3= xmax

#================================================
#for i in range(len(xmin)):
	#print (xmin[i] + xmax[i])/ 2.0
y0=func(plsq0,xrange)
y1=func(plsq1,xrange)
y3=func(plsq3,xrange)
print 'xrange = ', xrange
print 'len(xrange) =', len(xrange)
print 'xmin=', xmin 
print 'xmax=',xmax

#============ Save Results  =========================
bias0  =  np.empty(len(xrange)); bias0.fill(1.0)
kmax  =  np.empty(len(xrange)); kmax.fill(0.5)
err_z =np.empty(len(xrange)); err_z.fill(0.0)
volume =np.empty(len(xrange)); volume.fill(30000.0)

dataN = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(xmax,(len(xrange),1)),reshape(dNofdz_rms_0muJy*100000,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1))),axis=1)

dataN5000 = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(xmax,(len(xrange),1)),reshape(dNofdz_rms_0muJy*100000,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1))),axis=1)

data = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(xmax,(len(xrange),1)),reshape(dNofdz_rms_0muJy,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1))),axis=1)
data2 = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(xmax,(len(xrange),1)),reshape(dNofdz_rms_1muJy,(len(xrange),1)),reshape(y1,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1))),axis=1)

data3 = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(xmax,(len(xrange),1)),reshape( dNofdz_rms_3muJy,(len(xrange),1)),reshape(y3,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1))),axis=1)
#savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_0N_mJy_SKANew_S3_5000.txt' , dataN5000)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_0N_mJy_SKANew_S3.txt' , dataN)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_0_mJy_SKANew_S3.txt' , data)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_1_mJy_SKANew_S3.txt' , data2)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_3_mJy_SKANew_S3.txt' , data3)

#============== Plotting ============================

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_title("fitting  Mario's BIAS")
p0,  = ax.plot(x, bias_rms0, 'bo')
p1, = ax.plot(x, bias_rms1 , 'ro')
p1, = ax.plot(x, bias_rms3 , 'yo')
p01,  = plt.plot(xrange, y0, color='b')
p11, = plt.plot(xrange, y1, color='r')
p11, = plt.plot(xrange, y3, color='y')
plt.legend([p01, p11] ,['$ 0\mu$Jy', '$1\mu$Jy', 'S$_{rms} = 3\mu$Jy' ,'S$_{rms} = 23\mu$Jy' ,'S$_{rms} = 100\mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('fittingMario_bias_using_1muJy_0muJy_3muJy.eps')
#plt.show()
print '============ Program excuted successfully ==========='
print '==================thanks ======================='
