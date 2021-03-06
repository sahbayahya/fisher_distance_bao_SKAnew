import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
#===================Functions ==============================
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
   
  #==================Inputs=================================== 
(x73, bias_rs0, bias_rs1,rms3, rms5, rms6,  bias_rms7_3) =  np.loadtxt('bias_rms_73_2.txt', unpack= True)  
(x23, rms0 ,  rms1,  rms3, rms5,  rms6,  rms73,  rms10,  bias_rms23) =  np.loadtxt('bias_rms_23.txt', unpack= True)
#(xdn, dn0muJy,dn1muJy, dn3muJy, dn5muJy, dn6muJy, dn73muJy, dn10muJy, dn23muJy, dn40muJy, dn70muJy, dn100muJy, dn150muJy, dn200muJy) = loadtxt('n_table.txt', unpack=True)
(x, dNofdz_rms_0muJy,dNofdz_rms_1muJy 
,dNofdz_rms_3muJy,dNofdz_rms_73muJy , dNofdz_rms_23muJy, dNofdz_rms_70muJy ,dNofdz_rms_100muJy,dNofdz_rms_200muJy) = np.loadtxt ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)# np.loadtxt ('data_all_NOfz_SAX3_diff_14bin_new.txt', unpack=True)
#(xrange7, dNofdz_rms_73muJy , xrange23, dNofdz_rms_23muJy) = np.loadtxt ('data_all_dndz_SAX3_diff_14bin_new_73_23.txt', unpack=True)
p0=[0.3,1.]
p04=[5.74, 1.14]
#=================Fitting ===================================
bias_rms7_3 = sqrt(bias_rms7_3)
bias_rms23 = sqrt(bias_rms23)
plsq2 = opt.fmin(residuals, p0, args=(x73,bias_rms7_3), maxiter=10000, maxfun=10000)
plsq23 = opt.fmin(residuals, p0, args=(x23,bias_rms23), maxiter=10000, maxfun=10000)
print '====================='
print ' |   c1  | ',       '|         c2  |'
print '7.3muJy Bias = ',  plsq2[0],' |  ', plsq2[1]
print '23 muJy Bias = ',  plsq23[0], ' | ', plsq23[1]
print '====================='
x2 = x
#dNofdz_rms_73muJy = dNofdz_rms_73muJy[2:16]
#dNofdz_rms_23muJy = dNofdz_rms_23muJy[2:12]
#x= x[2:16]
#x2= x2[2:12]
xrange = array([ 0.02, 0.04, 0.06, 0.08,  0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8,  0.9 ,  1.0,  1.1 , 1.2, 1.3  , 1.4 , 1.5 , 1.6 ,1.7 , 1.8, 1.9, 2.0])
xmin = [ 0.01,  0.03 , 0.05,  0.07, 0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
xmax = [ 0.03,   0.05,  0.07,  0.09,  0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]
#=========== SAVE 7.3 output ====================================
xrange7= xrange
y2=func(plsq2,xrange7)
print '7.3 uJy range = ',  xrange7
#xmin =xrange7-0.1#
#xmax  =xrange7 + 0.1
#dNofdz_rms_73muJy = dNofdz_rms_73muJy[0:8]
kmax  =  np.empty(len(xrange7)); kmax.fill(0.5)
err_z =np.empty(len(xrange7)); err_z.fill(0.00)
volume =np.empty(len(xrange7)); volume.fill(30000.0)
data = concatenate((reshape(volume,(len(xrange7),1)),reshape(xmin,(len(xrange7),1)),reshape(xmax,(len(xrange7),1)),reshape( dNofdz_rms_73muJy,(len(xrange7),1)),reshape(y2,(len(xrange7),1)),reshape(kmax,(len(xrange7),1)),reshape(err_z,(len(xrange7),1))),axis=1)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_7point3_mJy_SKANew_S3.txt' , data)
#========= save 23 ouput =====================================
xrange23 = xrange
y23=func(plsq23,xrange23)
print '23 uJy range = ' , xrange23 
#xmin =xrange23-0.1 #[0.01, 0.031 ,0.054,  0.016, 0.044,0.10749 , 0.3078, 0.52360, 0.7277, 0.888 ,  1.2734, 1.2897,  1.5303, 1.6666, 1.969999]     # 
#xmax  =xrange23 + 0.1
print 'xmin = ', xmin 
print 'xmax = ', xmax
kmax  =  np.empty(len(xrange23)); kmax.fill(0.5)
err_z =np.empty(len(xrange23)); err_z.fill(0.00)
volume =np.empty(len(xrange23)); volume.fill(30000.0)
data23 = concatenate((reshape(volume,(len(xrange23),1)),reshape(xmin,(len(xrange23),1)),reshape(xmax,(len(xrange23),1)),reshape( dNofdz_rms_23muJy,(len(xrange23),1)),reshape(y23,(len(xrange23),1)),reshape(kmax,(len(xrange23),1)),reshape(err_z,(len(xrange23),1))),axis=1)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/w0_wa_marginlized_over_omega_K_Omega_m_h_Ob_SKANew/number_sax3_23_mJy_SKANew_S3.txt' , data23)

print '============ Program excuted successfully ==========='
print '===============Thanks! ========================='
#print data
# ============= Plotting======================================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#plt.xlim(0.,2.5,0.5)
#plt.ylim(1, 5)
p2, =  plt.plot(x73, bias_rms7_3 , 'ro')
p21, = plt.plot(xrange7, y2, color='r')
p23, =plt.plot(xrange23, y23, color='y')
plt.plot(x23, bias_rms23 , 'yo')
plt.legend([p21, p23] ,[' $ 7.3 \ \mu$Jy',  '$23  \ \mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('fittingMario_bias_using_ObreschkowFunc_7point3_23.eps')
plt.show()

