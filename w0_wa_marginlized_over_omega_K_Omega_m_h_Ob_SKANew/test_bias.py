'This program to fit the data points of the bias to a specific function  and plot the reults '
import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt

#==================Functions ===========================
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
   #==============Input data ==============================

(xt, bias_rs0t, bias_rs1t, rms3_7t) =  np.loadtxt('bias_rms_7point3.txt', unpack= True)  
(xtt, bias_rs0tt, bias_rs1tt,rms3_7tt,  bias_rms23tt) =  np.loadtxt('bias_all.txt', unpack= True)  
z23 =np.array( [0.5, 1.0]) ;  rms23= np.array([1.2, 3.5])

(x0, bias_rms0,bias_rms1,rms3) =  np.loadtxt('HIBias_0_1muJy.txt', unpack= True)  
(x23, rms0 ,  rms1,  rms3, rms5,  rms6,   rms73,  rms10,  bias_rms23) =  np.loadtxt('bias_rms_23_2.txt', unpack= True)  
(x73, bias_rs0, bias_rs1,rms3, rms5, rms6,  bias_rms7_3) =  np.loadtxt('bias_rms_73_2.txt', unpack= True)  
(x5,    rms5) =  np.loadtxt('HIBias_5muJy.txt', unpack= True)  
(z, rms70) =  np.loadtxt('HIBias_70muJy.txt', unpack= True)  
(z100, r, rms100) =  np.loadtxt('HIBias_100muJy.txt', unpack= True)  

bias_rms0 = sqrt(bias_rms0)
bias_rms1 = sqrt(bias_rms1)
bias_rms23 = sqrt(bias_rms23)
bias_rms7_3 = sqrt(bias_rms7_3)
rms5 = sqrt(rms5)
rms70 = sqrt(rms70)
rms100 = sqrt(rms100)

#==============Fitting ====================================
#==========The initial guess ==================
p0=[6.3,2.]
p04=[5.74, 1.14]
#=========================================

plsq1t = opt.fmin(residuals, p0, args=(xtt,bias_rs1tt), maxiter=10000, maxfun=10000)
plsq2t = opt.fmin(residuals, p0, args=(xt,rms3_7t), maxiter=10000, maxfun=10000)
plsq23t = opt.fmin(residuals, p0, args=(z23,rms23), maxiter=10000, maxfun=10000)

plsq0= opt.fmin(residuals, p0, args=(x0,bias_rms0), maxiter=10000, maxfun=10000)
plsq1 = opt.fmin(residuals, p0, args=(x0,bias_rms1), maxiter=10000, maxfun=10000)
plsq5 = opt.fmin(residuals, p0, args=(x5,rms5), maxiter=10000, maxfun=10000)
plsq2 = opt.fmin(residuals, p0, args=(x73,bias_rms7_3), maxiter=10000, maxfun=10000)
plsq3 = opt.fmin(residuals, p0, args=(x23,bias_rms23), maxiter=10000, maxfun=10000)
plsq70 =  opt.fmin(residuals, p0, args=(z,rms70), maxiter=10000, maxfun=10000)
plsq100 = opt.fmin(residuals, p0, args=(z100,rms100), maxiter=10000, maxfun=10000)

#=================print Output ==============================
print '__________________________________________________'
print 'Rms in muJy ', '     c1   '      ,       '       |    c2  '
print '__________________________________________________'
print'1 muJy Bias old = ',  plsq1t[0],  ' |  ', plsq1t[1]
print '7.3muJy Bias old= ',  plsq2t[0],' |  ', plsq2t[1]
print '23 muJy Bias old  = ',  plsq23t[0], ' | ', plsq23t[1]

print'0 muJy Bias = ' , plsq0[0], '  |  ', plsq0[1]
print'1 muJy Bias = ',  plsq1[0],  ' |  ', plsq1[1]
print '5 muJy Bias = ',  plsq5[0],' |  ', plsq5[1]
print '7.3muJy Bias = ',  plsq2[0],' |  ', plsq2[1]
print '23 muJy Bias = ',  plsq3[0], ' | ', plsq3[1]
print '70 muJy Bias = ',  plsq70[0], ' | ', plsq70[1]
print '100 muJy Bias = ',  plsq100[0], ' | ', plsq100[1]

#================Producing data from the fitting function ===========

xrange = np.linspace(0, 2.5, 100)
xmin =  xrange-0.1
xmax  = xrange
x4range = np.linspace(0, 1.0, 100)

y0t=func(plsq1t,xrange)
y1t=func(plsq2t,xrange)
y23t=func(plsq23t,xrange)

y0=func(plsq0,xrange)
y1=func(plsq1,xrange)
y5=func(plsq5,xrange)
y2=func(plsq2,xrange)
y3=func(plsq3,xrange)
y70=func(plsq70,xrange)
y100=func(plsq100,xrange)
#=================  PLOTING =============================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

p0t,  = ax.plot(xtt,bias_rs1tt, 'bo')
p1t, = ax.plot(xt,rms3_7t , 'ro')
p223, =  plt.plot(z23, rms23 , 'yo')

#p0,  = ax.plot(x0, bias_rms0, 'bo')
p1, = ax.plot(x0, bias_rms1 , 'co')
#p50, =  plt.plot(x5, rms5 , 'yo')
p2, =  plt.plot(x73, bias_rms7_3 , 'go')
p3, =  plt.plot(x23, bias_rms23 , 'ko')
#p4, =  plt.plot(z, rms70 , 'co')
#p5, =  plt.plot(z100, rms100 , 'mo')

p01,  = plt.plot(xrange, y0t, color='b')
p11, = plt.plot(xrange, y1t, color='r')
p23, = plt.plot(xrange, y23t, color='y')
#p01,  = plt.plot(xrange, y0, color='b')
p81, = plt.plot(xrange, y1, color='c')
#p55, = plt.plot(xrange, y5, color='y')
p21, = plt.plot(xrange, y2, color='g')
p31, = plt.plot(xrange, y3, color='k')
#p41, = plt.plot(xrange, y70, color='c')
#p51, = plt.plot(xrange, y100, color='m')
plt.legend([ p0t, p1t,p223, p1, p2, p3] ,['$ 1\mu$Jy', '$7. 3\mu$Jy','$23\mu$Jy' ,'$1\mu$Jy new', '$7.3\mu$Jy new ' ,'$23\mu$Jy new ' ,'$70 \mu$Jy', '$100\mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.xlim(0,2.1,0.5)
plt.ylim(0., 10)
plt.savefig('fitted_bias_test_no_lines.eps')
plt.show()
print '============Thanks! ======================'
