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

(x0, bias_rms0,bias_rms1,bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True) 
(x5, rs0, rs1, rs3,rms5, rms6, rs7_3) =  np.loadtxt('HIBias_5_6_73.txt', unpack= True)  
(x73, bias_rs0, bias_rs1,rs03, rs05, rs06,  bias_rms7_3) =  np.loadtxt('bias_rms_73_2.txt', unpack= True)  
(z10, rms10) =  np.loadtxt('HIBias_10muJy.txt', unpack= True)  
(x23, rs0 ,  rs1,  rms3, rs5,  rs6,   rs73,  rs10,  bias_rms23) =  np.loadtxt('bias_rms_23_2.txt', unpack= True)  
(x, rs00 ,  rs10,  rms30, rs50,  rs60,   rs730,  r10,  bias_rm23) =  np.loadtxt('bias_rms_23.txt', unpack= True)  
(z40, rms40) =  np.loadtxt('HIBias_40muJy.txt', unpack= True)  
(z70, rms70) =  np.loadtxt('HIBias_70muJy.txt', unpack= True)  
(z100, r, rms100) =  np.loadtxt('HIBias_100muJy.txt', unpack= True)  

bias_rms0 = sqrt(bias_rms0)
bias_rms1 = sqrt(bias_rms1)
bias_rms3 = sqrt(bias_rms3)
rms5 = sqrt(rms5)
rms6= sqrt(rms6)
bias_rms7_3 = sqrt(bias_rms7_3)
rs7_3 = sqrt(rs7_3)
rms10= sqrt(rms10)
rs10= sqrt(r10)
bias_rms23 = sqrt(bias_rms23)
bias_rm23 = sqrt(bias_rm23)
rms40= sqrt(rms40)
rms70 = sqrt(rms70)
rms100 = sqrt(rms100)


#==============Fitting ====================================
#==========The initial guess ==================
p0=[6.3,2.]
p04=[5.74, 1.14]
#=========================================

plsq0= opt.fmin(residuals, p0, args=(x0,bias_rms0), maxiter=10000, maxfun=10000)
plsq1 = opt.fmin(residuals, p0, args=(x0,bias_rms1), maxiter=10000, maxfun=10000)
plsq3= opt.fmin(residuals, p0, args=(x0,bias_rms3), maxiter=10000, maxfun=10000)
plsq5 = opt.fmin(residuals, p0, args=(x5,rms5), maxiter=10000, maxfun=10000)
plsq6 = opt.fmin(residuals, p0, args=(x5,rms6), maxiter=10000, maxfun=10000)
plsq7_3 = opt.fmin(residuals, p0, args=(x73,bias_rms7_3), maxiter=10000, maxfun=10000)
plsq10 = opt.fmin(residuals, p0, args=(x,rs10), maxiter=10000, maxfun=10000)
plsq23 = opt.fmin(residuals, p0, args=(x,bias_rm23), maxiter=10000, maxfun=10000)
plsq40 = opt.fmin(residuals, p0, args=(z40,rms40), maxiter=10000, maxfun=10000)
plsq70 =  opt.fmin(residuals, p0, args=(z70,rms70), maxiter=10000, maxfun=10000)
plsq100 = opt.fmin(residuals, p0, args=(z100,rms100), maxiter=10000, maxfun=10000)

#=================print Output ==============================
print '__________________________________________________'
print 'Rms in muJy ', '     c1   '      ,       '       |    c2  '
print '__________________________________________________'
print'0 muJy Bias = ' , plsq0[0], '  |  ', plsq0[1]
print'1 muJy Bias = ',  plsq1[0],  ' |  ', plsq1[1]
print '3 muJy Bias = ',  plsq3[0],' |  ', plsq3[1]
print '5 muJy Bias = ',  plsq5[0],' |  ', plsq5[1]
print '6 muJy Bias = ',  plsq6[0],' |  ', plsq6[1]
print '7.3muJy Bias = ',  plsq7_3[0],' |  ', plsq7_3[1]
print '10 muJy Bias = ',  plsq10[0],' |  ', plsq10[1]
print '23 muJy Bias = ',  plsq23[0], ' | ', plsq23[1]
print '40 muJy Bias = ',  plsq40[0],' |  ', plsq40[1]
print '70 muJy Bias = ',  plsq70[0], ' | ', plsq70[1]
print '100 muJy Bias = ',  plsq100[0], ' | ', plsq100[1]

#================Producing data from the fitting function ===========

xrange = np.linspace(0, 2.5, 100)
xmin =  xrange-0.1
xmax  = xrange
x4range = np.linspace(0, 1.0, 100)
#plsq0[0] = 0.764931E+00;  plsq0[1] =  0.453041E+00
y0=func(plsq0,xrange)
y1=func(plsq1,xrange)
y3=func(plsq3,xrange)
y5=func(plsq5,xrange)
y6=func(plsq6,xrange)
y7_3=func(plsq7_3,xrange)
y10=func(plsq10,xrange)
y23=func(plsq23,xrange)
y40=func(plsq40,xrange)
y70=func(plsq70,xrange)
y100=func(plsq100,xrange)
#=================  Plotting =============================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

p0,  = ax.plot(x0, bias_rms0, 'bo')
p1, = ax.plot(x0, bias_rms3 , 'ro')
ax.errorbar(x0, bias_rms1,  color = 'orange', fmt ='o')
p73, =  plt.plot(x5, rs7_3 , 'go')
plt.plot(x23, bias_rms23 , 'ko')
plt.errorbar(z10, rms10 , color='#cc0066', fmt='o')
plt.plot(z100, rms100 , 'mo')
plt.errorbar(z40, rms40 , color='#808000', fmt='o')
plt.errorbar(z70, rms70 , color='#008080', fmt='o')
plt.plot(x5, rms6 , 'yo')
plt.plot(x5, rms5 , 'co')
p01,  = ax.plot(xrange, y0, color='blue', linewidth=2.0, linestyle="-")
p31, =  ax.plot(xrange, y3, color='red', linewidth=2.0, linestyle="-")
p11, =  ax.plot(xrange, y1,color='orange', linewidth=2.0, linestyle="-")
p55, =  ax.plot(xrange, y5, color='cyan', linewidth=2.0, linestyle="-")
p6, =  ax.plot(xrange, y6, color='yellow', linewidth=2.0, linestyle="-")
p73, =  ax.plot(xrange, y7_3, color='green', linewidth=2.0, linestyle="-")
p10, =  ax.plot(xrange, y10, color='#cc0066', linewidth=2.0, linestyle="-")
p23, =  ax.plot(xrange, y23, color='black', linewidth=2.0, linestyle="-")
p40,  = ax.plot(xrange, y40, color='#808000', linewidth=2.0, linestyle="-")
p70,  = ax.plot(xrange, y70, color='#008080', linewidth=2.0, linestyle="-")
p100,  = ax.plot(xrange, y100, color='magenta', linewidth=2.0, linestyle="-")
#plt.legend([p31] ,['$7.3\mu$Jy', '$10\mu$Jy' , '$ 23\mu$Jy',  '$40\mu$Jy'  ,'$70\mu$Jy'  ,'$100\mu$Jy'], loc='best', borderaxespad=0.)
plt.legend([p01, p11, p31,p55, p6, p73, p10, p23, p40,  p70, p100] ,['$0 \mu$Jy','$1 \mu$Jy', '$ 3 \mu$Jy', '$ 5 \mu$Jy','$ 6 \mu$Jy','$ 7.3 \mu$Jy' ,'$ 10 \mu$Jy' ,  '$ 23\mu$Jy',  '$ 40 \mu$Jy' , '$70\mu$Jy'  ,'$100\mu$Jy'  ,'$150\mu$Jy',  '$200\mu$Jy'], loc='best', borderaxespad=0.)
plt.xlabel(r"$ {\rm redshift} (z)$", fontsize=15)
plt.ylabel(r"$b(z)$", fontsize=15)
plt.ylim(0., 10)
#======= x axis 
xticks = arange(0, 3.3, 0.3)
plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])
plt.xlim(0.,2.5)
plt.savefig('fitted_bias.eps')
plt.show()
print '============ Program excuted successfully ==========='
print '==================Thanks! ======================'
