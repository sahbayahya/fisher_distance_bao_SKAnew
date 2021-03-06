'This Program to fit dndz for different rms sensitvities, produce the fitting values and plot  '
import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
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
#(z0, bias_rms00, bias_rms11, bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True)  
#(x,rms0muJy,rms1muJy, rms73muJy, rms23muJy, rms100muJy) = np.loadtxt('dndz_modified.dat', unpack=True)
#(x1, total, rm0muJy,rm1muJy, rm3muJy, rms5muJy, rm6muJy, rm73muJy, rm10muJy, rm23muJy, rms40muJy, rms70muJy, rm100muJy) = np.loadtxt('HIdndz_modified.txt', unpack=True)
(x2, to, rm00muJy,rm01muJy, rm03muJy, rm05muJy, rm06muJy, rm073muJy, rm010muJy, rm023muJy, rm040muJy, rm070muJy, rm100muJy, rm150muJy, rm200muJy) = np.loadtxt('HIdndzb_modified_high.txt', unpack=True)
(x1, to, rm0muJy,rm1muJy, rm3muJy, rm5muJy, rm6muJy, rm73muJy, rm10muJy, rm23muJy, rm40muJy, rm70muJy, rm0100muJy, rm0150muJy, rm0200muJy) = np.loadtxt('HIdndzb_modified.txt', unpack=True)
(x, total, rms0muJy,rms1muJy, rms3muJy, rms5muJy, rms6muJy, rms73muJy, rms10muJy, rms23muJy, rms40muJy, rms70muJy, rms100muJy, rms150muJy, rms200muJy) = np.loadtxt('HIdndzb3.txt', unpack=True)
p0=[5.52,  0.6, 4.6]
p04=[5.74, 1.14, 3.95]
plsqtot= opt.fmin(residuals, p0, args=(x,total), maxiter=10000, maxfun=10000)
plsq0= opt.fmin(residuals, p0, args=(x,rms0muJy), maxiter=10000, maxfun=10000)
plsq1 = opt.fmin(residuals, p04, args=(x1,rm1muJy), maxiter=10000, maxfun=10000)
plsq3= opt.fmin(residuals, p0, args=(x1,rm3muJy), maxiter=10000, maxfun=10000)


plsq5= opt.fmin(residuals, p0, args=(x1,rm5muJy), maxiter=10000, maxfun=10000)
plsq6= opt.fmin(residuals, p0, args=(x1,rm6muJy), maxiter=10000, maxfun=10000)
plsq7point3= opt.fmin(residuals, p0, args=(x1,rm73muJy), maxiter=10000, maxfun=10000)
plsq10= opt.fmin(residuals, p0, args=(x1,rm10muJy), maxiter=10000, maxfun=10000)

plsq23= opt.fmin(residuals, p0, args=(x1,rm23muJy), maxiter=10000, maxfun=10000)
plsq40= opt.fmin(residuals, p0, args=(x1, rm40muJy), maxiter=10000, maxfun=10000)
plsq70= opt.fmin(residuals, p0, args=(x1,rm70muJy), maxiter=10000, maxfun=10000)

plsq100 = opt.fmin(residuals, p04, args=(x2,rm100muJy), maxiter=10000, maxfun=10000)
plsq150 = opt.fmin(residuals, p04, args=(x2,rm150muJy), maxiter=10000, maxfun=10000)
plsq200 = opt.fmin(residuals, p04, args=(x2,rm200muJy), maxiter=10000, maxfun=10000)
#data = concatenate((reshape(x,(len(x),1)),reshape(total,(len(x),1))),axis=1)
#savetxt('dndz_1mJy.dat' , data)
print ' |   c1  | ',       '|         c2  |',        '|         c3  |'
print  '0 muJy ',plsq0[0],   plsq0[1], plsq0[2]
print '1 muJy', plsq1[0], plsq1[1]  ,plsq1[2]
print '3 muJy', plsq3[0], plsq3[1], plsq3[2]
print '5 muJy', plsq5[0], plsq5[1], plsq5[2]
print '6 muJy', plsq6[0], plsq6[1], plsq6[2]
print '7.3 muJy',  plsq7point3[0], plsq7point3[1], plsq7point3[2]
print '10 muJy', plsq10[0], plsq10[1], plsq10[2]
print '23 muJy',  plsq23[0], plsq23[1], plsq23[2]
print '40 muJy', plsq40[0], plsq40[1], plsq40[2]
print '70 muJy', plsq70[0], plsq70[1], plsq70[2]
print '100 muJy',  plsq100[0], plsq100[1], plsq100[2]
print '150 muJy',  plsq150[0], plsq150[1], plsq150[2]
print '200 muJy',  plsq200[0], plsq200[1], plsq200[2]
xrange = np.linspace(0, 3.0, 200)
#xrange = array([ 0.02, 0.04, 0.06, 0.08,  0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8,  0.9 ,  1.0,  1.1 , 1.2, 1.3  , 1.4 , 1.5 , 1.6 ,1.7 , 1.8, 1.9, 2.0])
xrange4 = np.linspace(0.01,0.6,16)
print 'xrang4 =' ,xrange4
#x4range = np.linspace(0, 2.0, 200)
y0=func(plsq0,xrange)
y1=func(plsq1,xrange)
y3=func(plsq3,xrange)
y5 = func(plsq5,xrange)
y6 = func(plsq6,xrange)
y2=func(plsq7point3,xrange)
y10=func(plsq10,xrange)
y23=func(plsq23,xrange)
y70=func(plsq70,xrange)
y40=func(plsq40,xrange)
y100=func(plsq100,xrange)
y150=func(plsq150,xrange)
y200=func(plsq200,xrange)

#data= concatenate((reshape(xrange,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(y1,(len(xrange),1)),reshape(y3,(len(xrange),1)),reshape( y2,(len(xrange),1)),reshape(y23,(len(xrange),1)), reshape(y70,(len(xrange),1)), reshape(y100,(len(xrange),1)), reshape(y200,(len(xrange),1))),axis=1)
#savetxt('data_all_dndz_SAX3_diff_14bin_new.txt' , data)

#================ 23 and 7.3 ============================

xrange7 = np.linspace(0.1, 1.3 , 16)
xrange23 = np.linspace(0.1, 0.8 , 16)
print 'xrang7 =' ,xrange7
y23_2=func(plsq23,xrange23)
y73=func(plsq7point3,xrange7)
#data7= concatenate((reshape(xrange7,(len(xrange7),1)),reshape( y73,(len(xrange7),1)), reshape(xrange23,(len(xrange23),1)), reshape(y23_2,(len(xrange23),1))),axis=1)
#savetxt('data_all_dndz_SAX3_diff_14bin_new_73_23.txt' , data7)

#======================================================
y704=func(plsq70,xrange4)
y404=func(plsq40,xrange4)
y1004=func(plsq100,xrange4)
y1504=func(plsq150,xrange4)
y2004=func(plsq200,xrange4)
data4= concatenate((reshape(y704,(len(xrange4),1)), reshape(y1004,(len(xrange4),1)), reshape(y2004,(len(xrange4),1))),axis=1)
savetxt('data_all_dndz_SAX3_diff_14bin_new_100.txt' , data4)
print '============ Program excuted successfully ==========='
print '==================Thanks! ======================'
# =================== Plotting ========================

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')
ax.plot(x, rms0muJy, 'bo')
ax.plot(x, rms3muJy , 'ro')
ax.errorbar(x, rms1muJy,  color = 'orange', fmt ='o')
ax.plot(x, rms5muJy, 'co')
ax.plot(x,rms6muJy , 'yo')
ax.plot(x,rms73muJy , 'go')
ax.errorbar(x, rms10muJy,  color = '#cc0066', fmt ='o')
ax.plot(x,rms23muJy , 'ko')
ax.errorbar(x, rms40muJy,  color = '#808000', fmt ='o')
ax.errorbar(x, rms70muJy,  color = '#008080', fmt ='o')
ax.plot(x,rms100muJy , 'mo')
ax.errorbar(x, rms150muJy, color = 'skyblue', fmt ='o')
ax.errorbar(x, rms200muJy , color = '#00FF00',fmt= 'o')
p01,  = ax.plot(xrange, y0, color='blue', linewidth=2.0, linestyle="-")
p11, =  ax.plot(xrange, y1, color='orange', linewidth=2.0, linestyle="-")
p31, =  ax.plot(xrange, y3,color='red', linewidth=2.0, linestyle="-")
p55, =  ax.plot(xrange, y5, color='cyan', linewidth=2.0, linestyle="-")
p6, =  ax.plot(xrange, y6, color='yellow', linewidth=2.0, linestyle="-")
p73, =  ax.plot(xrange, y2, color='green', linewidth=2.0, linestyle="-")
p10, =  ax.plot(xrange, y10, color='#cc0066', linewidth=2.0, linestyle="-")
p23, =  ax.plot(xrange, y23, color='black', linewidth=2.0, linestyle="-")
p40,  = ax.plot(xrange, y40, color='#808000', linewidth=2.0, linestyle="-")
p70,  = ax.plot(xrange, y70, color='#008080', linewidth=2.0, linestyle="-")
p100,  = ax.plot(xrange, y100, color='magenta', linewidth=2.0, linestyle="-")
p150,  = ax.plot(xrange, y150, color='skyblue', linewidth=2.0, linestyle="-")
p200,  = ax.plot(xrange, y200, color='#00FF00', linewidth=2.0, linestyle="-")
plt.legend([p01, p11, p31,p55, p6, p73, p10, p23, p40,  p70, p100,p150 , p200] ,['$0 \mu$Jy','$1 \mu$Jy', '$ 3 \mu$Jy', '$ 5 \mu$Jy','$ 6 \mu$Jy','$ 7.3 \mu$Jy' ,'$ 10 \mu$Jy' ,  '$ 23\mu$Jy',  '$ 40 \mu$Jy' , '$70\mu$Jy'  ,'$100\mu$Jy'  ,'$150\mu$Jy',  '$200\mu$Jy'], loc=1, borderaxespad=0.)#p01, p11, p31, p73, p23,  p70, p100,
plt.xlim(0.1,3.0 ,0.5)
plt.ylim(1, 10**8)
#================ x axis 
xticks = arange(min(x), max(x)+1, 0.3)
plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])
plt.xlabel(r"${ \rm redshift} (z)$", fontsize=15)
#============= y axis
yticks = [1, 10,100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
plt.yticks(yticks,[r'$1$', r'$10$',r'$10^2$', r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$', r'$10^7$' ])
plt.ylabel(r'$\frac{d{\rm N}}{dz}(z) \ [ {\rm deg}^{-2} \ {\rm per} \ {\rm unit} \ z ]$', fontsize= 15)
#plt.ylabel(r"$\frac{{\rmd}N}{{\rm d}z} $/deg$^{-2}$")
#========= save fig#
plt.savefig('fittingMario_dNOverdz_using_ObreschkowFunc_diff_14bins_3.eps')
plt.show()

