from uncertainties import ufloat
from uncertainties.umath import *
print 
print '*************Values *********************'
print 
H = ufloat(67.4000000, 1.40000000)  # H = 67.4+/-1.4
h = H/100.0
print 'h = ', h
obh2 = ufloat(0.0220700000, 0.000330000000)
print 'obh2 = ' , obh2
ob = obh2/h**2
ErrorLnOb = 0.00124
LnOb = log(ob)
#*******************************

#print 'ErrorLnOb0 =' , ErrorLnOb
#print 'ob0 = ', ob
#print 'ErrorOb0 = ' , ob.s 

#*********************************
print 
print '**** check if ErrorLnOmega_b = ErrorOmega_b/  Omega_b******* '
print 
ErrorOb =  exp(ErrorLnOb)
print 'LnOb0 = ', LnOb 
print 'ErrorLnOb = ', (ErrorLnOb/LnOb)
print 'ErrorOb /Ob = ', ob.s / ob.n
print 'ErrorLnOb x Ob = ', LnOb.s * ob.n
