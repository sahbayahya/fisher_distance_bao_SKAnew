from scipy import *
import numpy
from scipy import linalg


#********************************************* 
#               FUNCTIONS
#*********************************************
def FoM(dx, dy, dxy, Delta_x):
	part1 = (dx**2 + dy**2)/ 2.0
	part2 =sqrt( (((dx**2 - dy**2 )**2)/ 4.0) + dxy**2)
	a = abs(part1 + part2)
	b = abs(part1 -  part2)
	FoM =pi/( pi * Delta_x * sqrt(a)* sqrt(b) )
	return FoM


def DET2x2(A):
	A = A[0:2,0:2]
	DET = linalg.det(A)
	return DET

def identity(n):
    return [[1 if i==j else 0 for j in range(n)] for i in range(n)]


Planks_prior = mat('[1.83928663e+05   5.16525685e+04  -7.42050738e+06   -3.98758357e+06  -1.11710442e+06   1.32438370e+06 ; 5.16525685e+04   1.45055577e+04  -2.08389634e+06   -1.11983054e+06  -3.13715719e+05   3.71925825e+05  ;  -7.42050738e+06  -2.08389634e+06   3.64943809e+08    1.58599621e+08   4.25932543e+07  -5.16878541e+07  ; -3.98758357e+06  -1.11983054e+06   1.58599621e+08    8.70535526e+07   2.48738854e+07  -2.91740427e+07 ;   -1.11710442e+06  -3.13715719e+05   4.25932543e+07    2.48738854e+07   7.49686718e+06  -8.54525588e+06  ;  1.32438370e+06   3.71925825e+05  -5.16878541e+07   -2.91740427e+07  -8.54525588e+06   9.88949015e+06]' )

A_SKA = mat('[58341.202505538284        14073.055169744477       -139637.73719383232        189553.84457727903        341626.91147019085        320162.24456380110      ;   14073.055169744479        3494.1175338113753       -39187.532991549851        50890.827061364827        88438.554420194880        75431.974105263085      ;  -139637.73719383232       -39187.532991549851        2294041.5937065133       -2063748.4996850546       -2114782.6498606685       -34821.066256911057      ;   189553.84457727903        50890.827061364827       -2063748.4996850546        1948166.8314751021        2194853.2874319619        442609.73203912599      ;   341626.91147019085        88438.554420194880       -2114782.6498606685        2194853.2874319619        2932159.6782722347        1405369.0909180786      ;   320162.24456380110        75431.974105263085       -34821.066256911057        442609.73203912599        1405369.0909180790        2043040.2571025942]')

Delta_x = 1.0
#============================================================
Added_prior_matrix= Planks_prior + A_SKA
SKA_plus_prior_inv = linalg.inv(Added_prior_matrix)
#print SKA_plus_prior_inv
#===========================================================
SKA_plus_prior_inv[0,4] =sqrt(SKA_plus_prior_inv[0,0]) * ( SKA_plus_prior_inv [4,4] + SKA_plus_prior_inv[3,3] + 2. * SKA_plus_prior_inv[3,4])
SKA_plus_prior_inv[4,0]  = SKA_plus_prior_inv[0,4]

SKA_plus_prior_inv[1,4] =sqrt(SKA_plus_prior_inv[1,1]) * ( SKA_plus_prior_inv [4,4] + SKA_plus_prior_inv[3,3] + 2. * SKA_plus_prior_inv[3,4])
SKA_plus_prior_inv[4,1]  = SKA_plus_prior_inv[1,4]

SKA_plus_prior_inv[2,4] =sqrt(SKA_plus_prior_inv[2,2]) * ( SKA_plus_prior_inv [4,4] + SKA_plus_prior_inv[3,3] + 2. * SKA_plus_prior_inv[3,4])
SKA_plus_prior_inv[4,2]  = SKA_plus_prior_inv[2,4]

SKA_plus_prior_inv[3,4] =sqrt(SKA_plus_prior_inv[3,3]) * ( SKA_plus_prior_inv [4,4] + SKA_plus_prior_inv[3,3] + 2. * SKA_plus_prior_inv[3,4])
SKA_plus_prior_inv[4,3]  = SKA_plus_prior_inv[3,4]

SKA_plus_prior_inv[4,4] = ( SKA_plus_prior_inv [4,4] + SKA_plus_prior_inv[3,3] )

SKA_plus_prior_inv[5,4] =sqrt(SKA_plus_prior_inv[5,5]) * ( SKA_plus_prior_inv [4,4] + SKA_plus_prior_inv[3,3] + 2. * SKA_plus_prior_inv[3,4])
SKA_plus_prior_inv[4,5]  = SKA_plus_prior_inv[5,4]
#==========================================================
#print 'SKA_plus_prior  Cov= ', SKA_plus_prior_inv
w0 = sqrt(SKA_plus_prior_inv[0,0]) 
#print 'sigma_w0 = ', w0
wa = sqrt(SKA_plus_prior_inv[1,1])
w0a = sqrt(SKA_plus_prior_inv[0,1])
wa0 = sqrt(SKA_plus_prior_inv[1,0])
#print' w0a = ', w0a, 'wa0 =', wa0
#print'sigma_wa =', wa
ob0 =sqrt(SKA_plus_prior_inv[2,2]) 
#print'sigma_ob0 =', ob0
om0 = sqrt(SKA_plus_prior_inv[4,4])
#print 'sigma_om0 = ', om0
h = sqrt(SKA_plus_prior_inv[5,5])
#print 'sigma_h = ',  h
ok0 = sqrt(SKA_plus_prior_inv[3,3])
#print 'sigma_ok0=', ok0
FoM2 =  1.0/sqrt(SKA_plus_prior_inv[0,0] * SKA_plus_prior_inv[1,1] - SKA_plus_prior_inv[1,0]* SKA_plus_prior_inv[0,1])
#print 'FoM DETF (Coe 2009) for 3 muJy SKA + Planck =  ', FoM(w0, wa, w0a**2, Delta_x)
print  '&', w0,  '&', wa,  '&', om0,  '&', ob0,  '&', ok0,  '&',  h,  '&', FoM2 , ' \ \ ' 
print  '&','%.2f' % w0,  '&', '%.2f' % wa,  '&','%.3f' % om0,  '&', '%.3e' % ob0,  '&', '%.3f' % ok0,  '&',  '%.3f'  % h,  '&', '%.2d' % FoM2 , ' \\\ ' 

