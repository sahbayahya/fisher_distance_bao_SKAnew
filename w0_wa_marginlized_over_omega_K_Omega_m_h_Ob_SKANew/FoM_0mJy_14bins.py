from scipy import *
from numpy import *
from scipy import linalg
#import warnings
#warnings.simplefilter("ignore", np.ComplexWarning)
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




def add_cmb_1(matrix):
	prior_inv1 = mat('[-3.20220084e+03   4.73166943e+04   2.13900236e-05  5.48234698e+02  -8.08642777e+02  -4.32094449e+02 ;  4.73166848e+04  -7.21954893e+05  -3.10228608e-04  -7.99336941e+03   1.38090118+04   9.16638205e+03 ; 2.13900192e-05  -3.10228604e-04   2.28286048e-08 -3.48191844e-06   6.23816361e-06   4.04255574e-06 ;  5.48234563e+02  -7.99336894e+03  -3.48191828e-06  -8.58740714e+01   1.54794800e+02   1.07622797e+02 ;  -8.08642535e+02   1.38090109e+04   6.23816332e-06  1.54794800e+02  -2.81726183e+02  -1.97826849e+02 ;  -4.32094280e+02   9.16638144e+03   4.04255552e-06  1.07622796e+02  -1.97826848e+02  -1.40315032e+02]')
	prior_inv2= linalg.inv(prior_inv1)
	full matrix = prior_inv2 + matrix
	w0 = sqrt(SKA_plus_prior[0,0])    #print 'sigma_w0 = ', w0
	wa = sqrt(SKA_plus_prior[1,1])    ; print 'sigma_wa =', wa
	w0a = (SKA_plus_prior[0,1])  ;  print 'sigma w0a = ', w0a   
	wa0 =((SKA_plus_prior[1,0]))   ; print' wa0 = ', wa0
	ob0 =sqrt(SKA_plus_prior[2,2])    #print'sigma_ob0 =', ob0
	ok0 = sqrt((SKA_plus_prior[3,3]))   #print 'sigma_ok0=', ok0
	om0 = sqrt(SKA_plus_prior[4,4])  #print 'sigma_om0 = ', om0
	h = sqrt(SKA_plus_prior[5,5])       #print 'sigma_h = ',  h 
	print  '&','%.2f' % w0,  '&', '%.2f' % wa,  '&','%.3f' % om0,  '&', '%.3e' % ob0,  '&', '%.3f' % ok0,  '&',  '%.3f'  % h,  '&', '%.3d' % FoM2 , ' \\\ ' 
	print'======================Thanks===================================================================='
	return 
	

#=================== SKA Fisher matrix ====================================================


Planks_prior = mat('[1.99579245e+05  -3.73667528e+04  -1.04936812e+04   1.39977603e+06    5.58643962e+05  -4.64225267e+04  -7.65181989e+04  -2.23806234e+03;  -3.73667528e+04   1.83928663e+05   5.16525685e+04  -7.42050738e+06   -3.98758357e+06  -1.11710442e+06   1.32438370e+06  -4.51559188e+02;  -1.04936812e+04   5.16525685e+04   1.45055577e+04  -2.08389634e+06   -1.11983054e+06  -3.13715719e+05   3.71925825e+05  -1.26811078e+02;  1.39977603e+06  -7.42050738e+06  -2.08389634e+06   3.64943809e+08    1.58599621e+08   4.25932543e+07  -5.16878541e+07   3.20338905e+04;   5.58643962e+05  -3.98758357e+06  -1.11983054e+06   1.58599621e+08    8.70535526e+07   2.48738854e+07  -2.91740427e+07   1.88438127e+04;  -4.64225267e+04  -1.11710442e+06  -3.13715719e+05   4.25932543e+07    2.48738854e+07   7.49686718e+06  -8.54525588e+06   1.25851649e+04;  -7.65181989e+04   1.32438370e+06   3.71925825e+05  -5.16878541e+07   -2.91740427e+07  -8.54525588e+06   9.88949015e+06  -1.01838183e+04; -2.23806234e+03  -4.51559188e+02  -1.26811078e+02   3.20338905e+04    1.88438127e+04   1.25851649e+04  -1.01838183e+04   1.51709659e+04]' )
print'========================================================'
Delta_x = 2.31
#===================Take the inverse of the prior matrix=========================================
#print Planks_prior
#prior_inv = linalg.inv(Planks_prior)
#print prior_inv
#================= take the ns and sigma8 out============================================

 #================invert again to get the Fisher matrix for w0, wa, ob0, ok0, 0de0, h ======================
#
#print prior_inv2
 #=================transform using the Jacobian rule from the parameters above to w0, wa, ob0, ok0, om0, h=============
print Planks_prior
M = mat('[1. 0.  0. 0. 0. 0. 0. 0.;  0. 1. 0. 0. 0. 0. 0. 0.  ; 0.  0. 1. 0. 0. 0. 0. 0. ;0.  0.  0. 1. 0. 0. 0. 0. ; 0. 0. 0. 0. 1. 0. 0.  0. ;0.  0. 0. 1. 1. 1. 0. 0. ;  0. 0. 0. 0. 0. 0. 1. 0.; 0. 0. 0. 0. 0. 0. 0. 1.]')
print 'M'
print M
print'========================================================'
#========== 
MT = M.T
print 'M^ T'
print MT
print'========================================================'
#======= Final stage is to (MT).F.M ===========================================

Final_prior_Fisher = dot( dot(MT, Planks_prior), M)
print 'final prior Fisher'
print Final_prior_Fisher 

##===========stack zeros on the matrix to match planks' one========================
A_SKA = mat('[302199.98663453473        71107.325974370295       -97261.582932366437        467298.79632609867        1354283.9790009155        1891870.1684092784      ;   71107.325974370295        17250.473439956440       -28791.958783779919        118665.45684335398        336860.14613860723        444483.19887715677      ;  -97261.582932366437       -28791.958783779919        2302012.3347055237       -2012059.2320107589       -1924306.9602475378        241525.60722084789      ;   467298.79632609867        118665.45684335398       -2012059.2320107589        2306377.6073531481        3430064.5315410281        2248113.0163574219      ;   1354283.9790009155        336860.14613860723       -1924306.9602475378        3430064.5315410281        7483930.7542739045        8009187.1311108023      ;   1891870.1684092784        444483.19887715677        241525.60722084789        2248113.0163574219        8009187.1311108023        12240539.325799888]')

newraw= linspace(0,0,6)
A_SKA = vstack((A_SKA,newraw))
A_SKA = vstack((newraw, A_SKA))
newcolumn = linspace(0., 0., 8)
A_SKA = column_stack((newcolumn, A_SKA))
A_SKA = column_stack(( A_SKA, newcolumn))
print A_SKA

print'========================================================'


#================Take the inverse again to get the fisher matrix to our parameters and add the SKA fisher matrix==========================================
SKA_plus_prior = A_SKA +  Final_prior_Fisher 
print 'SKA + priorFinal'
print SKA_plus_prior
print'========================================================'
SKA_plus_prior =  linalg.inv(SKA_plus_prior)
print 'The SKA Cov matrix + Plank'
print SKA_plus_prior 
print'========================================================'
#================Print final results for each parameter ==================
#print 'SKA_plus_prior = ', SKA_plus_prior_inv
w0 = sqrt(SKA_plus_prior[0,0])    #print 'sigma_w0 = ', w0
wa = sqrt(SKA_plus_prior[1,1])    ; print 'sigma_wa =', wa
w0a = (SKA_plus_prior[0,1])  ;  print 'sigma w0a = ', w0a   
wa0 =((SKA_plus_prior[1,0]))   ; print' wa0 = ', wa0
ob0 =sqrt(SKA_plus_prior[2,2])    #print'sigma_ob0 =', ob0
ok0 = sqrt((SKA_plus_prior[3,3]))   #print 'sigma_ok0=', ok0
om0 = sqrt(SKA_plus_prior[4,4])  #print 'sigma_om0 = ', om0
h = sqrt(SKA_plus_prior[5,5])       #print 'sigma_h = ',  h
#====================print all the results on one line ====================
FoM2 =  1.0/sqrt(SKA_plus_prior[0,0] * SKA_plus_prior[1,1] - SKA_plus_prior[1,0]* SKA_plus_prior[0,1])
print'========================================================'
print 'FoM DETF (Coe 2009) for 0 muJy SKA + Planck =  ', FoM(w0, wa, w0a, Delta_x)
print'===========================values of the sigmas of the cosmological parameters============================='
print  '&', w0,  '&', wa,  '&', om0,  '&', ob0,  '&', ok0,  '&',  h,  '&', FoM2 , ' \ \ ' 
print'================================================================================================'
print  '&','%.2f' % w0,  '&', '%.2f' % wa,  '&','%.3f' % om0,  '&', '%.3e' % ob0,  '&', '%.3f' % ok0,  '&',  '%.3f'  % h,  '&', '%.3d' % FoM2 , ' \\\ ' 
print'======================Thanks===================================================================='
