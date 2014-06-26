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


Planks_prior = mat('[1.99579245e+05  -3.73667528e+04  -1.04936812e+04   1.39977603e+06    5.58643962e+05  -4.64225267e+04  -7.65181989e+04  -2.23806234e+03;  -3.73667528e+04   1.83928663e+05   5.16525685e+04  -7.42050738e+06   -3.98758357e+06  -1.11710442e+06   1.32438370e+06  -4.51559188e+02;  -1.04936812e+04   5.16525685e+04   1.45055577e+04  -2.08389634e+06   -1.11983054e+06  -3.13715719e+05   3.71925825e+05  -1.26811078e+02;  1.39977603e+06  -7.42050738e+06  -2.08389634e+06   3.64943809e+08    1.58599621e+08   4.25932543e+07  -5.16878541e+07   3.20338905e+04;   5.58643962e+05  -3.98758357e+06  -1.11983054e+06   1.58599621e+08    8.70535526e+07   2.48738854e+07  -2.91740427e+07   1.88438127e+04;  -4.64225267e+04  -1.11710442e+06  -3.13715719e+05   4.25932543e+07    2.48738854e+07   7.49686718e+06  -8.54525588e+06   1.25851649e+04;  -7.65181989e+04   1.32438370e+06   3.71925825e+05  -5.16878541e+07   -2.91740427e+07  -8.54525588e+06   9.88949015e+06  -1.01838183e+04; -2.23806234e+03  -4.51559188e+02  -1.26811078e+02   3.20338905e+04    1.88438127e+04   1.25851649e+04  -1.01838183e+04   1.51709659e+04]' )
print'========================================================'
Delta_x = 1.0
#===================Take the inverse of the prior matrix=========================================
#print Planks_prior
prior_inv = linalg.inv(Planks_prior)
#print prior_inv
#================= take the ns and sigma8 out============================================
prior_inv1 = mat('[-3.20220084e+03   4.73166943e+04   2.13900236e-05  5.48234698e+02  -8.08642777e+02  -4.32094449e+02 ;  4.73166848e+04  -7.21954893e+05  -3.10228608e-04  -7.99336941e+03   1.38090118+04   9.16638205e+03 ; 2.13900192e-05  -3.10228604e-04   2.28286048e-08 -3.48191844e-06   6.23816361e-06   4.04255574e-06 ;  5.48234563e+02  -7.99336894e+03  -3.48191828e-06  -8.58740714e+01   1.54794800e+02   1.07622797e+02 ;  -8.08642535e+02   1.38090109e+04   6.23816332e-06  1.54794800e+02  -2.81726183e+02  -1.97826849e+02 ;  -4.32094280e+02   9.16638144e+03   4.04255552e-06  1.07622796e+02  -1.97826848e+02  -1.40315032e+02]')
 #================invert again to get the Fisher matrix for w0, wa, ob0, ok0, 0de0, h ======================
prior_inv2= linalg.inv(prior_inv1)
#print prior_inv2
 #=================transform using the Jacobian rule from the parameters above to w0, wa, ob0, ok0, om0, h=============
 
M = mat('[1. 0. 0. 0. 0. 0.  ; 0. 1. 0. 0. 0. 0. ; 0.  0. 1. 0. 0. 0. ; 0. 0. 0. 1. 0. 0. ; 0. 0. 1. 1. 1. 0. ; 0. 0. 0. 0. 0. 1.]')
print 'M'
print M
print'========================================================'
#========== 
MT = M.T
print 'M^ T'
print MT
print'========================================================'

A_SKA = mat('[55091.338971399244        11572.661536079018       -142188.65827027761        189868.37032801486        280667.90511872421        292747.08362111344      ;   11572.661536079018        2575.7847864760602       -40041.436512017390        48725.483669898262        68032.940695658937        56567.738925845348      ;  -142188.65827027761       -40041.436512017390        2293270.0484161945       -2066470.0070332121       -2133220.1604082338       -53411.538516925059      ;   189868.37032801486        48725.483669898262       -2066470.0070332121        1948668.4584894511        2129817.8046591738        428269.17793047766      ;   280667.90511872421        68032.940695658937       -2133220.1604082338        2129817.8046591738        2491561.0517790639        961115.15964626428      ;   292747.08362111344        56567.738925845348       -53411.538516925059        428269.17793047766        961115.15964626428        1847029.8602369505]')
#======= Final stage is to (MT).F.M ===========================================

Final_prior_Fisher = dot( dot(MT, prior_inv2), M)
print 'final prior Fisher'
print Final_prior_Fisher 
print'========================================================'
#================Take the inverse again to get the fisher matrix to our parameters and add the SKA fisher matrix==========================================
SKA_plus_prior = Final_prior_Fisher  + A_SKA
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
w0a = sqrt(SKA_plus_prior[0,1])  ;  print 'sigma w0a = ', w0a   
wa0 = sqrt(abs(SKA_plus_prior[1,0]))   ; print' wa0 = ', wa0
ob0 =sqrt(SKA_plus_prior[2,2])    #print'sigma_ob0 =', ob0
ok0 = sqrt(abs(SKA_plus_prior[3,3]))   #print 'sigma_ok0=', ok0
om0 = sqrt(SKA_plus_prior[4,4])  #print 'sigma_om0 = ', om0
h = sqrt(SKA_plus_prior[5,5])       #print 'sigma_h = ',  h
#====================print all the results on one line ====================
FoM2 =  1.0/sqrt(SKA_plus_prior[0,0] * SKA_plus_prior[1,1] - SKA_plus_prior[1,0]* SKA_plus_prior[0,1])
print'========================================================'
print 'FoM DETF (Coe 2009) for 0 muJy SKA + Planck =  ', FoM(w0, wa, w0a**2, Delta_x)
print'===========================values of the sigmas of the cosmological parameters============================='
print  '&', w0,  '&', wa,  '&', om0,  '&', ob0,  '&', ok0,  '&',  h,  '&', FoM2 , ' \ \ ' 
print'================================================================================================'
print  '&','%.2f' % w0,  '&', '%.2f' % wa,  '&','%.3f' % om0,  '&', '%.3e' % ob0,  '&', '%.3f' % ok0,  '&',  '%.3f'  % h,  '&', '%.3d' % FoM2 , ' \\\ ' 
print'======================Thanks===================================================================='
