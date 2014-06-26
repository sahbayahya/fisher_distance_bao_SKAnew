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

A_SKA = mat('[57653.329497105718        14238.042279285237       -139426.11496117630        191003.67339122304        346684.01854563423        315147.21042618755      ;   14238.042279285235        3602.9691334389890       -39081.815559567425        51704.527039432141        90964.869196133630        76609.079462818045      ;  -139426.11496117630       -39081.815559567425        2294142.4348854232       -2062987.3837260695       -2112372.8622265379       -33297.673006971083      ;   191003.67339122304        51704.527039432141       -2062987.3837260695        1953710.3370368842        2213041.5696984753        452206.56910674786      ;   346684.01854563423        90964.869196133630       -2112372.8622265379        2213041.5696984753        2989746.0383721478        1441773.4076397382      ;   315147.21042618755        76609.079462818045       -33297.673006971083        452206.56910674798        1441773.4076397382        2004530.7297580042]')

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

