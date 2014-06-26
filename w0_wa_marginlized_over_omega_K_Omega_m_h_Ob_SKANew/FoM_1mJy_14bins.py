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
	DET = A[0,0]*A[1,1] - A[0,1]*A[1,0]
	return DET
def Fom_DETF(A):
	inverse_of_A = linalg.inv(A)
	print inverse_of_A
	det_of_A = DET2x2(A)
	return 1./(sqrt(det_of_A))

A_SKA = mat('[166140.52587176536        39612.663537319502        227601.29440770994        10989.423167505282        744439.81222608045        744439.81222608045      ;   39612.663537319502        9802.3647956804944        57122.371298135753        2642.7452339608735        188207.14948274710        188207.14948274710      ;   227601.29440770994        57122.371298135753        351945.48912233161        15194.877122166072        1104192.2435014464        1104192.2435014464      ;   10989.423167505282        2642.7452339608735        15194.877122166072        741.83623269732641        50291.130145544543        50291.130145544543      ;   744439.81222608045        188207.14948274710        1104192.2435014464        50291.130145544543        3724973.4785104864        3684973.4785104864      ;   744439.81222608045        188207.14948274710        1104192.2435014464        50291.130145544543        3684973.4785104864        25551189.575014990]')

A_SKA_marginlized = mat('[0.02557      -0.00260  ; -0.00260       0.14397]')
print '************************ Coe 2009 FoM DETF***********************************'
w0 =   0.02557; 		 wa= 0.14397 ;			 w0a= -0.00260; 		Delta_x = 6.17
#w0_cmb = 0.04444;	 wa_cmb = 0.39337; 		w0a_cmb = -0.01550
#print 
print 'FOM DETF (Coe2009) for 1 muJy SKA = ', FoM(w0,wa,w0a,Delta_x)
#print 'FoM DETF (Coe 2009) for 23 muJy SKA + CMB =  ', FoM(w0_cmb, wa_cmb, w0a_cmb, Delta_x)
#print '************************* SKA det(sqrt(F)) FoM ***********************************************************'
#print    
M0mJy = linalg.sqrtm(A_SKA,disp=1)
#print  ' square root of  M0mJy= ' , M0mJy
#print  'FoM of M0mJy 5x5 = ' , DET5x5(M0mJy)/7.8
#print  'FoM of 1 M0mJy 2x2 /7.8= ' , DET2x2(M0mJy)/7.8
#print  'FoM of 23 M0mJy 2x2 = ' , DET2x2(M0mJy)

print '************************* SKA FoM using DETF defintion from fisher4cast***********************************************************'
print    
M0mJy = linalg.sqrtm(A_SKA,disp=1)
#print  ' square root of  M0mJy= ' , M0mJy
#print  'FoM of M0mJy 5x5 = ' , DET5x5(M0mJy)/7.8
print  'FoM of 1 muJy DETF = ' , Fom_DETF(A_SKA_marginlized)
#print  'FoM of 23 M0mJy 2x2 = ' , DET2x2(M0mJy)


#print '************************* SKA sqrt(det(F)) FoM ***********************************************************'
#print    
M0mJy = linalg.sqrtm(A_SKA,disp=1)
#print  ' square root of  M0mJy= ' , M0mJy
#print  'FoM of M0mJy 5x5 = ' , DET5x5(M0mJy)/7.8
#print  'FoM of 1 M0mJy 2x2 /7.8= ' , sqrt(DET2x2(A_SKA))/7.8
#print  'FoM of 23 M0mJy 2x2 = ' , sqrt(DET2x2(A_SKA))
#print '***************** CMB + SKA FoM**********************************'
#M0mJy_cmb = linalg.sqrtm(A_cmb,disp=1)
#print  ' square root of  M0mJy_cmb = ' , M0mJy_cmb
#print  'FoM of M0mJy_cmb 5x5 = ' , DET5x5(M0mJy_cmb)/7.8
#print  'FoM of 23 M0mJy_cmb  2x2 = ' , DET2x2(M0mJy_cmb)/7.8
#print 
#print '*********************************************************************************************'

