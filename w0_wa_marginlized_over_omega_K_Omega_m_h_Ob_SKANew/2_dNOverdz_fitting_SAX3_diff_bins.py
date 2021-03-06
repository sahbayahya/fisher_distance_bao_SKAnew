from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import special
#******************************************************
#FUNCTIONS 
#******************************************************
def dNOverdz(zmin,zmax, c1,c2,c3):
	z= 0.5 * (zmin+ zmax)
	return 10**c1 * z**c2 * exp(-c3 * z)
	
def NOfz(zmin, zmax, c1,c2,c3):
	return  quad(dNOverdz, zmin, zmax, args=(zmax,c1,c2,c3,))[0]

def D(z):
     return cp.fgrowth(z,omega_M_0, unnormed=False)	
	
if __name__ == "__main__":
#*****************Testing the Functions**************************

	import doctest
	doctest.testmod()
	
#*********************INPUT DATA*********************************
        (x, total, rms0muJy,rms1muJy, rms3muJy, rms5muJy, rms6muJy, rms73muJy, rms10muJy, rms23muJy, rms40muJy, rms70muJy, rms100muJy, rms150muJy, rms200muJy) = loadtxt('HIdndzb3.txt', unpack=True)
 	c11  =6.23000844621; c21  = 1.82232488932  ; c31  = 0.89608495919    #0 muJy
	c12 =  7.33473992074 ; c22  =  3.02002614773  ; c32  = 5.33842588543 # 1 muJy
	c5 = 6.91214734152 ;   c51  =     2.38186955537 ; c52  =  5.84080403921 # 3 muJy
	c13 = 6.75522312458   ; c23  = 2.13569693813  ; c33 = 7.35533385121   # 7.3 muJy
	c14 = 6.01593890751   ; c24 = 1.43281797508   ; c34  =9.03321046833 # 23 muJy
	c16 =5.6200078849      ; c26 = 1.11377237426     ; c36 =13.0257055979  # 70 muJy
	c15  =  5.6266673048   ; c25 =1.40867290563     ; c35 =15.4937962327 # 100 muJy
	c200 =  5.00377325178 ; c200_2 = 1.04281566255  ; c200_3 = 17.5261229518 #200 muJy
#****************************************************************************

	dn = [] ; dn_f1 = [] ; dn_f2 = []; dn_f3 = []; dn_f4 = []; dn_f5 = []; dn_f6 = [];  dn_f200 = []
	#print x
	
	xrange = array([ 0.02, 0.04, 0.06, 0.08,  0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8,  0.9 ,  1.0,  1.1 , 1.2, 1.3  , 1.4 , 1.5 , 1.6 ,1.7 , 1.8, 1.9, 2.0])
       #xrange = array([ 0.116,    0.144,    0.2075 ,  0.4079  , 0.6236  , 0.8277,  0.988 ,   1.1734  , 1.3897 ,  1.6303  , 1.7666 ,  2.07])

	#xmin = xrange -0.1
	#xmax  =xrange+ 0.1
	xmin = [ 0.01,  0.03 , 0.05,  0.07, 0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
	xmax = [ 0.03,   0.05,  0.07,  0.09,  0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]
	kmax  =  empty(len(xrange)); kmax.fill(0.2)
	err_z =empty(len(xrange)); err_z.fill(0.0)
	volume =empty(len(xrange)); volume.fill(30000.0)
	bias = empty(len(xrange));bias.fill(1.0)
	#z = xrange
	print 
	print '==============================RESULTS========================='
	print 
	print 'xrange = ', xrange
	print 'The Output dndz for 1, 0, 3, 7.3, 23, 70, 100 muJy has been saved in: data_all_NOfz_SAX3_diff_14bin_new.txt'
	print
	print '======================Thanks!========================================'
	for i in range(len(xmin)):
		dn.append(NOfz(xmin[i],xmax[i], c11,c21,c31)) ;dn_f1.append(NOfz(xmin[i],xmax[i], c12,c22,c32)); dn_f2.append(NOfz(xmin[i],xmax[i], c13,c23,c33)); dn_f3.append(NOfz(xmin[i],xmax[i], c14,c24,c34)); dn_f4.append(NOfz(xmin[i],xmax[i], c15,c25,c35)) ; dn_f5.append(NOfz(xmin[i],xmax[i], c16,c26,c36)); dn_f6.append(NOfz(xmin[i],xmax[i], c5,c51,c52)); dn_f200.append(NOfz(xmin[i],xmax[i],c200,c200_2,c200_3))
        data_all_NOfz= concatenate((reshape(xrange,(len(xrange),1)),reshape(dn,(len(xrange),1)),reshape(dn_f1,(len(xrange),1)),reshape(dn_f6,(len(xrange),1)),reshape( dn_f2,(len(xrange),1)),reshape(dn_f3,(len(xrange),1)), reshape(dn_f5,(len(xrange),1)), reshape(dn_f4,(len(xrange),1)), reshape(dn_f200,(len(xrange),1))),axis=1)
        savetxt('data_all_NOfz_SAX3_diff_14bin_new.txt' , data_all_NOfz)
        	
