from numpy import * 
from scipy import *


def S_rms(T_inst, A_eff,v , Dv, tp,D_dish, Nb):
	lamda = 3.e8 /( v )                            
	T_sky = 60. * (300./v)**2.55 		# MHz . K  
	part1= 368. *( (T_inst + T_sky) /20.) # uJy / K
	part2 = 25000. /A_eff 			# m^2
	part3 =sqrt( 0.01/Dv) 			# 1/ MHz
	#thetaB =pi/4. * sqrt(66.* lamda / D_dish) ; tp = tot * (thetaB**2/ S_area)
	part4 = sqrt(1./tp) 				# hour
	return (part1 * part2 * part3 * part4), 1.0# , (Nb * thetaB**2)

Area = 20000.0
Srms = [40.0, 70.0, 100.0, 200.0]
for i in range(len(Srms)):
	S_rms_new = Srms[i]*sqrt(Area/(5000.0))
	print S_rms_new


(ASKAP, S_ASKAP) = S_rms(50.0, 3257.0,  1000.0,  0.01, 60.0,  12.0 ,  36.0)
(Meerkat, S_Meerkat) = S_rms(20.0, 5955.0,  1420.0,  0.01,  3.34 ,  13.5,  1.0)
(SKA1_Sur, S_SKA1_Sur) = S_rms(30.0, 8482.0,  1150.0 ,  0.01, 50.0,  15.0,  1.0)
(SKA2, S_SKA2) = S_rms(20.0, 400000.0, 850.0 , 0.01, 10.0, 50.0,250.0 )
print 'ASKAP Srms = ',  ASKAP #, '  , S_area = ', S_ASKAP
print 'Meerkat Srms= ', Meerkat #, , '  , S_area = ', S_Meerkat
print 'SKA1_Sur Srms = ', SKA1_Sur # , ' , S_area = ',  S_SKA1_Sur
print 'SKA2_ Srms= ', SKA2
	 
