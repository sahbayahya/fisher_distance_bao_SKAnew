!-----------------------------------------------------------------------------
!
! Comoving volume per area (area in units of deg^2), returned 
! in units of h^-3 Mpc^3.  Note h^-3!
! 
! *** Radiation density is ignored in the calculation, so
! this routine cannot be used to calculate the volume to
! the last scattering surface! ***
!
! << sample code >>
!
! USE cosmo
! USE angular_distance
! double precision :: vol,redshift1,redshift2,area ! in deg^2
! ! Specify three cosmological parameters in double precision.
! ! The data type has been defined in MODULE cosmo.
! ode0=0.723d0
! om0=0.277d0
! w   = -1d0
! CALL setup_da
! redshift1=1.8d0
! redshift2=3.5d0
! area=426d0
! print*,'omega matter=',om0
! print*,'omega de=',ode0
! print*,'w=',w
! CALL volume(redshift1,redshift2,area,vol)
! print*,'Volume(z) per',area,' deg^2 between z=' &
!       ,redshift1,' and',redshift2,' is',vol,' h^-3 Mpc^3'
! end
!
! August 23, 2008
! E. Komatsu
!
!-----------------------------------------------------------------------------
SUBROUTINE volume(zmin,zmax,area,vol)
! vol is the comoving volume per area (where area is in units of deg^2),
! in units of h^-3 Mpc^3.  Note h^-3!
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: zmin,zmax,area ! area in units of deg^2
  DOUBLE PRECISION, intent(OUT):: vol
  DOUBLE PRECISION :: dvdz,rombint,tol=1d-7
  EXTERNAL dvdz
  vol=rombint(dvdz,zmin,zmax,tol)
  vol=area*(3.14159d0/180d0)**2.*vol
 print*, 'volume', vol
  return
END SUBROUTINE volume
!-----------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION dvdz(z)
! dV/dz is differential comoving volume per steradian, in units of h^-3 Mpc^3
! Note h^-3!
use cosmo
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: z
  DOUBLE PRECISION :: one_over_h,da
  dvdz = (1d0+z)**2d0*da(z)**2d0*((c/H_0)*one_over_h(z))
 ! print*, 'da', da(z)
  return
END FUNCTION dvdz
