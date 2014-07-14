PROGRAM Fisher_Distance
  USE cosmo
  USE growth
  USE linearpk
  USE angular_distance
  ! A sample program for computing the errors on ln(Da) and ln(H), 
  ! from the information of Baryon Acoustic Oscillations only.
  ! Ref: Seo & Eisenstein, ApJ, 665, 14 (2007)
  ! September 17, 2008: E.Komatsu
  ! Modified to report also on the errors on w, Omega_k, and Omega_m.
  ! November 2, 2008: E.Komatsu
  IMPLICIT none
  integer :: npara=2 ! # of parameters
  double precision, dimension(2)   :: der
  double precision, dimension(2,2) :: cov,fis,fistot,fistest
  double precision, dimension(6,6) :: fis3x3, fis6x6
  integer,          dimension(2)   :: work
  double precision :: sigma0=12.4d0*0.817d0/0.9d0 
  ! non-linearity scale in units of h Mpc^-1, rescaled by WMAP5's sigma8 value (0.817)
  double precision :: BAO_AMP=0.5817d0 ! A_0 in Seo&Eisenstein
  double precision :: kmax_ov_h,k_ov_h,k0_ov_h, kmin_ov_h, k0_new, k_new! h Mpc^-1
  double precision :: linear_pk,pk,p02,p01,Pb
 double precision :: wdamp_perp, wdamp_para, wdamp_silk, wdamp_factor, wdamp_sinterm
  double precision :: sigma_para,sigma_perp,sigma_silk=8.38d0 ! h^-1 Mpc
  double precision :: dlnk,mu,mu2,dmu,factor,wkmu,dummy,wdamp,Rmu
  double precision :: z,zin,beta,bias,g,dgdlna,sigma_z,delta_v,fz, w,dvdz
  double precision :: zmin,zmax,area,Vsurvey,Vtot,ngal,dndz,da
  character(len=128) :: filename
  integer :: n,i,j,ibin,nbins,ifile
  external linear_pk,dgdlna,g,dvdz,da
  ! Specify three cosmological parameters for computing the growth factor
  ! and the angular diameter distance. 
  ! The data type has been defined in MODULE cosmo.
  !############################################
  
  !******************* PLANCK ********************************
  H_0 = 67d0
  c = 300000d0
  om0=0.316d0
  ok0 = 0.0d0
  ode0= 0.686d0
  ob0 = 0.049
  ok0 = 0d0! 1d0 - om0 -ob0-ode0
  w0=-1d0
 ! w = -1d0
  w_a = 0d0
  
 ! ********* WMAP 3 ******************************************
!  H_0 = 73d0
!  c = 299999d0
!  om0=0.24d0
!  ok0 = 0.0d0
!  ob0 =(0.0223d0/(H_0/100d0)**2)
!  print'(1A20,1F9.5)', 'ob0 = ', ob0
!  ode0= 1d0 - ok0 -ob0 - om0 
!  print '(1A20,1F9.5)', 'ode0 = ', ode0
!  w0=-1d0
!  w_a = 0d0
  CALL setup_growth ! tabulate the growth factor
  CALL setup_da       ! tabulate the angular diameter distance
  ! ask for survey specific parameters
!==============Enter the File name ================================== 
!filename = 'number_sax3_7point3_mJy_SKANew_S3.txt'
!filename='number_7point3_z=1.txt'!
filename= 'number_EuclidmJy_ref.txt'
  open(2,file=filename,status='unknown')
  
!=============  Read in linear P(k) ===================================
  filename='linear_matterpower_1.dat' 
! filename= 'wmap5baosn_max_likelihood_matterpower_at_z=30.dat'
  n=127 ! no. of lines in the file
  zin=1d0 ! redshift of the input power spectrum
  CALL open_linearpk(filename,n) ! This fuction just to interpolate the values in the filename
!================================================================  
  ! loop over redshift bins
  fistot=0d0   ! total Fisher matrix integrated over all the redshift bins
  fis3x3=0d0 ! total Fisher matrix for w, Omega_k, and Omega_m
  fis6x6=0d0 ! total Fisher matrix for w, Omega_k, and Omega_m
  Vtot=0d0   ! total survey volume
  nbins=0    ! initialize the number of bins
  !=========== Loop over redshift bins =======================================
  do ibin=1,100
     read(2,*,end=10)area,zmin,zmax,dndz,bias,kmax_ov_h,kmin_ov_h,delta_v ! uncomment for Euclid
    ! read(2,*,end=10)area,zmin,zmax,dndz,bias,kmax_ov_h,delta_v ! Uncomment for SKA
     nbins=nbins+1
    CALL volume(zmin,zmax,area,Vsurvey)
     z=0.5d0*(zmin+zmax)
     beta=(1d0+dgdlna(z)/g(z))/bias
     w = w0 + w_a*(z/(1d0+z))
     ngal=dndz !  Uncomment for  Euclid
    !ngal=dndz/((3.14159d0/180d0)**2*dvdz(z)) ! Uncomment for SKA
     Vtot=Vtot+Vsurvey
    ! print*, 'da', da(z)
   ! Vsurvey =  13.38d0*10d0**9d0/(H_0/100d0)**3d0 ! Phil 
     fz =(1d0 + z)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(z/(1d0+ z)))
     sigma_z= (1d0+z)*(delta_v/2.998d5) &
          *((c/H_0)/dsqrt((om0+ob0)*(1d0+z)**3d0+ok0* (1d0+z)**2d0+ ode0* fz)) ! in units of h^-1 Mpc
     sigma_perp=sigma0*g(z)/(1d0+z) ! in units of h^-1 Mpc
     sigma_para=sigma0*g(z)/(1d0+z)*(2d0+dgdlna(z)/g(z)) ! in units of h^-1 Mpc
     !==================Test the growth =================================
    !    open(101,file='nz_fz_Dz_test_rms=73.txt' ,status='unknown')
   !	write(101,*) z,ngal, (1d0+dgdlna(z)/g(z)), g(z)/(1d0+z)
   !=================  computing the Fisher matrix... =================================
     open(11,file=filename,status='old')
     read(11,*)k0_ov_h,dummy
     k_ov_h=k0_ov_h
     fis=0d0
     fistest= 0d0
     !================ loop over k ===================================
     do while (k_ov_h<=kmax_ov_h)
        read(11,*)k_ov_h,dummy
        !dlnk=dlog(k_ov_h)-dlog(k0_ov_h) ! uncomment
        !print*, 'dlnk', dlnk 
        !===========================================================
         k0_new=6.8492999999999996E-005 
         k_new =9.5589999999999998E-005
         dlnk=dlog(k_new)-dlog(k0_new)                         !uncomment to use with Plank's parameters 
        factor=(k_ov_h)**3d0*dlnk/(2.d0 *2d0 * 2d0 *3.1415926535d0**2d0) ! h^3 Mpc^-3
        !=========== loop over mu..========================
         mu=0d0
        dmu=1d-3
        do while (mu<=1d0)
           mu2= mu*mu
           ! P02 is P(k) at k=0.2 h Mpc^-1 in unitys of h^-3 Mpc^3
           p02=linear_pk(0.2d0)*(g(z)/g(zin)*(1d0+zin)/(1d0+z))**2d0 &
              *bias**2d0*(1d0+beta*mu2)**2d0
           ! P(k) the galaxy power spectrum  in units of h^-3 Mpc^3  
           pk=linear_pk(k_ov_h)*(g(z)/g(zin)*(1d0+zin)/(1d0+z))**2d0 &
               *bias**2d0*(1d0+beta*mu2)**2d0
              ! print*, p02
            !P01  P(k) at k=0.1 h Mpc^-1  in units of h^-3 Mpc^3  
              p01=linear_pk(k_ov_h)*(g(z)/g(zin)*(1d0+zin)/(1d0+z))**2d0&
                *bias**2d0!*(1d0+beta*mu2)**2d0
           !==========================================     
           wkmu=Vsurvey*(ngal*p02/(1d0+ngal*pk))**2d0
           wdamp=exp(-(k_ov_h*sigma_z)**2d0*mu2) & 
                *exp(-(k_ov_h*sigma_perp)**2d0*(1d0-mu**2d0)) & 
                *exp(-(k_ov_h*sigma_para)**2d0*mu**2d0) &
                *exp(-2d0*(k_ov_h*sigma_silk)**1.4d0)
            wdamp_perp = exp(-(k_ov_h*sigma_perp)**2d0*(1d0-mu**2d0))
            wdamp_para =exp(-(k_ov_h*sigma_para)**2d0*mu**2d0)
            wdamp_silk = exp(-2d0*(k_ov_h*sigma_silk)**1.4d0)
            wdamp_factor = dsqrt(2d0 *3.1415926535d0**2d0) * BAO_AMP *  p02
            wdamp_sinterm = ( dsin(k_ov_h *100d0)/(k_ov_h*150d0))
           Pb= wdamp_factor *  wdamp_sinterm * wdamp_perp  &
           *  wdamp_para *  wdamp_silk
           der(1)=(1d0-mu2) ! dlnPkdlnDa
           der(2)=(-mu2) ! dlnPkdlnH
           do i=1,npara
               fistest=factor *  Vsurvey*(ngal*p01/(1d0+ngal*p01))**2d0
               fis(i,:)= fis(i,:) + factor*wkmu*der(i)*der(:)*dmu &
                *(8d0*3.1415926535d0**2d0 *BAO_AMP**2d0)*wdamp ! factor*wkmu!
           enddo
           mu=mu+dmu
        enddo
       !================= save files ====================================== 
      ! open(5,file='Pk_Euclid_testPlank.txt' ,status='unknown')
     ! write(5,*) k_ov_h, pk! , 1d0/ngal
     !  open(5,file='gz.txt' ,status='unknown')
     !  write(5,*)  z, g(z)
    !	open(100,file='cosmic_v_k_Euclid_fis_test.txt' ,status='unknown')
    !	write(100,*) k_ov_h,  dsqrt(1d0/fistest(1,1))!
    ! 	write(5,*) ((1d0/P02)*dvdz(z) ), (dndz*((3.14159d0/180d0)**2*dvdz(z)))
    ! 	print*, 1d0/P02, dvdz(z)
   ! 	write(5,*) z, ((1d0/P02)*dvdz(z) ), (dndz*((3.14159d0/180d0)**2*dvdz(z)))
   !  open(5,file='Pb_vs_k_testz=1.txt' ,status='unknown')
   !   write(5,*) k_ov_h , Pb
   !   open(150,file='wdamps_testz=1.txt' ,status='unknown')
   !  write(150,*) k_ov_h , wdamp_perp, wdamp_para, wdamp_silk, wdamp_factor, wdamp_sinterm
   !======================================================================
        k0_ov_h=k_ov_h
      	enddo
      	fistot=fistot+fis
     !====================================================
     ! open(3,file='no_galaxies_70_diff_14bins_S3.txt' ,status='unknown')
     ! write(3,'(6F20.5)') dndz*area, ngal
      print*,'=== redshift bin#=(',ibin,') ==='
      print'(1F6.3,1A7,1F5.3)',zmin,'< z < ',zmax
      print*,'Vsur =',Vsurvey,' h^-3 Mpc^3'
      print*,'ngal =',ngal,'h^3 Mpc^-3'
      print*,'1/ngal =',1/ngal,'h^-3 Mpc^3'
      print'(1A7,1F8.5)','Bias =',bias
      print'(1A7,1F8.5)','Beta =',beta
      print'(1A7,1F8.5)','g(z) =',g(z)
      print'(1A7,1F8.5)','f(z) =',1d0+dgdlna(z)/g(z)
      print'(1A7,1F8.5,1A9)','sigz =',sigma_z,' h^-1 Mpc'
      print'(1A7,1F8.5,1A9)','spara=',sigma_para,' h^-1 Mpc'
      print'(1A7,1F8.5,1A9)','sperp=',sigma_perp,' h^-1 Mpc'
      print'(1A7,1F8.5,1A9)','kmin =',kmin_ov_h,' h Mpc^-1'
      print'(1A7,1F8.5,1A9)','kmax =',kmax_ov_h,' h Mpc^-1'
      CALL report_result(z,bias,npara,fis)
      !print*,''
      close(11)
      CALL transform_fisher(z,fis,fis3x3)
 enddo
10 close(2)
  CALL close_linearpk
  if(nbins>1)then
     write(12,*),'=== combined ==='
     write(12,*),'Vsur =',Vtot*1d-9,' h^-3 Gpc^3'
     write(12,*), 'ngal =',ngal,' 10^-3 h^3 Mpc^-3'
    CALL report_result(z, bias,npara,fistot)
  endif
  print*,''
 CALL report_result3x3(fis3x3)
  print*,''
CALL add_cmb(fis3x3)
CALL report_result3x3(fis3x3)
END PROGRAM Fisher_Distance

!====================================================================
!===================! SUBROUTINES  ===================================

SUBROUTINE report_result(z,bias,npara,fis)
  IMPLICIT none
  integer, intent(IN) :: npara
  double precision, intent(IN) :: fis(npara,npara), z,bias
  double precision, allocatable, dimension(:,:) :: cov
  double precision, allocatable, dimension(:)   :: work
  double precision :: r12,err_lnda,err_lnh,err_lnR,beta, linear_pk,dgdlna,g
  integer :: i,j
  external linear_pk,dgdlna,g
  ALLOCATE(cov(2,2),work(2))
  cov=fis
  CALL DVERT(cov,2,2,work)
  beta=(1d0+dgdlna(z)/g(z))/bias
  r12=cov(1,2)/sqrt(cov(1,1)*cov(2,2))
  err_lnda=sqrt(cov(1,1))
  err_lnh=sqrt(cov(2,2))
  err_lnR=err_lnda*sqrt((1d0-r12**2d0) &
       /(1d0+2d0*r12*err_lnda/err_lnh+(err_lnda/err_lnh)**2d0))
 !==============save desired fisher files =======================
!open(12,file='Fisher_7point3mJy_diff_14bins_S3.txt' ,status='unknown')
!  write(12,'(4F18.5)') err_lnda,err_lnh,err_lnR
! open(13,file='output_7point3mJy_diff_14bins_S3.txt', status='unknown')
! write(13,'(6F18.5)') z, err_lnda*1d2,err_lnh*1d2,err_lnR*1d2, beta
!========================================================
  print'(1A15,1F9.5)','Err[lnDa](%) =',err_lnda*1d2
  print'(1A15,1F9.5)','Err[lnH](%)  =',err_lnh*1d2
  print'(1A15,1F9.5)','r(lnDa,lnH)  =',r12
  print'(1A15,1F9.5)','Err[lnR](%)  =',err_lnR*1d2
  DEALLOCATE(cov,work)
  return
END SUBROUTINE report_result

!===========================================================
SUBROUTINE report_result3x3(fis)
  USE cosmo
  IMPLICIT none
  double precision, intent(IN) :: fis(6,6)
  double precision :: work(6),cov(6,6), A(6,6), M55DET, DET5x5,DET2x2,A05(6,6)
  integer :: i, j 
  cov=fis
 !==============================Calculate DET 2x2 ===============
  CALL DVERT(cov,6,6,work)
  A = fis
   DET2x2 =  A(1,1)*A(2,2) - A(1,2)*A(2,1)  
  write(12,*)'*** Error on w, with  Om0, ob0, ok0 and  H0 marginalized over ***'
 write(12,*) 'The Fisher Matrix = ', '[', A(1,1) ,    A(1,2),  A(1,3) , A(1,4),  A(1,5), A(1,6),  ';'&
   , A(2,1) ,  A(2,2),  A(2,3), A(2,4),  A(2,5), A(2,6), ';' &
   , A(3,1),  A(3,2),  A(3,3), A(3,4),  A(3,5),  A(3,6),  ';'  &
   , A(4,1),  A(4,2),  A(4,3), A(4,4),  A(4,5),  A(4,6), ';'    &
   , A(5,1),  A(5,2),  A(5,3), A(5,4),  A(5,5),  A(5,6), ';'  &
   , A(6,1),  A(6,2),  A(6,3), A(6,4),  A(6,5),  A(6,6), ']'     
   
 !====================write the resutls to file 12 ========================
   write(12,*),'Figure of Merit  =' ,  1d0/dsqrt((cov(1,1) * cov(2,2) - cov(1,2)* cov(2,1)))    
  write(12,'(1A20,1F9.5)')'Err[w] =',dsqrt(cov(1,1))
  write(12,'(1A20,1F9.5)')'Err[wa] =',dsqrt(cov(2,2))
  print*,'Figure of Merit  =' , (1d0/ dsqrt((cov(1,1) * cov(2,2) - cov(1,2)* cov(2,1))) )!/7.8!(3.14d0*dsqrt(7.61d0))
  print'(1A20,1F9.5)', 'Err[w] =',dsqrt(cov(1,1))
  print'(1A20,1F9.5)', 'Err[wa] =',dsqrt(cov(2,2))
  write(12,'(1A20,1F9.5)')'Err[wa_w] =',(cov(1,2))
  write(12,'(1A20,1F9.5)')'Err[ Omega_m]=',(dsqrt(cov(5,5)))
  write(12,* )'Err[ Omega_b]=',(dsqrt(cov(6,6))* ob0)
  write(12, '(1A20,1F9.5)')'Err[Omega_m_w]=',(cov(1,5))
  write(12,'(1A20,1F9.5)')'Err[Omega_k]=',dsqrt(cov(3,3))
  write(12,'(1A20,1F9.5)')'Err[Omega_k_w]=',(cov(3,1))
 write(12,'(1A20,1F9.5)') 'Err[H_0] =',dsqrt(cov(4,4))
  return
END SUBROUTINE report_result3x3
!================================================================
SUBROUTINE transform_fisher(z,fisDH,fis3x3)
  USE cosmo
  IMPLICIT none
  integer :: a,b,i,j
  double precision, intent(IN)    :: z,fisDH(2,2)
  double precision, intent(INOUT) :: fis3x3(6,6)
  double precision :: dpdq(6,2)
  double precision :: chi,h2,func0,func1,func2,func3,func4,rombint,fz
  external h2,func0,func1,func2,func3,func4,rombint
  chi=rombint(func0,0d0,z,1d-7)
   fz =(1d0 + z)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(z/(1d0+ z)))
  dpdq(1,1)=-1.5d0*ode0*rombint(func1,0d0,z,1d-7)/chi         ! dlnDa/dw
  dpdq(1,2)= 1.5d0*ode0*dlog(1d0+z)*fz/h2(z)                       !dlnH/dw
  dpdq(2,1)=-1.5d0*ode0*rombint(func4,0d0,z,1d-7)/chi        !dlnDa/dwa
  dpdq(2,2)= 1.5d0*ode0*(dlog(1d0+z) - z/(1d0+ z))*fz/h2(z) !dlnH/dwa
  dpdq(4,1)=-0.5d0*rombint(func2,0d0,z,1d-7)/chi+chi**2d0/6d0 !dlnDa/dOmega_k
  dpdq(4,2)= 0.5d0*((1d0+z)**2d0 - fz)/h2(z)                         !dlnH/dOmega_k
  dpdq(6,1)= -100d0/H_0		                                      !dlnDa/dH_0
  dpdq(6,2)= 100d0/H_0		                                              !dlnH/dH_0
  dpdq(5,1)=-0.5d0*rombint(func3,0d0,z,1d-7)/chi                !dlnDa/d(Omega_m)
  dpdq(5,2)= 0.5d0*((1d0+z)**3d0 - fz)/h2(z)                        !dlnH/dOmega_m
  dpdq(3,1)=-0.5d0*ob0*rombint(func3,0d0,z,1d-7)/chi               !dlnDa/d(Omega_b)
  dpdq(3,2)= 0.5d0*ob0*((1d0+z)**3d0 - fz)/h2(z)                       !dlnH/dOmega_b
  do a=1,6
     do b=1,6
        do i=1,2
           do j=1,2
              ! transform and accumulate fis3x3
              fis3x3(a,b)=fis3x3(a,b)+dpdq(a,i)*dpdq(b,j)*fisDH(i,j) 
           enddo
        enddo
     enddo
  enddo
  return
END SUBROUTINE transform_fisher
!===============================================================
SUBROUTINE add_cmb(fis3x3)
  USE cosmo
  ! Add the distance information from CMB [Eq.(38) of Shoji et al. 2008]
  IMPLICIT none
  integer :: a,b,i,j
  double precision :: da_accuracy=0.1d0 ! Percent error in Da(zcmb) 
  double precision :: ob_accuracy, h_accuracy= 0.007,  om_accuracy!= 0.017
  double precision :: ob0h2= 0.02207d0, ob0h2_accuracy= 0.00017d0, Planck_h= 0.673d0, Planck_ob0 
   double precision :: om0h2= 0.12038d0, om0h2_accuracy=0.0013d0,  Planck_om0 
  double precision, intent(INOUT) :: fis3x3(6,6)
  double precision :: dpdq(6),zcmb=1090d0
  double precision :: chi,h2,func0,func1,func2,func3,func4, rombint
  external h2,func0,func1,func2,func3,func4,rombint
  chi=rombint(func0,0d0,zcmb,1d-7)
 dpdq(1)=-1.5d0*ode0*rombint(func1,0d0,zcmb,1d-7)/chi         ! dlnDa/dw
  dpdq(2)=-1.5d0*ode0*rombint(func4,0d0,zcmb,1d-7)/chi          ! dlnDa/dwa
  dpdq(4)=-0.5d0*rombint(func2,0d0,zcmb,1d-7)/chi+chi**2d0/6d0 ! dlnDa/dOmega_k
  dpdq(6)= -100d0/H_0	 						   ! dlnDa/dH0
  dpdq(5)=-0.5d0*rombint(func3,0d0,zcmb,1d-7)/chi          ! dlnDa/dlnOmega_m
  dpdq(3)=-0.5d0*ob0*rombint(func3,0d0,zcmb,1d-7)/chi  !dlnDa/dlnOmega_b
   do a=1,6
     do b=1,6
     ! add Da(z=1090)
       fis3x3(a,b)=fis3x3(a,b)+dpdq(a)*dpdq(b)*1d4/da_accuracy**2d0
     enddo
  enddo
  Planck_ob0  = ob0h2/(Planck_h**2d0) 
  ob_accuracy = dsqrt( (ob0h2_accuracy/ob0h2)**2d0  +  ((h_accuracy)/ Planck_h)**2d0  ) * Planck_ob0 
    Planck_om0  = om0h2/(Planck_h**2d0) 
  om_accuracy = dsqrt( (om0h2_accuracy/om0h2)**2d0  +  ((h_accuracy)/ Planck_h)**2d0  ) * Planck_om0 
  ! add Omega_matter
 ! fis3x3(1,1) = fis3x3(1,1) + .172276d6
 ! fis3x3(2,2)= fis3x3(2,2) + .139551d5
 !fis3x3(3,3) = fis3x3(3,3) + 1d0/(0.115717676721)**2d0 !.253489d8 
 !fis3x3(4,4)=fis3x3(4,4)+  1d0/(h_accuracy)**2d0
 ! fis3x3(5,5)=fis3x3(5,5)+   1d0/(om_accuracy)**2d0
 ! fis3x3(6,6)=fis3x3(6,6)+   1d0/(ob_accuracy)**2d0 
  print*,'==============Planck Priors Added==========================='
  print'(1A23,1F9.5)','Err[lnDa(z=1090)](%) =',da_accuracy
  print'(1A23,1F9.5)','Err[ h ]  =',h_accuracy
  print'(1A23,1F9.5)','Err[Omega_m]  =',om_accuracy
  print*,'Err[Omega_b]    =',ob_accuracy
  print*, '============================================================='
  return
END SUBROUTINE add_cmb
!========================================================
!==============! Functions ==================================
DOUBLE PRECISION FUNCTION h2(redshift)
  USE cosmo
  ! h2(z) = Omega_matter(1+z)^3+Omega_lambda
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: fz
    fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  h2 =  ( (om0+ob0)*(1d0+redshift)**3d0+ ok0 * (1d0+redshift)**2d0+ode0* fz)
  return
END FUNCTION h2
DOUBLE PRECISION FUNCTION func0(redshift)
  USE cosmo
  ! func0(z) = 1/[h2(z)]^0.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2
  external :: h2
  func0 = 1d0/dsqrt(h2(redshift))
  return
END FUNCTION func0
DOUBLE PRECISION FUNCTION func1(redshift)
  USE cosmo
  ! func1(z) = ln(1+z)/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func1 = dlog(1d0+redshift) * fz /h2(redshift)**1.5d0
  return
END FUNCTION func1
DOUBLE PRECISION FUNCTION func2(redshift)
  USE cosmo
  ! func2(z) = (1+z)^2/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2,fz
  external :: h2
   fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func2 =( (1d0+redshift)**2d0 - fz)/h2(redshift)**1.5d0
  return
END FUNCTION func2
DOUBLE PRECISION FUNCTION func3(redshift)
  USE cosmo
  ! func3(z) = (1+z)^3/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
   fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func3 =( (1d0+redshift)**3d0 - fz) /h2(redshift)**1.5d0
  return
END FUNCTION func3
DOUBLE PRECISION FUNCTION func4(redshift)
  USE cosmo
  ! func4(z) = ln(1+z) - z/1+z /[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func4 = fz * (dlog(1d0+redshift) - (redshift/(1d0+ redshift)))/h2(redshift)**1.5d0
  return
END FUNCTION func4
