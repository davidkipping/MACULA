MODULE like

use params
use maculamod

implicit none
      
contains
      
!=======================================================================

SUBROUTINE slikelihood(R,slhood)
         
 implicit none
      
 REAL(8), DIMENSION(nest_nPar) :: R
 REAL(8) :: slhood, loglike
 INTEGER :: i
 REAL(8), PARAMETER :: roottwo = 1.414213562373095D0

 ! Scaling of the parameters from hypercube
 DO i = 1, sdim
   IF( Rflag(i) .EQ. 0 ) THEN
     ! Delta function: Rsol = centering
     R(i) = Rsol(i)
   ELSE IF( Rflag(i) .EQ. 2 ) THEN
     ! Uniform: Rmax = max, Rmin = min
     R(i) = Rmin(i) + (Rmax(i)-Rmin(i))*R(i)
   ELSE IF( Rflag(i) .EQ. 3 ) THEN
     ! Gaussian: Rmax = mean, Rmin = stdev
     R(i) = Rmax(i) + roottwo*Rmin(i)*inverf(-1.0D0+2.0D0*R(i))
   ELSE IF( Rflag(i) .EQ. 4 ) THEN
     ! Jeffrey's: Rmax = max, Rmin = min
     R(i) = ( Rmax(i)**R(i) )*( Rmin(i)**(1.0D0-R(i)) )
   ELSE IF( Rflag(i) .EQ. 5 ) THEN
     ! Modified Jefrey's: Rmax = max, Rmin = inflection point
     R(i) = -( Rmin(i)**(1.0D0-R(i)) )*( Rmin(i)**R(i) - ( Rmin(i)+Rmax(i) )**R(i) )
   END IF
 END DO

 ! call models to get likelihood
 call models(R,loglike)
 slhood = loglike

END SUBROUTINE slikelihood
!=======================================================================

! ======================================================================
SUBROUTINE models(Rvec,loglike)

!implicit none
 REAL(8), DIMENSION(nest_nPar), INTENT(IN) :: Rvec   ! Fitted-parameter vector
 REAL(8), DIMENSION(12) :: Thetastar
 REAL(8), DIMENSION(1,mmax) :: Thetainst
 REAL(8), DIMENSION(8,Ns) :: Thetaspot
 REAL(8) :: loglike
 INTEGER :: j, k

 ! Default model shown here has:
 ! - infinitely long lived spots
 ! - quadratic limb darkening for star and spot with same and fixed coefficients
 ! - simple differential rotation law, with kappa_4=0 for simplicicity
 ! - all spots have common filling factor, fspot
 ! - each spot has unique longitude, latitude and size

 ! start_of_declarations
 ! Star parameters
 Thetastar(1)  = DACOS(Rvec(1))		! stellar inclination angle
 Thetastar(2)  = Rvec(2)		! P_EQ (equitorial rotation period)
 Thetastar(3)  = Rvec(3)		! kappa2 (differential rotation coeff)
 Thetastar(4)  = 0.0D0			! kappa4 (differential rotation coeff)
 Thetastar(5)  = 0.0D0			! c1 (star's limb darkening)
 Thetastar(6)  = 1.49070480644882844D0	! c2 (star's limb darkening)
 Thetastar(7)  = 0.0D0			! c3 (star's limb darkening)
 Thetastar(8)  = -0.67591259874875889D0	! c4 (star's limb darkening)
 Thetastar(9)  = Thetastar(5)		! d1 (spot's limb darkening)
 Thetastar(10) = Thetastar(6)		! d2 (spot's limb darkening)
 Thetastar(11) = Thetastar(7)		! d3 (spot's limb darkening)
 Thetastar(12) = Thetastar(8)		! d4 (spot's limb darkening)

 ! Instrument parameters
 Thetainst(1,:) = 1.0D0   ! Blend factor of each quarter

 ! Spot parameters
 DO k=1,Ns
   ! parameters for the k^th spot...
   Thetaspot(1,k) = Rvec(starparams+1+(k-1)*spotparams+1)     ! longitude
   Thetaspot(2,k) = Rvec(starparams+1+(k-1)*spotparams+2)     ! latitude
   Thetaspot(3,k) = Rvec(starparams+1+(k-1)*spotparams+3)     ! alpha_max
   Thetaspot(4,k) = Rvec(4)                                   ! fspot
   Thetaspot(5,k) = T_start(1)+0.5D0*(T_end(mmax)-T_start(1)) ! tmax
   Thetaspot(6,k) = 2.0D0*(T_end(mmax)-T_start(1))            ! lifetime
   Thetaspot(7,k) = 0.02D0                                    ! ingress
   Thetaspot(8,k) = Thetaspot(7,k)                            ! egress
 END DO
 ! end_of_declarations

 ! call macula to get likelihood
 call macula(tp,fp,sigfp,fpwei,sigfpwei,mval,&
             nplen,Ns,mmax,.FALSE.,.FALSE.,.FALSE.,&
             Thetastar,Thetaspot,Thetainst,&
             loglike)

END SUBROUTINE models
! ======================================================================

! ======================================================================
FUNCTION inverf(x)

 implicit none

 REAL(8) :: x
 REAL(8), PARAMETER :: awil = 0.14001228868666646D0
 REAL(8), PARAMETER :: bwil = 4.546884979448289D0
 REAL(8) :: factor, xsq, inverf

 IF( x .LT. 0.0D0 ) THEN
  factor = -1.0D0
 ELSE
  factor = 1.0D0
 END IF

 xsq = 1.0D0 - x**2
 x = bwil + 0.5D0*DLOG(xsq)
 x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
 inverf = factor*DSQRT(x)

END FUNCTION
! ======================================================================

END MODULE like
