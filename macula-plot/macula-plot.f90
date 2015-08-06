!===============================================================================
! macula-plot
!===============================================================================

PROGRAM maculaplot
use maculamod
implicit none

 ! Variables you need to choose...
 INTEGER, PARAMETER :: Ns = 3
 INTEGER, PARAMETER :: mmax = 1
 REAL(8), PARAMETER :: Rstar = 1.0D0
 REAL(8), PARAMETER :: Vkappa = 2000.0D0
 ! Variables you don't need to change...
 INTEGER :: nplen
 REAL(8), DIMENSION(:), ALLOCATABLE :: tp, filling, Fmod, Fdot
 REAL(8), DIMENSION(:), ALLOCATABLE :: deltaratio
 REAL(8), DIMENSION(mmax) :: T_start, T_end
 REAL(8), DIMENSION(12) :: Thetastar
 REAL(8), DIMENSION(1,mmax) :: Thetainst
 REAL(8), DIMENSION(8,Ns) :: Thetaspot, Thetaspot_temp
 INTEGER :: j, k, i
 REAL(8) :: tfirst, tlast, cadence
 REAL(8), DIMENSION(4) :: ranvec
 REAL(8), DIMENSION(:,:), ALLOCATABLE :: Fmod_array, Fdot_array, deltaratio_array
 REAL(8), DIMENSION(:,:), ALLOCATABLE :: RVrot, RVcon, RVtot, cosalpha
 REAL(8), PARAMETER :: deg = 1.745329251994330D-2
 REAL(8), PARAMETER :: RSun = 6.955D8

 ! times to simulate
 tfirst = 0.0D0; tlast = 1000.0D0; cadence = 0.1D0

 ! start_of_declaration
 ! Star parameters
 Thetastar(1)  = 65.0D0*deg		! stellar inclination angle
 Thetastar(2)  = 27.0D0			! P_EQ (equitorial rotation period)
 Thetastar(3)  = 0.20D0			! kappa2 (differential rotation coeff)
 Thetastar(4)  = 0.0D0			! kappa4 (differential rotation coeff)
 Thetastar(5)  = 0.0D0			! c1 (star's limb darkening)
 Thetastar(6)  = 1.49070480644882844D0	! c2 (star's limb darkening)
 Thetastar(7)  = 0.0D0			! c3 (star's limb darkening)
 Thetastar(8)  = -0.67591259874875889D0	! c4 (star's limb darkening)
 Thetastar(9)  = Thetastar(5)		! d1
 Thetastar(10) = Thetastar(6)		! d2
 Thetastar(11) = Thetastar(7)		! d3
 Thetastar(12) = Thetastar(8)		! d4 

 ! Instrument parameters: blend factors
 Thetainst(1,:) = 1.0D0   ! B(j)

 ! Spot parameters
 DO k=1,Ns
   call random_flat(4,ranvec)
   ! parameters for the k^th spot...
   Thetaspot(1,k) = ranvec(1)*360.0D0*deg		! longitude
   Thetaspot(2,k) = (ranvec(2)-0.5D0)*180.0D0*deg	! latitude
   Thetaspot(3,k) = 0.05D0+ranvec(3)*0.10D0		! alpha_max
   Thetaspot(4,k) = 0.25D0 + 0.75D0*ranvec(4)		! fspot
   Thetaspot(5,k) = 500.0D0				! tmax
   Thetaspot(6,k) = 2000.0D0				! lifetime
   Thetaspot(7,k) = 0.02D0				! ingress
   Thetaspot(8,k) = Thetaspot(7,k)			! egress
 END DO
 ! end_of_declarations

 ! generate times
 nplen = 1 + NINT( (tlast - tfirst)/cadence )
 ALLOCATE(tp(nplen))
 ALLOCATE(Fmod(nplen))
 ALLOCATE(Fdot(nplen))
 ALLOCATE(deltaratio(nplen))
 ALLOCATE(cosalpha(Ns,nplen))
 DO i=1,nplen
  tp(i) = tfirst + (i-1)*cadence
 END DO
 T_start(1) = tfirst-cadence; T_end(1) = tlast+cadence

 ! call macula
 call macula(tp,&
             nplen,Ns,mmax,.FALSE.,.FALSE.,.TRUE.,.TRUE.,&
             Thetastar,Thetaspot,Thetainst,&
             T_start,T_end,Fmod,Fdot,cosalpha,deltaratio)

 ! basic information
 open(unit=12,file='macula-basics.dat')
 write(12,*) Thetastar(1),Thetaspot(4,:)
 close(12)

 ALLOCATE(filling(nplen))
 ALLOCATE(Fmod_array(Ns,nplen))
 ALLOCATE(Fdot_array(Ns,nplen))
 ALLOCATE(deltaratio_array(Ns,nplen))
 ALLOCATE(RVrot(Ns,nplen))
 ALLOCATE(RVcon(Ns,nplen))
 ALLOCATE(RVtot(Ns,nplen))

 DO k=1,Ns
   Thetaspot_temp = Thetaspot		! copy the original array
   Thetaspot_temp(4,:) = 1.0D0		! but set all contrasts to unity
   Thetaspot_temp(4,k) = Thetaspot(4,k)	! except for the k^th spot
   call macula(tp,&
               nplen,Ns,mmax,.TRUE.,.TRUE.,.TRUE.,.FALSE.,&
               Thetastar,Thetaspot_temp,Thetainst,&
               T_start,T_end,Fmod,Fdot,cosalpha,deltaratio)
   Fmod_array(k,:) = Fmod(:)
   Fdot_array(k,:) = Fdot(:)
   deltaratio_array(k,:) = deltaratio(:)
   DO i=1,nplen
     filling(i) = 2.0D0*(1.0D0-Thetaspot(4,k))*(1.0D0-cosalpha(k,i))
     RVrot(k,i) = -(1.0D0-Fmod(i))*Fdot(i)*Rstar*RSun/filling(i)/86400.0D0
     RVcon(k,i) = -(1.0D0-Fmod(i))**2*Vkappa/filling(i)
     RVtot(k,i) = RVrot(k,i) + RVcon(k,i)
   END DO
 END DO

 ! sum the effects
 open(unit=14,file='macula-RV.dat')
 DO i = 1,nplen
   write(14,*) tp(i),RVtot(:,i)
 END DO
 close(14)
 
 ! sum the effects
 open(unit=15,file='macula-lc_components.dat')
 DO i = 1,nplen
   write(15,*) tp(i),Fmod_array(:,i)
 END DO
 close(15)

 ! sum the effects
 open(unit=15,file='macula-TdeltaV.dat')
 DO i = 1,nplen
   write(15,*) tp(i),deltaratio_array(:,i)
 END DO
 close(15)

CONTAINS

! =============================================================
 SUBROUTINE random_flat(n,r)

 implicit none

 INTEGER :: i, n
 REAL(8) :: seeda, r(n)

 call random_seed()
 DO i=1,n
   call random_number(seeda)
   r(i)=seeda
 END DO

 END SUBROUTINE random_flat
! =============================================================

! =============================================================
 SUBROUTINE random_norm(n,r)

 implicit none

 INTEGER :: i, n
 REAL(8) :: seeda, seedb, r(n)
 REAL(8), PARAMETER :: pi = 3.141592653589793D0

 call random_seed()

 DO i=1,n
   call random_number(seeda)
   call random_number(seedb)
   r(i) = SQRT(-2.0D0*LOG(seeda))*DCOS(2.0D0*pi*seedb)
 END DO

 END SUBROUTINE random_norm
! =============================================================

END PROGRAM maculaplot
