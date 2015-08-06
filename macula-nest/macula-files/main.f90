PROGRAM main

 use params
 use nestwrapper
      
 implicit none
	
 INTEGER :: i, j, m

 ! # of parameters to wrap around
 nest_pWrap(:) = 0
 DO j=1,Ns
   ! longitude of j^th spot
   nest_pWrap(starparams+1+(j-1)*spotparams+1) = 1
 END DO

 ! read-in priors
 open(1,FILE='macula-files/priors.txt',FORM='FORMATTED',STATUS='UNKNOWN')
 DO j=1,starparams+spotparams+1
   read(1,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
 END DO
 ! if mulitple spots, then set priors of the other spots to be the same as those 
 ! for the first spot
 IF( Ns .GE. 2 ) THEN
 DO j=2,Ns
   Rsol(starparams+1+(j-1)*spotparams+1) = Rsol(starparams+1+1)
   Rsol(starparams+1+(j-1)*spotparams+2) = Rsol(starparams+1+2)
   Rsol(starparams+1+(j-1)*spotparams+3) = Rsol(starparams+1+3)
   Rdel(starparams+1+(j-1)*spotparams+1) = Rdel(starparams+1+1)
   Rdel(starparams+1+(j-1)*spotparams+2) = Rdel(starparams+1+2)
   Rdel(starparams+1+(j-1)*spotparams+3) = Rdel(starparams+1+3)
   Rmin(starparams+1+(j-1)*spotparams+1) = Rmin(starparams+1+1)
   Rmin(starparams+1+(j-1)*spotparams+2) = Rmin(starparams+1+2)
   Rmin(starparams+1+(j-1)*spotparams+3) = Rmin(starparams+1+3)
   Rmax(starparams+1+(j-1)*spotparams+1) = Rmax(starparams+1+1)
   Rmax(starparams+1+(j-1)*spotparams+2) = Rmax(starparams+1+2)
   Rmax(starparams+1+(j-1)*spotparams+3) = Rmax(starparams+1+3)
   Rflag(starparams+1+(j-1)*spotparams+1) = Rflag(starparams+1+1)
   Rflag(starparams+1+(j-1)*spotparams+2) = Rflag(starparams+1+2)
   Rflag(starparams+1+(j-1)*spotparams+3) = Rflag(starparams+1+3)
 END DO
 END IF
 close(1)
        
 ! read-in light curve data
 write(*,*) 'Reading in primary data...'
 open(unit=30301,file='macula-files/data.dat')
 DO i=1,nplen
   read(30301,*) tp(i),fp(i),sigfp(i)
   fpwei(i) = fp(i)/(sigfp(i)**2)
   sigfpwei(i) = 1.0D0/(sigfp(i)**2)
 END DO
 close(30301)

 ! assign Tstart/Tend points
 DO m=1,mmax
   T_start(m) = QXstart( Qused(m) )
   T_end(m) = QXend( Qused(m) )
 END DO

 ! box-car function (labelled as Pi_m in the paper)
 DO i=1,nplen
   DO m=1,mmax
     IF( tp(i) .GT. T_start(m) .AND. tp(i) .LT. T_end(m) ) THEN
       mval(i) = m
     END IF
   END DO
 END DO

 call nest_Sample
 
END
