! Include file for example MultiNest program 'Gaussian' (see arXiv:1001.0719)

module params
implicit none

! ==============================================================================
! Model Parameters

 ! parameters to be explored by MaculaNest
 INTEGER, PARAMETER :: Ns = 1		! # of star spots
 INTEGER, PARAMETER :: starparams = 3	! # of star parameters
 INTEGER, PARAMETER :: spotparams = 3	! # of spot parameters
 INTEGER, PARAMETER :: sdim = starparams+Ns*spotparams+1

 ! data.dat controls
 INTEGER, PARAMETER :: nplen = 3386
 REAL(8), PARAMETER :: facP = 1.0D0

 ! quarters used
 INTEGER, PARAMETER :: mmax = 1		! # of quarters/instruments
 INTEGER, DIMENSION(mmax), PARAMETER :: Qused = (/ 5 /)

! ==============================================================================
! Model Variables

 ! input data
 REAL(8), DIMENSION(nplen) :: tp, fp, sigfp, fpmod, fpwei, sigfpwei
 INTEGER, DIMENSION(nplen) :: mval ! quarter number of each datum

 ! hypercube terms
 REAL(8), DIMENSION(sdim) :: Rmin, Rmax, Rsol, Rdel
 INTEGER, DIMENSION(sdim) :: Rflag

 ! Tstart(m) = Time stamp @ start of m^th data series 
 ! Tend(m)   = Time stamp @ end of m^th data series
 REAL(8), DIMENSION(mmax) :: T_start, T_end

! ==============================================================================
! Useful constants

 ! constants
 REAL(8), PARAMETER :: radian = 0.017453292519943295D0
 REAL(8), PARAMETER :: third = 0.333333333333333333D0
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: Grv = 6.67384D-11

 ! Quarter-to-quarter start/end points
 REAL(8), DIMENSION(17), PARAMETER :: QXstart = (/ 54953.52862037D0, &	! Q1
                                                55002.50923614D0, &	! Q2
						55093.21414610D0, &	! Q3
						55185.36716152D0, &	! Q4
						55276.48039581D0, &	! Q5
						55372.43908262D0, &	! Q6
						55463.16447070D0, &	! Q7
						55568.35386174D0, &	! Q8
						55641.50607599D0, &	! Q9
						55739.83509146D0, &	! Q10
						55834.19777420D0, &	! Q11
						55932.39901141D0, &	! Q12
						56015.72702732D0, &	! Q13
						56107.12890281D0, &	! Q14
						56206.47750326D0, &	! Q15
						56305.08740097D0, &	! Q16
						56390.97969746D0 /)	! Q17

 REAL(8), DIMENSION(17), PARAMETER :: QXend = (/ 54997.99334075D0, &	! Q1
						55091.47732973D0, &	! Q2
						55182.50650574D0, &	! Q3
						55275.21348059D0, &	! Q4
						55371.17221239D0, &	! Q5
						55462.30628212D0, &	! Q6
						55552.55903682D0, &	! Q7
						55635.35547074D0, &	! Q8
						55738.93604628D0, &	! Q9
						55833.27831544D0, &	! Q10
						55931.33647373D0, &	! Q11
						56015.03229749D0, &	! Q12
						56106.06638110D0, &	! Q13
						56204.33204943D0, &	! Q14
						56304.14747395D0, &	! Q15
						56390.96969746D0, &	! Q16
						66390.96969746D0 /)	! Q17

! ==============================================================================

! Parameters for Nested Sampler
	
      	!whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.false.)
	
      	!sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.false.)
	
      	!max no. of live points
      	integer nest_nlive
	parameter(nest_nlive=4000)
      
      	!tot no. of parameters, should be sdim in most cases but if you need to
      	!store some additional parameters with the actual parameters then
      	!you need to pass them through the likelihood routine
	integer nest_nPar 
	parameter(nest_nPar=sdim)
      
      	!seed for nested sampler, -ve means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-1)
      
      	!evidence tolerance factor
      	double precision nest_tol 
      	parameter(nest_tol=1.0D0)
      
      	!enlargement factor reduction parameter
      	double precision nest_efr
      	parameter(nest_efr=0.01D0)
      
      	!root for saving posterior files
      	character*100 nest_root
	parameter(nest_root='chains/macula-nest-')
	
	!after how many iterations feedback is required & the output files should be updated
	!note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=1000)
	
	!null evidence (set it to very high negative no. if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.0D5)
      
      	!max modes expected, for memory allocation
      	integer nest_maxModes 
      	parameter(nest_maxModes=100)
      
      	!no. of parameters to cluster (for mode detection)
      	integer nest_nClsPar
      	parameter(nest_nClsPar=sdim)
      
      	!whether to resume from a previous run
      	logical nest_resume
      	parameter(nest_resume=.true.)
      
      	!whether to write output files
      	logical nest_outfile
      	parameter(nest_outfile=.true.)
      
      	!initialize MPI routines?, relevant only if compiling with MPI
	!set it to F if you want your main program to handle MPI initialization
      	logical nest_initMPI
      	parameter(nest_initMPI=.true.)
      
      	!points with loglike < nest_logZero will be ignored by MultiNest
      	double precision nest_logZero
      	parameter(nest_logZero=-huge(1.0D0))

      	!max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
	!has done max no. of iterations or convergence criterion (defined through nest_tol) has been satisfied
      	integer nest_maxIter
      	parameter(nest_maxIter=0)
	
	!parameters to wrap around (0 is F & non-zero T)
	integer nest_pWrap(sdim)
	
      	!feedback on the sampling progress?
      	logical nest_fb 
      	parameter(nest_fb=.true.)
!=======================================================================

end module params
