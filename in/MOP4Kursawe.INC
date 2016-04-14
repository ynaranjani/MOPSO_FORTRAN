! Setup your problem here:
! The following parameters must be set by the user:
!	f			(~line 83) the solver name (refer to ProblemBank.f90)
!	NF			# of objectives (fitness functions)
!	NP			# of parameters on design space (variables)
!	pops		# of individuals in each generation
!	gens		# of generations or flight cycles
!	maxArchive	the capacity of archive (must be greater than pops)
!	N			array of the divisions of hyperspace in objective space 
!				(recommended: 30 for each objective)
!	ub			upper bound on design space
!	lb			lower bound on design space
!	probName	the name of the problem!
!	generateReport	
!				a flag indicating if you need a report file to be generated
!	intermResults	
!				a flag indicating if you need to collect intermediate results
!	intermResultsInterval	
!				the generations interval for the intermediate results to be collected	


!		The archive related parameters
!	maxArchive	the capacity of archive
!	indArchive	index of the last filled spot on archive
!	popArchive	archive population
!	fitArchive	archive corresponding fitness values
!	violArchive	constraint violations for each member of archive
!	partPosArchive	hypercube location of archive members
!	ubFit		upper bound on objective space
!	lbFit		lower bound on objective space
!	h			the division size (ubFit-lbFit)/N
!	hyperLen	length of the hyperspace

INTEGER, PARAMETER:: 	maxArchive = 1000, &
						gens = 1000, &	
						pops = 100, &
						NP = 3, &
						NF = 2
INTEGER, PARAMETER, DIMENSION(NF):: N = (/100, 100/)
REAL::					lb(NP) = (/-5.0	,-5.0, -5.0/), &
						ub(NP) = (/5.0, 5.0, 5.0/)
LOGICAL, PARAMETER::	generateReport = .FALSE., &
						intermResults = .TRUE.
INTEGER, PARAMETER:: 	intermResultsInterval = 100
CHARACTER(len = 30), PARAMETER::	probName = "MOP4Kursawe                   "

! No need to change anything beyond this point except the SOLVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL::					lbFit(NF) = 0.0, &
						ubFit(NF) = 0.0
INTEGER::				indArchive = 0, &
						hyperLen = 1
INTEGER*2, DIMENSION(:), ALLOCATABLE :: hyperSpace
REAL:: 					mutRate = 0.05

REAL::					popArchive(NP*maxArchive) = 0.0, &
						fitArchive(NF*maxArchive) = 0.0, &
						pop(NP*pops) = 0.0, &
						fit(NF*pops) = 0.0
INTEGER:: 				violArchive(maxArchive) = 0, &
						violations(pops) = 0, &
						partPosArchive(maxArchive) = 0
						
REAL::					popPbest(NP*pops) = 0.0, &
						fitPbest(NF*pops) = 0.0
INTEGER::				violPbest(pops) = 0

REAL::					velocity(NP*pops) = 0.0
INTEGER::				lead(pops)
REAL::					t00, t01, t10, t11, t20, t21 		! the timers		
ABSTRACT INTERFACE
	SUBROUTINE objective_function(fit, indv, viol, NF, NP)
		IMPLICIT NONE
		INTEGER, INTENT(IN)::NF, NP
		REAL, INTENT(IN):: indv(NP)
		REAL, INTENT(OUT):: fit(NF)
		INTEGER, INTENT(OUT)::viol
	END SUBROUTINE objective_function
END INTERFACE				

PROCEDURE (objective_function), pointer :: f => null()
!!!!!!!!!!!!!!!!!!!!!!!!!! SET THE SOLVER NAME HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
	f => MOP4_Kursawe

	DO i = 1 , NF
		hyperLen = hyperLen * N(i)
	END DO
	ALLOCATE(hyperSpace(hyperLen))
	hyperSpace = 0

