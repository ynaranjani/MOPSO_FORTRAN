!	Copyright (c) 2016 Yousef Naranjani
!	Email: yo.na1118@gmail.com
!	Please read the Readme.txt file to get started
PROGRAM MOPSOtools
	USE ProblemBank
	
	IMPLICIT NONE
	
	INTEGER:: i, j, gen, seed, time(8), funCnt = 0
	INTEGER::fpLog, fpPF, fpPS, fpProgress	! output file references
	REAL::RAND	
	
	!!!!!!!!!!!!!!! Choose the proper problem settings here  !!!!!!!!!!!!!!!!!!!
	!INCLUDE 'in/MOP1circles2D.INC'
	INCLUDE 'in/MOP2Fonseca.INC'
	!INCLUDE 'in/MOP3Poloni.INC'
	!INCLUDE 'in/MOP4Kursawe.INC'
	!INCLUDE 'in/MOP5Viennet.INC'
	!INCLUDE 'in/ZDT3_10D.INC'
	!INCLUDE 'in/MOP_C4_Tanaka.INC'
	
	PRINT*,'Starting up!'
	! IDs for the output files
	fpLog = 10
	fpPF = 11
	fpPS = 12
	fpProgress = 13
	
	CALL CPU_TIME(t00)
	! for OpenMP: seconds = omp_get_wtime ( )
	
!	CALL execute_command_line('mkdir -p out/' // trim( probName ) )
	CALL system('mkdir -p out/' // trim( probName ))
	IF (generateReport) OPEN(UNIT = fpLog, 		FILE = ('out/'//trim( probName )//'/log.dat'))
	OPEN(UNIT = fpPF, 		FILE = ('out/'//trim( probName )//'/PF.dat'))
	OPEN(UNIT = fpPS, 		FILE = ('out/'//trim( probName )//'/PS.dat'))
	IF (intermResults) OPEN(UNIT = fpProgress,	FILE = ('out/'//trim( probName )//'/progress.dat'))

	IF (generateReport) CALL DATE_AND_TIME(values=time)
	seed = 1000*time(7)+time(8)
	CALL SRAND(seed)

!	initialize and evaluate the population
	DO i = 1, pops
		DO j = 1, NP
			pop((i-1)*NP+j) = RAND(0) * (ub(j)-lb(j)) + lb(j)
		ENDDO
	ENDDO
!	function evaluations
	CALL evaluate(fit, violations, pop, NF, NP, pops, funCnt, f)

!	initializing the archive
	CALL gBestUpdate(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
		pop, fit, violations, NP, NF, pops, maxArchive, &
		hyperSpace, hyperLen, lbFit, ubFit, N)
!	initializing the local best
	popPbest(1:NP*pops) = pop(1:NP*pops)
	fitPbest(1:NF*pops) = fit(1:NF*pops)
	violPbest(1:pops) = violations(1:pops)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  THE FLIGHT CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO gen = 1, gens
	
		IF (generateReport) THEN
			WRITE(fpLog,"('Generation: ',I5,' began. Archive length: ', I5)") gen,indArchive
		END IF
	
		! choosing the leaders from Archive for velocity update	
		CALL rouletteWheelSelection(lead, hyperSpace, hyperLen, indArchive, pops,&
					partPosArchive)
		
		! update the velocity here
		CALL velocityUpdate(velocity, pop, popPbest, popArchive, maxArchive, lead, pops, NP)
		
		! Updating the population
		pop(1:NP*pops) = pop(1:NP*pops) + velocity(1:NP*pops)
		
		! Correcting the individuals that are out of bound	
		CALL keepIn(pop, velocity, lb, ub, NP, pops)
	
		! include mutation operator here
		CALL mutate(pop, pops, NP, gen, gens, mutRate, lb, ub)
		
		IF (generateReport) CALL CPU_TIME(t20)	
		!	function evaluations
		CALL evaluate(fit, violations, pop, NF, NP, pops, funCnt, f)
		IF (generateReport) THEN
			CALL CPU_TIME(t21)
			WRITE(fpLog,"('Function Eval: ',F10.2,' s')") t21-t20
		END IF
		
		!	updating the archive
		
		CALL gBestUpdate(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
			pop, fit, violations, NP, NF, pops, maxArchive, &
			hyperSpace, hyperLen, lbFit, ubFit, N)
				
		!	updating the local best
		CALL pBestUpdate(popPbest, fitPbest, violPbest, pop, fit, violations, NP, NF, pops)
		
		!	write to progress
		IF ((intermResults.EQV. .TRUE.).AND.(MOD(gen,intermResultsInterval)==0)) THEN
			CALL writeTOfileProgress(fpProgress,	popArchive, fitArchive, &
									violArchive, gen, indArchive, NP, NF)
		END IF
		
		IF (generateReport) THEN
			CALL CPU_TIME(t11)
			WRITE(fpLog,"('Generation: ',I5,' finished in ',F10.2,' s')") gen,t11-t10
		END IF
		
	END DO
	
	CALL CPU_TIME(t01)
	IF (generateReport) WRITE(fpLog,"('Total time: ',F10.2,' s')")t01-t00
	
	CALL writeTOfile(fpPF, fitArchive, NF, indArchive)
	CALL writeTOfile(fpPS, popArchive, NP, indArchive)
	
	IF (generateReport) CLOSE(fpLog)
	CLOSE(fpPF)
	CLOSE(fpPS)
	IF (intermResults) CLOSE(fpProgress)
	
	PRINT*,'Done! Time=',t01-t00,' ,',indArchive, ' solutions found'
	
END PROGRAM MOPSOtools
