
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE evaluate(fit, viol, indv, NF, NP, length, funCnt, f)
!------------------------------------------------------------------
!	evaluate: calculates the objective values corresponding to each individual
!	Input:
!			indv	the individual
!			NF		dim of fitness space
!			NP		dim of parameter space
!			length	# of individuales to be evaluated in the batch
!			f		the objective function (solver) of the form:
!					f(fit, indv, viol, NF, NP)
!	Output:
!			fit		the fitness values
!			viol	# of constraint violations for each indv
!	Input/Output:
!			funCnt	the number of function evaluations
!	By: Yousef Naranjani-Apr 09 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP, length
	INTEGER, INTENT(INOUT)::funCnt
	REAL, INTENT(IN):: indv(NP*length)
	REAL, INTENT(OUT):: fit(NF*length)
	INTEGER, INTENT(OUT)::viol(length)
!	POINTER, INTENT(IN)::f
	
	INTEGER:: i
	
	DO i = 1, length ! Recommended for parallelization for big solvers
		CALL f(fit((i-1)*NF+1:i*NF), indv((i-1)*NP+1:i*NP), viol(i), NF, NP)
	END DO
	funCnt = funCnt + length
END SUBROUTINE evaluate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE xtoz(z, x, h, lb, ub, nd)
!------------------------------------------------------------------
!	xtoz converts x coordinates to z coordinates on the discretized space
!	Input:
!			x		the coordinates of the point to be converted
!			h		the space grid size
!			lb		the lower bound of the space
!			ub		the upper bound of the space
!			nd		the dimension of apace
!	Output:
!			z		the z coordinates of the point
!	By: Yousef Naranjani-Mar 30 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN)::nd
	REAL, DIMENSION (1 : nd), INTENT(IN):: x, h, lb, ub
	INTEGER, DIMENSION (1 : nd), INTENT(OUT):: z
	INTEGER:: i
	DO i = 1, nd
		IF ((x(i) < lb(i)).OR.(x(i) > ub(i))) THEN
!			ERROR STOP "xtoz: the x coordinate was out of bound"
			z = -1
			EXIT
		ELSE IF (x(i) == ub(i)) THEN
			z(i) = FLOOR((ub(i)-lb(i))/h(i))
		ELSE
			z(i) = FLOOR((x(i)-lb(i))/h(i)) + 1
		END IF
	END DO
END SUBROUTINE xtoz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION ztocell(z, N, nd)
!------------------------------------------------------------------
!	ztocell converts z coordinates to cell number
!	Input:
!			z		z coordinate of the cell
!			N		the partitioning of spcae
!			nd		the dimension of space
!	Output:
!			cell	the cell number
!	By: Yousef Naranjani-Mar 30 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN)::nd
	INTEGER, DIMENSION (1 : nd), INTENT(IN):: z, N
	INTEGER:: cell
	INTEGER:: i, b
	cell = z(1)
	b = N(1)
	DO i = 2, nd
		cell = cell + (z(i)-1) * b
		b = b * N(i)
	END DO
	ztocell = cell
	RETURN
END FUNCTION ztocell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION domChkSingle(x,vx,y,vy,NF)
!------------------------------------------------------------------
!	domChkSingle performs a dominancy test
!	Input:
!			x,y		Objective function values for two points (f(x),f(y))
!			vx,vy	constraint violations of x and y
!			NF		the length of x
!	Output:
!			-1 if the input1 is dominated by input2: f(x)>f(y) (or y has less violations)
!			0 if the are non-dominant (and equal violations)
!			1 if input one dominates input2 f(x)<f(y) (or x has less violations)
!			2 if equal  
!	By: Yousef Naranjani-Mar 30 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER::NF, vx, vy
	REAL::x(NF), y(NF)
	LOGICAL::dominant, dominated, equal
	INTEGER i
	dominant = .TRUE.
	dominated = .TRUE.
	equal = .TRUE.
	IF (vx > vy) THEN
!		input 2 has violated less constraints (better)
		domChkSingle = -1
	ELSE IF (vx < vy) THEN
!		input 1 has violated less constraints (better)
		domChkSingle = 1
	ELSE
		DO i = 1, NF
			IF (y(i)-x(i)>1e-5) THEN
				dominated = .FALSE.
				equal = .FALSE.
			ELSE IF (x(i)-y(i)>1e-5) THEN
				dominant = .FALSE.
				equal = .FALSE.
			END IF
		END DO
		IF ((dominated.EQV. .TRUE.).AND.(equal.EQV. .FALSE.)) THEN
!			the cell is dominated by input2 (input2 is better)
			domChkSingle = -1
		ELSE IF ((dominant.EQV. .TRUE.).AND.(equal.EQV. .FALSE.)) THEN
!			input1 is dominant to input2 (input1 is better)
			domChkSingle = 1
		ELSE IF (equal.EQV. .TRUE.) THEN
			domChkSingle = 2
		ELSE
			domChkSingle = 0
		END IF
	END IF
	RETURN
END FUNCTION domChkSingle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE pBestUpdate(popPbest, fitPbest, violPbest, &
						pop, fit, violations, NP, NF, pops)
!------------------------------------------------------------------
!	pBestUpdate updates the particles best position of all time
!
!
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN):: NP, NF, pops
	REAL, INTENT(IN):: pop(pops*NP), fit(pops*NF)
	INTEGER, INTENT(IN):: violations(pops)
	REAL, INTENT(INOUT):: popPbest(pops*NP), fitPbest(pops*NF)
	INTEGER, INTENT(INOUT):: violPbest(pops)
	
	INTEGER:: i, domChkSingle, domVal
	REAL:: RAND
	DO i = 1, pops
		domVal = domChkSingle(fit((i-1)*NF+1:i*NF),violations(i),fitPbest((i-1)*NF+1:i*NF),violPbest(i),NF)
		IF (domVal == 1) THEN ! current value is dominant so the Pbest should be updated
			popPbest((i-1)*NP+1:i*NP) = pop((i-1)*NP+1:i*NP)
			fitPbest((i-1)*NF+1:i*NF) = fit((i-1)*NF+1:i*NF)
			violPbest(i) = violations(i)
		ELSE IF (domVal == 0) THEN ! if they are nondominant we would exchange them with %50 probability
			IF (RAND(0) > 0.5) THEN
				popPbest((i-1)*NP+1:i*NP) = pop((i-1)*NP+1:i*NP)
				fitPbest((i-1)*NF+1:i*NF) = fit((i-1)*NF+1:i*NF)
				violPbest(i) = violations(i)
			END IF
		END IF
	END DO
END SUBROUTINE pBestUpdate	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gBestUpdate(popArchive, fitArchive, violArchive, partPosArchive, &
		indArchive, pop, fit, violations, &
		NP, NF, pops, maxArchive, &
		hyperSpace, hyperLen, lbFit, ubFit, N)
!------------------------------------------------------------------
!	gBestUpdate inserts the population of current cycle (pop, fit, violations)
!		into archive while keeping archive always checked with dominancy.
!	
!	Input and output:
!			popArchive		the archived global best population
!			fitArchive		the corresponding global best fitness
!			violArchive		the corresponding constraint violations
!			partPosArchive	the hyperSpace location of each archive point
!			indArchive		length of archive
!	Input:
!			NP			# of parameters
!			NF			# of fitness functions
!			pops		# of individuals in each population
!			maxArchive	archive size (should be bigger than pops)
!			pop			current population
!			fit			current fitness
!			violations	current violations
!			
!	By: Yousef Naranjani-Apr 02 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN):: NP, NF, pops, maxArchive, hyperLen
	REAL, INTENT(IN):: pop(pops*NP), fit(pops*NF)
	INTEGER, INTENT(IN):: violations(pops), N(NF)
	INTEGER, INTENT(INOUT):: violArchive(maxArchive), indArchive,&
							 partPosArchive(maxArchive)
	REAL, INTENT(INOUT):: popArchive(maxArchive*NP), fitArchive(maxArchive*NF)
	INTEGER*2, INTENT(INOUT)::hyperSpace(hyperLen)
	REAL, INTENT(INOUT):: lbFit(NF), ubFit(NF)

	LOGICAL:: flag(pops), flag2(maxArchive), updateFlag
	REAL:: popTmp(pops*NP), fitTmp(pops*NF)
	INTEGER:: violTmp(pops), lenTmp, partPosTmp(pops), indTmpStart
	INTEGER::i, j, k, dom, domChkSingle

	popTmp = pop
	fitTmp = fit
	violTmp = violations
	flag = .FALSE.
	lenTmp = 0
	DO i = 1, pops-1
		IF (flag(i).EQV. .TRUE.) THEN
			CONTINUE
		ELSE	! if we have not decided to remove it yet
			DO j = i+1, pops
				dom = domChkSingle(fit((i-1)*NF+1:i*NF),violations(i), &
					fit((j-1)*NF+1:j*NF),violations(j),NF)
				IF (dom == 1) THEN	! i dominates j
					flag(j) = .TRUE.
				ELSE IF (dom == -1) THEN
					flag(i) = .TRUE.
				END IF
			END DO
		END IF
	END DO
!	Eliminating the dominated ones out of the list	
	CALL shrink(popTmp, fitTmp, violTmp, lenTmp, pop, fit, violations, &
				pops, flag, NP, NF, pops)
				
!	Comparing the shrinked list with archive
	flag = .FALSE.
	flag2 = .FALSE.
	DO i = 1, indArchive 	! on archive
		DO j = 1, lenTmp	! on list
			dom = domChkSingle(fitArchive((i-1)*NF+1:i*NF),violArchive(i), &
				fitTmp((j-1)*NF+1:j*NF),violTmp(j),NF)
			IF (dom == 1) THEN	! i dominates j
				flag(j) = .TRUE.
			ELSE IF (dom == -1) THEN
				flag2(i) = .TRUE.
			END IF
		END DO
	END DO
	!debugging
	k = 0
	DO i = 1, indArchive
		IF(flag2(i).EQV..TRUE.) k = k + 1
	END DO
!	Clean up the dominated ones from archive and the list
	CALL shrink(popTmp, fitTmp, violTmp, lenTmp, popTmp, fitTmp, violTmp, &
				lenTmp, flag, NP, NF, pops)

	CALL shrinkArchive(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
				flag2, NP, NF, maxArchive, hyperSpace, hyperLen)

!	initialize the lbFit and ubFit if they are empty
	IF ((lbFit(1)==0.0).AND.(ubFit(1)==0.0).AND.(lenTmp>0)) THEN
		DO i = 1, NF
			lbFit(i) = fitTmp(i)
			lbFit(i) = fitTmp(i)
		END DO
	END IF
	IF (lenTmp > 0) THEN ! update might be necessary
		CALL makeHyper(partPosTmp, fitTmp, lenTmp, lbFit, ubFit, updateFlag, NF, N)
		!	Updating the archive
		IF (updateFlag.EQV. .TRUE.) THEN
			CALL updateHyper(hyperSpace, partPosArchive, fitArchive, indArchive,&
							hyperLen, lbFit, ubFit,  NF, N)
		END IF
	
		IF (indArchive+lenTmp < maxArchive) THEN	! the archive has plenty of vacancies
			popArchive(indArchive*NP+1:(indArchive+lenTmp)*NP) = popTmp(1:lenTmp*NP)
			fitArchive(indArchive*NF+1:(indArchive+lenTmp)*NF) = fitTmp(1:lenTmp*NF)
			violArchive(indArchive+1 : indArchive+lenTmp) = violTmp(1:lenTmp)
			partPosArchive(indArchive+1 : indArchive+lenTmp) = partPosTmp(1:lenTmp)
			!	update the hyperSpace
			DO i = 1, lenTmp
				hyperSpace(partPosTmp(i)) = hyperSpace(partPosTmp(i)) + 1
			END DO
			
			indArchive = indArchive+lenTmp
		ELSE ! 2 cases: 1)there are some space avalable 2)archive is full 
			IF (indArchive < maxArchive) THEN
			 	! 1)there are some apaces, lets fill the vacant parts first
			
				popArchive(indArchive*NP+1:(maxArchive)*NP) = popTmp(1:(maxArchive-indArchive)*NP)
				fitArchive(indArchive*NF+1:(maxArchive)*NF) = fitTmp(1:(maxArchive-indArchive)*NF)
				violArchive(indArchive+1 : maxArchive) = violTmp(1:(maxArchive-indArchive))
				partPosArchive(indArchive+1 : maxArchive) = partPosTmp(1:(maxArchive-indArchive))
				!	update the hyperSpace
				DO i = 1, (maxArchive-indArchive)
					hyperSpace(partPosTmp(i)) = hyperSpace(partPosTmp(i)) + 1
				END DO
				
				indTmpStart = maxArchive-indArchive
				indArchive = maxArchive
			ELSE ! 2) archive is fill
				indTmpStart = 1
			END IF
			! points from populated areas should be replaced
			DO i = indTmpStart, lenTmp
				DO j = 1, maxArchive ! that is the total length of hyperCubeCells
					IF (hyperSpace(partPosArchive(j)) > hyperSpace(partPosTmp(i))) THEN
						! They have to be replaced
						popArchive((j-1)*NP+1 : j*NP) = popTmp((i-1)*NP+1 : i*NP)
						fitArchive((j-1)*NF+1 : j*NF) = fitTmp((i-1)*NF+1 : i*NF)
						violArchive(j) = violTmp(i)
						hyperSpace(partPosTmp(i)) = hyperSpace(partPosTmp(i)) + 1
						hyperSpace(partPosArchive(j)) = hyperSpace(partPosArchive(j)) - 1
						partPosArchive(j) = partPosTmp(i)
						
					END IF
				END DO
			END DO
		END IF
	END IF
END SUBROUTINE gBestUpdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE writeTOfile(fileID, info, dims, length)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::fileID, dims, length
	REAL, INTENT(IN):: info(dims*length)
	
	INTEGER:: i, j
	
	DO i = 1, length
		DO j = 1, dims
			WRITE(fileID,"(E14.6E3,1X)",advance='no') info((i-1)*dims+j)
		END DO
		IF (i.NE.length) WRITE(fileID,*)
	END DO
	
END SUBROUTINE writeTOfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE writeTOfileProgress(fileID, popArchive, fitArchive, &
								violArchive, gen, length, NP, NF)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::fileID, NP, NF, length, gen
	INTEGER, INTENT(IN)::violArchive(length)
	REAL, INTENT(IN):: popArchive(NP*length), fitArchive(NF*length)
	
	INTEGER:: i, j
	
	DO i = 1, length
		WRITE(fileID,"(I4,1X)",advance='no') gen
		DO j = 1, NP
			WRITE(fileID,"(E14.6E3,1X)",advance='no') popArchive((i-1)*NP+j)
		END DO
		DO j = 1, NF
			WRITE(fileID,"(E14.6E3,1X)",advance='no') fitArchive((i-1)*NF+j)
		END DO
		WRITE(fileID,"(I3,1X)") violArchive(i)
	END DO
	
END SUBROUTINE writeTOfileProgress

