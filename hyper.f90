!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE makeHyper(partPos, fit, length, lb, ub, updateFlag, NF, N)
!------------------------------------------------------------------
!	makeHyper calculated the hyperCube location for the inputs
!	
!	Input and output:
!			lb, ub			lower and upper bound on fit space
!	Input:
!			fit				input points 
!			NF				# of fitness functions
!			N				divisions in fit space
!			length			# of individuals in fit
!	Output:
!			partPos			hyperLoc of points in fit
!			updateFlag		indicates if the bounds have changed or not
!
!	By: Yousef Naranjani-Apr 04 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	LOGICAL, INTENT(INOUT)::updateFlag
	INTEGER, INTENT(IN)::NF
	INTEGER, INTENT(IN)::length, N(NF)
	INTEGER, INTENT(INOUT):: partPos(length)
	REAL, INTENT(IN)::fit(length*NF)
	REAL, INTENT(INOUT):: lb(NF), ub(NF)
	
	INTEGER:: i, j
	INTEGER:: z(NF), ztocell
	REAL::h(NF)
	updateFlag = .FALSE.
	DO i = 1, length
		DO j = 1, NF
			IF (fit((i-1)*NF+j) < lb(j)) THEN
				updateFlag = .TRUE. ! archive should be updated according to new bound
				lb(j) = fit((i-1)*NF+j)
			ELSE IF (fit((i-1)*NF+j) > ub(j)) THEN
				updateFlag = .TRUE.
				ub(j) = fit((i-1)*NF+j)
			END IF
		END DO
	END DO
	
	DO i = 1, NF
		 h(i) = (ub(i)-lb(i)) / N(i)
	END DO
	
!	Now the bounds are uptodate, we can find the hypercubes, parallelizable loop
	DO i = 1, length
		CALL xtoz(z, fit((i-1)*NF+1:i*NF), h, lb, ub, NF)
		IF (z(1) == -1) THEN
			ERROR STOP "makeHyper: the x coordinate was out of bound"
		ELSE
			partPos(i) = ztocell(z, N, NF)
		END IF
	END DO
	
END SUBROUTINE makeHyper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE updateHyper(hyperSpace, partPosArchive, fitArchive, indArchive,&
						hyperLen, lbFit, ubFit,  NF, N)
!------------------------------------------------------------------
!	updateHyper updates the hyperCubes based on new ubFit and lbFit
!	
!	Input:
!			fitArchive		archive fitnesses 
!			NF				# of fitness functions
!			N				divisions in fit space
!			indArchive		last filled index of archive		
!			hyperLen		length of hyperSpace
!			lbFit, ubFit	lower and upper bound on fit space
!	Output:
!			hyperSpace		hyperSpace array (counts of points in each cubic)
!			partPosArchive	each archive particle's hyperCube number	
!
!	By: Yousef Naranjani-Apr 04 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, hyperLen, indArchive
	INTEGER, INTENT(IN)::N(NF)
	INTEGER, INTENT(INOUT)::partPosArchive(indArchive) ! or maxArchive
	INTEGER*2, INTENT(INOUT)::hyperSpace(hyperLen)
	REAL, INTENT(IN)::fitArchive(indArchive*NF)
	REAL, INTENT(INOUT):: lbFit(NF), ubFit(NF)
	
	INTEGER:: i, j
	INTEGER:: z(NF), ztocell
	REAL::h(NF)
	DO i = 1, NF
		 h(i) = (ubFit(i)-lbFit(i)) / N(i)
	END DO
	
	hyperSpace = 0 ! might not be efficiet (check if it is time consuming)
	
	DO i = 1, indArchive
		CALL xtoz(z, fitArchive((i-1)*NF+1:i*NF), h, lbFit, ubFit, NF)
		IF (z(1) == -1) THEN
			ERROR STOP "makeHyper: the x coordinate was out of bound"
		ELSE
			partPosArchive(i) = ztocell(z, N, NF)
			hyperSpace(partPosArchive(i)) = hyperSpace(partPosArchive(i)) + 1
		END IF
	END DO
	
	
END SUBROUTINE updateHyper
