!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rouletteWheelSelection(lead, hyperSpace, hyperLen, indArchive, pops,&
				partPosArchive)
!------------------------------------------------------------------
!	rouletteWheelSelection selects a member of the archive as the best leader
!		for updating the velocity of the particles
!	Input:
!			hyperSpace
!			hyperLen
!			indArchive		
!			pops	
!			partPosArchive		
!	Output:
!			lead	the lead points chosen from archive for velocity update		
!	By: Yousef Naranjani-Apr 10 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN):: hyperLen, indArchive, pops
	INTEGER*2, INTENT(IN):: hyperSpace(hyperLen)
	INTEGER, INTENT(IN):: partPosArchive(indArchive)
	
	INTEGER, INTENT(INOUT):: lead(pops)
	
	INTEGER:: cubeList(indArchive), cube, which
	REAL:: cubeListCumSumProb(indArchive)
	INTEGER:: i, j, c
	INTEGER:: findArchive
	REAL:: RAND, cubeChance
	INTEGER:: IRAND
	
	c = 0
	cubeList = 0
	cubeListCumSumProb = 0.0
!	calculating the cummulative sum probability for the cubes to be chosen	
	DO i = 1, hyperLen
		IF (hyperSpace(i).NE.0) THEN
			c = c + 1
			cubeList(c) = i
			cubeListCumSumProb(c) = 1.0 / hyperSpace(i) 
			IF (c>1) THEN
				cubeListCumSumProb(c) = cubeListCumSumProb(c) + cubeListCumSumProb(c-1)
			END IF
			!print*,'c=',c,'cube:', i, 'Probability:', cubeListCumSumProb(c), 'hyperSpace:',hyperSpace(i)
		END IF
	END DO
	
	IF (c==0) ERROR STOP "rouletteWheelSelection: is hyperSpace empty?"
	
	!print*,'the indArchive  = ', indArchive
	!DO i = 1, indArchive
	!	print *, i, partPosArchive(i)
	!END DO
	
	DO i = 1, pops !paralellizable loop
		cubeChance = RAND(0)*cubeListCumSumProb(c)
!		finding the corresponding cube
		DO j = 1, c
			IF (cubeChance <= cubeListCumSumProb(j)) EXIT
		END DO
		cube = cubeList(j)
!		now that the cube is found, one of the points in that cube should be randomly chosen
!		which = FLOOR(RAND(0)*hyperSpace(cube)+1.0)
		which = MOD(IRAND(0),INT(hyperSpace(cube)))+1
		!print*,'we should find cube=',cube,'which=',which
		lead(i) = findArchive(partPosArchive, cube, which, indArchive)
		!print*, 'lead is updated!, index = ', lead(i)
	END DO
	
END SUBROUTINE rouletteWheelSelection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION findArchive(partPosArchive, cube, which, indArchive)
	IMPLICIT NONE
	INTEGER, INTENT(IN):: indArchive, cube, which
	INTEGER, INTENT(IN):: partPosArchive(indArchive)
	
	INTEGER:: i, c
	
	
	c = 0
	DO i = 1, indArchive
		!print*,partPosArchive(i)
		IF (partPosArchive(i)==cube) c = c + 1
		IF (c == which) EXIT
	END DO
	!print*,'cube=',cube,'which=',which, 'index=', i , ',c=',c
	IF (c.NE.which) ERROR STOP "findArchive: the cube cannot be found in archive"
	findArchive = i
	RETURN
END FUNCTION findArchive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE velocityUpdate(velocity, pop, popPbest, popArchive, indArchive, lead, pops, NP)
!------------------------------------------------------------------
!	velocityUpdate updates the particle velocity.		
!	By: Yousef Naranjani-Apr 10 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN)::pops, NP, indArchive
	INTEGER, INTENT(IN)::lead(pops)
	REAL, INTENT(IN):: pop(pops*NP), popPbest(pops*NP), popArchive(indArchive*NP)
	REAL, INTENT(INOUT)::velocity(NP * pops)
	
	INTEGER:: i, j
	REAL::RAND
	
	DO i = 1, pops
		DO j = 1, NP
			velocity((i-1)*NP+j) = 0.4 * velocity((i-1)*NP+j) + &
								RAND(0) * (popPbest((i-1)*NP+j)-pop((i-1)*NP+j)) + &
								RAND(0) * (popArchive((lead(i)-1)*NP+j)-pop((i-1)*NP+j))
		END DO
	END DO
	
END SUBROUTINE velocityUpdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE keepIn(pop, velocity, lb, ub, NP, pops)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NP, pops
	REAL, INTENT(IN):: lb(NP), ub(NP)
	REAL, INTENT(INOUT):: pop(NP*pops), velocity(NP*pops)
	
	INTEGER:: i, j
	
	DO i = 1, pops
		DO j = 1, NP
			IF (pop((i-1)*NP+j) > ub(j)) THEN
				pop((i-1)*NP+j) = ub(j)
				velocity((i-1)*NP+j) = - velocity((i-1)*NP+j)
			ELSE IF (pop((i-1)*NP+j) < lb(j)) THEN
				pop((i-1)*NP+j) = lb(j)
				velocity((i-1)*NP+j) = - velocity((i-1)*NP+j)
			END IF
		END DO
	END DO
	
END SUBROUTINE keepIn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







