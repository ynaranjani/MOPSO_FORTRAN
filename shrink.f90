!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE shrink(list1Out, list2Out, list3OutInt, listLenOut, list1In, list2In, list3InInt, listLenIn, flags, nd1, nd2, listLen)
!------------------------------------------------------------------
!	shrink removes certain flagged indices from 3 arrays and returns a shorter
!	array (the 3 arrays could be: the variables, function evals, and the 
!			contraint violations for a certain point).
!	Input:
!			list1In		1st input list to be shrinked
!			list2In		2nd input list to be shrinked
!			list3InInt	3rd input list (integer)
!			listLenIn	# of list entries (total entries)
!			flags		flags for removing the members (T=remove, F=keep)
!			nd1			the dimention of the array 1
!			nd2			the dimention of the array 2
!			listLen		total list capacity (total length = listLen*(nd1 or nd2)) 
!	Output:
!			list1Out	updated list1 (supposed to be shorter)
!			list2Out	updated list2 (supposed to be shorter)
!			list3OutInt updated list3 (integer)
!			listLenOut	length of the output list (total legth = listLenOut*(nd1 or nd2))
!			
!	By: Yousef Naranjani-Apr 02 2016
!------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER, INTENT(IN):: listLen, listLenIn, nd1, nd2
	LOGICAL, INTENT(IN):: flags(listLen)
	REAL, INTENT(IN):: list1In(listLen*nd1), list2In(listLen*nd2)
	INTEGER, INTENT(IN):: list3InInt(listLen)
	INTEGER, INTENT(OUT):: listLenOut, list3OutInt(listLen)
	REAL, INTENT(OUT):: list1Out(listLen*nd1),  list2Out(listLen*nd2)
	
	INTEGER::i, flagInd, loopVal
	
	flagInd = 0
	loopVal = listLenIn
	listLenOut = 0
	DO i = 1 , loopVal
		IF (flags(i).EQV. .FALSE.) THEN ! we want to keep it
			flagInd = flagInd + 1
			list1Out(nd1*(flagInd-1)+1 : nd1*flagInd) = list1In(nd1*(i-1)+1 : nd1*i)
			list2Out(nd2*(flagInd-1)+1 : nd2*flagInd) = list2In(nd2*(i-1)+1 : nd2*i)
			list3OutInt(flagInd) = list3InInt(i)
		END IF
	END DO
	listLenOut = flagInd
	list1Out(nd1*flagInd+1 : nd1 * listLen) = 0.0
	list2Out(nd2*flagInd+1 : nd2 * listLen) = 0.0
	list3OutInt(flagInd + 1 : listLen) = 0
END SUBROUTINE shrink
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE shrinkArchive(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
						flags, NP, NF, maxArchive, hyperSpace, hyperLen)
!------------------------------------------------------------------
!	shrinkArchive removes certain flagged indices from 3 arrays and returns a shorter
!	array (the 3 arrays could be: the variables, function evals, and the 
!			contraint violations for a certain point).
!	Input:
!			list1In		1st input list to be shrinked
!			list2In		2nd input list to be shrinked
!			list3InInt	3rd input list (integer)
!			list4InInt	4th input list (integer)
!			listLenIn	# of list entries (total entries)
!			flags		flags for removing the members (T=remove, F=keep)
!			nd1			the dimention of the array 1
!			nd2			the dimention of the array 2
!			listLen		total list capacity (total length = listLen*(nd1 or nd2)) 
!	Output:
!			list1Out	updated list1 (supposed to be shorter)
!			list2Out	updated list2 (supposed to be shorter)
!			list3OutInt updated list3 (integer)
!			list4OutInt updated list4 (integer)
!			listLenOut	length of the output list (total legth = listLenOut*(nd1 or nd2))
!			
!	By: Yousef Naranjani-Apr 04 2016
!------------------------------------------------------------------
	IMPLICIT NONE

	INTEGER, INTENT(IN):: maxArchive, hyperLen
	INTEGER, INTENT(INOUT):: indArchive, NP, NF
	LOGICAL, INTENT(INOUT):: flags(maxArchive)
	REAL, INTENT(INOUT):: popArchive(maxArchive*NP), fitArchive(maxArchive*NF)
	INTEGER, INTENT(INOUT):: violArchive(maxArchive), partPosArchive(maxArchive)
	INTEGER*2, INTENT(INOUT)::hyperSpace(hyperLen)

	
	INTEGER::i, flagInd
	
	flagInd = 0
	DO i = 1 , indArchive
		IF (flags(i).EQV. .FALSE.) THEN ! we want to keep it
			flagInd = flagInd + 1
			popArchive(NP*(flagInd-1)+1 : NP*flagInd) = popArchive(NP*(i-1)+1 : NP*i)
			fitArchive(NF*(flagInd-1)+1 : NF*flagInd) = fitArchive(NF*(i-1)+1 : NF*i)
			violArchive(flagInd) = violArchive(i)
			partPosArchive(flagInd) = partPosArchive(i)
		ELSE ! it is removed. Updating the hyper space:
			hyperSpace(partPosArchive(i)) = hyperSpace(partPosArchive(i))-1
			IF (hyperSpace(partPosArchive(i))==-1) ERROR STOP "shrinkArchive: hyperSpace value cannot be -1"
		END IF
	END DO

	indArchive = flagInd
	popArchive(NP*flagInd+1 : NP * maxArchive) = 0.0
	fitArchive(NF*flagInd+1 : NF * maxArchive) = 0.0
	violArchive(flagInd + 1 : maxArchive) = 0
	partPosArchive(flagInd + 1 : maxArchive) = 0
END SUBROUTINE shrinkArchive
