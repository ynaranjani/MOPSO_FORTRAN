SUBROUTINE mutate(pop, pops, NP, gen, gens, mutRate, lb, ub)
	IMPLICIT NONE
	INTEGER, INTENT(IN):: NP, gen, pops, gens
	REAL, INTENT(IN):: mutRate, lb(NP), ub(NP)
	REAL, INTENT(INOUT):: pop(pops*NP)
	
	INTEGER:: i, dm
	REAL::rang, chance, inf, sup
	REAL::RAND
	INTEGER:: IRAND
	
	DO i = 1, pops
		chance = (1.0 - gen/(gens*mutRate))**1.5
		IF (RAND(0) <= chance) THEN
			dm = MOD(IRAND(0),NP)+1
			rang = (ub(dm)-lb(dm))*chance/2
			inf  = pop((i-1)*NP+dm) - rang
			IF (inf < lb(dm)) inf = lb(dm)
			sup = pop((i-1)*NP+dm) + rang
			IF (sup > ub(dm)) sup = ub(dm)
			pop((i-1)*NP+dm) = RAND(0)*(sup-inf)+sup
		END IF
	END DO
	
END SUBROUTINE mutate
