
C --------------------------------------------------------------------
      SUBROUTINE SUBCONDPROB(UBAR,CONDPROB,SNNDPROB,NOS)
C     former:      PROGRAM SNND_RIS
C --------------------------------------------------------------------
C
C  SUBROUTINE TO CALCULATE THE CONDITIONAL PROBABILITES FOR THE
C  2NND_RIS MODEL FOR PE (NOS=1) AND POE (NOS=3)
C 
C  FOR DETAILS SEE PROGRAM SNND_RIS

      INTEGER   NU
      PARAMETER (NU=9)

      INTEGER   I, II, J, O, P, K
      INTEGER   NOS, INDEX, IFIX, POWER

      REAL      DUMMY
      REAL      NORMFACTOR, NORMFACTORPROD
      REAL      ZC
      REAL      SNNDPROB(4,4)
      REAL      CONDPROB(4,4,4)
      REAL      UTEMP(NU), UTEMP2(NU)
      REAL      U1(NU), UN(NU)
      REAL      PCELEMENT(NU,NU)
      REAL      QCELEMENT(NU,NU)
      REAL      UZETAETA(NU,NU)
      REAL      UBAR(NU,NU,4)

C  UBAR  : COARSE GRAINED STATISTICAL WEIGHT MATRIX WITH DIMENSIONS
C          (9X9X4) FOR THE FOUR DIFFERENT LOCAL CHEMISTRY 


C  `````````` BORDER MATRICES U1  AND  UN ```````
      DO I = 1, 9
        U1(I) = 0.0
        UN(I) = 1.0
      ENDDO
      U1(1)  = 1.0
C  `````````` BORDER MATRICES U1  AND  UN  ```````

C **********CALCULATIONS FOR POLYETHYLENE ***************
      IF(NOS.EQ.1) THEN

C  `````````` PARTITION FUNCTION OF MODEL C :  - ZC - ```````````
C
C                - 28
C    ZC = U(1) . U   . U(n)
C
C      -------------------------------------
       CALL SUBPARTFUNCT(U1,UBAR,UN,NOS,28,ZC)
C      -------------------------------------

C  `````````` PARTITION FUNCTION OF MODEL C :  - ZC - ```````````

C  ``````` PROBABILITIES OF TWO BONDS PCELEMENT(ZETA,ETA,30) ``````````````
      DO O = 1, 9
      DO P = 1, 9
        NORMFACTORPROD = 1.0

C  FIRST PART:  U(1) * UBAR^13
C       ....................................
        CALL SUBVECMAT2(U1, UBAR, UTEMP, 9, 9, 1)
C       ````````````````````````````````````
        NORMFACTOR = MAX(UTEMP(1),UTEMP(2),UTEMP(3),
     >                   UTEMP(4),UTEMP(5),UTEMP(6),
     >                   UTEMP(7),UTEMP(8),UTEMP(9))
        NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
        DO I = 1, 9
          UTEMP(I) = UTEMP(I) / NORMFACTOR
        ENDDO
        DO I = 1, 12
C         ..........................................
          CALL SUBVECMAT2(UTEMP, UBAR, UTEMP2, 9, 9, 1)
C         ``````````````````````````````````````````
          NORMFACTOR = MAX(UTEMP2(1),UTEMP2(2),UTEMP2(3),
     >                     UTEMP2(4),UTEMP2(5),UTEMP2(6),
     >                     UTEMP2(7),UTEMP2(8),UTEMP2(9))
          NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
          DO J = 1, 9
            UTEMP(J) = UTEMP2(J) / NORMFACTOR
          ENDDO
        ENDDO
C  FIRST PART:  U(1) * UBAR^13

C  INSERT THE MATRIX WITH THE SPECIFIED STATE
        DO I = 1, 9
          DO J = 1, 9
            UZETAETA(J,I) = 0.0
          ENDDO
        ENDDO
        UZETAETA(P,O) = UBAR(P,O,1)
C       .............................................
        CALL SUBVECMAT(UTEMP, UZETAETA, UTEMP2, 9, 9)
C       `````````````````````````````````````````````
        NORMFACTOR = MAX(UTEMP2(1),UTEMP2(2),UTEMP2(3),
     >                   UTEMP2(4),UTEMP2(5),UTEMP2(6),
     >                   UTEMP2(7),UTEMP2(8),UTEMP2(9))
        NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
        DO J = 1, 9
          UTEMP(J) = UTEMP2(J) / NORMFACTOR
        ENDDO
C  INSERT THE MATRIX WITH THE SPECIFIED STATE

C  LAST PART:  UTEMP * U^14 * U(N)
        DO I = 1, 14
C         ..........................................
          CALL SUBVECMAT2(UTEMP, UBAR, UTEMP2, 9, 9, 1)
C         ``````````````````````````````````````````
          NORMFACTOR = MAX(UTEMP2(1),UTEMP2(2),UTEMP2(3),
     >                     UTEMP2(4),UTEMP2(5),UTEMP2(6),
     >                     UTEMP2(7),UTEMP2(8),UTEMP2(9))
          NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
          DO J = 1, 9
            UTEMP(J) = UTEMP2(J) / NORMFACTOR
          ENDDO
        ENDDO
        PCELEMENT(P,O) = UTEMP(1)*UN(1) + UTEMP(2)*UN(2) +
     >                   UTEMP(3)*UN(3) +
     >                   UTEMP(4)*UN(4) + UTEMP(5)*UN(5) +
     >                   UTEMP(6)*UN(6) +
     >                   UTEMP(7)*UN(7) + UTEMP(8)*UN(8) +
     >                   UTEMP(9)*UN(9)
C  LAST PART:  UTEMP * U^14 * U(N)

C  NORMALIZE AND DIVIDE BY THE PARTITION FUNCTION
        PCELEMENT(P,O) = PCELEMENT(P,O) * NORMFACTORPROD
        PCELEMENT(P,O) = PCELEMENT(P,O) / ZC
C  NORMALIZE AND DIVIDE BY THE PARTITION FUNCTION

      ENDDO
      ENDDO
C  ``````` PROBABILITIES OF TWO BONDS PCELEMENT(ZETA,ETA,30) ``````````````

C  ````````` CONDITIONAL PROBABILITIES QCELEMENT(ETA,30) ``````````````````
      DO I = 1, 9
        DUMMY = 0.0
        DO J = 1, 9
          DUMMY = DUMMY + PCELEMENT(I,J)
        ENDDO
        DO J = 1, 9
          QCELEMENT(I,J) = PCELEMENT(I,J) / DUMMY
        ENDDO
      ENDDO
C  ````````` CONDITIONAL PROBABILITIES QCELEMENT(ETA,30) ``````````````````


   12 FORMAT(4F16.11)
      WRITE(7,*) ' '
      WRITE(7,*) '##############  POLYETHYLENE ##############'
      WRITE(7,*) ' '
      WRITE(7,*) '****** BOND PROBABILITIES ******'
      WRITE(7,*) ' '
      DUMMY = PCELEMENT(1,1) + PCELEMENT(1,2) + PCELEMENT(1,6) +
     >        PCELEMENT(1,8) +
     >        PCELEMENT(2,1) + PCELEMENT(2,2) + PCELEMENT(2,6) +
     >        PCELEMENT(2,8) +
     >        PCELEMENT(6,1) + PCELEMENT(6,2) + PCELEMENT(6,6) +
     >        PCELEMENT(6,8) +
     >        PCELEMENT(8,1) + PCELEMENT(8,2) + PCELEMENT(8,6) +
     >        PCELEMENT(8,8)
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITHOUT NORMALIZATION:'
      DO I = 1, 4
        II = I
        IF( I .EQ. 3 ) THEN
          II = 6
        ENDIF
        IF( I .EQ. 4 ) THEN
          II = 8
        ENDIF
        WRITE(7,12) PCELEMENT(II,1), PCELEMENT(II,2),
     >              PCELEMENT(II,6), PCELEMENT(II,8)
      ENDDO

      DO I = 1, 4
        SNNDPROB(I,1) = 0.0
      ENDDO

      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITH NORMALIZATION:'
      DO I = 1, 4
        II = I
        IF( I .EQ. 3 ) THEN
          II = 6
        ENDIF
        IF( I .EQ. 4 ) THEN
          II = 8
        ENDIF
        WRITE(7,12) PCELEMENT(II,1)/DUMMY, PCELEMENT(II,2)/DUMMY,
     >              PCELEMENT(II,6)/DUMMY, PCELEMENT(II,8)/DUMMY
        SNNDPROB(1,1) = SNNDPROB(1,1) + PCELEMENT(II,1)/DUMMY
        SNNDPROB(2,1) = SNNDPROB(2,1) + PCELEMENT(II,2)/DUMMY
        SNNDPROB(3,1) = SNNDPROB(3,1) + PCELEMENT(II,6)/DUMMY
        SNNDPROB(4,1) = SNNDPROB(4,1) + PCELEMENT(II,8)/DUMMY
      ENDDO

      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) '1ST VECTOR:'
      WRITE(7,12) SNNDPROB(1,1), SNNDPROB(2,1), 
     >		  SNNDPROB(3,1), SNNDPROB(4,1)

      WRITE(7,*) ' '
      WRITE(7,*) ' '
      WRITE(7,*) '****** CONDITIONAL PROBABILITIES ******'
      WRITE(7,*) ' '
      CONDPROB(1,1,1) = QCELEMENT(1,1)
      CONDPROB(1,2,1) = QCELEMENT(1,2)
      CONDPROB(1,3,1) = QCELEMENT(1,6)
      CONDPROB(1,4,1) = QCELEMENT(1,8)
      CONDPROB(2,1,1) = QCELEMENT(2,1)
      CONDPROB(2,2,1) = QCELEMENT(2,2)
      CONDPROB(2,3,1) = QCELEMENT(2,6)
      CONDPROB(2,4,1) = QCELEMENT(2,8)
      CONDPROB(3,1,1) = QCELEMENT(6,1)
      CONDPROB(3,2,1) = QCELEMENT(6,2)
      CONDPROB(3,3,1) = QCELEMENT(6,6)
      CONDPROB(3,4,1) = QCELEMENT(6,8)
      CONDPROB(4,1,1) = QCELEMENT(8,1)
      CONDPROB(4,2,1) = QCELEMENT(8,2)
      CONDPROB(4,3,1) = QCELEMENT(8,6)
      CONDPROB(4,4,1) = QCELEMENT(8,8)
C  THE CONDITIONAL PROBABILITIES HAVE TO BE RESCALED BECAUSE
C  THEY ARE A SELECTION CONDPROB(4X4) OF A BIGGER MATRIX 
C  QCELEMENT(9X9)
      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITHOUT NORMALIZATION:'
      DO I = 1, 4
        WRITE(7,12) CONDPROB(I,1,1), CONDPROB(I,2,1),
     >              CONDPROB(I,3,1), CONDPROB(I,4,1)
      ENDDO
      DO II = 1, 4
        DUMMY = CONDPROB(II,1,1) + CONDPROB(II,2,1) +
     >          CONDPROB(II,3,1) + CONDPROB(II,4,1)
        DO I = 1, 4
          CONDPROB(II,I,1) = CONDPROB(II,I,1) / DUMMY
        ENDDO
      ENDDO
      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITH NORMALIZATION:'
      DO I = 1, 4
        WRITE(7,12) CONDPROB(I,1,1), CONDPROB(I,2,1),
     >              CONDPROB(I,3,1), CONDPROB(I,4,1)
      ENDDO
      WRITE(7,*) ' '
      WRITE(7,*)  'sum of first row condprob  ',CONDPROB(1,1,1) +
     >          CONDPROB(1,2,1) + CONDPROB(1,3,1) + CONDPROB(1,4,1)
      WRITE(7,*)  'sum of second row condprob ',CONDPROB(2,1,1) +
     >          CONDPROB(2,2,1) + CONDPROB(2,3,1) + CONDPROB(2,4,1)
      WRITE(7,*)  'sum of third row condprob  ',CONDPROB(3,1,1) +
     >          CONDPROB(3,2,1) + CONDPROB(3,3,1) + CONDPROB(3,4,1)
      WRITE(7,*)  'sum of fourth row condprob ',CONDPROB(4,1,1) +
     >          CONDPROB(4,2,1) + CONDPROB(4,3,1) + CONDPROB(4,4,1)

      ENDIF
C **********CALCULATIONS FOR POLYETHYLENE ***************

C **********CALCULATIONS FOR POLYOXYETHYLENE ***************
      IF(NOS.EQ.3) THEN
      POWER = 13
C  `````````` PARTITION FUNCTION OF MODEL C :  - ZC - ```````````
C
C                - 30
C    ZC = U(1) . U   . U(n)
C
C      -------------------------------------
       CALL SUBPARTFUNCT(U1,UBAR,UN,NOS,30,ZC)
C      -------------------------------------

C  `````````` PARTITION FUNCTION OF MODEL C :  - ZC - ```````````

C  ``````` PROBABILITIES OF TWO BONDS PCELEMENT(ZETA,ETA,30) ``````````````
      DO K = 1, 3
      POWER = POWER + 1
      DO O = 1, 9
      DO P = 1, 9
        NORMFACTORPROD = 1.0

C  FIRST PART:  U(1) * UBAR^14(15,16)
C       ....................................
	INDEX = 2
        CALL SUBVECMAT2(U1, UBAR, UTEMP, 9, 9, INDEX)
C       ````````````````````````````````````
        NORMFACTOR = MAX(UTEMP(1),UTEMP(2),UTEMP(3),
     >                   UTEMP(4),UTEMP(5),UTEMP(6),
     >                   UTEMP(7),UTEMP(8),UTEMP(9))
	IF(NORMFACTOR.GT.1E-24) THEN
        NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
        DO J = 1, 9
          UTEMP(J) = UTEMP(J) / NORMFACTOR
        ENDDO
	ELSE
        DO J = 1, 9
          UTEMP(J) = UTEMP(J)
        ENDDO
	ENDIF
        DO I = 1, POWER-1 
C         ..........................................
	  IF (INDEX.EQ.4) INDEX = 1
	  INDEX = INDEX + 1
          CALL SUBVECMAT2(UTEMP, UBAR, UTEMP2, 9, 9, INDEX)
C         ``````````````````````````````````````````
          NORMFACTOR = MAX(UTEMP2(1),UTEMP2(2),UTEMP2(3),
     >                     UTEMP2(4),UTEMP2(5),UTEMP2(6),
     >                     UTEMP2(7),UTEMP2(8),UTEMP2(9))
	  IF(NORMFACTOR.GT.1E-24) THEN
          NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
          DO J = 1, 9
            UTEMP(J) = UTEMP2(J) / NORMFACTOR
          ENDDO
	  ELSE
          DO J = 1, 9
            UTEMP(J) = UTEMP2(J)
          ENDDO
	  ENDIF
        ENDDO
C  FIRST PART:  U(1) * UBAR^14(15,16)

C  INSERT THE MATRIX WITH THE SPECIFIED STATE
        DO I = 1, 9
          DO J = 1, 9
            UZETAETA(J,I) = 0.0
          ENDDO
        ENDDO
	IF (INDEX.EQ.4) INDEX = 1
	INDEX = INDEX + 1
	IFIX = INDEX
        UZETAETA(P,O) = UBAR(P,O,INDEX)
C       .............................................
        CALL SUBVECMAT(UTEMP, UZETAETA, UTEMP2, 9, 9)
C       `````````````````````````````````````````````
        NORMFACTOR = MAX(UTEMP2(1),UTEMP2(2),UTEMP2(3),
     >                   UTEMP2(4),UTEMP2(5),UTEMP2(6),
     >                   UTEMP2(7),UTEMP2(8),UTEMP2(9))
	IF(NORMFACTOR.GT.1E-24) THEN
        NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
        DO J = 1, 9
          UTEMP(J) = UTEMP2(J) / NORMFACTOR
        ENDDO
	ELSE
        DO J = 1, 9
          UTEMP(J) = UTEMP2(J)
        ENDDO
	ENDIF
C  INSERT THE MATRIX WITH THE SPECIFIED STATE
C  LAST PART:  UTEMP * U^15(14,13) * U(N)
        DO I = 1, 30-POWER-1
C         ..........................................
	  IF (INDEX.EQ.4) INDEX = 1
	  INDEX = INDEX + 1
          CALL SUBVECMAT2(UTEMP, UBAR, UTEMP2, 9, 9, INDEX)
C         ``````````````````````````````````````````
          NORMFACTOR = MAX(UTEMP2(1),UTEMP2(2),UTEMP2(3),
     >                     UTEMP2(4),UTEMP2(5),UTEMP2(6),
     >                     UTEMP2(7),UTEMP2(8),UTEMP2(9))
	  IF(NORMFACTOR.GT.1E-24) THEN
          NORMFACTORPROD = NORMFACTORPROD * NORMFACTOR
          DO J = 1, 9
            UTEMP(J) = UTEMP2(J) / NORMFACTOR
          ENDDO
	  ELSE
          DO J = 1, 9
          UTEMP(J) = UTEMP2(J)
          ENDDO
	  ENDIF
        ENDDO
        PCELEMENT(P,O) = UTEMP(1)*UN(1) + UTEMP(2)*UN(2) +
     >                   UTEMP(3)*UN(3) +
     >                   UTEMP(4)*UN(4) + UTEMP(5)*UN(5) +
     >                   UTEMP(6)*UN(6) +
     >                   UTEMP(7)*UN(7) + UTEMP(8)*UN(8) +
     >                   UTEMP(9)*UN(9)
C  LAST PART:  UTEMP * U^15(14,13) * U(N)

C  NORMALIZE AND DIVIDE BY THE PARTITION FUNCTION
        PCELEMENT(P,O) = PCELEMENT(P,O) * NORMFACTORPROD
        PCELEMENT(P,O) = PCELEMENT(P,O) / ZC
C  NORMALIZE AND DIVIDE BY THE PARTITION FUNCTION

      ENDDO
      ENDDO
C  ``````` PROBABILITIES OF TWO BONDS PCELEMENT(ZETA,ETA,30) ``````````````

C  ````````` CONDITIONAL PROBABILITIES QCELEMENT(ETA,30) ``````````````````
      DO I = 1, 9
        DUMMY = 0.0
        DO J = 1, 9
          DUMMY = DUMMY + PCELEMENT(I,J)
        ENDDO
	IF(DUMMY.GT.1E-24) THEN
        DO J = 1, 9
          QCELEMENT(I,J) = PCELEMENT(I,J) / DUMMY
        ENDDO
	ELSE
        DO J = 1, 9
          QCELEMENT(I,J) = 0. 
        ENDDO
	ENDIF
      ENDDO
C  ````````` CONDITIONAL PROBABILITIES QCELEMENT(ETA,30) ``````````````````

      WRITE(7,*) ' '
      WRITE(7,*) '###########  POLYOXYETHYLENE_',IFIX,' ############'
      WRITE(7,*) ' '
      WRITE(7,*) '****** BOND PROBABILITIES ******'
      WRITE(7,*) ' '
      DUMMY = PCELEMENT(1,1) + PCELEMENT(1,2) + PCELEMENT(1,6) +
     >        PCELEMENT(1,8) +
     >        PCELEMENT(2,1) + PCELEMENT(2,2) + PCELEMENT(2,6) +
     >        PCELEMENT(2,8) +
     >        PCELEMENT(6,1) + PCELEMENT(6,2) + PCELEMENT(6,6) +
     >        PCELEMENT(6,8) +
     >        PCELEMENT(8,1) + PCELEMENT(8,2) + PCELEMENT(8,6) +
     >        PCELEMENT(8,8)
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITHOUT NORMALIZATION:'
      DO I = 1, 4
        II = I
        IF( I .EQ. 3 ) THEN
          II = 6
        ENDIF
        IF( I .EQ. 4 ) THEN
          II = 8
        ENDIF
        WRITE(7,12) PCELEMENT(II,1), PCELEMENT(II,2),
     >              PCELEMENT(II,6), PCELEMENT(II,8)
      ENDDO

      DO I = 1, 4
        SNNDPROB(I,IFIX) = 0.0
      ENDDO

      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITH NORMALIZATION:'
      DO I = 1, 4
        II = I
        IF( I .EQ. 3 ) THEN
          II = 6
        ENDIF
        IF( I .EQ. 4 ) THEN
          II = 8
        ENDIF
        WRITE(7,12) PCELEMENT(II,1)/DUMMY, PCELEMENT(II,2)/DUMMY,
     >              PCELEMENT(II,6)/DUMMY, PCELEMENT(II,8)/DUMMY
        SNNDPROB(1,IFIX) = SNNDPROB(1,IFIX) + PCELEMENT(II,1)/DUMMY
        SNNDPROB(2,IFIX) = SNNDPROB(2,IFIX) + PCELEMENT(II,2)/DUMMY
        SNNDPROB(3,IFIX) = SNNDPROB(3,IFIX) + PCELEMENT(II,6)/DUMMY
        SNNDPROB(4,IFIX) = SNNDPROB(4,IFIX) + PCELEMENT(II,8)/DUMMY
      ENDDO

      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) '1ST VECTOR:'
      WRITE(7,12) SNNDPROB(1,IFIX), SNNDPROB(2,IFIX), 
     >		  SNNDPROB(3,IFIX), SNNDPROB(4,IFIX)

      WRITE(7,*) ' '
      WRITE(7,*) ' '
      WRITE(7,*) '****** CONDITIONAL PROBABILITIES ******'
      WRITE(7,*) ' '
      CONDPROB(1,1,IFIX) = QCELEMENT(1,1)
      CONDPROB(1,2,IFIX) = QCELEMENT(1,2)
      CONDPROB(1,3,IFIX) = QCELEMENT(1,6)
      CONDPROB(1,4,IFIX) = QCELEMENT(1,8)
      CONDPROB(2,1,IFIX) = QCELEMENT(2,1)
      CONDPROB(2,2,IFIX) = QCELEMENT(2,2)
      CONDPROB(2,3,IFIX) = QCELEMENT(2,6)
      CONDPROB(2,4,IFIX) = QCELEMENT(2,8)
      CONDPROB(3,1,IFIX) = QCELEMENT(6,1)
      CONDPROB(3,2,IFIX) = QCELEMENT(6,2)
      CONDPROB(3,3,IFIX) = QCELEMENT(6,6)
      CONDPROB(3,4,IFIX) = QCELEMENT(6,8)
      CONDPROB(4,1,IFIX) = QCELEMENT(8,1)
      CONDPROB(4,2,IFIX) = QCELEMENT(8,2)
      CONDPROB(4,3,IFIX) = QCELEMENT(8,6)
      CONDPROB(4,4,IFIX) = QCELEMENT(8,8)
C  THE CONDITIONAL PROBABILITIES HAVE TO BE RESCALED BECAUSE
C  THEY ARE A SELECTION CONDPROB(4X4) OF A BIGGER MATRIX 
C  QCELEMENT(9X9)
      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITHOUT NORMALIZATION:'
      DO I = 1, 4
        WRITE(7,12) CONDPROB(I,1,IFIX), CONDPROB(I,2,IFIX),
     >              CONDPROB(I,3,IFIX), CONDPROB(I,4,IFIX)
      ENDDO
      DO II = 1, 4
        DUMMY = CONDPROB(II,1,IFIX) + CONDPROB(II,2,IFIX) +
     >          CONDPROB(II,3,IFIX) + CONDPROB(II,4,IFIX)
	IF (DUMMY.GT.1E-24) THEN
        DO I = 1, 4
          CONDPROB(II,I,IFIX) = CONDPROB(II,I,IFIX) / DUMMY
        ENDDO
	ENDIF
      ENDDO
      WRITE(7,*) ' '
      WRITE(7,*) 'tt ; t+,t-,+t,-t ; ++,-- ; +-,-+'
      WRITE(7,*) 'WITH NORMALIZATION:'
      DO I = 1, 4
        WRITE(7,12) CONDPROB(I,1,IFIX), CONDPROB(I,2,IFIX),
     >              CONDPROB(I,3,IFIX), CONDPROB(I,4,IFIX)
      ENDDO
      WRITE(7,*) ' '
      WRITE(7,*)  'sum of first row condprob  ',CONDPROB(1,1,IFIX) +
     >  CONDPROB(1,2,IFIX) + CONDPROB(1,3,IFIX) + CONDPROB(1,4,IFIX)
      WRITE(7,*)  'sum of second row condprob ',CONDPROB(2,1,IFIX) +
     >  CONDPROB(2,2,IFIX) + CONDPROB(2,3,IFIX) + CONDPROB(2,4,IFIX)
      WRITE(7,*)  'sum of third row condprob  ',CONDPROB(3,1,IFIX) +
     >  CONDPROB(3,2,IFIX) + CONDPROB(3,3,IFIX) + CONDPROB(3,4,IFIX)
      WRITE(7,*)  'sum of fourth row condprob ',CONDPROB(4,1,IFIX) +
     >  CONDPROB(4,2,IFIX) + CONDPROB(4,3,IFIX) + CONDPROB(4,4,IFIX)

      ENDDO
      ENDIF
C **********CALCULATIONS FOR POLYOXYETHYLENE ***************

      CLOSE(7)

      RETURN
      END
