
C  *******************************************************************
C --------------------------------------------------------------------
      SUBROUTINE SUBSTATWEIGHTS(TEMP, RT, UBAR)
C --------------------------------------------------------------------
C
C  SUBROUTINE TO CALCULATE THE STATISTICAL WEIGHTS FOR THE COARSE  
C  GRAINED (9x9) MATRIX, UBAR, IN THE 2NND_RIS MODEL
C 
C  GEOMETRIC AVERAGING IS USED FOR PE. 
C  ARITHMETIC AVERAGING IS USED FOR POE.
C  FOR DETAILS SEE PROGRAM SNND_RIS

      INTEGER   NU
      PARAMETER (NU=9)

      INTEGER   I, J, II
      INTEGER   TEMP

      REAL      ESIGMA, EOMEGA, RT, RTCONFIRM
      REAL      ESIGMAA, ESIGMBB, EOMEGAA, EOMEGAB
      REAL      TEMPERATURE
      REAL      A, B, C, D
      REAL      DUMMY
      REAL      SIGMA, OMEGA
      REAL      SIGMAA, SIGMBB, OMEGAA, OMEGAB
      REAL      SIG1, SIG2, OME1, OME2
      REAL	UBAR(NU,NU,4)


          TEMPERATURE = REAL(TEMP)

C  RT IN  (kJ / mol)
      RTCONFIRM = 8.3145 / 1000.0 * TEMPERATURE
      IF( ( RT - RTCONFIRM ) .GT. 0.00001 ) THEN
        WRITE(6,*) 'PROBLEMS WITH -RT-'
        WRITE(6,*) 'IN SUBROUTINE -SUBSTATWEIGHTS-'
        STOP 65
      ENDIF

       
C  `````````````````` PARAMETERS TO BE EVALUATED ````````````````
C  set the values for the matrices_polyethylene
      ESIGMA = 2.10
      EOMEGA = 8.40
      SIGMA = EXP(-ESIGMA/RT)
      OMEGA = EXP(-EOMEGA/RT)

C  set the values for the matrices_polyoxyethylene
      ESIGMAA =  3.20
      ESIGMBB = -2.20
C     EOMEGAA =  infinite
      EOMEGAB =  2.20
      SIGMAA = EXP(-ESIGMAA/RT)
      SIGMBB = EXP(-ESIGMBB/RT)
C     OMEGAA = EXP(-EOMEGAA/RT)
      OMEGAA = 0.0
      OMEGAB = EXP(-EOMEGAB/RT)
C  `````````````````` PARAMETERS TO BE EVALUATED ````````````````


   11 FORMAT(F16.11,A30)      
      OPEN(UNIT=7,FILE='INPUT_BOND_COND_PROB')
      WRITE(7,11) ESIGMA,         ' ; E_sigma[kJ/mol]    (f16.11)'
      WRITE(7,11) EOMEGA,         ' ; E_omega[kJ/mol]    (f16.11)'
      WRITE(7,11) SIGMA,          ' ; Sigma              (f16.11)'
      WRITE(7,11) OMEGA,          ' ; Omega              (f16.11)'
      WRITE(7,11) ESIGMAA,        ' ; E_sigma_aa[kJ/mol] (f16.11)'
      WRITE(7,11) ESIGMBB,        ' ; E_sigma_bb[kJ/mol] (f16.11)'
      WRITE(7,11) EOMEGAA,        ' ; E_omega_aa[kJ/mol] (f16.11)'
      WRITE(7,11) EOMEGAB,        ' ; E_omega_ab[kJ/mol] (f16.11)'
      WRITE(7,11) SIGMAA,         ' ; Sigma_aa           (f16.11)'
      WRITE(7,11) SIGMBB,         ' ; Sigma_bb           (f16.11)'
      WRITE(7,11) OMEGAA,         ' ; Omega_aa           (f16.11)'
      WRITE(7,11) OMEGAB,         ' ; Omega_ab           (f16.11)'
      WRITE(7,11) TEMPERATURE,    ' ; Temperature        (f16.11)'
      WRITE(7,*) ' '


C                                                   -
C  `````````````````` STATISTICAL WEIGHT MATRICES - U - ```````````

C  **CALCULATIONS FOR POLYETHYLENE**

C  GEOMETRIC AVERAGES:  - A - , - B- , - C - , - D -
        D = (SIGMA*SIGMA) **(1.0/2.0)
        A = D * OMEGA **(1.0/8.0)
        B = D * OMEGA **(1.0/4.0)
        C = D * D * OMEGA **(1.0/2.0)
C  GEOMETRIC AVERAGES:  - A - , - B- , - C - , - D -

        DO I = 1, 9
          UBAR(I,1,1) = 1.0
        ENDDO
        I = 1
        DO J = 2, 5
          UBAR(I,J,1) = D
        ENDDO
        DO J = 6, 7
          UBAR(I,J,1) = D * D
        ENDDO
        DO J = 8, 9
          UBAR(I,J,1) = D * D * OMEGA
        ENDDO
        DO I = 2, 5
          DO J = 2, 5
            UBAR(I,J,1) = A
          ENDDO
          DO J = 6, 7
            UBAR(I,J,1) = D * B
          ENDDO
          DO J = 8, 9
            UBAR(I,J,1) = D * B * OMEGA
          ENDDO
        ENDDO
        DO I = 6, 9
          DO J = 2, 5
            UBAR(I,J,1) = B
          ENDDO
          DO J = 6, 7
            UBAR(I,J,1) = C
          ENDDO
          DO J = 8, 9
            UBAR(I,J,1) = C * OMEGA
          ENDDO
        ENDDO

C  **CALCULATIONS FOR POLYETHYLENE**
	DO I = 1, 9
	WRITE(7,13)(UBAR(I,J,1), J = 1, 9)
	ENDDO

13	FORMAT(9F8.4)

C  **CALCULATIONS FOR POLYOXYETHYLENE**
      
      DO II = 2, 4
	SIG1 = SIGMAA
	SIG2 = SIGMBB
	OME1 = OMEGAB
	OME2 = OMEGAB
        IF (II.EQ.3) THEN
	  SIG1 = SIGMBB
	  SIG2 = SIGMAA
	  OME1 = OMEGAA
        ENDIF
        IF (II.EQ.4) THEN
	  SIG2 = SIGMAA
	  OME2 = OMEGAA
        ENDIF
	
C  ARITHMETIC AVERAGES:  - A - , - B- , - C - , - D -
      DUMMY = SIG1*SIG2
      A = (4.*SIG1 + SIG2*(3.+OME1))/8.
      B = DUMMY * (3.+OME1) / 4.
      C = (2.*SIG1 + SIG2*(1.+OME1))/4.
      D = DUMMY * (1.+OME1) / 2.
C  ARITHMETIC AVERAGES:  - A - , - B- , - C - , - D -

      DO I = 1, 9
        UBAR(I,1,II) = 1.0
      ENDDO
      I = 1
      DO J = 2, 5
        UBAR(I,J,II) = (SIG1+SIG2)/2.
      ENDDO
      DO J = 6, 7
        UBAR(I,J,II) = DUMMY
      ENDDO
      DO J = 8, 9
        UBAR(I,J,II) = DUMMY * OME2
      ENDDO
      DO I = 2, 5
        DO J = 2, 5
          UBAR(I,J,II) = A
        ENDDO
        DO J = 6, 7
          UBAR(I,J,II) = B
        ENDDO
        DO J = 8, 9
          UBAR(I,J,II) = B * OME2
        ENDDO
      ENDDO
      DO I = 6, 9
        DO J = 2, 5
          UBAR(I,J,II) = C
        ENDDO
        DO J = 6, 7
          UBAR(I,J,II) = D
        ENDDO
        DO J = 8, 9
          UBAR(I,J,II) = D * OME2
        ENDDO
      ENDDO

      ENDDO

C  **CALCULATIONS FOR POLYETHYLENE**
C                                                   -
C  `````````````````` STATISTICAL WEIGHT MATRICES - U - ```````````

	

      RETURN
      END
