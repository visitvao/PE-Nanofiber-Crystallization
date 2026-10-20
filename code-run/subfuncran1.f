C  *************************
      FUNCTION RAN1(IDUM)
C  *************************
C  RANDOM NUMBER GENERATOR FROM
C  'NUMERICAL RECIPES IN FORTRAN' SECOND EDITION, p.271

c  'minimal' random number generator returns a uniform random deviate
c  between 0.0 and 1.0 (exclusive of the endpoint values).
c  Call with IDUM a negative integer to initialize; thereafter, do not
c  alter IDUM between successive deviates in a sequence. RNMX should
c  approximate the largest floating value that is less than 1.

      INTEGER    IDUM, IA, IM, IQ, IR, NTAB, NDIV
      REAL       RAN1, AM, EPS, RNMX

      PARAMETER (IA=16807, IM=2147483647, AM=1./IM, IQ=127773, IR=2836,
     >           NTAB=32, NDIV=1+(IM-1)/NTAB, EPS=1.2E-7, RNMX=1.-EPS)

      INTEGER    J, K, IV(NTAB), IY
      SAVE       IV, IY
      DATA       IV /NTAB*0/,  IY /0/

C  INITIALIZE:
      IF(IDUM.LE.0.OR.IY.EQ.0) THEN
	IDUM = MAX(-IDUM,1)
	DO J = NTAB + 8, 1, -1
	  K = IDUM / IQ
	  IDUM = IA * (IDUM-K*IQ) - IR * K
	  IF(IDUM.LT.0) IDUM = IDUM + IM
	  IF(J.LE.NTAB) IV(J) = IDUM
        ENDDO
	IY = IV(1)
      ENDIF

C  START WHEN NOT INITIALIZING
      K = IDUM / IQ
      IDUM = IA * (IDUM - K * IQ) - IR * K
      IF(IDUM.LT.0) IDUM = IDUM + IM
      J = 1 + IY/NDIV
      IY = IV(J)
      IV(J) = IDUM
      RAN1 = MIN(AM*IY, RNMX)
      RETURN
      END
