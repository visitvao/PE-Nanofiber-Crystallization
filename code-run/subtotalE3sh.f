
C  *******************************************************************
      SUBROUTINE SUBTOTALENERGY(BONDVECTOR, BONDVECTORSECONDSHELL,
     >  RT, IFRING, ENERGYINTERFIRSTSHELL, ENERGYINTERSECONDSHELL, 
     >  CONDPROB, TOTALENERGY, SNNDPROB, LC, MCS, ACCEPTRATE,
     >  bondvectorthirdshell, energyinterthirdshell)

      INTEGER  CHAINL, LATTICES, CHAINN

      PARAMETER (CHAINL=100,LATTICES=100,CHAINN=200)

      INTEGER  I, II, J
      INTEGER  A, B, c
      integer  AONEBEAD, APREVIOUS, BONEBEAD, BPREVIOUS
      integer  conebead, cprevious
      INTEGER  IFRING
      INTEGER  IFIRSTVECTOR, ISECONDVECTOR
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  XTEST, YTEST, ZTEST
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      integer  bondvectorthirdshell(92,3)
      INTEGER  X(CHAINN,-2:CHAINL+3), Y(CHAINN,-2:CHAINL+3),
     >         Z(CHAINN,-2:CHAINL+3)
      INTEGER  XP(CHAINN,-2:CHAINL+3), YP(CHAINN,-2:CHAINL+3),
     >         ZP(CHAINN,-2:CHAINL+3)
      INTEGER  BEADOCCUPANCY(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)
      INTEGER  LC(-2:CHAINL+3)
      INTEGER  MCS

      REAL     RT
      REAL     ENERGYINTERFIRSTSHELL, ENERGYINTERSECONDSHELL
      real     energyinterthirdshell
      REAL     TOTALENERGY, TOTALINTRA, TOTALINTER
      REAL     SNNDPROB(4,4)
      REAL     CONDPROB(4,4,4)
      REAL     CHAINENERGY(CHAINN)
      REAL     ACCEPTRATE, DUMMY1, DUMMY2, DUMMY3
      REAL     AOVRALLMEAN, AMEAN, BOVRALLMEAN, BMEAN
      REAL     COVRALLMEAN, CMEAN

      CHARACTER CHTAB

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER, 
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      COMMON  /BEADOCCUPANCY/   BEADOCCUPANCY
      COMMON  /COORDINATES/     X, Y, Z, XP, YP, ZP

   20 FORMAT(3F12.3)

      CHTAB = CHAR(9)          !INSERT TAB BETWEEN DATA POINTS

      A = 0	! COUNT THE OCCUPANCY IN THE FIRST SHELL
      AOVRALLMEAN = 0.
      DO II = 1, CHAINNUMBER
	AONEBEAD = 0
	AMEAN = 0.
        DO I = 1, CHAINLENGTH
          APREVIOUS = A
          DO J = 1, 12
            XTEST = XP(II,I) + BONDVECTOR(J,1)
            YTEST = YP(II,I) + BONDVECTOR(J,2)
            ZTEST = ZP(II,I) + BONDVECTOR(J,3)
            IF( I .GE. 3 ) THEN
              IF( (XTEST .EQ. XP(II,I-2)) .AND.
     >            (YTEST .EQ. YP(II,I-2)) .AND.
     >            (ZTEST .EQ. ZP(II,I-2)) ) THEN
                A = A - 1
              ENDIF
            ENDIF
            IF( I .LE. CHAINLENGTH-2 ) THEN
              IF( (XTEST .EQ. XP(II,I+2)) .AND.
     >            (YTEST .EQ. YP(II,I+2)) .AND.
     >            (ZTEST .EQ. ZP(II,I+2)) ) THEN
                A = A - 1
              ENDIF
            ENDIF
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST) .EQ. 1) A = A + 1
          ENDDO
	  AONEBEAD = A - APREVIOUS
C          WRITE(6,*) 'A (1st neighbor) = ', AONEBEAD
          AMEAN = AMEAN + AONEBEAD  !NOT AN INTEGER
        ENDDO
        AMEAN = AMEAN/CHAINLENGTH
        AOVRALLMEAN = AOVRALLMEAN + AMEAN  !NOT AN INTEGER
C	WRITE(6,*) '<A> = ', AMEAN
      ENDDO
      AOVRALLMEAN = AOVRALLMEAN/CHAINNUMBER
C      WRITE(6,*) 'Overall <A> = ', AOVRALLMEAN

      B = 0	!COUNT THE OCCUPANCY IN THE SECOND SHELL
      BOVRALLMEAN = 0.
      DO II = 1, CHAINNUMBER
        BONEBEAD = 0
	BMEAN = 0.
        DO I = 1, CHAINLENGTH
          BPREVIOUS = B
          DO J = 1, 42
            XTEST = XP(II,I) + BONDVECTORSECONDSHELL(J,1)
            YTEST = YP(II,I) + BONDVECTORSECONDSHELL(J,2)
            ZTEST = ZP(II,I) + BONDVECTORSECONDSHELL(J,3)
            IF( I .GE. 3 ) THEN
              IF( (XTEST .EQ. XP(II,I-2)) .AND.
     >            (YTEST .EQ. YP(II,I-2)) .AND.
     >            (ZTEST .EQ. ZP(II,I-2)) ) THEN
                B = B - 1
              ENDIF
            ENDIF
            IF( I .LE. CHAINLENGTH-2 ) THEN
              IF( (XTEST .EQ. XP(II,I+2)) .AND.
     >            (YTEST .EQ. YP(II,I+2)) .AND.
     >            (ZTEST .EQ. ZP(II,I+2)) ) THEN
                B = B - 1
              ENDIF
            ENDIF
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST) .EQ. 1) B = B + 1
          ENDDO
          BONEBEAD = B - BPREVIOUS
          BMEAN = BMEAN + BONEBEAD  !NOT AN INTEGER
        ENDDO
        BMEAN = BMEAN/CHAINLENGTH
        BOVRALLMEAN = BOVRALLMEAN + BMEAN  !NOT AN INTEGER
      ENDDO
      BOVRALLMEAN = BOVRALLMEAN/CHAINNUMBER


      c = 0	!COUNT THE OCCUPANCY IN THE third SHELL
      covrallmean = 0.
      DO II = 1, CHAINNUMBER
        conebead = 0
	cmean = 0.
        DO I = 1, CHAINLENGTH
          cprevious = c
          DO J = 1, 92
            XTEST = XP(II,I) + bondvectorthirdshell(J,1)
            YTEST = YP(II,I) + bondvectorthirdshell(J,2)
            ZTEST = ZP(II,I) + bondvectorthirdshell(J,3)
C            IF( I .GE. 3 ) THEN
C              IF( (XTEST .EQ. XP(II,I-2)) .AND.
C     >            (YTEST .EQ. YP(II,I-2)) .AND.
C     >            (ZTEST .EQ. ZP(II,I-2)) ) THEN
C                c = c - 1
C              ENDIF
C            ENDIF
C            IF( I .LE. CHAINLENGTH-2 ) THEN
C              IF( (XTEST .EQ. XP(II,I+2)) .AND.
C     >            (YTEST .EQ. YP(II,I+2)) .AND.
C     >            (ZTEST .EQ. ZP(II,I+2)) ) THEN
C                c = c - 1
C              ENDIF
C            ENDIF
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST) .EQ. 1) c = c + 1
          ENDDO
          conebead = c - cprevious
          cmean = cmean + conebead  !NOT AN INTEGER
        ENDDO
        cmean = cmean/CHAINLENGTH
        covrallmean = covrallmean + cmean  !NOT AN INTEGER
      ENDDO
      covrallmean = covrallmean/CHAINNUMBER

      A = A / 2
      B = B / 2
      c = c / 2



      A = A - CHAINNUMBER * ( CHAINLENGTH - 1 )


      IF( A .LT. 0 ) THEN
        WRITE(56,*) 'PROBLEM WITH NEAREST NEIGHBOUR INTERACTIONS'
        WRITE(56,*) 'IN SUBROUTINE SUBTOTALENERGY'
        STOP 56
      ENDIF





      DO II = 1, CHAINNUMBER

        CALL SUBSNNDRISSTATE(XP(II,1),YP(II,1),ZP(II,1),
     >                       XP(II,3),YP(II,3),ZP(II,3),
     >                       IFIRSTVECTOR)

        IF( IFRING .EQ. 1 ) THEN
          CHAINENERGY(II) = 0.0
          J = CHAINLENGTH + 3
        ELSE
          CHAINENERGY(II) = - ( LOG(SNNDPROB(IFIRSTVECTOR,LC(1))) )* RT
          J = CHAINLENGTH
        ENDIF
        DO I = 4, J
          CALL SUBSNNDRISSTATE(XP(II,I-2),YP(II,I-2),ZP(II,I-2),
     >                         XP(II,I),YP(II,I),ZP(II,I),
     >                         ISECONDVECTOR)
          CHAINENERGY(II) = CHAINENERGY(II) - 
     >        ( LOG(CONDPROB(IFIRSTVECTOR,ISECONDVECTOR,LC(I-2))) )* RT
          IFIRSTVECTOR = ISECONDVECTOR
        ENDDO
      ENDDO




      TOTALINTER = FLOAT(A) * ENERGYINTERFIRSTSHELL
     >           + FLOAT(B) * ENERGYINTERSECONDSHELL
     >           + float(c) * energyinterthirdshell

      TOTALINTRA = 0.0
      DO II = 1, CHAINNUMBER
        TOTALINTRA  = TOTALINTRA  + CHAINENERGY(II)
      ENDDO 

      TOTALENERGY = TOTALINTRA + TOTALINTER

      WRITE(14,20) TOTALENERGY, TOTALINTRA, TOTALINTER

C  SUBTRACT THE BONDED MERS
      AOVRALLMEAN = AOVRALLMEAN - 2.*(1.-1./CHAINLENGTH)
C      WRITE(18,25) AOVRALLMEAN,BOVRALLMEAN
C25    FORMAT(F10.5,6X,F10.5)
C      WRITE(18,*) AOVRALLMEAN, CHTAB, BOVRALLMEAN

      DUMMY1 = AOVRALLMEAN/12.
      DUMMY2 = BOVRALLMEAN/42.
c      DUMMY3 = DUMMY1 * DUMMY2 * ACCEPTRATE
      DUMMY3 = covrallmean/92.
      WRITE(18,*) MCS, CHTAB, DUMMY1, CHTAB, DUMMY2, CHTAB, DUMMY3


      RETURN
      END

