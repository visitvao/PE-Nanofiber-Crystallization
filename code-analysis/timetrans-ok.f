
C --------------------------------------------------------------------
      PROGRAM SNNDANALYSE_2NND
C --------------------------------------------------------------------
C --------------------------------------------------------------------
c   This program is used to calculate trans conformation
C   as a function of MCS time
C --------------------------------------------------------------------
C  VARIABLES:
C  ----------
C  INTERVAL1            =  INTERVAL AT WHICH THE VARIABLES ARE WRITTEN TO
C                          THE FILES 'DIMENSIONS'
C  INTERVAL2            =  BIGGER INTERVAL WRITTEN TO 'AVERAGE_DIMENSIONS'
C			   AND	'TOTAL_AVERAGES'
C  SIMULATIONLENGTH     =  TOTAL NUMBER OF MCM IN SIMULATION
C  STOREINTERVAL        =  INTERVAL AT WHICH THE CONFORMATIONS ARE STORED IN THE
C                          INPUT FILE 'CONFORMATIONS'
C  X(CHAINN,CHAINL)     =  FOLDED COORDINATES; INTEGER; SNND-LATTICE
C  XP(CHAINN,CHAINL)    =  PARENT COORDINATES; INTEGER; SNND-LATTICE
C  XPD(CHAINN,CHAINL)   =  PARENT COORDINATES; REAL; CARTHESIAN SPACE

C  LIST OF FILES:
C  --------------
C      UNIT=11-19  --> FILE=CONFORMATIONS1(-9)
c        OPEN(35,FILE='t_transnum.OUT')
c        OPEN(36,FILE='t_transper.OUT')
c        OPEN(37,FILE='t_transper2.OUT')
C --------------------------------------------------------------------

      INTEGER  CHAINL, CHAINN, NSNAPSHOT,MAXTTTT
      PARAMETER (CHAINL=100,CHAINN=200,NSNAPSHOT=8000,MAXTTTT=30)

      INTEGER  REALLENGTH
      INTEGER  NFILE, NDATA
      INTEGER  I, II, J, K, KK, JJ
      INTEGER  TEMP, IFRING,IFILM, IFCOLLAPSE, IFEXCLUDEDVOLUME
      INTEGER  INTERVAL1, INTERVAL2
      INTEGER  EQUILIBRATIONTIME, CURRENTMCM, STARTPOINT
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  CHAINMOVES, TOTALMOVES
      INTEGER  ICHAINNUMBER, ISEQUENCE, ISEED
      INTEGER  SIMULATIONLENGTH, STOREINTERVAL
      INTEGER  X(CHAINN,CHAINL), Y(CHAINN,CHAINL), Z(CHAINN,CHAINL)
      INTEGER  XP(CHAINN,CHAINL),YP(CHAINN,CHAINL),ZP(CHAINN,CHAINL)
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      INTEGER  bondvectorthirdshell(92,3)
      integer  bondvectorfourthshell(162,3)
      integer  bondvectorfifthshell(252,3)
      INTEGER  BONDDIRECTION(CHAINN,CHAINL)
      INTEGER  NOS

      REAL     BONDLENGTH, SNNDLENGTH
      REAL     ENERGYINTERFIRSTSHELL, ENERGYINTERSECONDSHELL
      REAL     ENERGYINTERTHIRDSHELL

      INTEGER  TRANG(CHAINN,CHAINL*2)
      REAL     TRANSGAUCH(3),TRANS4(MAXTTTT)
      REAL     OVERALLTG(3),OVERALLT4(MAXTTTT)
      INTEGER  TTTT,MCSTRANSGAUCH
      REAL     TOTALTG,TOTALTTT,TOTALT,aaa,bbb

      INTEGER REALBONDSTATEL(12,12)
      INTEGER REALBONDSTATER(12,12)

      CHARACTER  GARB1*50
      CHARACTER  POLYMER*18, LITERATURE*47
      CHARACTER  FILENAME19*12,FILENAME20*12,FILENAME21*19

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      IDATA2=0
      IDATA=0

      DO I=1,3
       OVERALLTG(I)=0
      END DO

      DO I=1,MAXTTTT
	OVERALLT4(I)=0
      END DO

C  format

    6 FORMAT(A18)
    7 FORMAT(A47)
    8 FORMAT(A50)
    9 FORMAT(I8)
   10 FORMAT(4I7)
   13 FORMAT(10I2)
   14 FORMAT(F8.3)
   15 FORMAT(6F10.4)
   16 FORMAT(I1)
   17 FORMAT(I2.2)
   18 FORMAT(I8,6F10.4)

      CALL SUBREALBONDSTATES2(REALBONDSTATEL,REALBONDSTATER)

C  INITIALIZE


      DO I = 1, CHAINL
        DO II = 1, CHAINN
          X(II,I) = 0
          Y(II,I) = 0
          Z(II,I) = 0
          XP(II,I) = 0
          YP(II,I) = 0
          ZP(II,I) = 0
        ENDDO
      ENDDO


      BONDLENGTH = 1.53
      SNNDLENGTH = BONDLENGTH * SQRT(8./3.)

      WRITE(6,*) 'BONDLENGTH = ',BONDLENGTH
      WRITE(6,*) 'SNNDLENGTH = ',SNNDLENGTH
      WRITE(6,*) ' '

      WRITE(6,*) 'NUMBER OF CONFORMATIONS FILES TO BE ANALYZED ?'
      READ(5,*)  NFILE
      WRITE(6,*) 'NFILE = ',NFILE
      WRITE(6,*) ' '

      WRITE(6,*) 'NUMBER OF CHAINMOVES IN EACH CONFORMATIONS FILE ?'
      READ(5,*)  CHAINMOVES
      WRITE(6,*) 'CHAINMOVES = ',CHAINMOVES
      WRITE(6,*) ' '
      TOTALMOVES = NFILE*CHAINMOVES
      WRITE(6,*) 'TOTAL SIMULATIONLENGTH = ',TOTALMOVES,' MCS'
      WRITE(6,*) ' '

c      WRITE(6,*) 'EQUILIBRATION TIME (NOT USED FOR STATISTICS) ?'
c      WRITE(6,*) ' IN MONTE-CARLO STEPS'
c      READ(5,9)  EQUILIBRATIONTIME
      EQUILIBRATIONTIME=0
c      WRITE(6,*) ' '


C  ----------------------
      CALL SUBGETBONDVECTORS(BONDVECTOR)

c      CALL SUBGETBONDVECTORSECONDSHELL(BONDVECTOR,BONDVECTORSECONDSHELL)

c      CALL SUBGETBONDVECTORTHIRDSHELL(BONDVECTOR,
c     >       BONDVECTORSECONDSHELL, bondvectorthirdshell)
C     ...................................
C     NO COMMON BLOCK
C  ----------------------

       FILENAME21 = 'conpe2001.dat'
c      FILENAME21 = 'CONFORMATIONS1'

      DO II = 1, NFILE
	J = 10+II
        WRITE(FILENAME21(8:9),17) II
        OPEN(UNIT=J,FILE=FILENAME21)
      ENDDO
     
        OPEN(35,FILE='t_transnum.OUT')
        OPEN(36,FILE='t_transper.OUT')
        OPEN(37,FILE='t_transper2.OUT')


C  ========================================
C  ///  BIG LOOP FOR DYNAMIC-MC MOVES  \\\
C  ----------------------------------------

      DO KK = 1, NFILE

      READ(10+KK,6)  POLYMER
      READ(10+KK,7)  LITERATURE
      READ(10+KK,8)  GARB1
      READ(10+KK,9)  CHAINLENGTH
      READ(10+KK,9)  LATTICESIZEX
      READ(10+KK,9)  LATTICESIZEY
      READ(10+KK,9)  LATTICESIZEZ
      READ(10+KK,9)  CHAINNUMBER
      READ(10+KK,9)  CHAINMOVES
      READ(10+KK,9)  STOREINTERVAL
      READ(10+KK,9)  ISEQUENCE
      READ(10+KK,9)  ISEED
      READ(10+KK,9)  TEMP
      READ(10+KK,9)  IFEXCLUDEDVOLUME
      READ(10+KK,9)  IFRING
      READ(10+KK,9)  NOS
      READ(10+KK,8)  IFLIM
      READ(10+KK,14) ENERGYINTERFIRSTSHELL
      READ(10+KK,14) ENERGYINTERSECONDSHELL
      READ(10+KK,14) ENERGYINTERTHIRDSHELL
      READ(10+KK,8)  GARB1
c      READ(10+KK,8)  GARB1
	
      READ(10+KK,9) CURRENTMCM


c *******
c	  CHAINMOVES=780000

      DO II = 1, CHAINNUMBER
        READ(10+KK,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
        READ(10+KK,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)
      ENDDO

      SIMULATIONLENGTH = CURRENTMCM

C ///// initialize ////////////////////
       IF (KK.EQ.1) THEN 

        STARTPOINT = CURRENTMCM
	REALLENGTH = CHAINLENGTH*2-1

       WRITE(6,*) 'CHAINLENGTH = ',CHAINLENGTH
       WRITE(6,*) ' '
       WRITE(6,*) 'REALLENGTH = ',REALLENGTH
       WRITE(6,*) ' '
       WRITE(6,*) '# OF CHAINS = ',CHAINNUMBER
       WRITE(6,*) ' '

       WRITE(6,*) 'INTERVALS AT WHICH THE CONFORMATIONS ARE RECORDED'
       WRITE(6,*) ' -> SMALLEST UNIT POSSIBLE', STOREINTERVAL
       WRITE(6,*) ' '

        NDATA = CHAINMOVES/STOREINTERVAL

	print *,NDATA

       ENDIF
C ///// initialize ////////////////////


      DO JJ = 1, NDATA

	READ(10+KK,9) CURRENTMCM
        DO II = 1, CHAINNUMBER
          READ(10+KK,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
          READ(10+KK,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)
        ENDDO
        IF((EQUILIBRATIONTIME+STARTPOINT).GE.CURRENTMCM) THEN
          SIMULATIONLENGTH = CURRENTMCM
 	  GOTO 60
	ENDIF


C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      DO II = 1, CHAINNUMBER
        DO I = 2, CHAINLENGTH
          XP(II,I) = XP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),1)
          YP(II,I) = YP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),2)
          ZP(II,I) = ZP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),3)
        ENDDO
34    ENDDO

      IDATA = IDATA + 1

      IF(STOREINTERVAL .NE. CURRENTMCM-SIMULATIONLENGTH) THEN
	WRITE(4,*) 'PROBLEM WITH STOREINTERVAL'
      	STOP 3
      ENDIF

      SIMULATIONLENGTH = CURRENTMCM


C ******** begin to calculate trans and gauche ********************
C------(1)---T
C------(2)--G+
C------(3)--G-


	IDATA2=IDATA2+1

       DO I = 1, CHAINNUMBER
	 NNNN=3
         DO J = 2, CHAINLENGTH - 1
	   TRANG(I,NNNN)=REALBONDSTATEL(BONDDIRECTION(I,J),
     +       BONDDIRECTION(I,J+1))
	   TRANG(I,NNNN+1)=REALBONDSTATER(BONDDIRECTION(I,J),
     +       BONDDIRECTION(I,J+1))
	  NNNN=NNNN+2
	  END DO
	END DO


	DO I=1,3
	 TRANSGAUCH(I)=0
	END DO


       DO I = 1, CHAINNUMBER
         DO J = 3, REALLENGTH-2
c	 IF (TRANG(I,J).GT.3.OR.TRANG(I,J).LE.0) print *,'error',TRANG(I,J)
	   TRANSGAUCH(TRANG(I,J))=TRANSGAUCH(TRANG(I,J))+1
	 END DO
       END DO


       WRITE(35,*) CURRENTMCM,TRANSGAUCH(1),TRANSGAUCH(2),TRANSGAUCH(3)


	overall=TRANSGAUCH(1)+TRANSGAUCH(2)+TRANSGAUCH(3)

          WRITE(36,*) CURRENTMCM,TRANSGAUCH(1)/overall
     &    ,TRANSGAUCH(2)/overall,TRANSGAUCH(3)/overall

          WRITE(37,*) CURRENTMCM-STARTPOINT,TRANSGAUCH(1)/overall
     &    ,TRANSGAUCH(2)/overall,TRANSGAUCH(3)/overall


C ********************************************************************

C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   60 ENDDO
      ENDDO

	CLOSE(35)
	CLOSE(36)
	CLOSE(37)



       DO I=1,MAXTTTT
	  TRANS4(I)=0
       END DO

       DO I = 1, CHAINNUMBER
	   TTTT=0
         DO J = 3, REALLENGTH-2
	  IF (TRANG(I,J).EQ.1) THEN
	    TTTT=TTTT+1
	  ELSE
	    IF (TTTT.GT.MAXTTTT) THEN
		print *,'TTTT>MAXTTTT'
	        TTTT=0
	    ELSE
	      TRANS4(TTTT)=TRANS4(TTTT)+1
	      TTTT=0
	    ENDIF
	  ENDIF
	END DO
      END DO


	  TOTALTTT=0.0
	  DO I=1,MAXTTTT
	    TOTALTTT=TOTALTTT+I*TRANS4(I)
	  END DO

	
	OPEN(32,FILE='TTTTPERCENT.OUT')
	  TOTALT=0
	  DO I=1,MAXTTTT
	   bbb=I*TRANS4(I)/real(TOTALTTT)
	    WRITE(32,*)	I,bbb
	    TOTALT=TOTALT+bbb
	  END DO
	  WRITE(32,*) TEMP,TOTALT
	CLOSE(32)


C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      END




C  *******************************************************************
C  *                                                                 *
C  *   B O N D V E C T O R F I F T H S H E L L                       *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           07/22/96                     *
C  *                    last modified   07/22/96                     *
C  *                                                                 *
C  * (C) copyright  Junhan Cho                                       *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: jcho@polymer.uakron.edu                  *
C  *                                                                 *
C  *******************************************************************
      subroutine subgetbondvectorfifthshell(BONDVECTOR,
     >       BONDVECTORSECONDSHELL, bondvectorthirdshell,
     >       bondvectorfourthshell, bondvectorfifthshell)


      INTEGER I, II, J, K
      INTEGER BONDVECTOR(12,3)
      INTEGER BONDVECTORSECONDSHELL(42,3)
      INTEGER bondvectorthirdshell(92,3)
      integer bondvectorfourthshell(162,3)
      integer bondvectorfifthshell(252,3)
      integer bondvectorfifthshell2(350,3)
c      INTEGER COUNTSQDISTANCE(18)

C----INITIALIZE THE COUNTSQDISTANCE()
c      DO I = 1,18
c        COUNTSQDISTANCE(I) = 0
c      ENDDO


      II = 0
      do I = 1, 162
c       DO I = 1, 42
        DO J = 1, 12
            BONDVECTORfifthSHELL2(II+1,1) = 
     >               BONDVECTORfourthSHELL(I,1) + BONDVECTOR(J,1)
            BONDVECTORfifthSHELL2(II+1,2) = 
     >               BONDVECTORfourthSHELL(I,2) + BONDVECTOR(J,2)
            BONDVECTORfifthSHELL2(II+1,3) = 
     >               BONDVECTORfourthSHELL(I,3) + BONDVECTOR(J,3)

C----WE START WITH 2ND SHELL VECTORS!
            DO K = 1, 12
              IF((BONDVECTORfifthSHELL2(II+1,1) .EQ. BONDVECTOR(K,1))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,2) .EQ. BONDVECTOR(K,2))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,3) .EQ. BONDVECTOR(K,3)))
     >           THEN
               GO TO 28
              ENDIF
            ENDDO

            DO K = 1, 42
             IF((BONDVECTORfifthSHELL2(II+1,1) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,1))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,2) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,2))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,3) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,3)))
     >           THEN
               GO TO 28
              ENDIF
             ENDDO
             
             DO K = 1, 92
             IF((BONDVECTORfifthSHELL2(II+1,1) .EQ. 
     >                         BONDVECTORthirdSHELL(K,1))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,2) .EQ. 
     >                         BONDVECTORthirdSHELL(K,2))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,3) .EQ. 
     >                         BONDVECTORthirdSHELL(K,3)))
     >           THEN
               GO TO 28
              ENDIF
             ENDDO

             DO K = 1, 162
             IF((BONDVECTORfifthSHELL2(II+1,1) .EQ. 
     >                         BONDVECTORfourthSHELL(K,1))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,2) .EQ. 
     >                         BONDVECTORfourthSHELL(K,2))
     >           .AND.
     >         (BONDVECTORfifthSHELL2(II+1,3) .EQ. 
     >                         BONDVECTORfourthSHELL(K,3)))
     >           THEN
               GO TO 28
              ENDIF
             ENDDO

             DO K = 1, II
              IF(( BONDVECTORfifthSHELL2(II+1,1) .EQ.
     >           BONDVECTORfifthSHELL2(K,1) ) .AND.
     >         ( BONDVECTORfifthSHELL2(II+1,2) .EQ.
     >           BONDVECTORfifthSHELL2(K,2) ) .AND.
     >         ( BONDVECTORfifthSHELL2(II+1,3) .EQ.
     >           BONDVECTORfifthSHELL2(K,3) )) THEN
                GO TO 28
               ENDIF
              ENDDO

             IF( ( BONDVECTORfifthSHELL2(II+1,1) .EQ. 0 ) .AND.
     >         ( BONDVECTORfifthSHELL2(II+1,2) .EQ. 0 ) .AND.
     >         ( BONDVECTORfifthSHELL2(II+1,3) .EQ. 0 ) ) THEN
               GO TO 28
             ENDIF

             II = II + 1
28       ENDDO
c       ENDDO
      enddo
C      IF(II .NE. 42) THEN
C        WRITE(15,*) 'PROBLEM WITH -BONDVECTORTHIRDSHELL2-'
C        WRITE(15,*) ' -- NOT ? VECTORS -- '
C        STOP 111
C     ENDIF

      DO J = 1, II
        DO I = 1, 3
          bondvectorfifthshell(J,I) = BONDVECTORfifthSHELL2(J,I)
        ENDDO
c        write(444,*) J,(bondvectorfifthshell(J,k),k=1,3)
      ENDDO

     
      return

      END
C  *******************************************************************
C  *                                                                 *
C  *   B O N D V E C T O R F O U R T H S H E L L                     *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           07/22/96                     *
C  *                    last modified   07/22/96                     *
C  *                                                                 *
C  * (C) copyright  Junhan Cho                                       *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: jcho@polymer.uakron.edu                  *
C  *                                                                 *
C  *******************************************************************
      subroutine subgetbondvectorfourthshell(BONDVECTOR,
     >       BONDVECTORSECONDSHELL, bondvectorthirdshell,
     >       bondvectorfourthshell)


      INTEGER I, II, J, K
      INTEGER BONDVECTOR(12,3)
      INTEGER BONDVECTORSECONDSHELL(42,3)
      INTEGER bondvectorthirdshell(92,3)
      integer bondvectorfourthshell(162,3)
      integer bondvectorfourthshell2(200,3)
c      INTEGER COUNTSQDISTANCE(18)

C----INITIALIZE THE COUNTSQDISTANCE()
c      DO I = 1,18
c        COUNTSQDISTANCE(I) = 0
c      ENDDO


      II = 0
      do I = 1, 92
c       DO I = 1, 42
        DO J = 1, 12
            BONDVECTORfourthSHELL2(II+1,1) = 
     >               BONDVECTORthirdSHELL(I,1) + BONDVECTOR(J,1)
            BONDVECTORfourthSHELL2(II+1,2) = 
     >               BONDVECTORthirdSHELL(I,2) + BONDVECTOR(J,2)
            BONDVECTORfourthSHELL2(II+1,3) = 
     >               BONDVECTORthirdSHELL(I,3) + BONDVECTOR(J,3)

C----WE START WITH 2ND SHELL VECTORS!
            DO K = 1, 12
              IF((BONDVECTORfourthSHELL2(II+1,1) .EQ. BONDVECTOR(K,1))
     >           .AND.
     >         (BONDVECTORfourthSHELL2(II+1,2) .EQ. BONDVECTOR(K,2))
     >           .AND.
     >         (BONDVECTORfourthSHELL2(II+1,3) .EQ. BONDVECTOR(K,3)))
     >           THEN
               GO TO 28
              ENDIF
            ENDDO

            DO K = 1, 42
             IF((BONDVECTORfourthSHELL2(II+1,1) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,1))
     >           .AND.
     >         (BONDVECTORfourthSHELL2(II+1,2) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,2))
     >           .AND.
     >         (BONDVECTORfourthSHELL2(II+1,3) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,3)))
     >           THEN
               GO TO 28
              ENDIF
             ENDDO
             
             DO K = 1, 92
             IF((BONDVECTORfourthSHELL2(II+1,1) .EQ. 
     >                         BONDVECTORthirdSHELL(K,1))
     >           .AND.
     >         (BONDVECTORfourthSHELL2(II+1,2) .EQ. 
     >                         BONDVECTORthirdSHELL(K,2))
     >           .AND.
     >         (BONDVECTORfourthSHELL2(II+1,3) .EQ. 
     >                         BONDVECTORthirdSHELL(K,3)))
     >           THEN
               GO TO 28
              ENDIF
             ENDDO
             DO K = 1, II
              IF(( BONDVECTORfourthSHELL2(II+1,1) .EQ.
     >           BONDVECTORfourthSHELL2(K,1) ) .AND.
     >         ( BONDVECTORfourthSHELL2(II+1,2) .EQ.
     >           BONDVECTORfourthSHELL2(K,2) ) .AND.
     >         ( BONDVECTORfourthSHELL2(II+1,3) .EQ.
     >           BONDVECTORfourthSHELL2(K,3) )) THEN
                GO TO 28
               ENDIF
              ENDDO

             IF( ( BONDVECTORfourthSHELL2(II+1,1) .EQ. 0 ) .AND.
     >         ( BONDVECTORfourthSHELL2(II+1,2) .EQ. 0 ) .AND.
     >         ( BONDVECTORfourthSHELL2(II+1,3) .EQ. 0 ) ) THEN
               GO TO 28
             ENDIF

             II = II + 1
28       ENDDO
c       ENDDO
      enddo
C      IF(II .NE. 42) THEN
C        WRITE(15,*) 'PROBLEM WITH -BONDVECTORTHIRDSHELL2-'
C        WRITE(15,*) ' -- NOT ? VECTORS -- '
C        STOP 111
C     ENDIF

      DO J = 1, II
        DO I = 1, 3
          bondvectorfourthshell(J,I) = BONDVECTORfourthSHELL2(J,I)
        ENDDO
      ENDDO

      return

      END
C  *******************************************************************
C  *                                                                 *
C  *   S U B G E T B O N D V E C T O R S                             *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           04/27/94                     *
C  *                    last modified   04/27/94                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
      SUBROUTINE SUBGETBONDVECTORS(BONDVECTOR)

      INTEGER BONDVECTOR(12,3)

      BONDVECTOR(1,1) = 1
      BONDVECTOR(1,2) = 0
      BONDVECTOR(1,3) = 0

      BONDVECTOR(2,1) = -1
      BONDVECTOR(2,2) = 0
      BONDVECTOR(2,3) = 0

      BONDVECTOR(3,1) = 0
      BONDVECTOR(3,2) = 1
      BONDVECTOR(3,3) = 0

      BONDVECTOR(4,1) = 0
      BONDVECTOR(4,2) = -1
      BONDVECTOR(4,3) = 0

      BONDVECTOR(5,1) = 0
      BONDVECTOR(5,2) = 0
      BONDVECTOR(5,3) = 1

      BONDVECTOR(6,1) = 0
      BONDVECTOR(6,2) = 0
      BONDVECTOR(6,3) = -1

      BONDVECTOR(7,1) = 1
      BONDVECTOR(7,2) = -1
      BONDVECTOR(7,3) = 0

      BONDVECTOR(8,1) = -1
      BONDVECTOR(8,2) = 1
      BONDVECTOR(8,3) = 0

      BONDVECTOR(9,1) = 1
      BONDVECTOR(9,2) = 0
      BONDVECTOR(9,3) = -1

      BONDVECTOR(10,1) = -1
      BONDVECTOR(10,2) = 0
      BONDVECTOR(10,3) = 1

      BONDVECTOR(11,1) = 0
      BONDVECTOR(11,2) = 1
      BONDVECTOR(11,3) = -1

      BONDVECTOR(12,1) = 0
      BONDVECTOR(12,2) = -1
      BONDVECTOR(12,3) = 1

      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *   S U B G E T B O N D V E C T O R S E C O N D S H E L L         *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           08/18/94                     *
C  *                    last modified   08/18/94                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
      SUBROUTINE SUBGETBONDVECTORSECONDSHELL(BONDVECTOR,
     >                                       BONDVECTORSECONDSHELL)

      INTEGER I, II, J, K
      INTEGER BONDVECTOR(12,3)
      INTEGER BONDVECTORSECONDSHELL(42,3)
      INTEGER BONDVECTORSECONDSHELL2(50,3)

      II = 0
      DO I = 1, 12
        DO J = 1, 12
          BONDVECTORSECONDSHELL2(II+1,1) = BONDVECTOR(I,1) +
     >                                     BONDVECTOR(J,1)
          BONDVECTORSECONDSHELL2(II+1,2) = BONDVECTOR(I,2) +
     >                                     BONDVECTOR(J,2)
          BONDVECTORSECONDSHELL2(II+1,3) = BONDVECTOR(I,3) +
     >                                     BONDVECTOR(J,3)

          DO K = 1, 12
            IF((BONDVECTORSECONDSHELL2(II+1,1) .EQ. BONDVECTOR(K,1))
     >           .AND.
     >         (BONDVECTORSECONDSHELL2(II+1,2) .EQ. BONDVECTOR(K,2))
     >           .AND.
     >         (BONDVECTORSECONDSHELL2(II+1,3) .EQ. BONDVECTOR(K,3)))
     >           THEN
              GO TO 28
            ENDIF
          ENDDO

          DO K = 1, II
            IF(( BONDVECTORSECONDSHELL2(II+1,1) .EQ.
     >           BONDVECTORSECONDSHELL2(K,1) ) .AND.
     >         ( BONDVECTORSECONDSHELL2(II+1,2) .EQ.
     >           BONDVECTORSECONDSHELL2(K,2) ) .AND.
     >         ( BONDVECTORSECONDSHELL2(II+1,3) .EQ.
     >           BONDVECTORSECONDSHELL2(K,3) )) THEN
              GO TO 28
            ENDIF
          ENDDO

          IF( ( BONDVECTORSECONDSHELL2(II+1,1) .EQ. 0 ) .AND.
     >        ( BONDVECTORSECONDSHELL2(II+1,2) .EQ. 0 ) .AND.
     >        ( BONDVECTORSECONDSHELL2(II+1,3) .EQ. 0 ) ) THEN
            GO TO 28
          ENDIF

          II = II + 1
   28   ENDDO
      ENDDO
      IF(II .NE. 42) THEN
        WRITE(15,*) 'PROBLEM WITH -BONDVECTORSECONDSHELL2-'
        WRITE(15,*) ' -- NOT 42 VECTORS -- '
        STOP 111
      ENDIF

      DO I = 1, 3
        DO J = 1, 42
          BONDVECTORSECONDSHELL(J,I) = BONDVECTORSECONDSHELL2(J,I)
        ENDDO
      ENDDO

      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *   B O N D V E C T O R T H I R D S H E L L                       *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           11/22/95                     *
C  *                    last modified   11/22/95                     *
C  *                                                                 *
C  * (C) copyright  Junhan Cho                                       *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: jcho@polymer.uakron.edu                  *
C  *                                                                 *
C  *******************************************************************
      subroutine subgetbondvectorthirdshell(BONDVECTOR,
     >       BONDVECTORSECONDSHELL, bondvectorthirdshell)

      INTEGER I, II, J, K
      INTEGER BONDVECTOR(12,3)
      INTEGER BONDVECTORSECONDSHELL(42,3)
      INTEGER bondvectorthirdshell(92,3)
      INTEGER BONDVECTORTHIRDSHELL2(100,3)
c      INTEGER COUNTSQDISTANCE(18)

C----INITIALIZE THE COUNTSQDISTANCE()
c      DO I = 1,18
c        COUNTSQDISTANCE(I) = 0
c      ENDDO


      II = 0
      DO I = 1, 42
        DO J = 1, 12
            BONDVECTORTHIRDSHELL2(II+1,1) = 
     >               BONDVECTORSECONDSHELL(I,1) + BONDVECTOR(J,1)
            BONDVECTORTHIRDSHELL2(II+1,2) = 
     >               BONDVECTORSECONDSHELL(I,2) + BONDVECTOR(J,2)
            BONDVECTORTHIRDSHELL2(II+1,3) = 
     >               BONDVECTORSECONDSHELL(I,3) + BONDVECTOR(J,3)

C----WE START WITH 2ND SHELL VECTORS!
            DO K = 1, 12
              IF((BONDVECTORTHIRDSHELL2(II+1,1) .EQ. BONDVECTOR(K,1))
     >           .AND.
     >         (BONDVECTORTHIRDSHELL2(II+1,2) .EQ. BONDVECTOR(K,2))
     >           .AND.
     >         (BONDVECTORTHIRDSHELL2(II+1,3) .EQ. BONDVECTOR(K,3)))
     >           THEN
               GO TO 28
              ENDIF
            ENDDO

            DO K = 1, 42
             IF((BONDVECTORTHIRDSHELL2(II+1,1) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,1))
     >           .AND.
     >         (BONDVECTORTHIRDSHELL2(II+1,2) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,2))
     >           .AND.
     >         (BONDVECTORTHIRDSHELL2(II+1,3) .EQ. 
     >                         BONDVECTORSECONDSHELL(K,3)))
     >           THEN
               GO TO 28
              ENDIF
             ENDDO
             DO K = 1, II
              IF(( BONDVECTORTHIRDSHELL2(II+1,1) .EQ.
     >           BONDVECTORTHIRDSHELL2(K,1) ) .AND.
     >         ( BONDVECTORTHIRDSHELL2(II+1,2) .EQ.
     >           BONDVECTORTHIRDSHELL2(K,2) ) .AND.
     >         ( BONDVECTORTHIRDSHELL2(II+1,3) .EQ.
     >           BONDVECTORTHIRDSHELL2(K,3) )) THEN
                GO TO 28
               ENDIF
              ENDDO

             IF( ( BONDVECTORTHIRDSHELL2(II+1,1) .EQ. 0 ) .AND.
     >         ( BONDVECTORTHIRDSHELL2(II+1,2) .EQ. 0 ) .AND.
     >         ( BONDVECTORTHIRDSHELL2(II+1,3) .EQ. 0 ) ) THEN
               GO TO 28
             ENDIF

             II = II + 1
28       ENDDO
       ENDDO
C      IF(II .NE. 42) THEN
C        WRITE(15,*) 'PROBLEM WITH -BONDVECTORTHIRDSHELL2-'
C        WRITE(15,*) ' -- NOT ? VECTORS -- '
C        STOP 111
C     ENDIF

      DO I = 1, 3
        DO J = 1, II
          bondvectorthirdshell(J,I) = BONDVECTORTHIRDSHELL2(J,I)
        ENDDO
      ENDDO


      return
      END


C --------------------------------------------------------------------
      SUBROUTINE SUBMATMAT(A,B,AB,MFORM,NFORM,LFORM)
C --------------------------------------------------------------------
C  MATRIX MULTIPLICATION A.B=AB
C
      INTEGER     MFORM, NFORM, LFORM
      INTEGER     I, J, K
      REAL        A(MFORM,NFORM), B(NFORM,LFORM), AB(MFORM,LFORM)

      DO 5 I = 1, MFORM
        DO 6 J = 1, LFORM
          AB(I,J) = 0.0
    6   CONTINUE
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, LFORM
          DO 30 K = 1, NFORM
            AB(I,J) = AB(I,J) + ( A(I,K) * B(K,J) )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      RETURN
      END

C --------------------------------------------------------------------
      SUBROUTINE SUBMATVEC(A, V, AV, MFORM, NFORM)
C --------------------------------------------------------------------
C  MATRIX- VECTOR MULTIPLICATION, A.V=AV
C

      INTEGER      MFORM, NFORM
      INTEGER      I, J
      REAL         A(MFORM,NFORM), V(NFORM), AV(MFORM)

      DO 5 I = 1, MFORM
	AV(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, NFORM
          AV(I) = AV(I) + ( A(I,J) * V(J) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END



C  *******************************************************************
C  *                                                                 *
C  *   S U B R E A L B O N D S T A T E S                             *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * VERSION 0.1        DATE:           7/27/97                      *
C  *                                       			     *
C  *                                                                 *
C  * (C) COPYRIGHT  TURKAN HALILOGLU                                 *
C  *                INSTITUTE OF POLYMER SCIENCE                     *
C  *                UNIVERSITY OF AKRON, OH 44325-3909               *
C  *                E-MAIL: HALILO@POLYMER.UAKRON.EDU                *
C  *                                                                 *
C  *******************************************************************
      SUBROUTINE SUBREALBONDSTATES2(REALBONDSTATEL,REALBONDSTATER)
C ***** ACCORDING TABLE II OF J PHYS. CHEM. 6989 (1998) *********
C------(1)---T
C------(2)--G+
C------(3)--G-
C------(5)---COL
C------(4)---REV

       INTEGER REALBONDSTATEL(12,12)
       INTEGER REALBONDSTATER(12,12)

	REALBONDSTATEL(5,5)=1
	REALBONDSTATER(5,5)=1
	REALBONDSTATEL(5,1)=2
	REALBONDSTATER(5,1)=1
	REALBONDSTATEL(5,3)=3
	REALBONDSTATER(5,3)=1
	REALBONDSTATEL(5,9)=2
	REALBONDSTATER(5,9)=3
	REALBONDSTATEL(5,12)=1
	REALBONDSTATER(5,12)=3
	REALBONDSTATEL(5,8)=3
	REALBONDSTATER(5,8)=3
	REALBONDSTATEL(5,7)=2
	REALBONDSTATER(5,7)=2
	REALBONDSTATEL(5,11)=3
	REALBONDSTATER(5,11)=2
	REALBONDSTATEL(5,10)=1
	REALBONDSTATER(5,10)=2
	REALBONDSTATEL(5,4)=5
	REALBONDSTATER(5,4)=5
	REALBONDSTATEL(5,2)=5
	REALBONDSTATER(5,2)=5
	REALBONDSTATEL(5,6)=4
	REALBONDSTATER(5,6)=4

	REALBONDSTATEL(1,5)=3
	REALBONDSTATER(1,5)=1
	REALBONDSTATEL(1,1)=1
	REALBONDSTATER(1,1)=1
	REALBONDSTATEL(1,3)=2
	REALBONDSTATER(1,3)=1
	REALBONDSTATEL(1,9)=1
	REALBONDSTATER(1,9)=3
	REALBONDSTATEL(1,12)=3
	REALBONDSTATER(1,12)=3
	REALBONDSTATEL(1,8)=2
	REALBONDSTATER(1,8)=3
	REALBONDSTATEL(1,7)=1
	REALBONDSTATER(1,7)=2
	REALBONDSTATEL(1,11)=2
	REALBONDSTATER(1,11)=2
	REALBONDSTATEL(1,10)=3
	REALBONDSTATER(1,10)=2
	REALBONDSTATEL(1,4)=5
	REALBONDSTATER(1,4)=5
	REALBONDSTATEL(1,2)=4
	REALBONDSTATER(1,2)=4
	REALBONDSTATEL(1,6)=5
	REALBONDSTATER(1,6)=5

	REALBONDSTATEL(3,5)=2
	REALBONDSTATER(3,5)=1
	REALBONDSTATEL(3,1)=3
	REALBONDSTATER(3,1)=1
	REALBONDSTATEL(3,3)=1
	REALBONDSTATER(3,3)=1
	REALBONDSTATEL(3,9)=3
	REALBONDSTATER(3,9)=3
	REALBONDSTATEL(3,12)=2
	REALBONDSTATER(3,12)=3
	REALBONDSTATEL(3,8)=1
	REALBONDSTATER(3,8)=3
	REALBONDSTATEL(3,7)=3
	REALBONDSTATER(3,7)=2
	REALBONDSTATEL(3,11)=1
	REALBONDSTATER(3,11)=2
	REALBONDSTATEL(3,10)=2
	REALBONDSTATER(3,10)=2
	REALBONDSTATEL(3,4)=4
	REALBONDSTATER(3,4)=4
	REALBONDSTATEL(3,2)=5
	REALBONDSTATER(3,2)=5
	REALBONDSTATEL(3,6)=5
	REALBONDSTATER(3,6)=5

	REALBONDSTATEL(9,5)=5
	REALBONDSTATER(9,5)=5
	REALBONDSTATEL(9,1)=1
	REALBONDSTATER(9,1)=2
	REALBONDSTATEL(9,3)=3
	REALBONDSTATER(9,3)=3
	REALBONDSTATEL(9,9)=1
	REALBONDSTATER(9,9)=1
	REALBONDSTATEL(9,12)=5
	REALBONDSTATER(9,12)=5
	REALBONDSTATEL(9,8)=3
	REALBONDSTATER(9,8)=2
	REALBONDSTATEL(9,7)=1
	REALBONDSTATER(9,7)=3
	REALBONDSTATEL(9,11)=3
	REALBONDSTATER(9,11)=1
	REALBONDSTATEL(9,10)=4
	REALBONDSTATER(9,10)=4
	REALBONDSTATEL(9,4)=2
	REALBONDSTATER(9,4)=2
	REALBONDSTATEL(9,2)=2
	REALBONDSTATER(9,2)=3
	REALBONDSTATEL(9,6)=2
	REALBONDSTATER(9,6)=1

	REALBONDSTATEL(12,5)=1
	REALBONDSTATER(12,5)=2
	REALBONDSTATEL(12,1)=3
	REALBONDSTATER(12,1)=3
	REALBONDSTATEL(12,3)=5
	REALBONDSTATER(12,3)=5
	REALBONDSTATEL(12,9)=3
	REALBONDSTATER(12,9)=2
	REALBONDSTATEL(12,12)=1
	REALBONDSTATER(12,12)=1
	REALBONDSTATEL(12,8)=5
	REALBONDSTATER(12,8)=5
	REALBONDSTATEL(12,7)=3
	REALBONDSTATER(12,7)=1
	REALBONDSTATEL(12,11)=4
	REALBONDSTATER(12,11)=4
	REALBONDSTATEL(12,10)=1
	REALBONDSTATER(12,10)=3
	REALBONDSTATEL(12,4)=2
	REALBONDSTATER(12,4)=1
	REALBONDSTATEL(12,2)=2
	REALBONDSTATER(12,2)=2
	REALBONDSTATEL(12,6)=2
	REALBONDSTATER(12,6)=3

	REALBONDSTATEL(8,5)=3
	REALBONDSTATER(8,5)=3
	REALBONDSTATEL(8,1)=5
	REALBONDSTATER(8,1)=5
	REALBONDSTATEL(8,3)=1
	REALBONDSTATER(8,3)=2
	REALBONDSTATEL(8,9)=5
	REALBONDSTATER(8,9)=5
	REALBONDSTATEL(8,12)=3
	REALBONDSTATER(8,12)=2
	REALBONDSTATEL(8,8)=1
	REALBONDSTATER(8,8)=1
	REALBONDSTATEL(8,7)=4
	REALBONDSTATER(8,7)=4
	REALBONDSTATEL(8,11)=1
	REALBONDSTATER(8,11)=3
	REALBONDSTATEL(8,10)=3
	REALBONDSTATER(8,10)=1
	REALBONDSTATEL(8,4)=2
	REALBONDSTATER(8,4)=3
	REALBONDSTATEL(8,2)=2
	REALBONDSTATER(8,2)=1
	REALBONDSTATEL(8,6)=2
	REALBONDSTATER(8,6)=2

	REALBONDSTATEL(7,5)=2
	REALBONDSTATER(7,5)=2
	REALBONDSTATEL(7,1)=1
	REALBONDSTATER(7,1)=3
	REALBONDSTATEL(7,3)=5
	REALBONDSTATER(7,3)=5
	REALBONDSTATEL(7,9)=1
	REALBONDSTATER(7,9)=2
	REALBONDSTATEL(7,12)=2
	REALBONDSTATER(7,12)=1
	REALBONDSTATEL(7,8)=4
	REALBONDSTATER(7,8)=4
	REALBONDSTATEL(7,7)=1
	REALBONDSTATER(7,7)=1
	REALBONDSTATEL(7,11)=5
	REALBONDSTATER(7,11)=5
	REALBONDSTATEL(7,10)=2
	REALBONDSTATER(7,10)=3
	REALBONDSTATEL(7,4)=3
	REALBONDSTATER(7,4)=1
	REALBONDSTATEL(7,2)=3
	REALBONDSTATER(7,2)=2
	REALBONDSTATEL(7,6)=3
	REALBONDSTATER(7,6)=3

	REALBONDSTATEL(11,5)=5
	REALBONDSTATER(11,5)=5
	REALBONDSTATEL(11,1)=2
	REALBONDSTATER(11,1)=2
	REALBONDSTATEL(11,3)=1
	REALBONDSTATER(11,3)=3
	REALBONDSTATEL(11,9)=2
	REALBONDSTATER(11,9)=1
	REALBONDSTATEL(11,12)=4
	REALBONDSTATER(11,12)=4
	REALBONDSTATEL(11,8)=1
	REALBONDSTATER(11,8)=2
	REALBONDSTATEL(11,7)=2
	REALBONDSTATER(11,7)=3
	REALBONDSTATEL(11,11)=1
	REALBONDSTATER(11,11)=1
	REALBONDSTATEL(11,10)=5
	REALBONDSTATER(11,10)=5
	REALBONDSTATEL(11,4)=3
	REALBONDSTATER(11,4)=2
	REALBONDSTATEL(11,2)=3
	REALBONDSTATER(11,2)=3
	REALBONDSTATEL(11,6)=3
	REALBONDSTATER(11,6)=1

	REALBONDSTATEL(10,5)=1
	REALBONDSTATER(10,5)=3
	REALBONDSTATEL(10,1)=5
	REALBONDSTATER(10,1)=5
	REALBONDSTATEL(10,3)=2
	REALBONDSTATER(10,3)=2
	REALBONDSTATEL(10,9)=4
	REALBONDSTATER(10,9)=4
	REALBONDSTATEL(10,12)=1
	REALBONDSTATER(10,12)=2
	REALBONDSTATEL(10,8)=2
	REALBONDSTATER(10,8)=1
	REALBONDSTATEL(10,7)=5
	REALBONDSTATER(10,7)=5
	REALBONDSTATEL(10,11)=2
	REALBONDSTATER(10,11)=3
	REALBONDSTATEL(10,10)=1
	REALBONDSTATER(10,10)=1
	REALBONDSTATEL(10,4)=3
	REALBONDSTATER(10,4)=3
	REALBONDSTATEL(10,2)=3
	REALBONDSTATER(10,2)=1
	REALBONDSTATEL(10,6)=3
	REALBONDSTATER(10,6)=2

	REALBONDSTATEL(4,5)=3
	REALBONDSTATER(4,5)=2
	REALBONDSTATEL(4,1)=2
	REALBONDSTATER(4,1)=3
	REALBONDSTATEL(4,3)=4
	REALBONDSTATER(4,3)=4
	REALBONDSTATEL(4,9)=2
	REALBONDSTATER(4,9)=2
	REALBONDSTATEL(4,12)=3
	REALBONDSTATER(4,12)=1
	REALBONDSTATEL(4,8)=5
	REALBONDSTATER(4,8)=5
	REALBONDSTATEL(4,7)=2
	REALBONDSTATER(4,7)=1
	REALBONDSTATEL(4,11)=5
	REALBONDSTATER(4,11)=5
	REALBONDSTATEL(4,10)=3
	REALBONDSTATER(4,10)=3
	REALBONDSTATEL(4,4)=1
	REALBONDSTATER(4,4)=1
	REALBONDSTATEL(4,2)=1
	REALBONDSTATER(4,2)=2
	REALBONDSTATEL(4,6)=1
	REALBONDSTATER(4,6)=3

	REALBONDSTATEL(2,5)=2
	REALBONDSTATER(2,5)=3
	REALBONDSTATEL(2,1)=4
	REALBONDSTATER(2,1)=4
	REALBONDSTATEL(2,3)=3
	REALBONDSTATER(2,3)=2
	REALBONDSTATEL(2,9)=5
	REALBONDSTATER(2,9)=5
	REALBONDSTATEL(2,12)=2
	REALBONDSTATER(2,12)=2
	REALBONDSTATEL(2,8)=3
	REALBONDSTATER(2,8)=1
	REALBONDSTATEL(2,7)=5
	REALBONDSTATER(2,7)=5
	REALBONDSTATEL(2,11)=3
	REALBONDSTATER(2,11)=3
	REALBONDSTATEL(2,10)=2
	REALBONDSTATER(2,10)=1
	REALBONDSTATEL(2,4)=1
	REALBONDSTATER(2,4)=3
	REALBONDSTATEL(2,2)=1
	REALBONDSTATER(2,2)=1
	REALBONDSTATEL(2,6)=1
	REALBONDSTATER(2,6)=2

	REALBONDSTATEL(6,5)=4
	REALBONDSTATER(6,5)=4
	REALBONDSTATEL(6,1)=3
	REALBONDSTATER(6,1)=2
	REALBONDSTATEL(6,3)=2
	REALBONDSTATER(6,3)=3
	REALBONDSTATEL(6,9)=3
	REALBONDSTATER(6,9)=1
	REALBONDSTATEL(6,12)=5
	REALBONDSTATER(6,12)=5
	REALBONDSTATEL(6,8)=2
	REALBONDSTATER(6,8)=2
	REALBONDSTATEL(6,7)=3
	REALBONDSTATER(6,7)=3
	REALBONDSTATEL(6,11)=2
	REALBONDSTATER(6,11)=1
	REALBONDSTATEL(6,10)=5
	REALBONDSTATER(6,10)=5
	REALBONDSTATEL(6,4)=1
	REALBONDSTATER(6,4)=2
	REALBONDSTATEL(6,2)=1
	REALBONDSTATER(6,2)=3
	REALBONDSTATEL(6,6)=1
	REALBONDSTATER(6,6)=1

     
      RETURN
      END
