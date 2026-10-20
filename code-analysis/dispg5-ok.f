C
C --------------------------------------------------------------------
      PROGRAM SNNDANALYSE_2NND
C --------------------------------------------------------------------
C
C   5 kinds of mean square displacement
C   g1(t)=<[rN/2(t)-rN/2(0)]2>
C   g2(t)=<[rN/2(t)-rCM(t)-(rN/2(0)-rCM(0))]2>
C   g3(t)=<[rCM(t)-rCM(0)]2>
C   g4(t)=<[rend(t)-rend(0)]2>
C   g5(t)=<[rend(t)-rCM(t)-(rend(0)-rCM(0))]2>
C
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
C      UNIT=21 --> FILE='G1MIDDL.OUT'
C      UNIT=22 --> FILE='G2MIDCM.OUT'
C      UNIT=23 --> FILE='G3CM.OUT'
C      UNIT=24 --> FILE='G4END.OUT'
C      UNIT=25 --> FILE='G5ENDCM.OUT'
C --------------------------------------------------------------------

      INTEGER  CHAINL, CHAINN, NSNAPSHOT
      PARAMETER (CHAINL=50,CHAINN=100,NSNAPSHOT=20000)

      INTEGER  REALLENGTH
      INTEGER  NFILE, NDATA
      INTEGER  I, II, J, K, KK, JJ, I1, I2
      INTEGER  IDATA, IDATA2, IMAX, IMAX2
      INTEGER  TEMP
      INTEGER  EQUILIBRATIONTIME
      INTEGER  CURRENTMCM
      INTEGER  STARTPOINT
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  CHAINMOVES, TOTALMOVES
      INTEGER  SIMULATIONLENGTH, STOREINTERVAL, INTERVAL2
      INTEGER  X(CHAINN,CHAINL), Y(CHAINN,CHAINL), Z(CHAINN,CHAINL)
      INTEGER  XP(CHAINN,CHAINL),YP(CHAINN,CHAINL),ZP(CHAINN,CHAINL)
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      INTEGER  bondvectorthirdshell(92,3)
      INTEGER  BONDDIRECTION(CHAINN,CHAINL+1)
      INTEGER  INTERTYPE(12)
      INTEGER  TYPE(CHAINN,CHAINL)
      REAL     TMATRIXKD(3,3)
      REAL     DIAMONDVECTOR(4,3)

      REAL     DUMMY, AA
      REAL     BONDLENGTH, SNNDLENGTH

      REAL     COMX(CHAINN,NSNAPSHOT),COMY(CHAINN,NSNAPSHOT),
     >	       COMZ(CHAINN,NSNAPSHOT)

      REAL     XPD(CHAINN,CHAINL), YPD(CHAINN,CHAINL),
     >         ZPD(CHAINN,CHAINL)

c 1 head, 2 middle bead, 3 end /////////
      REAL     XPS(CHAINN,3,NSNAPSHOT),
     >	       YPS(CHAINN,3,NSNAPSHOT),
     >	       ZPS(CHAINN,3,NSNAPSHOT)

      REAL     g1x(NSNAPSHOT/2),
     >	       g1y(NSNAPSHOT/2),
     >	       g1z(NSNAPSHOT/2)

      REAL     g2x(NSNAPSHOT/2),
     >	       g2y(NSNAPSHOT/2),
     >	       g2z(NSNAPSHOT/2)

      REAL     g3x(NSNAPSHOT/2),
     >	       g3y(NSNAPSHOT/2),
     >	       g3z(NSNAPSHOT/2)

      REAL     g4x(NSNAPSHOT/2),
     >	       g4y(NSNAPSHOT/2),
     >	       g4z(NSNAPSHOT/2)

      REAL     g5x(NSNAPSHOT/2),
     >	       g5y(NSNAPSHOT/2),
     >	       g5z(NSNAPSHOT/2)


      CHARACTER  GARB1*50
      CHARACTER  POLYMER*18, LITERATURE*47 ,FILENAME22*14
      CHARACTER  FILENAME19*12,FILENAME20*12,FILENAME21*13

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD
      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c  *** the confromation from it begin to analysize *****
       nbegin=1


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

    6 FORMAT(A18)
    7 FORMAT(A47)
    8 FORMAT(A50)
    9 FORMAT(I8)
   10 FORMAT(4I7)
   13 FORMAT(10I3)
   14 FORMAT(F8.3)
   15 FORMAT(6F10.4)
   16 FORMAT(I1)
   17 FORMAT(I2.2)
   18 FORMAT(I8,6F10.4)


      BONDLENGTH = 1.53
      SNNDLENGTH = BONDLENGTH * SQRT(8./3.)

      WRITE(6,*) 'BONDLENGTH = ',BONDLENGTH
      WRITE(6,*) 'SNNDLENGTH = ',SNNDLENGTH
      WRITE(6,*) ' '

      WRITE(6,*) 'file begins to be analyze: nbegin'
      READ(5,*)  nbegin
      WRITE(6,*) 'nbegin = ',nbegin
      WRITE(6,*) ' '

      WRITE(6,*) 'NUMBER OF CONFORMATIONS FILES TO BE ANALYZED ?'
      READ(5,*)  NFILE
      WRITE(6,*) 'NFILE = ',NFILE-nbegin+1
      WRITE(6,*) ' '

      WRITE(6,*) 'NUMBER OF CHAINMOVES IN EACH CONFORMATIONS FILE ?'
      READ(5,*)  CHAINMOVES
      WRITE(6,*) 'CHAINMOVES = ',CHAINMOVES
      WRITE(6,*) ' '
      TOTALMOVES = (NFILE-nbegin+1)*CHAINMOVES
      WRITE(6,*) 'TOTAL SIMULATIONLENGTH = ', 
     >		  TOTALMOVES, ' MCS'
      WRITE(6,*) ' '

c      WRITE(6,*) 'EQUILIBRATION TIME (NOT USED FOR STATISTICS) ?'
c      WRITE(6,*) ' IN MONTE-CARLO STEPS'
c      READ(5,9)  EQUILIBRATIONTIME
      EQUILIBRATIONTIME=0
      WRITE(6,*) ' '

c      WRITE(6,*) 'INTERVAL LENGTH FOR STATISTICS'
c      WRITE(6,*) '-> FOR AVERAGE DIMENSION TRAJECTORIES, (1/100)'
c      INTERVAL1 = (TOTALMOVES-EQUILIBRATIONTIME) / 100
c      WRITE(6,*) INTERVAL1
c      WRITE(6,*) ' '

c      WRITE(6,*) 'INTERVAL LENGTH FOR STATISTICS'
c      WRITE(6,*) '-> FOR TOTAL AVERAGES, (1/10)'
c      INTERVAL2 = (TOTALMOVES-EQUILIBRATIONTIME) / 10
c      WRITE(6,*) INTERVAL2
c      WRITE(6,*) ' '

C  ----------------------
      CALL SUBTRANSFORMMATRIX
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C  ----------------------
C  ----------------------

      CALL SUBGETBONDVECTORS(BONDVECTOR)

c      CALL SUBGETBONDVECTORSECONDSHELL(BONDVECTOR,
c     &                           BONDVECTORSECONDSHELL)

c      CALL SUBGETBONDVECTORTHIRDSHELL(BONDVECTOR,
c     >       BONDVECTORSECONDSHELL, bondvectorthirdshell)

c      CALL subgetbondvectorfourthshell(BONDVECTOR,
c     >       BONDVECTORSECONDSHELL, bondvectorthirdshell,
c     >       bondvectorfourthshell)

c      CALL subgetbondvectorfifthshell(BONDVECTOR,
c     >       BONDVECTORSECONDSHELL, bondvectorthirdshell,
c     >       bondvectorfourthshell, bondvectorfifthshell)
C     ...................................
C     NO COMMON BLOCK
C  ----------------------

C  ----------------------
      CALL SUBINTERBEADVECTORS(BONDLENGTH,INTERTYPE,DIAMONDVECTOR) 
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C  ----------------------

      FILENAME21 = 'conpe2000.dat'
      DO II = nbegin, NFILE
	J = 10+II
        WRITE(FILENAME21(9:9),16) II
        OPEN(UNIT=J,FILE=FILENAME21)
      ENDDO
c        OPEN(UNIT=11,FILE=FILENAME21)

C  ========================================
C  ///  BIG LOOP FOR DYNAMIC-MC MOVES  \\\
C  ----------------------------------------

      IDATA = 0
      IDATA2 = 0

      DO KK = nbegin, NFILE

c	print *,KK

      READ(10+KK,6)  POLYMER
      READ(10+KK,7)  LITERATURE
      READ(10+KK,8)  GARB1
c      READ(10+KK,8)  GARB1
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
      READ(10+KK,9)  IFILM
      READ(10+KK,14) ENERGYINTERFIRSTSHELL
      READ(10+KK,14) ENERGYINTERSECONDSHELL
      READ(10+KK,14) ENERGYINTERTHIRDSHELL
c      READ(10+KK,14)  probmeso
      READ(10+KK,9)  idum

      READ(10+KK,*) CURRENTMCM

c ////////
c	  CHAINMOVES=4000000

c	print *,' '
c	print *,'CURRENTMCM =',CURRENTMCM

      DO II = 1, CHAINNUMBER
        READ(10+KK,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
        READ(10+KK,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)
      ENDDO

      SIMULATIONLENGTH = CURRENTMCM

      IF(KK.EQ.nbegin) THEN 

      STARTPOINT = CURRENTMCM
      REALLENGTH = CHAINLENGTH*2-1


c/// INTERVAL2 /////
c      INTERVAL2= STOREINTERVAL*10

      INTERVAL2= STOREINTERVAL

      WRITE(6,*) 'INTERVAL LENGTH FOR STATISTICS'
      WRITE(6,*) '-> FOR TOTAL AVERAGES, (1/10)'
      WRITE(6,*) INTERVAL2
      WRITE(6,*) ' '
c/// INTERVAL2 /////


      WRITE(6,*) 'CHAINLENGTH = ',CHAINLENGTH
      WRITE(6,*) ' '
      WRITE(6,*) 'REALLENGTH = ',REALLENGTH
      WRITE(6,*) ' '
      WRITE(6,*) '# OF CHAINS = ',CHAINNUMBER
      WRITE(6,*) ' '

      WRITE(6,*) 'INTERVALS AT WHICH THE CONFORMATIONS ARE RECORDED'
      WRITE(6,*) ' -> SMALLEST UNIT POSSIBLE', STOREINTERVAL
      WRITE(6,*) ' '

      IF((TOTALMOVES-EQUILIBRATIONTIME)/STOREINTERVAL.GT.NSNAPSHOT)
     >	THEN
	WRITE(6,*) 'PROBLEM WITH PARAMETER NSNAPSHOT_MAKE IT BIGGER'
      	STOP 5
      ENDIF

      NDATA = CHAINMOVES/STOREINTERVAL

	print *,'NDATA=',NDATA

      ENDIF


      DO JJ = 1, NDATA

	READ(10+KK,*) CURRENTMCM
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
	  TYPE(II,I-1) = INTERTYPE(BONDDIRECTION(II,I))
        ENDDO
      ENDDO

      IF(STOREINTERVAL.NE.CURRENTMCM-SIMULATIONLENGTH) THEN
	print *,STOREINTERVAL,CURRENTMCM,SIMULATIONLENGTH
	WRITE(4,*) 'PROBLEM WITH STOREINTERVAL'
      	STOP 3
      ENDIF

      SIMULATIONLENGTH = CURRENTMCM


C  TRANSFORM BEAD COORDINATES TO THE CARTHESIAN FRAME
C  ----------------------
      CALL SUBTRANSFORMCHAIN(XP,YP,ZP,XPD,YPD,ZPD,SNNDLENGTH)
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C     COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
C			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
C  ----------------------


c //// save the all coordinates //////////
c 1 head, 2 middle bead, 3 end /////////


c      IF(MOD((CURRENTMCM-EQUILIBRATIONTIME-STARTPOINT),
c     >	(STOREINTERVAL*10)).EQ.0) THEN

        IDATA2 = IDATA2 + 1

	DO II=1,CHAINNUMBER
	  XPS(II,1,IDATA2) = XPD(II,1)
	  YPS(II,1,IDATA2) = YPD(II,1)
	  ZPS(II,1,IDATA2) = ZPD(II,1)

	  XPS(II,2,IDATA2) = XPD(II,CHAINLENGTH/2)
	  YPS(II,2,IDATA2) = YPD(II,CHAINLENGTH/2)
	  ZPS(II,2,IDATA2) = ZPD(II,CHAINLENGTH/2)

	  XPS(II,3,IDATA2) = XPD(II,CHAINLENGTH)
	  YPS(II,3,IDATA2) = YPD(II,CHAINLENGTH)
	  ZPS(II,3,IDATA2) = ZPD(II,CHAINLENGTH)
	END DO

c	ENDIF

c //// save the all coordinates //////////


c //// save the coordinates of mass center //////////
      IDATA = IDATA + 1

      DO II = 1, CHAINNUMBER

	COMX(II,IDATA) = 0.0
	COMY(II,IDATA) = 0.0
	COMZ(II,IDATA) = 0.0

        DO I = 1, CHAINLENGTH
          COMX(II,IDATA) = COMX(II,IDATA) + XPD(II,I)
          COMY(II,IDATA) = COMY(II,IDATA) + YPD(II,I)
          COMZ(II,IDATA) = COMZ(II,IDATA) + ZPD(II,I)
        ENDDO

        COMX(II,IDATA) = COMX(II,IDATA) / FLOAT(CHAINLENGTH)
        COMY(II,IDATA) = COMY(II,IDATA) / FLOAT(CHAINLENGTH)
        COMZ(II,IDATA) = COMZ(II,IDATA) / FLOAT(CHAINLENGTH)

      ENDDO
c //// save the coordinates of mass center //////////

   60 ENDDO
      ENDDO

C  ----------------------------------------
C  \\\  BIG LOOP FOR DYNAMIC-MC MOVES  ///
C  ========================================

C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      IF ((SIMULATIONLENGTH-STARTPOINT).NE.(TOTALMOVES)) THEN
	WRITE(4,*)'PROBLEM WITH TOTALMOVES'
	WRITE(4,*)'TOTALMOVES = ', TOTALMOVES
	WRITE(4,*)'RECORDED DATA POINTS = ', 
     >		SIMULATIONLENGTH-STARTPOINT
	STOP 6
      ENDIF


C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C  CALCULATE THE CENTER OF MASS DISPLACEMENT OF THE CHAINS 
C  AS A FUNTION OF MCS (TIME)

      IMAX2 = IDATA2 / 2

      MID=CHAINLENGTH/2


c ************************** for g1 *************************************
      OPEN(UNIT=21,FILE='NG1MIDDL.OUT')

  	DO I = 0, IMAX2
	  g1x(I) = 0.0
	  g1y(I) = 0.0
	  g1z(I) = 0.0

 	  DO J = 1, (IDATA2-I)
	     K = J+I
             DO II = 1, CHAINNUMBER
	      g1x(I) = g1x(I) + (XPS(II,2,K)-XPS(II,2,J))**2.0
              g1y(I) = g1y(I) + (YPS(II,2,K)-YPS(II,2,J))**2.0
              g1z(I) = g1z(I) + (ZPS(II,2,K)-ZPS(II,2,J))**2.0
	     ENDDO

	  ENDDO
	    g1x(I) = g1x(I) / (IDATA2-I) / real(CHAINNUMBER)
	    g1y(I) = g1y(I) / (IDATA2-I) / real(CHAINNUMBER)
	    g1z(I) = g1z(I) / (IDATA2-I) / real(CHAINNUMBER)

	  WRITE(21,*)I*INTERVAL2,g1x(i)+g1y(i)+g1z(i)
     &     ,g1x(i),g1y(i),g1z(i)
        ENDDO

      CLOSE(21)
c ************************** for g1 *************************************


c ************************** for g2 *************************************
	OPEN(22,FILE='NG2MIDCM.OUT')

  	DO I = 0, IMAX2
	  g2x(I) = 0.0
	  g2y(I) = 0.0
	  g2z(I) = 0.0

 	  DO J = 1, (IDATA2-I)
	     K = J+I
             DO II = 1, CHAINNUMBER
	      g2x(I) = g2x(I) + (XPS(II,2,K)-XPS(II,2,J)
     &                           -COMX(II,K)+COMX(II,J))**2.0
	      g2y(I) = g2y(I) + (YPS(II,2,K)-YPS(II,2,J)
     &                           -COMY(II,K)+COMY(II,J))**2.0
	      g2z(I) = g2z(I) + (ZPS(II,2,K)-ZPS(II,2,J)
     &                           -COMZ(II,K)+COMZ(II,J))**2.0
	     ENDDO
	  ENDDO

	    g2x(I) = g2x(I) / (IDATA2-I) / real(CHAINNUMBER)
	    g2y(I) = g2y(I) / (IDATA2-I) / real(CHAINNUMBER)
	    g2z(I) = g2z(I) / (IDATA2-I) / real(CHAINNUMBER)

	WRITE(22,*)I*INTERVAL2,g2x(i)+g2y(i)+g2z(i)
     &    ,g2x(i),g2y(i),g2z(i)
        ENDDO

      CLOSE(22)
c ************************** for g2 *************************************




c ************************** for g3 *************************************
      IMAX = IDATA/2

      OPEN(UNIT=23,FILE='NG3CM.OUT')

  	DO I = 0, IMAX

	  g3x(I) = 0.0
	  g3y(I) = 0.0
	  g3z(I) = 0.0

 	  DO J = 1, (IDATA-I)
	    K = J+I

             DO II = 1, CHAINNUMBER
	      g3x(I) = g3x(I) + (COMX(II,K)-COMX(II,J))**2.0
	      g3y(I) = g3y(I) + (COMY(II,K)-COMY(II,J))**2.0
	      g3z(I) = g3z(I) + (COMZ(II,K)-COMZ(II,J))**2.0
	     ENDDO

	  ENDDO

	    g3x(I) = g3x(I) / (IDATA-I) / real(CHAINNUMBER)
	    g3y(I) = g3y(I) / (IDATA-I) / real(CHAINNUMBER)
	    g3z(I) = g3z(I) / (IDATA-I) / real(CHAINNUMBER)

	  WRITE(23,*)I*STOREINTERVAL,g3x(i)+g3y(i)+g3z(i)
     &    ,g3x(i),g3y(i),g3z(i)

	END DO
	CLOSE(23)
c ************************** for g3 *************************************


c ************************** for g4 *************************************
      OPEN(UNIT=24,FILE='NG4END.OUT')

  	DO I = 0, IMAX2
	  g4x(I) = 0.0
	  g4y(I) = 0.0
	  g4z(I) = 0.0

 	  DO J = 1, (IDATA2-I)
	     K = J+I

             DO II = 1, CHAINNUMBER
	      g4x(I) = g4x(I) + (XPS(II,1,K)-XPS(II,1,J))**2.0
     &               + (XPS(II,3,K)-XPS(II,3,J))**2.0

	      g4y(I) = g4y(I) + (YPS(II,1,K)-YPS(II,1,J))**2.0
     &               + (YPS(II,3,K)-YPS(II,3,J))**2.0

	      g4z(I) = g4z(I) + (ZPS(II,1,K)-ZPS(II,1,J))**2.0
     &               + (ZPS(II,3,K)-ZPS(II,3,J))**2.0

	     ENDDO
	  ENDDO

	    g4x(I) = g4x(I) / (IDATA2-I) / real(CHAINNUMBER) / 2.0
	    g4y(I) = g4y(I) / (IDATA2-I) / real(CHAINNUMBER) / 2.0
	    g4z(I) = g4z(I) / (IDATA2-I) / real(CHAINNUMBER) / 2.0

	  WRITE(24,*)I*INTERVAL2,g4x(i)+g4y(i)+g4z(i)
     &     ,g4x(i),g4y(i),g4z(i)
        ENDDO

      CLOSE(24)

c ************************** for g4 *************************************


c ************************** for g5 *************************************
        OPEN(UNIT=25,FILE='NG5ENDCM.OUT')

  	DO I = 0, IMAX2
	  g5x(I) = 0.0
	  g5y(I) = 0.0
	  g5z(I) = 0.0

 	  DO J = 1, (IDATA2-I)
	     K = J+I

             DO II = 1, CHAINNUMBER

	      g5x(I) = g5x(I) 
     &    + (XPS(II,1,K)-XPS(II,1,J)-COMX(II,K)+COMX(II,J))**2.0
     &    + (XPS(II,3,K)-XPS(II,3,J)-COMX(II,K)+COMX(II,J))**2.0

	      g5y(I) = g5y(I)
     &    + (YPS(II,1,K)-YPS(II,1,J)-COMY(II,K)+COMY(II,J))**2.0
     &    + (YPS(II,3,K)-YPS(II,3,J)-COMY(II,K)+COMY(II,J))**2.0

	      g5z(I) = g5z(I) 
     &    + (ZPS(II,1,K)-ZPS(II,1,J)-COMZ(II,K)+COMZ(II,J))**2.0
     &    + (ZPS(II,3,K)-ZPS(II,3,J)-COMZ(II,K)+COMZ(II,J))**2.0

	     ENDDO
	  ENDDO

	    g5x(I) = g5x(I) / (IDATA2-I) / real(CHAINNUMBER) / 2.0
	    g5y(I) = g5y(I) / (IDATA2-I) / real(CHAINNUMBER) / 2.0
	    g5z(I) = g5z(I) / (IDATA2-I) / real(CHAINNUMBER) / 2.0

	  WRITE(25,*)I*INTERVAL2,g5x(i)+g5y(i)+g5z(i)
     &     ,g5x(i),g5y(i),g5z(i)
        ENDDO

      CLOSE(25)
c ************************** for g5 *************************************

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


C23456789012345678901234567890123456789012345678901234567890123456789012
C --------------------------------------------------------------------
      SUBROUTINE SUBINTERBEADVECTORS(BONDLENGTH,INTERTYPE,DIAMONDVECTOR) 
C --------------------------------------------------------------------
C
C  SUBROUTINE TO DETERMINE THE FOUR UNIT VECTORS THAT LOCATE   
C  THE INTERMEDIATE BEADS ON THE DIAMOND LATTICE

      INTEGER I
      INTEGER INTERTYPE(12)
      INTEGER INTERVECTOR(4,3)

      REAL     BONDLENGTH
      REAL     MAGNITUDE
      REAL     DUMMYVECTOR(3), DUMMYVECTOR2(3)
      REAL     DIAMONDVECTOR(4,3)
      REAL     TMATRIXKD(3,3)

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD

      INTERTYPE(1) = 1
      INTERTYPE(2) = 3
      INTERTYPE(3) = 1
      INTERTYPE(4) = 4
      INTERTYPE(5) = 1
      INTERTYPE(6) = 2
      INTERTYPE(7) = 4
      INTERTYPE(8) = 3
      INTERTYPE(9) = 2
      INTERTYPE(10) = 3
      INTERTYPE(11) = 2
      INTERTYPE(12) = 4

      INTERVECTOR(1,1) = 1
      INTERVECTOR(1,2) = 1
      INTERVECTOR(1,3) = 1

      INTERVECTOR(2,1) = 1
      INTERVECTOR(2,2) = 1
      INTERVECTOR(2,3) = -3

      INTERVECTOR(3,1) = -3
      INTERVECTOR(3,2) = 1
      INTERVECTOR(3,3) = 1

      INTERVECTOR(4,1) = 1
      INTERVECTOR(4,2) = -3
      INTERVECTOR(4,3) = 1

      DO I = 1, 4
	DUMMYVECTOR(1) = FLOAT(INTERVECTOR(I,1))
        DUMMYVECTOR(2) = FLOAT(INTERVECTOR(I,2))
        DUMMYVECTOR(3) = FLOAT(INTERVECTOR(I,3))

        CALL SUBMATVEC(TMATRIXKD,DUMMYVECTOR,DUMMYVECTOR2,3,3)
	MAGNITUDE = DUMMYVECTOR2(1)*DUMMYVECTOR2(1) +
     >	     	    DUMMYVECTOR2(2)*DUMMYVECTOR2(2) +
     >	     	    DUMMYVECTOR2(3)*DUMMYVECTOR2(3)

	MAGNITUDE = SQRT(MAGNITUDE) 

	DIAMONDVECTOR(I,1) = DUMMYVECTOR2(1) / MAGNITUDE * BONDLENGTH
	DIAMONDVECTOR(I,2) = DUMMYVECTOR2(2) / MAGNITUDE * BONDLENGTH
	DIAMONDVECTOR(I,3) = DUMMYVECTOR2(3) / MAGNITUDE * BONDLENGTH
      ENDDO

      RETURN
      END


C --------------------------------------------------------------------
      SUBROUTINE SUBMATVEC(A, V, AV, MFORM, NFORM)
C --------------------------------------------------------------------
C
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
C  *   S U B T R A N S F O R M C H A I N                             *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           04/29/94                     *
C  *                    last modified   09/08/94                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
C
C --------------------------------------------------------------------
      SUBROUTINE SUBTRANSFORMCHAIN(X,Y,Z,XD,YD,ZD,SNNDLENGTH)
C --------------------------------------------------------------------

C  - TRANSFORM BEAD COORDINATES TO THE CARTHESIAN FRAME

      INTEGER  CHAINL, CHAINN

      PARAMETER (CHAINL=50,CHAINN=100)

      INTEGER  I, II
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  X(CHAINN,CHAINL), Y(CHAINN,CHAINL), Z(CHAINN,CHAINL)

      REAL     SNNDLENGTH  
      REAL     DUMMYVECTOR(3), DUMMYVECTOR2(3)
      REAL     XD(CHAINN,CHAINL), YD(CHAINN,CHAINL),
     >         ZD(CHAINN,CHAINL)
      REAL     TMATRIXKD(3,3)

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD
      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ


      DO I = 1, CHAINLENGTH
        DO II = 1, CHAINNUMBER	
          DUMMYVECTOR(1) = FLOAT(X(II,I))
          DUMMYVECTOR(2) = FLOAT(Y(II,I))
          DUMMYVECTOR(3) = FLOAT(Z(II,I))
C  ----------------------
          CALL SUBMATVEC(TMATRIXKD,DUMMYVECTOR,DUMMYVECTOR2,3,3)
C         ...................................
C         NO COMMON BLOCKS
C  ----------------------
C  BRING NEW COORDINATES TO SYSTEM OF THE ATOMISTIC CHAIN
C  ( => ONE UNIT ON THE SNND-LATTICE IS (SNNDLENGTH) ANGSTROM )
          XD(II,I) = DUMMYVECTOR2(1) * SNNDLENGTH
          YD(II,I) = DUMMYVECTOR2(2) * SNNDLENGTH
          ZD(II,I) = DUMMYVECTOR2(3) * SNNDLENGTH

        ENDDO
      ENDDO


      RETURN
      END


C  *******************************************************************
C  *                                                                 *
C  *   S U B T R A N S F O R M M A T R I X                           *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           05/02/94                     *
C  *                    last modified   05/02/94                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
C
C --------------------------------------------------------------------
      SUBROUTINE SUBTRANSFORMMATRIX
C --------------------------------------------------------------------

C  - SET UP TRANSFORMATION MATRIX TO CALCULATE THE CARTHESIAN COORDINATES

      INTEGER  I

      REAL     R
      REAL     DEGRAD
      REAL     TMATRIXKD(3,3)

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD

      DEGRAD=ASIN(1.0)/90.0

      TMATRIXKD(1,1) = 1.0
      TMATRIXKD(2,1) = 0.0
      TMATRIXKD(3,1) = 0.0
      TMATRIXKD(1,2) = COS(60.0*DEGRAD)
      TMATRIXKD(2,2) = SIN(60.0*DEGRAD)
      TMATRIXKD(3,2) = 0.0
      TMATRIXKD(1,3) = COS(30.0*DEGRAD) / SQRT(3.0)
      TMATRIXKD(2,3) = SIN(30.0*DEGRAD) / SQRT(3.0)
      TMATRIXKD(3,3) = SQRT(2.0/3.0)
      DO I = 1, 3
        R = TMATRIXKD(1,I)**2.0 + TMATRIXKD(2,I)**2.0 +
     >      TMATRIXKD(3,I)**2.0
        TMATRIXKD(1,I) = TMATRIXKD(1,I) / SQRT(R)
        TMATRIXKD(2,I) = TMATRIXKD(2,I) / SQRT(R)
        TMATRIXKD(3,I) = TMATRIXKD(3,I) / SQRT(R)
      ENDDO

      RETURN
      END


C  *******************************************************************
C  *                                                                 *
C  *   S U B T R A N S R E A L C H A I N                             *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           08/17/95                     *
C  *                    last modified   11/01/96                     *
C  *                                                                 *
C  * (C) copyright  Pemra Doruker                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: pemra@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
C
C --------------------------------------------------------------------
      SUBROUTINE SUBTRANSREALCHAIN(X,Y,Z,XD,YD,ZD,TYPE,
     >		DIAMONDVECTOR,IFRING,REALLENGTH) 
C --------------------------------------------------------------------
C
C  SUBROUTINE TO CONVERT THE 2NND CHAIN TO THE REAL DIAMOND CHAIN   
C  ON CARTESIAN COORDINATES BY LOCATING THE INTERMEDIATE BEADS

      INTEGER  CHAINL, CHAINN
      PARAMETER (CHAINL=50,CHAINN=100)

      INTEGER  REALLENGTH
      INTEGER  IFRING
      INTEGER  I, II, ICOUNT
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  TYPE(CHAINN,CHAINL)

      REAL     DIAMONDVECTOR(4,3)
      REAL     X(CHAINN,CHAINL), Y(CHAINN,CHAINL), Z(CHAINN,CHAINL)
      REAL     XD(CHAINN,CHAINL*2), YD(CHAINN,CHAINL*2), 
     >	       ZD(CHAINN,CHAINL*2)

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ

      DO II = 1, CHAINNUMBER
          XD(II,1) = X(II,1) 
          YD(II,1) = Y(II,1) 
          ZD(II,1) = Z(II,1)
      ENDDO

      DO II = 1, CHAINNUMBER
        ICOUNT = 1
        DO I = 1, CHAINLENGTH-1
	  ICOUNT = ICOUNT + 1
          XD(II,ICOUNT) = XD(II,ICOUNT-1) + DIAMONDVECTOR(TYPE(II,I),1)
          YD(II,ICOUNT) = YD(II,ICOUNT-1) + DIAMONDVECTOR(TYPE(II,I),2)
          ZD(II,ICOUNT) = ZD(II,ICOUNT-1) + DIAMONDVECTOR(TYPE(II,I),3)
	  ICOUNT = ICOUNT + 1
          XD(II,ICOUNT) = X(II,I+1) 
          YD(II,ICOUNT) = Y(II,I+1) 
          ZD(II,ICOUNT) = Z(II,I+1) 
        ENDDO

	IF(IFRING.EQ.1) THEN
	  ICOUNT = ICOUNT + 1
          XD(II,ICOUNT) = XD(II,ICOUNT-1) + 
     >			DIAMONDVECTOR(TYPE(II,CHAINLENGTH),1)
          YD(II,ICOUNT) = YD(II,ICOUNT-1) + 
     >			DIAMONDVECTOR(TYPE(II,CHAINLENGTH),2)
          ZD(II,ICOUNT) = ZD(II,ICOUNT-1) + 
     >			DIAMONDVECTOR(TYPE(II,CHAINLENGTH),3)
	ENDIF
      ENDDO

      IF (ICOUNT.NE.REALLENGTH) THEN
	WRITE(4,*)'PROBLEM WITH REALLENGTH IN SUBTRANSREALCHAIN'
	STOP 8
      ENDIF


      RETURN
      END
