
C --------------------------------------------------------------------
      PROGRAM SNNDANALYSE_2NND
C --------------------------------------------------------------------
C
C  PROGRAM TO ANALYSE THE TRAJECTORY OF A SIMULATION ON THE
C  'SECOND NEAREST NEIGHBOR DIAMOND' LATTICE.
C
C  NO EXCLUDED VOLUME OF THE INTERMEDIATE (DIAMOND LATTICE) BEADS
C  --> NO PHYSICAL REPRESENTATION ON THE DIAMOND LATTICE
C  --> TORSION ANGLES ARE NOT VALID
C
C  THE 2NND LATTICE IS A DISTORTED SIMPLE CUBIC LATTICE BUT WITH
C  COORDINATION NUMBER 12 INSTEAD OF 6.
C  TO SEE THIS LATTICE AS A SIMPLE CUBIC LATTICE ALLOWS FOR USING
C  THE FAST INTEGER NUMBERS.
C  A TRANSFORMATION ALWAYS ALLOWS THE TRANSITION BETWEEN THESE TWO
C  REPRESENTATIONS.

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
c      UNIT=28,FILE='TT_COOR_PAR_2NND.DAT')
c      UNIT=29,FILE='TT_COOR_2NND.DAT')
c      UNIT=30,FILE='TT_TRANGAUCH.DAT')
c      UNIT=31,FILE='TTGPERCENT.OUT')
c      UNIT=32,FILE='TTTTPERCENT.OUT')
C --------------------------------------------------------------------

      INTEGER  CHAINL, CHAINN, MAXBIN, NSNAPSHOT,MAXTTTT
      PARAMETER (CHAINL=100,CHAINN=50,NSNAPSHOT=8000,MAXTTTT=30)

      INTEGER  REALLENGTH, MCSCONFIG
      INTEGER  NFILE, NDATA
      INTEGER  INO, ICONF, INOCOUNT, ICOLLAPSE
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
      INTEGER  BONDDIRECTION(CHAINN,CHAINL+1)
      INTEGER  NOS
      INTEGER  LC(CHAINL), LCR(CHAINL)
      INTEGER  INTERTYPE(12),TYPE(CHAINN,CHAINL)

      REAL     DUMMY, AA
      REAL     BONDLENGTH, SNNDLENGTH
      REAL     ENERGYINTERFIRSTSHELL, ENERGYINTERSECONDSHELL
      REAL     ENERGYINTERTHIRDSHELL
      REAL     XPD(CHAINN,CHAINL), YPD(CHAINN,CHAINL),
     >         ZPD(CHAINN,CHAINL)
      REAL     XR(CHAINN,CHAINL*2), YR(CHAINN,CHAINL*2), 
     >	       ZR(CHAINN,CHAINL*2)
      REAL     XPS(CHAINN,CHAINL,NSNAPSHOT/10),
     >	       YPS(CHAINN,CHAINL,NSNAPSHOT/10),
     >	       ZPS(CHAINN,CHAINL,NSNAPSHOT/10)
      REAL     TMATRIXKD(3,3)
      REAL     DIAMONDVECTOR(4,3)

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

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD
      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
c      COMMON  /HISTOGRAM/       HISZ, HIST, HISMAG, HISE 
c      COMMON  /HISTOGRAM2/      HIS1ST, HIS2ND
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
c      EQUILIBRATIONTIME=0
c      WRITE(6,*) ' '


      WRITE(6,*) 'TIME  BEGIN TO CALCULATE TRANS AND GAUCHE?'
      WRITE(6,*) ' IN MONTE-CARLO STEPS'
      READ(5,*)  MCSTRANSGAUCH
      WRITE(6,*) ' '


      WRITE(6,*) 'MONTE-CARLO STEP OF THE SPECIFIC CONFIGURATION ?'
      READ(5,*)  MCSCONFIG
      WRITE(6,*) 'MCSCONFIG = ',MCSCONFIG
      WRITE(6,*) ' '

c      WRITE(6,*) 'DO YOU WANT TO TRACK COLLAPSED CONFIGURATIONS?'
c     >	,'IF YES->TYPE 1, IF NO->TYPE 0'
c      READ(5,*)  IFCOLLAPSE
c       IFCOLLAPSE=1
c       WRITE(6,*) 'IFCOLLAPSE = ',IFCOLLAPSE
c      WRITE(6,*) ' '
c      IF (IFCOLLAPSE.EQ.1) OPEN(UNIT=40,FILE='COLLAPSE')

C  ----------------------
      CALL SUBTRANSFORMMATRIX
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C  ----------------------
C  ----------------------
      CALL SUBGETBONDVECTORS(BONDVECTOR)

c      CALL SUBGETBONDVECTORSECONDSHELL(BONDVECTOR,BONDVECTORSECONDSHELL)

c      CALL SUBGETBONDVECTORTHIRDSHELL(BONDVECTOR,
c     >       BONDVECTORSECONDSHELL, bondvectorthirdshell)
C     ...................................
C     NO COMMON BLOCK
C  ----------------------

C  ----------------------
      CALL SUBINTERBEADVECTORS(BONDLENGTH,INTERTYPE,DIAMONDVECTOR) 
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C  ----------------------

      FILENAME21 = 'conpe2001.dat'
      DO II = 1, NFILE
	J = 10+II
        WRITE(FILENAME21(8:9),17) II
        OPEN(UNIT=J,FILE=FILENAME21)
      ENDDO
      

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
c      READ(10+KK,8)  GARB1
      READ(10+KK,8)  GARB1
	
      READ(10+KK,*) CURRENTMCM

      DO II = 1, CHAINNUMBER
        READ(10+KK,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
        READ(10+KK,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)
      ENDDO

      SIMULATIONLENGTH = CURRENTMCM

C ///// initialize ////////////////////
       IF (KK.EQ.1) THEN 

        STARTPOINT = CURRENTMCM
	REALLENGTH = CHAINLENGTH*2-1

c      IF (IFRING.EQ.0) THEN
c	REALLENGTH = CHAINLENGTH*2-1
c      ELSE
c	REALLENGTH = CHAINLENGTH*2
c      ENDIF 

       WRITE(6,*) 'CHAINLENGTH = ',CHAINLENGTH
       WRITE(6,*) ' '
       WRITE(6,*) 'REALLENGTH = ',REALLENGTH
       WRITE(6,*) ' '
       WRITE(6,*) '# OF CHAINS = ',CHAINNUMBER
       WRITE(6,*) ' '

c      WRITE(6,*) 'NUMBER OF CHAINS FOR WHICH OVERALL DIMENSION'
c      WRITE(6,*) 'TRAJECTORIES WILL BE PRINTED (MAX=30)?'
c      READ(5,*)  MAXCHAIN
c      WRITE(6,*) 'MAXCHAIN = ',MAXCHAIN
c      WRITE(6,*) ' '
c      IF (MAXCHAIN.GT.CHAINNUMBER) MAXCHAIN = CHAINNUMBER

c      WRITE(6,*) 'INTERVALS AT WHICH THE CONFORMATIONS ARE RECORDED'
c      WRITE(6,*) ' -> SMALLEST UNIT POSSIBLE', STOREINTERVAL
c      WRITE(6,*) ' '

c      IF((TOTALMOVES-EQUILIBRATIONTIME)/STOREINTERVAL.GT.NSNAPSHOT)
c     >	THEN
c	WRITE(6,*) 'PROBLEM WITH PARAMETER NSNAPSHOT_MAKE IT BIGGER'
c      	STOP 5
c      ENDIF

        NDATA = CHAINMOVES/STOREINTERVAL

       ENDIF
C ///// initialize ////////////////////


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
34    ENDDO

      IDATA = IDATA + 1

      IF(STOREINTERVAL .NE. CURRENTMCM-SIMULATIONLENGTH) THEN
	WRITE(4,*) 'PROBLEM WITH STOREINTERVAL'
      	STOP 3
      ENDIF

      SIMULATIONLENGTH = CURRENTMCM


C ******** begin to calculate trans and gauche ********************
	IF (IDATA*STOREINTERVAL.GT.MCSTRANSGAUCH) THEN


C  **** save the conformation at MCSCONFIG time ********************

       IF(CURRENTMCM.EQ.MCSCONFIG) THEN

C  TRANSFORM BEAD COORDINATES TO THE CARTHESIAN FRAME
C  ----------------------
       CALL SUBTRANSFORMCHAIN(XP,YP,ZP,XPD,YPD,ZPD,SNNDLENGTH)
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C     COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
C			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
C  ----------------------

C  TRANSFORM 2NND CHAIN TO REAL CHAIN ON THE CARTHESIAN FRAME
C  ----------------------
      CALL SUBTRANSREALCHAIN(XPD,YPD,ZPD,XR,YR,ZR,TYPE,
     >		DIAMONDVECTOR,IFRING,REALLENGTH) 
C     ...................................
C     COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
C			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
C  ----------------------


c        OPEN(28,FILE='TT_COOR_PAR_2NND.DAT')
c        WRITE(28,*) MCSCONFIG
c        DO II = 1, CHAINNUMBER
c         DO I = 1, CHAINLENGTH
c	  WRITE(28,*)i,ii,XPD(II,I),YPD(II,I),ZPD(II,I)
c         ENDDO
c         WRITE(28,*)
c        ENDDO
c        CLOSE(28)

c        OPEN(29,FILE='TT_COOR_2NND.DAT')
c        WRITE(29,*) MCSCONFIG
c        DO II = 1, CHAINNUMBER
c          DO I = 1, REALLENGTH
c	   WRITE(29,*)i,ii,XR(II,I),YR(II,I),ZR(II,I)
c          ENDDO
c        WRITE(29,*)
c        ENDDO
c        CLOSE(29)


      ENDIF
C  **** save the conformation at MCSCINFIG time ********************

C **************************** OLD ****************************************
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

C //// SAVE CERTAIN TRANG ////////////////////////
      IF(CURRENTMCM.EQ.MCSCONFIG) THEN
	OPEN(30,FILE='TT_TRANGAUCH.DAT')
        DO I = 1, CHAINNUMBER
          DO J = 1, REALLENGTH
   	    write(30,*)I,J,TRANG(I,J)
	  END DO
        END DO
	CLOSE(30)
       ENDIF
C //// SAVE CERTAIN TRANG ////////////////////////


	DO I=1,3
	 TRANSGAUCH(I)=0
	END DO


       DO I = 1, CHAINNUMBER
         DO J = 3, REALLENGTH-2
	 IF (TRANG(I,J).GT.3.OR.TRANG(I,J).LE.0) print *,'error',TRANG(I,J)
	   TRANSGAUCH(TRANG(I,J))=TRANSGAUCH(TRANG(I,J))+1
	 END DO
       END DO

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


      DO I=1,3
       OVERALLTG(I)=(OVERALLTG(I)*(IDATA2-1)+TRANSGAUCH(I))/real(IDATA2)
      END DO


      DO I=1,MAXTTTT
	OVERALLT4(I)=(OVERALLT4(I)*(IDATA2-1)+TRANS4(I)
     &     /REAL(REALLENGTH-4)/REAL(CHAINNUMBER))/real(IDATA2)
      END DO
   

      ENDIF		
C ********************************************************************

C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   60 ENDDO
      ENDDO

C  ----------------------------------------
C  \\\  BIG LOOP FOR DYNAMIC-MC MOVES  ///
C  ========================================

	TOTALTG=0.0
        DO I=1,3
         TOTALTG=TOTALTG+OVERALLTG(I)
        END DO

        DO I=1,3
         OVERALLTG(I)=OVERALLTG(I)/TOTALTG
        END DO

	OPEN(31,FILE='TTGPERCENT.OUT')
	  WRITE(31,*) '     t            g+              g-'
	  WRITE(31,*) TEMP,OVERALLTG(1),OVERALLTG(2),OVERALLTG(3)
	CLOSE(31)


	  TOTALTTT=0.0
	  DO I=1,MAXTTTT
	    TOTALTTT=TOTALTTT+I*OVERALLT4(I)
	  END DO

	

	OPEN(32,FILE='TTTTPERCENT.OUT')
	  TOTALT=0
	  DO I=1,MAXTTTT
	   aaa=OVERALLT4(I)/TOTALTTT
	   bbb=I*OVERALLT4(I)/TOTALTTT
	    WRITE(32,*)	I,bbb,aaa
	    TOTALT=TOTALT+bbb
	  END DO
	  WRITE(32,*) TEMP,TOTALT
	CLOSE(32)


C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      END



C  *******************************************************************
C  *                                                                 *
C  *   S U B C A L C G E O M E T R Y                                 *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           05/17/94                     *
C  *                    last modified   01/19/95                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
C
C --------------------------------------------------------------------
      SUBROUTINE SUBCALCGEOMETRY(XPD,YPD,ZPD,BONDLENGTH,RVEC,
     >                           COM,RSQUARE,SSQUARE,CHARR,CHARS)
C --------------------------------------------------------------------
C
C  SUBROUTINE TO TEST THE GEOMETRY OF THE 2NND-LATTICE CHAIN.
C
C  HERE FOR THE ANALYSE OF THE MC RUNS
C
C  THE LATTICE SPACING IS _ONE_ (1.0): 
C  -> THE TETRAHEDRAL ANGLE IS 2.0 * ARCSIN( 0.5 / SQRT(3/8) )
C
C  - TRANSFORM BEAD COORDINATES TO THE CARTHESIAN FRAME
C    -> THE COORDINATEST IN THE CARTHESIAN FRAME ARE CORRESPONDING
C       TO THE BONDLENGTH OF A 'REAL' CHAIN = (BONDLENGTH) ANGSTROM
C       I.E. THE SPACING OF THE LATTICE IS  =****** ANGSTROM
C       THE COORDINATES EXPRESSED IN XPD,YPD,ZPD ARE ALREADY IN ANGSTROM !!!!!
C
C  VARIABLES:
C  ----------
C  XPD(CHAINN,CHAINL)        : PARENT COORDINATES; REAL; CARTHESIAN SPACE

      INTEGER  CHAINL, CHAINN
      PARAMETER (CHAINL=50, CHAINN=50)
      INTEGER  I, II, J
      INTEGER  CHAINLENGTH
      INTEGER  CHAINNUMBER
      INTEGER  LATTICESIZEX,LATTICESIZEY,LATTICESIZEZ

      REAL     BONDLENGTH
      REAL     R
      REAL     RSQUARE(CHAINN), SSQUARE(CHAINN)
      REAL     CHARR(CHAINN), CHARS(CHAINN)
      REAL     COM(CHAINN,3), RVEC(CHAINN,3)
      REAL     XPD(CHAINN,CHAINL),YPD(CHAINN,CHAINL),ZPD(CHAINN,CHAINL)
      REAL     TMATRIXKD(3,3)

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD
      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ

      DO I = 1, 3
        DO II = 1, CHAINN
          COM(II,I) = 0.0
        ENDDO
      ENDDO

C  *****************************************

C  CALCULATE THE CENTER OF MASS  - COM -
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
          COM(II,1) = COM(II,1) + XPD(II,I)
          COM(II,2) = COM(II,2) + YPD(II,I)
          COM(II,3) = COM(II,3) + ZPD(II,I)
        ENDDO
      COM(II,1) = COM(II,1) / FLOAT(CHAINLENGTH)
      COM(II,2) = COM(II,2) / FLOAT(CHAINLENGTH)
      COM(II,3) = COM(II,3) / FLOAT(CHAINLENGTH)
      ENDDO
C  CALCULATE THE CENTER OF MASS  - COM -

C  *****************************************

C  CALCULATE THE END-TO-END DISTANCE  - RSQUARE -
      DO II = 1, CHAINNUMBER

        RVEC(II,1) = XPD(II,CHAINLENGTH) - XPD(II,1)
        RVEC(II,2) = YPD(II,CHAINLENGTH) - YPD(II,1)
        RVEC(II,3) = ZPD(II,CHAINLENGTH) - ZPD(II,1)

        RSQUARE(II) = RVEC(II,1)**2.0 + RVEC(II,2)**2.0 +
     >                RVEC(II,3)**2.0
     >               
      ENDDO
C  CALCULATE THE END-TO-END DISTANCE  - RSQUARE -

C  *****************************************

C  CALCULATE THE RADIUS OF GYRATION
      DO II = 1, CHAINNUMBER
        SSQUARE(II) = 0.0
        DO I = 1, CHAINLENGTH-1
          DO J = I+1, CHAINLENGTH
            R = ( ( XPD(II,J) - XPD(II,I) ) ) ** 2.0 +
     >          ( ( YPD(II,J) - YPD(II,I) ) ) ** 2.0 +
     >          ( ( ZPD(II,J) - ZPD(II,I) ) ) ** 2.0
            SSQUARE(II) = SSQUARE(II) + R
          ENDDO
        ENDDO
        SSQUARE(II) = SSQUARE(II) / FLOAT(CHAINLENGTH*CHAINLENGTH)
      ENDDO
C  CALCULATE THE RADIUS OF GYRATION

C  *****************************************

C  CALCULATE THE CHARACTERISTIC RATIO
C  NUMBER OF 'REAL' BONDS = 2.00 * (CHAINLENGTH - 1)
C  LENGTH OF 'REAL' BONDS = (BONDLENGTH) ANGSTROM
      DO II = 1, CHAINNUMBER
        CHARR(II) = RSQUARE(II) /
     >           ( 2.0 * FLOAT(CHAINLENGTH-1) * BONDLENGTH**2.0 )
        CHARS(II) = SSQUARE(II) /
     >           ( 2.0 * FLOAT(CHAINLENGTH-1) * BONDLENGTH**2.0 )
      ENDDO
C  CALCULATE THE CHARACTERISTIC RATIO

C  *****************************************

      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *   S U B F O L D E D C H A I N                                   *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           05/09/94                     *
C  *                    last modified   05/26/95                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
      SUBROUTINE SUBFOLDEDCHAIN(XPP,YPP,ZPP,SNNDLENGTH,CURRENTMCM,
     >	REALLENGTH,DELR,MCSCONFIG,INDEX,IFILM,SLABVOL,IFCOLLAPSE)

      INTEGER  CHAINL, CHAINN, MAXBIN

      PARAMETER (CHAINL=50,CHAINN=50,MAXBIN=50)

      INTEGER  I, II, IBIN
      INTEGER  CHAINLENGTH
      INTEGER  LATTICESIZEX,LATTICESIZEY,LATTICESIZEZ
      INTEGER  CHAINNUMBER
      INTEGER  NCOLLAPSE,IFCOLLAPSE
      INTEGER  CURRENTMCM
      INTEGER  REALLENGTH
      INTEGER  MCSCONFIG
      INTEGER  INDEX, IFILM

      REAL     DELR,SLABVOL
      REAL     RCMZ,RZNEW
      REAL     SNNDLENGTH,DISTANCE  
      REAL     LATTICELENGTHX,LATTICELENGTHY,LATTICELENGTHZ
      REAL     AX(CHAINN,CHAINL*2)
      REAL     AY(CHAINN,CHAINL*2)
      REAL     AZ(CHAINN,CHAINL*2)
      REAL     XPP(CHAINN,CHAINL*2)
      REAL     YPP(CHAINN,CHAINL*2)
      REAL     ZPP(CHAINN,CHAINL*2)
      REAL     HISZ(10,MAXBIN),HISE(10,MAXBIN)
      REAL     XX,YY,ZZ,MAG,COSTH
      REAL     HIST(10,MAXBIN),HISMAG(10,MAXBIN)

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      COMMON  /HISTOGRAM/       HISZ, HIST, HISMAG, HISE 


	LATTICELENGTHX = SNNDLENGTH * LATTICESIZEX
	LATTICELENGTHY = SNNDLENGTH * LATTICESIZEY * 0.86603
	LATTICELENGTHZ = SNNDLENGTH * LATTICESIZEZ * 0.81594
      SLABVOL = 2 * LATTICELENGTHX * LATTICELENGTHY * DELR
      SLABVOL = SLABVOL / 140. * 6.022



      DO II = 1, CHAINNUMBER
      DO I = 1, REALLENGTH
        AX(II,I) = XPP(II,I)
        AY(II,I) = YPP(II,I)
        AZ(II,I) = ZPP(II,I)

  103   IF( AX(II,I) .GE. LATTICELENGTHX ) THEN
          AX(II,I) = AX(II,I) - LATTICELENGTHX
          GO TO 103
        END IF
  104   IF( AX(II,I) .LT. 0. ) THEN
          AX(II,I) = AX(II,I) + LATTICELENGTHX
          GO TO 104
        END IF

  105   IF( AY(II,I) .GE. LATTICELENGTHY ) THEN
          AY(II,I) = AY(II,I) - LATTICELENGTHY
          GO TO 105
        END IF
  106   IF( AY(II,I) .LT. 0. ) THEN
          AY(II,I) = AY(II,I) + LATTICELENGTHY
          GO TO 106
        END IF

  107   IF( AZ(II,I) .GE. LATTICELENGTHZ ) THEN
          AZ(II,I) = AZ(II,I) - LATTICELENGTHZ
          GO TO 107
        END IF
  108   IF( AZ(II,I) .LT. 0. ) THEN
          AZ(II,I) = AZ(II,I) + LATTICELENGTHZ
          GO TO 108
        END IF
      ENDDO
      ENDDO

      IF (IFCOLLAPSE.EQ.1) THEN

      NCOLLAPSE = 0	
      DO II = 1, CHAINNUMBER
      DO I = 1, REALLENGTH-2

	DISTANCE = (AX(II,I)-AX(II,I+2))**2. + (AY(II,I)-AY(II,I+2))**2.
     >			+ (AZ(II,I)-AZ(II,I+2))**2.
	
	IF(DISTANCE.LT.1E-2) THEN
c	write(40,*)ii,i,I+2,sqrt(distance)
	NCOLLAPSE = NCOLLAPSE + 1
	ENDIF

      ENDDO
      ENDDO
	write(40,*)CURRENTMCM,NCOLLAPSE

      ENDIF


c      IF(IFILM.EQ.1) THEN

c	DO II = 1, CHAINNUMBER
c        DO I = 1, REALLENGTH
c	  RCMZ = RCMZ + ZPP(II,I)
c	ENDDO
c	ENDDO
c 	RCMZ = RCMZ / CHAINNUMBER / REALLENGTH


c        DO II = 1, CHAINNUMBER
c        DO I = 1, REALLENGTH
c	  RZNEW = ZPP(II,I) - RCMZ
c	  IBIN = INT(ABS(RZNEW)/DELR) + 1
c	  HISZ(INDEX,IBIN) = HISZ(INDEX,IBIN) + 1.
c	  IF(I.EQ.1.OR.I.EQ.REALLENGTH) THEN
c	    HISE(INDEX,IBIN) = HISE(INDEX,IBIN) + 1.
c	  ENDIF
c        ENDDO
c        ENDDO

c	DO II = 1, CHAINNUMBER
c        DO I = 2, REALLENGTH
c	  XX = XPP(II,I) - XPP(II,I-1)  
c	  YY = YPP(II,I) - YPP(II,I-1)  
c	  ZZ = ZPP(II,I) - ZPP(II,I-1)  
c	  MAG = XX*XX + YY*YY + ZZ*ZZ
c	  COSTH = ZZ*ZZ / MAG
c
c	  ZZ = ZPP(II,I-1) - RCMZ + ZZ/2.
c	  IBIN = INT(ABS(ZZ)/DELR) + 1
c	  HIST(INDEX,IBIN) = HIST(INDEX,IBIN) + 1.
c	  HISMAG(INDEX,IBIN) = HISMAG(INDEX,IBIN) + COSTH
c        ENDDO
c        ENDDO
c
c      ENDIF


C      IF(CURRENTMCM.EQ.MCSCONFIG) THEN
C      WRITE(30,*) MCSCONFIG
C      DO II = 1, CHAINNUMBER
C        DO I = 1, REALLENGTH
C	  WRITE(30,*)i,ii,AX(II,I),AY(II,I),AZ(II,I)
C        ENDDO
C      WRITE(30,*)
C      ENDDO
C      ENDIF


      RETURN
      END


C  *******************************************************************
C  *                                                                 *
C  *   S U B F O L D E D C H A I N 2                                 *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           05/09/94                     *
C  *                    last modified   05/26/95                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
      SUBROUTINE SUBFOLDEDCHAIN2(XP,YP,ZP,X,Y,Z,BONDVECTOR,
     >	BONDVECTORSECONDSHELL, bondvectorthirdshell,
     >	bondvectorfourthshell,bondvectorfifthshell,
     >	atotal,btotal,ctotal,dtotal,etotal,
     >	atotal2,btotal2,ctotal2,dtotal2,etotal2,
     >	atotal3,btotal3,ctotal3,dtotal3,etotal3)

      INTEGER  CHAINL, LATTICES, CHAINN

      PARAMETER(CHAINL=50, LATTICES=25, CHAINN=50)

      INTEGER  I, II, J, K, A, B, C, d, e
      INTEGER  CHAINLENGTH
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  CHAINNUMBER
      INTEGER  XTEST,YTEST,ZTEST
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      INTEGER  bondvectorthirdshell(92,3)
      integer bondvectorfourthshell(162,3)
      integer bondvectorfifthshell(252,3)
      INTEGER  X(CHAINN,CHAINL)
      INTEGER  Y(CHAINN,CHAINL)
      INTEGER  Z(CHAINN,CHAINL)
      INTEGER  XP(CHAINN,CHAINL)
      INTEGER  YP(CHAINN,CHAINL)
      INTEGER  ZP(CHAINN,CHAINL)
      INTEGER  HIS1ST(0:12)
      INTEGER  HIS2ND(0:42)
      INTEGER  BEADOCCUPANCY(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)
      REAL     ATOTAL, BTOTAL, CTOTAL, dtotal, etotal 
      REAL     atotal2, btotal2, ctotal2, dtotal2, etotal2
      REAL     atotal3, btotal3, ctotal3, dtotal3, etotal3

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      COMMON  /HISTOGRAM2/      HIS1ST, HIS2ND


      DO I = 0, LATTICESIZEZ-1
        DO J = 0, LATTICESIZEY-1
          DO K = 0, LATTICESIZEX-1
            BEADOCCUPANCY(K,J,I) = 0
          ENDDO
        ENDDO
      ENDDO

      DO II = 1, CHAINNUMBER
      DO I = 1, CHAINLENGTH

        X(II,I) = XP(II,I)
        Y(II,I) = YP(II,I)
        Z(II,I) = ZP(II,I)

  103   IF( X(II,I) .GE. LATTICESIZEX ) THEN
          X(II,I) = X(II,I) - LATTICESIZEX
          GO TO 103
        END IF
  104   IF( X(II,I) .LT. 0 ) THEN
          X(II,I) = X(II,I) + LATTICESIZEX
          GO TO 104
        END IF

  105   IF( Y(II,I) .GE. LATTICESIZEY ) THEN
          Y(II,I) = Y(II,I) - LATTICESIZEY
          GO TO 105
        END IF
  106   IF( Y(II,I) .LT. 0 ) THEN
          Y(II,I) = Y(II,I) + LATTICESIZEY
          GO TO 106
        END IF

  107   IF( Z(II,I) .GE. LATTICESIZEZ ) THEN
          Z(II,I) = Z(II,I) - LATTICESIZEZ
          GO TO 107
        END IF
  108   IF( Z(II,I) .LT. 0 ) THEN
          Z(II,I) = Z(II,I) + LATTICESIZEZ
          GO TO 108
        END IF

      ENDDO
      ENDDO

       DO II = 1, CHAINNUMBER
         DO I = 1, CHAINLENGTH
           BEADOCCUPANCY(X(II,I),Y(II,I),Z(II,I)) = 1
         ENDDO
       ENDDO

	atotal=0. 
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
	  A = 0
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
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1) A = A + 1
	  ENDDO
	  IF(I.EQ.1.OR.I.EQ.CHAINLENGTH) THEN
		A = A - 1
	  ELSE
		A = A - 2
	  ENDIF

	  HIS1ST(A) = HIS1ST(A) + 1
	atotal=atotal+real(A)
	ENDDO
      ENDDO

      atotal=atotal/real(CHAINNUMBER)/real(CHAINLENGTH)

	btotal=0. 
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
	  B = 0
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
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)B = B + 1
          ENDDO
	  HIS2ND(B) = HIS2ND(B) + 1
	  btotal=btotal+real(B)
        ENDDO
      ENDDO
      btotal=btotal/real(CHAINNUMBER)/real(CHAINLENGTH)

      c = 0	!COUNT THE OCCUPANCY IN THE third SHELL
	ctotal=0. 
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
	  C = 0
          DO J = 1, 92
            XTEST = XP(II,I) + bondvectorthirdshell(J,1)
            YTEST = YP(II,I) + bondvectorthirdshell(J,2)
            ZTEST = ZP(II,I) + bondvectorthirdshell(J,3)
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)c = c+1
          ENDDO
	  ctotal=ctotal+real(C)
        ENDDO
      ENDDO
      ctotal=ctotal/real(CHAINNUMBER)/real(CHAINLENGTH)

	dtotal=0. 
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
	  d = 0
          DO J = 1, 162
            XTEST = XP(II,I) + bondvectorfourthshell(J,1)
            YTEST = YP(II,I) + bondvectorfourthshell(J,2)
            ZTEST = ZP(II,I) + bondvectorfourthshell(J,3)
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)d = d+1
          ENDDO
	  dtotal=dtotal+real(d)
        ENDDO
      ENDDO
      dtotal=dtotal/real(CHAINNUMBER)/real(CHAINLENGTH)

	etotal=0. 
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
	  e = 0
          DO J = 1, 252
            XTEST = XP(II,I) + bondvectorfifthshell(J,1)
            YTEST = YP(II,I) + bondvectorfifthshell(J,2)
            ZTEST = ZP(II,I) + bondvectorfifthshell(J,3)
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)e = e+1
          ENDDO
	  etotal=etotal+real(e)
        ENDDO
      ENDDO
      etotal=etotal/real(CHAINNUMBER)/real(CHAINLENGTH)

	atotal2=0. 
	btotal2=0. 
	ctotal2=0.
	dtotal2=0.
	etotal2=0.
 
      DO II = 1, CHAINNUMBER

      DO I = 0, LATTICESIZEZ-1
        DO J = 0, LATTICESIZEY-1
          DO K = 0, LATTICESIZEX-1
            BEADOCCUPANCY(K,J,I) = 0
          ENDDO
        ENDDO
      ENDDO

      DO I = 1, CHAINLENGTH

        X(II,I) = XP(II,I)
        Y(II,I) = YP(II,I)
        Z(II,I) = ZP(II,I)

  203   IF( X(II,I) .GE. LATTICESIZEX ) THEN
          X(II,I) = X(II,I) - LATTICESIZEX
          GO TO 203
        END IF
  204   IF( X(II,I) .LT. 0 ) THEN
          X(II,I) = X(II,I) + LATTICESIZEX
          GO TO 204
        END IF

  205   IF( Y(II,I) .GE. LATTICESIZEY ) THEN
          Y(II,I) = Y(II,I) - LATTICESIZEY
          GO TO 205
        END IF
  206   IF( Y(II,I) .LT. 0 ) THEN
          Y(II,I) = Y(II,I) + LATTICESIZEY
          GO TO 206
        END IF

  207   IF( Z(II,I) .GE. LATTICESIZEZ ) THEN
          Z(II,I) = Z(II,I) - LATTICESIZEZ
          GO TO 207
        END IF
  208   IF( Z(II,I) .LT. 0 ) THEN
          Z(II,I) = Z(II,I) + LATTICESIZEZ
          GO TO 208
        END IF

      ENDDO

         DO I = 1, CHAINLENGTH
           BEADOCCUPANCY(X(II,I),Y(II,I),Z(II,I)) = 1
         ENDDO

        DO I = 1, CHAINLENGTH
	  A = 0
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
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1) A = A + 1
	  ENDDO
	  IF(I.EQ.1.OR.I.EQ.CHAINLENGTH) THEN
		A = A - 1
	  ELSE
		A = A - 2
	  ENDIF

	atotal2=atotal2+real(A)
	ENDDO


        DO I = 1, CHAINLENGTH
	  B = 0
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
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)B = B + 1
          ENDDO
	  btotal2=btotal2+real(B)
        ENDDO

      c = 0	!COUNT THE OCCUPANCY IN THE third SHELL
        DO I = 1, CHAINLENGTH
	  C = 0
          DO J = 1, 92
            XTEST = XP(II,I) + bondvectorthirdshell(J,1)
            YTEST = YP(II,I) + bondvectorthirdshell(J,2)
            ZTEST = ZP(II,I) + bondvectorthirdshell(J,3)
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)c = c+1
          ENDDO
	  ctotal2=ctotal2+real(C)
        ENDDO

        DO I = 1, CHAINLENGTH
	  d = 0
          DO J = 1, 162
            XTEST = XP(II,I) + bondvectorfourthshell(J,1)
            YTEST = YP(II,I) + bondvectorfourthshell(J,2)
            ZTEST = ZP(II,I) + bondvectorfourthshell(J,3)
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)d = d+1
          ENDDO
	  dtotal2=dtotal2+real(d)
        ENDDO

        DO I = 1, CHAINLENGTH
	  e = 0
          DO J = 1, 252
            XTEST = XP(II,I) + bondvectorfifthshell(J,1)
            YTEST = YP(II,I) + bondvectorfifthshell(J,2)
            ZTEST = ZP(II,I) + bondvectorfifthshell(J,3)
            CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
            IF(BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)e = e+1
          ENDDO
	  etotal2=etotal2+real(e)
        ENDDO


      ENDDO

      atotal2=atotal2/real(CHAINNUMBER)/real(CHAINLENGTH)
      btotal2=btotal2/real(CHAINNUMBER)/real(CHAINLENGTH)
      ctotal2=ctotal2/real(CHAINNUMBER)/real(CHAINLENGTH)
      dtotal2=dtotal2/real(CHAINNUMBER)/real(CHAINLENGTH)
      etotal2=etotal2/real(CHAINNUMBER)/real(CHAINLENGTH)

	atotal3=0. 
	btotal3=0. 
	ctotal3=0.
	dtotal3=0.
	etotal3=0.
 
      DO II = 1, CHAINNUMBER

        DO I = 1, CHAINLENGTH
	  A = 0
          DO J = 1, 12
            XTEST = XP(II,I) + BONDVECTOR(J,1)
            YTEST = YP(II,I) + BONDVECTOR(J,2)
            ZTEST = ZP(II,I) + BONDVECTOR(J,3)

	  DO K = 1, CHAINLENGTH
            IF(K.LE.(I+2).AND.K.GE.(I-2)) GOTO 88
              IF( (XTEST .EQ. XP(II,K)) .AND.
     >            (YTEST .EQ. YP(II,K)) .AND.
     >            (ZTEST .EQ. ZP(II,K)) ) THEN
                A = A + 1
              ENDIF
88	  ENDDO
	  ENDDO
	atotal3=atotal3+real(A)
	ENDDO


        DO I = 1, CHAINLENGTH
	  B = 0
          DO J = 1, 42
            XTEST = XP(II,I) + BONDVECTORSECONDSHELL(J,1)
            YTEST = YP(II,I) + BONDVECTORSECONDSHELL(J,2)
            ZTEST = ZP(II,I) + BONDVECTORSECONDSHELL(J,3)

	  DO K = 1, CHAINLENGTH
            IF(K.LE.(I+2).AND.K.GE.(I-2)) GOTO 89
              IF( (XTEST .EQ. XP(II,K)) .AND.
     >            (YTEST .EQ. YP(II,K)) .AND.
     >            (ZTEST .EQ. ZP(II,K)) ) THEN
                B = B + 1
              ENDIF
89	  ENDDO
	  ENDDO
	  btotal3=btotal3+real(B)
        ENDDO


        DO I = 1, CHAINLENGTH
	  C = 0
          DO J = 1, 92
            XTEST = XP(II,I) + bondvectorthirdshell(J,1)
            YTEST = YP(II,I) + bondvectorthirdshell(J,2)
            ZTEST = ZP(II,I) + bondvectorthirdshell(J,3)

	  DO K = 1, CHAINLENGTH
              IF( (XTEST .EQ. XP(II,K)) .AND.
     >            (YTEST .EQ. YP(II,K)) .AND.
     >            (ZTEST .EQ. ZP(II,K)) ) THEN
                C = C + 1
              ENDIF
          ENDDO
	  ENDDO
	  ctotal3=ctotal3+real(C)
        ENDDO

        DO I = 1, CHAINLENGTH
	  d = 0
          DO J = 1, 162
            XTEST = XP(II,I) + bondvectorfourthshell(J,1)
            YTEST = YP(II,I) + bondvectorfourthshell(J,2)
            ZTEST = ZP(II,I) + bondvectorfourthshell(J,3)
	  DO K = 1, CHAINLENGTH
              IF( (XTEST .EQ. XP(II,K)) .AND.
     >            (YTEST .EQ. YP(II,K)) .AND.
     >            (ZTEST .EQ. ZP(II,K)) ) THEN
                D = D + 1
              ENDIF
          ENDDO
          ENDDO
	  dtotal3=dtotal3+real(d)
        ENDDO

        DO I = 1, CHAINLENGTH
	  e = 0
          DO J = 1, 252
            XTEST = XP(II,I) + bondvectorfifthshell(J,1)
            YTEST = YP(II,I) + bondvectorfifthshell(J,2)
            ZTEST = ZP(II,I) + bondvectorfifthshell(J,3)
	  DO K = 1, CHAINLENGTH
              IF( (XTEST .EQ. XP(II,K)) .AND.
     >            (YTEST .EQ. YP(II,K)) .AND.
     >            (ZTEST .EQ. ZP(II,K)) ) THEN
                E = E + 1
              ENDIF
          ENDDO
          ENDDO
	  etotal3=etotal3+real(e)
        ENDDO


      ENDDO

      atotal3=atotal3/real(CHAINNUMBER)/real(CHAINLENGTH)
      btotal3=btotal3/real(CHAINNUMBER)/real(CHAINLENGTH)
      ctotal3=ctotal3/real(CHAINNUMBER)/real(CHAINLENGTH)
      dtotal3=dtotal3/real(CHAINNUMBER)/real(CHAINLENGTH)
      etotal3=etotal3/real(CHAINNUMBER)/real(CHAINLENGTH)

      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *   S U B F O L D E D P O I N T                                   *
C  *                                                                 *
C  *******************************************************************
C  *                                                                 *
C  * Version 0.1        Date:           08/10/94                     *
C  *                    last modified   08/10/94                     *
C  *                                                                 *
C  * (C) copyright  Roland Rapold                                    *
C  *                Institute of Polymer Science                     *
C  *                University of Akron, OH 44325-3909               *
C  *                e-mail: rapold@polymer.uakron.edu                *
C  *                                                                 *
C  *******************************************************************
      SUBROUTINE SUBFOLDEDPOINT(X,Y,Z)

      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  X, Y, Z

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ


  103 IF( X .GE. LATTICESIZEX ) THEN
        X = X - LATTICESIZEX
        GO TO 103
      END IF
  104 IF( X .LT. 0 ) THEN
        X = X + LATTICESIZEX
        GO TO 104
      END IF

  105 IF( Y .GE. LATTICESIZEY ) THEN
        Y = Y - LATTICESIZEY
        GO TO 105
      END IF
  106 IF( Y .LT. 0 ) THEN
        Y = Y + LATTICESIZEY
        GO TO 106
      END IF

  107 IF( Z .GE. LATTICESIZEZ ) THEN
        Z = Z - LATTICESIZEZ
        GO TO 107
      END IF
  108 IF( Z .LT. 0 ) THEN
        Z = Z + LATTICESIZEZ
        GO TO 108
      END IF

      RETURN
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
      INTEGER INTERTYPE(12),INTERVECTOR(4,3)
      REAL     BONDLENGTH, MAGNITUDE
      REAL     DUMMYVECTOR(3), DUMMYVECTOR2(3)
      REAL     DIAMONDVECTOR(4,3), TMATRIXKD(3,3)

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
C23456789012345678901234567890123456789012345678901234567890123456789012
C --------------------------------------------------------------------
      SUBROUTINE SUBLOCALCHEMISTRY(CHAINLENGTH, NOS, LC, LCR) 
C --------------------------------------------------------------------
C
C  SUBROUTINE TO DETERMINE THE LOCAL CHEMISTRY OF THE CHAIN   
C  IN THE 2NND_RIS MODEL
C  *** NEEDS MODIFICATION FOR RING CHAINS, POE AND COPOLYMERS***

      INTEGER   I
      INTEGER   CHAINL
      PARAMETER (CHAINL=50)

      INTEGER   NOS
      INTEGER   LC(CHAINL)
      INTEGER	LCR(CHAINL)
      INTEGER   CHAINLENGTH

      IF( NOS.EQ.1 ) THEN
	DO I = 1, CHAINLENGTH
	  LC(I) = 1
	  LCR(I) = 1
	ENDDO
      ENDIF

      IF( NOS.EQ.3 ) THEN
	DO I = 1, CHAINLENGTH, 3
	  LC(I) = 4
	  LCR(I) = 4
	  IF((I+1).LE.CHAINLENGTH) THEN
		LC(I+1) = 2
		LCR(I+1) = 3
	  ENDIF
	  IF((I+2).LE.CHAINLENGTH) THEN 
		LC(I+2) = 3
		LCR(I+2) = 2
	  ENDIF
	ENDDO
      ENDIF
	
      RETURN
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
      PARAMETER (CHAINL=50,CHAINN=50)
      INTEGER  I, II
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  X(CHAINN,CHAINL), Y(CHAINN,CHAINL), Z(CHAINN,CHAINL)
      REAL     SNNDLENGTH
      REAL     DUMMYVECTOR(3), DUMMYVECTOR2(3)
      REAL     XD(CHAINN,CHAINL),YD(CHAINN,CHAINL),ZD(CHAINN,CHAINL)
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
      PARAMETER (CHAINL=50,CHAINN=50)
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
