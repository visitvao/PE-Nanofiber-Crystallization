
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
C      UNIT=33,FILE='SCORRELATION.OUT')

C --------------------------------------------------------------------

      INTEGER  CHAINL, CHAINN, MAXBIN, NSNAPSHOT
      PARAMETER (CHAINL=100,CHAINN=200,LATTICES=100,NSNAPSHOT=100)

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
      INTEGER  X(CHAINN,CHAINL)
      INTEGER  Y(CHAINN,CHAINL)
      INTEGER  Z(CHAINN,CHAINL)
      INTEGER  XP(CHAINN,CHAINL)
      INTEGER  YP(CHAINN,CHAINL)
      INTEGER  ZP(CHAINN,CHAINL)
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      INTEGER  bondvectorthirdshell(92,3)
      integer  bondvectorfourthshell(162,3)
      integer  bondvectorfifthshell(252,3)
      INTEGER  BONDDIRECTION(CHAINN,CHAINL)
      REAL     REALBONDVECTOR(12,3)
      INTEGER  NOS
      INTEGER  LC(CHAINL), LCR(CHAINL)
      INTEGER  INTERTYPE(12),TYPE(CHAINN,CHAINL)
      INTEGER  SITECHAIN(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)
      INTEGER  SITEBEAD(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)
      INTEGER  XPOINT,YPOINT,ZPOINT,XPOINT2,YPOINT2,ZPOINT2,NCOS
      REAL     COSTHETA,ALLCOSTHETA,COSTHETAMEAN(CHAINN,CHAINL),ALLCOS
      REAL     OVERALLCOSTHETA(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)
      REAL     REALSITE(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)

      REAL     aside2,bside2,cside2,COSSITE,aII,aJJ,aKK
      INTEGER  BONDDIRECTION1,BONDDIRECTION2,NCOSITE


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
      REAL     TMATRIXKD(3,3),matrix2nnd(3),matrixcart(3)
      REAL     DIAMONDVECTOR(4,3)

      INTEGER REALBONDSTATEL(12,12)
      INTEGER REALBONDSTATER(12,12)

      CHARACTER  GARB1*50
      CHARACTER  POLYMER*18, LITERATURE*47
      CHARACTER  FILENAME19*12,FILENAME20*12,FILENAME21*19

      COMMON  /TRANSFORMMATRIX/ TMATRIXKD
      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      IDATA2=0
      IDATA=0

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


	DO I=0,LATTICES-1
	DO J=0,LATTICES-1
	DO K=0,LATTICES-1
	   OVERALLCOSTHETA(I,J,K)=0.0
	END DO
	END DO
	END DO


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
      CALL SUBTRANSFORMMATRIX
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C  ----------------------
C  ----------------------
      CALL SUBGETBONDVECTORS(BONDVECTOR)

      CALL SUBGETBONDVECTORSECONDSHELL(BONDVECTOR,BONDVECTORSECONDSHELL)

      CALL SUBGETBONDVECTORTHIRDSHELL(BONDVECTOR,
     >       BONDVECTORSECONDSHELL, bondvectorthirdshell)

C     ...................................
C     NO COMMON BLOCK
C  ----------------------

C  ----------------------
c      CALL SUBINTERBEADVECTORS(BONDLENGTH,INTERTYPE,DIAMONDVECTOR) 
C     ...................................
C     COMMON  /TRANSFORMMATRIX/ TMATRIXKD
C  ----------------------


	DO I=1,12
	   matrix2nnd(1)=BONDVECTOR(I,1)
	   matrix2nnd(2)=BONDVECTOR(I,2)
	   matrix2nnd(3)=BONDVECTOR(I,3)
          CALL SUBMATVEC(TMATRIXKD,matrix2nnd,matrixcart,3,3)
           REALBONDVECTOR(I,1)=matrixcart(1)*SNNDLENGTH
	   REALBONDVECTOR(I,2)=matrixcart(2)*SNNDLENGTH
	   REALBONDVECTOR(I,3)=matrixcart(3)*SNNDLENGTH

	 WRITE (41,*) REALBONDVECTOR(I,1),REALBONDVECTOR(I,2)
     &    ,REALBONDVECTOR(I,3)
	END DO


      FILENAME21 = 'conpe2001.dat'
      DO II = 1, NFILE
	J = 10+II
        WRITE(FILENAME21(8:9),17) II
        OPEN(UNIT=J,FILE=FILENAME21)
      ENDDO


        OPEN(40,FILE='tn_order.OUT')
        OPEN(41,FILE='tn_order2.OUT')

     
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


      DO II = 1, CHAINNUMBER
        DO I = 2, CHAINLENGTH
          XP(II,I) = XP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),1)
          YP(II,I) = YP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),2)
          ZP(II,I) = ZP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),3)
        ENDDO
      ENDDO

      CALL SUBFOLDEDCHAIN(XP,YP,ZP,X,Y,Z)

      DO I = 0, LATTICESIZEX - 1
       DO J = 0, LATTICESIZEY - 1
        DO K = 0, LATTICESIZEZ - 1
	 SITECHAIN(I,J,K)=0
	 SITEBEAD(I,J,K)=0
	END DO
       END DO
      END DO

      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
          SITECHAIN(X(II,I),Y(II,I),Z(II,I)) = II
          SITEBEAD(X(II,I),Y(II,I),Z(II,I)) = I
        ENDDO
      ENDDO

	
      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH
	   XPOINT=XP(II,I)
	   YPOINT=YP(II,I)
	   ZPOINT=ZP(II,I)

	   IF (I.EQ.1) THEN
		BONDDIRECTION1=BONDDIRECTION(II,2)
	   ELSE
		BONDDIRECTION1=BONDDIRECTION(II,I)
	   ENDIF

	   NCOS=0
	   ALLCOSTHETA=0.0
	   COSTHETAMEAN(II,I)=0.0


C ********************* SECOND shell *****************************
	   DO J=1,42

	     XPOINT2=XPOINT+BONDVECTORSECONDSHELL(J,1)
	     YPOINT2=YPOINT+BONDVECTORSECONDSHELL(J,2)
	     ZPOINT2=ZPOINT+BONDVECTORSECONDSHELL(J,3)
	
	     CALL SUBFOLDEDPOINT(XPOINT2,YPOINT2,ZPOINT2)

             IF (SITECHAIN(XPOINT2,YPOINT2,ZPOINT2).NE.0
     &         .AND.SITECHAIN(XPOINT2,YPOINT2,ZPOINT2).NE.II) THEN

	       IF (SITEBEAD(XPOINT2,YPOINT2,ZPOINT2).NE.1) THEN
		BONDDIRECTION2=
     &		  BONDDIRECTION(SITECHAIN(XPOINT2,YPOINT2,ZPOINT2)
     &                          ,SITEBEAD(XPOINT2,YPOINT2,ZPOINT2))
	       ELSE 
        	BONDDIRECTION2=BONDDIRECTION
     &    	   (SITECHAIN(XPOINT2,YPOINT2,ZPOINT2),2)
	       ENDIF


		IF (BONDDIRECTION1.EQ.BONDDIRECTION2) THEN
		   COSTHETA=1.0
	        ELSE

	          aside2=REALBONDVECTOR(BONDDIRECTION1,1)**2
     &                  +REALBONDVECTOR(BONDDIRECTION1,2)**2
     &                  +REALBONDVECTOR(BONDDIRECTION1,3)**2

	          bside2=REALBONDVECTOR(BONDDIRECTION2,1)**2
     &                  +REALBONDVECTOR(BONDDIRECTION2,2)**2
     &                  +REALBONDVECTOR(BONDDIRECTION2,3)**2

                  cside2=(REALBONDVECTOR(BONDDIRECTION1,1)-
     &                    REALBONDVECTOR(BONDDIRECTION2,1))**2
     &                  +(REALBONDVECTOR(BONDDIRECTION1,2)-
     &                    REALBONDVECTOR(BONDDIRECTION2,2))**2
     &                  +(REALBONDVECTOR(BONDDIRECTION1,3)-
     &                    REALBONDVECTOR(BONDDIRECTION2,3))**2

	   COSTHETA=(aside2+bside2-cside2)/2.0/SQRT(aside2*bside2)

c	     print *,2,aside2,bside2,cside2

	        ENDIF

		NCOS=NCOS+1
	        ALLCOSTHETA=ALLCOSTHETA+COSTHETA**2

	     ENDIF
	  END DO
C ********************* SECOND shell *****************************

C ********************* Third shell ******************************
	   DO J=1,92

	     XPOINT2=XPOINT+BONDVECTORTHIRDSHELL(J,1)
	     YPOINT2=YPOINT+BONDVECTORTHIRDSHELL(J,2)
	     ZPOINT2=ZPOINT+BONDVECTORTHIRDSHELL(J,3)
	
	     CALL SUBFOLDEDPOINT(XPOINT2,YPOINT2,ZPOINT2)

             IF (SITECHAIN(XPOINT2,YPOINT2,ZPOINT2).NE.0
     &         .AND.SITECHAIN(XPOINT2,YPOINT2,ZPOINT2).NE.II) THEN

	       IF (SITEBEAD(XPOINT2,YPOINT2,ZPOINT2).NE.1) THEN
		BONDDIRECTION2=
     &		  BONDDIRECTION(SITECHAIN(XPOINT2,YPOINT2,ZPOINT2)
     &                          ,SITEBEAD(XPOINT2,YPOINT2,ZPOINT2))
	       ELSE 
	        BONDDIRECTION2=
     &   	  BONDDIRECTION(SITECHAIN(XPOINT2,YPOINT2,ZPOINT2),2)
	       ENDIF

		IF (BONDDIRECTION1.EQ.BONDDIRECTION2) THEN
		   COSTHETA=1.0
	        ELSE
	          aside2=REALBONDVECTOR(BONDDIRECTION1,1)**2
     &                  +REALBONDVECTOR(BONDDIRECTION1,2)**2
     &                  +REALBONDVECTOR(BONDDIRECTION1,3)**2

	          bside2=REALBONDVECTOR(BONDDIRECTION2,1)**2
     &                  +REALBONDVECTOR(BONDDIRECTION2,2)**2
     &                  +REALBONDVECTOR(BONDDIRECTION2,3)**2

                  cside2=(REALBONDVECTOR(BONDDIRECTION1,1)-
     &                    REALBONDVECTOR(BONDDIRECTION2,1))**2
     &                  +(REALBONDVECTOR(BONDDIRECTION1,2)-
     &                    REALBONDVECTOR(BONDDIRECTION2,2))**2
     &                  +(REALBONDVECTOR(BONDDIRECTION1,3)-
     &                    REALBONDVECTOR(BONDDIRECTION2,3))**2

	   COSTHETA=(aside2+bside2-cside2)/2.0/SQRT(aside2*bside2)

c	     print *,3,aside2,bside2,cside2

	        ENDIF

		NCOS=NCOS+1
	        ALLCOSTHETA=ALLCOSTHETA+COSTHETA**2

	     ENDIF
	  END DO
C ********************* Third shell ******************************

	     IF (NCOS.NE.0) THEN
	       COSTHETAMEAN(II,I)=ALLCOSTHETA/REAL(NCOS)*1.5-0.5
	     ELSE
	       COSTHETAMEAN(II,I)=0.0
	     ENDIF

        ENDDO
      ENDDO



C  ========================================
c **** average for all beads ***********


	 NCOSITE=0
	 COSSITE=0.0

         DO II = 1, CHAINNUMBER
           DO I = 1, CHAINLENGTH
	    NCOSITE=NCOSITE+1
	    COSSITE=COSSITE+COSTHETAMEAN(II,I)
	   END DO
         END DO

	  IF (NCOSITE.NE.0) THEN
	    COSSITE=COSSITE/NCOSITE
	  ENDIF


 	 WRITE(40,*) CURRENTMCM,COSSITE
 	 WRITE(41,*) CURRENTMCM-STARTPOINT,COSSITE


c **** average for every lattice point ***********
C  ========================================

       END DO
       END DO

	CLOSE(40)
	CLOSE(41)


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
      PARAMETER (CHAINL=100, CHAINN=200)
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

C --------------------------------------------------------------------
      SUBROUTINE SUBLOCALCHEMISTRY(CHAINLENGTH, NOS, LC, LCR) 
C --------------------------------------------------------------------
C
C  SUBROUTINE TO DETERMINE THE LOCAL CHEMISTRY OF THE CHAIN   
C  IN THE 2NND_RIS MODEL
C  *** NEEDS MODIFICATION FOR RING CHAINS, POE AND COPOLYMERS***

      INTEGER   I
      INTEGER   CHAINL
      PARAMETER (CHAINL=100)

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
      PARAMETER (CHAINL=100,CHAINN=200)
      INTEGER  I, II
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  X(CHAINN,CHAINL)
      INTEGER  Y(CHAINN,CHAINL)
      INTEGER  Z(CHAINN,CHAINL)
      REAL  XD(CHAINN,CHAINL)
      REAL  YD(CHAINN,CHAINL)
      REAL  ZD(CHAINN,CHAINL)
      REAL     SNNDLENGTH
      REAL     DUMMYVECTOR(3), DUMMYVECTOR2(3)
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
      SUBROUTINE SUBFOLDEDCHAIN(XP,YP,ZP,X,Y,Z)
      INTEGER  CHAINL, LATTICES, CHAINN
      PARAMETER(CHAINL=100, LATTICES=100, CHAINN=200)

      INTEGER  I, II
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  CHAINNUMBER, CHAINLENGTH
      INTEGER  X(CHAINN,CHAINL)
      INTEGER  Y(CHAINN,CHAINL)
      INTEGER  Z(CHAINN,CHAINL)
      INTEGER  XP(CHAINN,CHAINL)
      INTEGER  YP(CHAINN,CHAINL)
      INTEGER  ZP(CHAINN,CHAINL)

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER, 
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ

      DO II = 1, CHAINNUMBER
       DO I = 1, CHAINLENGTH

        X(II,I) = XP(II,I)
        Y(II,I) = YP(II,I)
        Z(II,I) = ZP(II,I)

        CALL SUBFOLDEDPOINT(X(II,I),Y(II,I),Z(II,I))

       ENDDO
      ENDDO

      RETURN
      END
C  *******************************************************************


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
C  **********************************************************************
