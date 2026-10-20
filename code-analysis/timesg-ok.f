
C --------------------------------------------------------------------
      PROGRAM SNNDANALYSE_2NND
C --------------------------------------------------------------------
C --------------------------------------------------------------------
C   This program calculate the orientation of longest principle axes 
C      between two chains.  SG ~ t
C      AND the intermolecular orientational correlation function
C      P2(r)~r
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
c      UNIT=40,FILE='TTTTPERCENT.OUT')
C      UNIT=41,FILE='SCORRELATION.OUT')

C --------------------------------------------------------------------

      INTEGER  CHAINL, CHAINN, MAXBIN, NSNAPSHOT
      PARAMETER (CHAINL=100,CHAINN=200,LATTICES=100,NSNAPSHOT=1000)

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
      INTEGER  NOS
      INTEGER  INTERTYPE(12),TYPE(CHAINN,CHAINL)

      REAL     BONDLENGTH, SNNDLENGTH
      REAL     ENERGYINTERFIRSTSHELL, ENERGYINTERSECONDSHELL
      REAL     ENERGYINTERTHIRDSHELL
      REAL     XPD(CHAINN,CHAINL), YPD(CHAINN,CHAINL),
     >         ZPD(CHAINN,CHAINL)

      REAL     TMATRIXKD(3,3),matrix2nnd(3),matrixcart(3)
      REAL     DIAMONDVECTOR(4,3)

      REAL    LONGAXIS(CHAINN,3), LONGSCALAR(CHAINN), a2, b2, c2
      REAL     RCMZ,RZNEW,RXNEW,RYNEW
      REAL     xyz(3,3),XX,YY,ZZ,DSCALAR(3),DVECTOR(3,3)
      REAL     sglobal,nsglobal

	
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


        OPEN(40,FILE='SGLOABL.OUT')
        OPEN(41,FILE='SGLOABL2.OUT')

     
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


      DO II = 1, CHAINNUMBER
        DO I = 2, CHAINLENGTH
          XP(II,I) = XP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),1)
          YP(II,I) = YP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),2)
          ZP(II,I) = ZP(II,I-1) + BONDVECTOR(BONDDIRECTION(II,I),3)
        ENDDO
      ENDDO


	DO II=1,CHAINNUMBER
        DO I = 1, CHAINLENGTH
	   matrix2nnd(1)=XP(II,I)
	   matrix2nnd(2)=YP(II,I)
	   matrix2nnd(3)=ZP(II,I)
          CALL SUBMATVEC(TMATRIXKD,matrix2nnd,matrixcart,3,3)
           XPD(II,I)=matrixcart(1)*SNNDLENGTH
	   YPD(II,I)=matrixcart(2)*SNNDLENGTH
	   ZPD(II,I)=matrixcart(3)*SNNDLENGTH
c	 write(999,*) XPD(II,I), YPD(II,I), ZPD(II,I)
	END DO
	END DO


        DO II = 1, CHAINNUMBER

	 DO I=1,3
	  LONGAXIS(II,I)=0.0
	 END DO
	 LONGSCALAR(II)=0.0

	  RXNEW = 0.0
	  RYNEW = 0.0
	  RZNEW = 0.0
	  DO I = 1, CHAINLENGTH
	    RXNEW = RXNEW + XPD(II,I)
	    RYNEW = RYNEW + YPD(II,I)
	    RZNEW = RZNEW + ZPD(II,I)
	  END DO
	  RXNEW = RXNEW / REAL (CHAINLENGTH)
	  RYNEW = RYNEW / REAL (CHAINLENGTH)
	  RZNEW = RZNEW / REAL (CHAINLENGTH)

	  xyz(1,1)=0.0
	  xyz(1,2)=0.0
	  xyz(1,3)=0.0
	  xyz(2,2)=0.0
	  xyz(2,3)=0.0
	  xyz(3,3)=0.0

	  DO I = 1, CHAINLENGTH
	    XX= XPD(II,I)-RXNEW
	    YY= YPD(II,I)-RYNEW
	    ZZ= ZPD(II,I)-RZNEW
	    xyz(1,1)=xyz(1,1)+XX*XX
	    xyz(1,2)=xyz(1,2)+XX*YY
	    xyz(1,3)=xyz(1,3)+XX*ZZ
	    xyz(2,2)=xyz(2,2)+YY*YY
	    xyz(2,3)=xyz(2,3)+YY*ZZ
	    xyz(3,3)=xyz(3,3)+ZZ*ZZ
	  END DO
	  xyz(1,1)=xyz(1,1)/REAL(CHAINLENGTH)
	  xyz(1,2)=xyz(1,2)/REAL(CHAINLENGTH)
	  xyz(1,3)=xyz(1,3)/REAL(CHAINLENGTH)
	  xyz(2,2)=xyz(2,2)/REAL(CHAINLENGTH)
	  xyz(2,3)=xyz(2,3)/REAL(CHAINLENGTH)
	  xyz(3,3)=xyz(3,3)/REAL(CHAINLENGTH)

	  xyz(2,1)=xyz(1,2)
	  xyz(3,1)=xyz(1,3)
	  xyz(3,2)=xyz(2,3)

	  CALL DIAGONALIZATION(XYZ,DSCALAR,DVECTOR)
	 
	 DO I=1,3
	  LONGAXIS(II,I)=DVECTOR(1,I)
	 END DO
	 LONGSCALAR(II)=DSCALAR(1)

	END DO


	  sglobal=0.0
	  nsglobal=0.0

        DO II = 1, CHAINNUMBER
        DO II1 = 1, CHAINNUMBER
	 IF (II.LT.II1) THEN
	  c2=LONGAXIS(II,1)*LONGAXIS(II1,1)
     & +LONGAXIS(II,2)*LONGAXIS(II1,2)+LONGAXIS(II,3)*LONGAXIS(II1,3)

	  sglobal=sglobal+(c2/LONGSCALAR(II)/LONGSCALAR(II1))**2.0
	  nsglobal=nsglobal+1.0
	 ENDIF
	END DO
	END DO

 	 WRITE(40,*) CURRENTMCM,1.5*sglobal/nsglobal-0.5
 	 WRITE(41,*) CURRENTMCM-STARTPOINT,1.5*sglobal/nsglobal-0.5








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



	DO II=1,CHAINNUMBER
        DO I = 1, CHAINLENGTH
	   matrix2nnd(1)=XP(II,I)
	   matrix2nnd(2)=YP(II,I)
	   matrix2nnd(3)=ZP(II,I)
          CALL SUBMATVEC(TMATRIXKD,matrix2nnd,matrixcart,3,3)
           XPD(II,I)=matrixcart(1)*SNNDLENGTH
	   YPD(II,I)=matrixcart(2)*SNNDLENGTH
	   ZPD(II,I)=matrixcart(3)*SNNDLENGTH
	END DO
	END DO



        DO II = 1, CHAINNUMBER

	 DO I=1,3
	  LONGAXIS(II,I)=0.0
	 END DO
	 LONGSCALAR(II)=0.0

	  RXNEW = 0.0
	  RYNEW = 0.0
	  RZNEW = 0.0
	  DO I = 1, CHAINLENGTH
	    RXNEW = RXNEW + XPD(II,I)
	    RYNEW = RYNEW + YPD(II,I)
	    RZNEW = RZNEW + ZPD(II,I)
	  END DO
	  RXNEW = RXNEW / REAL (CHAINLENGTH)
	  RYNEW = RYNEW / REAL (CHAINLENGTH)
	  RZNEW = RZNEW / REAL (CHAINLENGTH)

	  xyz(1,1)=0.0
	  xyz(1,2)=0.0
	  xyz(1,3)=0.0
	  xyz(2,2)=0.0
	  xyz(2,3)=0.0
	  xyz(3,3)=0.0

	  DO I = 1, CHAINLENGTH
	    XX= XPD(II,I)-RXNEW
	    YY= YPD(II,I)-RYNEW
	    ZZ= ZPD(II,I)-RZNEW
	    xyz(1,1)=xyz(1,1)+XX*XX
	    xyz(1,2)=xyz(1,2)+XX*YY
	    xyz(1,3)=xyz(1,3)+XX*ZZ
	    xyz(2,2)=xyz(2,2)+YY*YY
	    xyz(2,3)=xyz(2,3)+YY*ZZ
	    xyz(3,3)=xyz(3,3)+ZZ*ZZ
	  END DO
	  xyz(1,1)=xyz(1,1)/REAL(CHAINLENGTH)
	  xyz(1,2)=xyz(1,2)/REAL(CHAINLENGTH)
	  xyz(1,3)=xyz(1,3)/REAL(CHAINLENGTH)
	  xyz(2,2)=xyz(2,2)/REAL(CHAINLENGTH)
	  xyz(2,3)=xyz(2,3)/REAL(CHAINLENGTH)
	  xyz(3,3)=xyz(3,3)/REAL(CHAINLENGTH)

	  xyz(2,1)=xyz(1,2)
	  xyz(3,1)=xyz(1,3)
	  xyz(3,2)=xyz(2,3)

	  CALL DIAGONALIZATION(XYZ,DSCALAR,DVECTOR)
	 
	 DO I=1,3
	  LONGAXIS(II,I)=DVECTOR(1,I)
	 END DO
	 LONGSCALAR(II)=DSCALAR(1)

	END DO


	  sglobal=0.0
	  nsglobal=0.0
        DO II = 1, CHAINNUMBER
        DO II1 = 1, CHAINNUMBER
	 IF (II.LT.II1) THEN
	  c2=LONGAXIS(II,1)*LONGAXIS(II1,1)
     & +LONGAXIS(II,2)*LONGAXIS(II1,2)+LONGAXIS(II,3)*LONGAXIS(II1,3)
	  
	  sglobal=sglobal+(c2/LONGSCALAR(II)/LONGSCALAR(II1))**2.0

	  nsglobal=nsglobal+1.0

	 ENDIF
	END DO
	END DO

 	 WRITE(40,*) CURRENTMCM,1.5*sglobal/nsglobal-0.5
 	 WRITE(41,*) CURRENTMCM-STARTPOINT,1.5*sglobal/nsglobal-0.5


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


C*****************************************************************************
C  Eigenvalues/Matrix Diagonalization
C
C The following routines (based on the original EISPACK library) perform a 
C diagonalization of a real symmetric matrix based on the QL algorithm. The
C first step is to reduce the matrix to tridiagonal form (TRED2) and then a 
C further routine (TQLI) finds the eigenvalues and eigenvectors of a tridiagonal
C matrix: 
C  (Fortran 77 version) 
C      PROGRAM D11R4
C     Driver for routine TQLI
C ***************************************************************************
C *** C(NP,NP) original value, F eigenvector, D scale eigenvalue *****

      SUBROUTINE DIAGONALIZATION(C,DD,FF)
      PARAMETER(NP=3,TINY=1.0E-6)
      DIMENSION A(NP,NP),C(NP,NP),D(NP),E(NP),F(NP,NP),DD(NP),FF(NP,NP)
      INTEGER LARGE(NP)
      REAL LARGEST,SMALLEST

      DO 12 I=1,NP
            DO 11 J=1,NP
                  A(I,J)=C(I,J)
11          CONTINUE
12    CONTINUE

      CALL TRED2(A,NP,NP,D,E)
      CALL TQLI(D,E,NP,NP,A)

      DO 16 I=1,NP
            DO 14 J=1,NP
                  F(I,J)=0.0
                  DO 13 K=1,NP
                        F(I,J)=F(I,J)+C(J,K)*A(K,I)
13                CONTINUE
14          CONTINUE
16    CONTINUE


C *** arrange the last value *******
	LARGEST=ABS(D(1))
	SMALLEST=ABS(D(1))
	LARGE(1)=1
	LARGE(3)=1
	DO I = 2,NP
	   IF (LARGEST.LT.ABS(D(I))) THEN
		LARGEST=ABS(D(I))
		LARGE(1)=I
	   ENDIF
	   IF (SMALLEST.GT.ABS(D(I))) THEN
		SMALLEST=ABS(D(I))
		LARGE(3)=I
	   ENDIF
	END DO
	LARGE(2)=6-LARGE(3)-LARGE(1)

	DO I=1,NP
	  DD(I)=D(LARGE(I))
	  DO J=1,NP
	   FF(I,J)=F(LARGE(I),J)
	  END DO
	END DO
C FF j is x, y or z *************** 
C *** arrange the last value *******


      RETURN
      END

      SUBROUTINE TRED2(A,N,NP,D,E)
      DIMENSION A(NP,NP),D(NP),E(NP)
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1  
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END

      SUBROUTINE TQLI(D,E,N,NP,Z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
C ************ END OF DIAGONALIZATION OF MATRIX **************************

