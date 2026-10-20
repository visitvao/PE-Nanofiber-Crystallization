
C --------------------------------------------------------------------
      PROGRAM SNNDANALYSE_2NND
C --------------------------------------------------------------------
C --------------------------------------------------------------------
c      calculate the number and percentage of the continuous bondvector
C --------------------------------------------------------------------
      INTEGER  CHAINL, CHAINN, MAXCONBOND
      PARAMETER (CHAINL=100,CHAINN=50,MAXCONBOND=20)

      INTEGER  REALLENGTH
      INTEGER  NFILE, NDATA
      INTEGER  I, II, J, K, KK, JJ
      INTEGER  TEMP, IFRING,IFILM, IFCOLLAPSE, IFEXCLUDEDVOLUME
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
      INTEGER  BONDDIRECTION(CHAINN,CHAINL)
      INTEGER  NOS
      REAL     BONDLENGTH, SNNDLENGTH
      REAL     ENERGYINTERFIRSTSHELL, ENERGYINTERSECONDSHELL
      REAL     ENERGYINTERTHIRDSHELL
      INTEGER  MCSTRANSGAUCH
      INTEGER  nbcont,BONDPREV
      real bondcon(MAXCONBOND),avgbondcon(MAXCONBOND),TOTALBOND

      CHARACTER  GARB1*50
      CHARACTER  POLYMER*18, LITERATURE*47
      CHARACTER  FILENAME19*12,FILENAME20*12,FILENAME21*19


      IDATA=0

      DO I=1,MAXCONBOND
        bondcon(I)=0.0
	avgbondcon(I)=0.0
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

C  INITIALIZE
      BONDLENGTH = 1.53
      SNNDLENGTH = BONDLENGTH * SQRT(8./3.)

      WRITE(6,*) 'BONDLENGTH = ',BONDLENGTH
      WRITE(6,*) 'SNNDLENGTH = ',SNNDLENGTH
      WRITE(6,*) ' '

      WRITE(6,*) 'number of conformation BEGIN to be analyzed?'
      READ(5,*)  nbegin
      WRITE(6,*) 'nbegin = ',nbegin
      WRITE(6,*) ' '

      WRITE(6,*) 'NUMBER OF CONFORMATIONS FILES TO BE ANALYZED ?'
      READ(5,*)  NFILE
      WRITE(6,*) 'NFILE = ',NFILE
      WRITE(6,*) ' '

      WRITE(6,*) 'NUMBER OF CHAINMOVES IN EACH CONFORMATIONS FILE ?'
      READ(5,*)  CHAINMOVES
      WRITE(6,*) 'CHAINMOVES = ',CHAINMOVES
      WRITE(6,*) ' '
      TOTALMOVES = (NFILE-nbegin+1)*CHAINMOVES
      WRITE(6,*) 'TOTAL SIMULATIONLENGTH = ',TOTALMOVES,' MCS'
      WRITE(6,*) ' '

      WRITE(6,*) 'EQUILIBRATION TIME (NOT USED FOR STATISTICS) ?'
      WRITE(6,*) ' IN MONTE-CARLO STEPS'
      READ(5,9)  EQUILIBRATIONTIME
       EQUILIBRATIONTIME=0
      WRITE(6,*) ' '


      FILENAME21 = 'conpe2001.dat'
      DO II = nbegin, NFILE
	J = 10+II
        WRITE(FILENAME21(8:9),17) II
        OPEN(UNIT=J,FILE=FILENAME21)
      ENDDO
      

C  ========================================
C  ///  BIG LOOP FOR DYNAMIC-MC MOVES  \\\
C  ----------------------------------------

      DO KK = nbegin, NFILE

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
c     READ(10+KK,8)  GARB1
	
      READ(10+KK,*) CURRENTMCM

      DO II = 1, CHAINNUMBER
        READ(10+KK,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
        READ(10+KK,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)
      ENDDO

      SIMULATIONLENGTH = CURRENTMCM

C ///// initialize ////////////////////
       IF (KK.EQ.nbegin) THEN 

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

       ENDIF
C ///// initialize ////////////////////


      DO JJ = 1, NDATA

	READ(10+KK,*) CURRENTMCM
c	write(77,*) CURRENTMCM
        DO II = 1, CHAINNUMBER
          READ(10+KK,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
          READ(10+KK,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)
c        write(77,10) ICHAINNUMBER, XP(II,1), YP(II,1), ZP(II,1)
c        write(77,13) (BONDDIRECTION(II,I), I=2,CHAINLENGTH)

        ENDDO


        IF((EQUILIBRATIONTIME+STARTPOINT).GE.CURRENTMCM) THEN
          SIMULATIONLENGTH = CURRENTMCM
c 	  GOTO 60
	ENDIF


C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      IDATA = IDATA + 1

      IF(STOREINTERVAL .NE. CURRENTMCM-SIMULATIONLENGTH) THEN
	WRITE(4,*) 'PROBLEM WITH STOREINTERVAL'
      	STOP 3
      ENDIF

      SIMULATIONLENGTH = CURRENTMCM


       DO I=1,MAXCONBOND
	  bondcon(I)=0
         avgbondcon(I)=0.0 
       END DO

       DO I = 1, CHAINNUMBER
	   nbcont=1
	   BONDPREV=BONDDIRECTION(I,2)
         DO J = 3, CHAINLENGTH
	  IF (BONDDIRECTION(I,J).EQ.BONDPREV) THEN
	    nbcont=nbcont+1
	  ELSE
	    IF (nbcont.GT.MAXCONBOND) THEN
		print *,'nbcont>MAXCONBOND'
	    ELSE
	      bondcon(nbcont)=bondcon(nbcont)+1
	    ENDIF

           nbcont=1
	   BONDPREV=BONDDIRECTION(I,J)

	  ENDIF
	END DO
        bondcon(nbcont)=bondcon(nbcont)+1
      END DO

      DO I=1,MAXCONBOND
       avgbondcon(I)= (avgbondcon(I)*(IDATA-1)+bondcon(I))/real(IDATA)
      END DO

C ********************************************************************

C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   60 ENDDO
      ENDDO

C  ----------------------------------------
C  \\\  BIG LOOP FOR DYNAMIC-MC MOVES  ///
C  ========================================

      TOTALBOND=0.0
      DO I=1,MAXCONBOND
	TOTALBOND=TOTALBOND+I*avgbondcon(I)
      END DO

      OPEN(51,FILE='CONBOND.OUT')
      DO I=1,MAXCONBOND
       WRITE(51,*)I,avgbondcon(I),I*avgbondcon(I)
     &    ,I*avgbondcon(I)/TOTALBOND
      END DO

C  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      END
