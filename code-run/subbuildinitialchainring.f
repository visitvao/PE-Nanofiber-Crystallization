
      SUBROUTINE SUBBUILDINITIALCHAINRING(BONDVECTOR,ISEED)

      INTEGER  CHAINL, LATTICES, CHAINN

      PARAMETER (CHAINL=100,LATTICES=100,CHAINN=200)

      INTEGER  I, II, III, K, KMINUS, L
      INTEGER  ISEED
      INTEGER  XTEST, YTEST, ZTEST
      integer  xtest1, ytest1, ztest1
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  CHAINLENGTHHALF
      INTEGER  X(CHAINN,-2:CHAINL+3), Y(CHAINN,-2:CHAINL+3),
     >         Z(CHAINN,-2:CHAINL+3)
      INTEGER  XP(CHAINN,-2:CHAINL+3), YP(CHAINN,-2:CHAINL+3),
     >         ZP(CHAINN,-2:CHAINL+3)
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BEADOCCUPANCY(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)

      REAL     RAN1, DUMMY

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER, 
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      COMMON  /BEADOCCUPANCY/   BEADOCCUPANCY
      COMMON  /COORDINATES/     X, Y, Z, XP, YP, ZP

      IF( MOD(CHAINLENGTH,2) .EQ. 0 ) THEN
        CHAINLENGTHHALF = CHAINLENGTH / 2
      ELSE
        CHAINLENGTHHALF = (CHAINLENGTH - 1) / 2
      ENDIF

      DO II = 1, CHAINNUMBER
      
 
   10 DUMMY = RAN1(ISEED)*FLOAT(LATTICESIZEX)
      X(II,1) = INT(DUMMY)
      DUMMY = RAN1(ISEED)*FLOAT(LATTICESIZEY)
      Y(II,1) = INT(DUMMY)
      DUMMY = RAN1(ISEED)*FLOAT(LATTICESIZEZ)/4.0
      Z(II,1) = INT(DUMMY)
      IF(BEADOCCUPANCY(X(II,1),Y(II,1),Z(II,1)).EQ.1) then
       
        GO TO 10
      endif
      
      XP(II,1) = X(II,1)
      YP(II,1) = Y(II,1)
      ZP(II,1) = Z(II,1)

      BEADOCCUPANCY(X(II,1),Y(II,1),Z(II,1)) = 1
 
      DUMMY = RAN1(ISEED)*12.0
      K = INT(DUMMY) + 1

      DO I = 1, 12
        IF( ( BONDVECTOR(I,1) .EQ. - BONDVECTOR(K,1) ) .AND.
     >      ( BONDVECTOR(I,2) .EQ. - BONDVECTOR(K,2) ) .AND.
     >      ( BONDVECTOR(I,3) .EQ. - BONDVECTOR(K,3) ) ) THEN
          KMINUS = I
          GO TO 89
        ENDIF
      ENDDO
      WRITE(4,*) 'PROBLEM WITH REVERSE DIRECTION OF RING'
      
   89 DO I = 2, CHAINLENGTHHALF
        XTEST = XP(II,I-1) + BONDVECTOR(K,1)
        YTEST = YP(II,I-1) + BONDVECTOR(K,2)
        ZTEST = ZP(II,I-1) + BONDVECTOR(K,3)
        XP(II,I) = XTEST
        YP(II,I) = YTEST
        ZP(II,I) = ZTEST
                
        CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
c        IF((BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1).or.
c     >		(zp(ii,i).lt.0)) THEN

        IF((BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)) THEN
           iyes =  0
c          DO III = 1, I
   	   do iii = 1, chainlength
            BEADOCCUPANCY(X(II,III),Y(II,III),Z(II,III)) = 0
           ENDDO
          GO TO 10
        ENDIF
        X(II,I) = XTEST
        Y(II,I) = YTEST
        Z(II,I) = ZTEST
        BEADOCCUPANCY(XTEST,YTEST,ZTEST) = 1
      ENDDO
	
 	
   88 DUMMY = RAN1(ISEED)*12.0
      L = INT(DUMMY) + 1
      IF( ( L .NE. K ) .AND. ( L .NE. KMINUS ) ) THEN
        XTEST = XP(II,CHAINLENGTHHALF) + BONDVECTOR(L,1)
        YTEST = YP(II,CHAINLENGTHHALF) + BONDVECTOR(L,2)
        ZTEST = ZP(II,CHAINLENGTHHALF) + BONDVECTOR(L,3)
      ELSE
        GO TO 88
      ENDIF
      XP(II,CHAINLENGTHHALF+1) = XTEST
      YP(II,CHAINLENGTHHALF+1) = YTEST
      ZP(II,CHAINLENGTHHALF+1) = ZTEST
       
      CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
c      IF((BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)
c     >		.or.(zp(ii,chainlengthhalf+1).lt.0)) THEN

      IF((BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1))  THEN
c      	beadoccupancy(xtest,ytest,ztest) = 0
         GO TO 88
      ENDIF
      X(II,I) = XTEST
      Y(II,I) = YTEST
      Z(II,I) = ZTEST
      BEADOCCUPANCY(XTEST,YTEST,ZTEST) = 1

      DO I = CHAINLENGTHHALF+2, CHAINLENGTH
        XTEST = XP(II,I-1) + BONDVECTOR(KMINUS,1)
        YTEST = YP(II,I-1) + BONDVECTOR(KMINUS,2)
        ZTEST = ZP(II,I-1) + BONDVECTOR(KMINUS,3)
        XP(II,I) = XTEST
        YP(II,I) = YTEST
        ZP(II,I) = ZTEST
        CALL SUBFOLDEDPOINT(XTEST,YTEST,ZTEST)
        IF((BEADOCCUPANCY(XTEST,YTEST,ZTEST).EQ.1)
     >		.or.(zp(ii,i).lt.0)) THEN
c          DO III = 1, I
	   do iii = 1, chainlength
            BEADOCCUPANCY(X(II,III),Y(II,III),Z(II,III)) = 0
           ENDDO
          GO TO 10
        ENDIF
        X(II,I) = XTEST
        Y(II,I) = YTEST
        Z(II,I) = ZTEST
        BEADOCCUPANCY(XTEST,YTEST,ZTEST) = 1
      ENDDO


      DO I = 1, 12
        IF((( X(II,CHAINLENGTH) - X(II,1)) .EQ. BONDVECTOR(I,1)) .AND.
     >     (( Y(II,CHAINLENGTH) - Y(II,1)) .EQ. BONDVECTOR(I,2)) .AND.
     >     (( Z(II,CHAINLENGTH) - Z(II,1)) .EQ. BONDVECTOR(I,3))) THEN
          GO TO 200
        ENDIF
      ENDDO

ckk I add this segment; otherwise, there will always errors for ring      
      do iii=1,chainlength
  	beadoccupancy(x(ii,iii),y(ii,iii),z(ii,iii)) = 0       
      enddo
ckk I add

      GO TO 10      
	
  200 ENDDO
	 
	
      WRITE(6,*) 'CHAIN SUCCESSFULLY BUILT'
      
      count = 0
      do i = 0,latticesizez-1
        do j = 0, latticesizey-1
          do k = 0, latticesizex-1
	  if (beadoccupancy(k,j,i).eq.1) then
	     count = count + 1
	  endif          
          enddo
        enddo
      enddo
      print*,count
      
      count = 0
      do ii=1,chainnumber
      do i=1,chainlength
         if (beadoccupancy(x(ii,i),y(ii,i),z(ii,i)).eq.1) then
c            print*,ii,i,x(ii,i),y(ii,i),z(ii,i)
            count = count + 1
         endif
      enddo
      enddo
      print*,count
      RETURN
      END
