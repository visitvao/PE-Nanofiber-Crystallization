
      SUBROUTINE SUBFORMDROPLET(X,XP)

      INTEGER  CHAINL, LATTICES, CHAINN

      PARAMETER (CHAINL=100,LATTICES=100,CHAINN=200)

      INTEGER  I, II, ICOUNTP, ICOUNTM, ICOUNTLEFT
      INTEGER  CHAINLENGTH
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  CHAINNUMBER
      INTEGER  X(CHAINN,-2:CHAINL+3)
      INTEGER  XP(CHAINN,-2:CHAINL+3)
      INTEGER  ILABEL(CHAINN), XDUM

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER, 
     >                  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ

c-- Label for every chain as 1 [0 means the chain is selected]
c-- start with parent chain at bead position 1
c-- ZP(II,I) is the parent coordinate
c-- Z(II,I)  is the foled chain

      DO II = 1, CHAINNUMBER
        ILABEL(II) = 1
        XP(II,1) = X(II,1)
      DO I = 2, CHAINLENGTH

        XDUM = X(II,I)

88      IF(XDUM-XP(II,I-1).LT.-1) THEN 
                XDUM = XDUM + LATTICESIZEX
                GO TO 88
        ENDIF

99      IF(XDUM-XP(II,I-1).GT.1) THEN 
                XDUM = XDUM - LATTICESIZEX
                GO TO 99
        ENDIF

        XP(II,I) = XDUM

      ENDDO
      ENDDO

        
      DO II = 1, CHAINNUMBER
        ICOUNTP = 0
        ICOUNTM = 0
      DO I = 1, CHAINLENGTH
        IF(XP(II,I).GE.LATTICESIZEX) ICOUNTP = ICOUNTP + 1
        IF(XP(II,I).LT.0) ICOUNTM = ICOUNTM + 1
      ENDDO
        ICOUNTLEFT = CHAINLENGTH - ICOUNTP - ICOUNTM
        IF(ICOUNTLEFT.GE.(CHAINLENGTH/2)) ILABEL(II) = 0
        IF(ICOUNTLEFT.GT.ICOUNTM.AND.ICOUNTLEFT.GT.ICOUNTP)
     >          ILABEL(II) = 0
      ENDDO

      DO II = 1, CHAINNUMBER
        IF(ILABEL(II).EQ.1) THEN
          XP(II,1) = X(II,1) - LATTICESIZEX
        ELSE 
          GOTO 77
        ENDIF
      DO I = 2, CHAINLENGTH

        XDUM = X(II,I)

188     IF(XDUM-XP(II,I-1).LT.-1) THEN 
                XDUM = XDUM + LATTICESIZEX
                GO TO 188
        ENDIF

199     IF(XDUM-XP(II,I-1).GT.1) THEN 
                XDUM = XDUM - LATTICESIZEX
                GO TO 199
        ENDIF

        XP(II,I) = XDUM

      ENDDO
77    ENDDO

      DO II = 1, CHAINNUMBER
        IF(ILABEL(II).EQ.1) THEN
        ICOUNTP = 0
        ICOUNTM = 0
        DO I = 1, CHAINLENGTH
          IF(XP(II,I).GE.LATTICESIZEX) ICOUNTP = ICOUNTP + 1
          IF(XP(II,I).LT.0) ICOUNTM = ICOUNTM + 1
        ENDDO
        ICOUNTLEFT = CHAINLENGTH - ICOUNTP - ICOUNTM
        IF(ICOUNTLEFT.GE.(CHAINLENGTH/2)) ILABEL(II) = 0
        IF(ICOUNTLEFT.GT.ICOUNTM.AND.ICOUNTLEFT.GT.ICOUNTP)
     >          ILABEL(II) = 0
        ENDIF
      ENDDO

      DO II = 1, CHAINNUMBER
        IF(ILABEL(II).EQ.1) THEN
          XP(II,1) = X(II,1) + LATTICESIZEX
        ELSE 
          GOTO 177
        ENDIF
      DO I = 2, CHAINLENGTH

        XDUM = X(II,I)

888     IF(XDUM-XP(II,I-1).LT.-1) THEN 
                XDUM = XDUM + LATTICESIZEX
                GO TO 888
        ENDIF

999     IF(XDUM-XP(II,I-1).GT.1) THEN 
                XDUM = XDUM - LATTICESIZEX
                GO TO 999
        ENDIF

        XP(II,I) = XDUM

      ENDDO
177    ENDDO

      DO II = 1, CHAINNUMBER
        IF(ILABEL(II).EQ.1) THEN
        ICOUNTP = 0
        ICOUNTM = 0
        DO I = 1, CHAINLENGTH
          IF(XP(II,I).GE.LATTICESIZEX) ICOUNTP = ICOUNTP + 1
          IF(XP(II,I).LT.0) ICOUNTM = ICOUNTM + 1
        ENDDO
        ICOUNTLEFT = CHAINLENGTH - ICOUNTP - ICOUNTM
        IF(ICOUNTLEFT.GE.(CHAINLENGTH/2)) ILABEL(II) = 0
        IF(ICOUNTLEFT.GT.ICOUNTM.AND.ICOUNTLEFT.GT.ICOUNTP)
     >          ILABEL(II) = 0
        ENDIF
      ENDDO

      DO II = 1, CHAINNUMBER
        IF(ILABEL(II).NE.0) THEN
          WRITE(15,*)'PROBLEM WITH FORMING THE FILM'
          STOP 11
        ENDIF
      ENDDO
        
      RETURN
      END
