
      SUBROUTINE SUBFOLDEDCHAIN_FIBER(XP,YP,ZP,AX,AY,AZ,BONDVECTOR,
     >	BONDVECTORSECONDSHELL,bondvectorthirdshell,
     >	FIRSTOCCUP,SECONDOCCUP,THIRDOCCUP,FOURTHOCCUP,
     >	bondvectorfourthshell,IFOCC4TH,IFRING)
      INTEGER  CHAINL, CHAINN

      PARAMETER (CHAINL=100,CHAINN=200)

      INTEGER  I, II, J, K, IFOCC4TH, IFRING
      INTEGER  CHAINLENGTH,CHAINNUMBER
      INTEGER  LATTICESIZEX,LATTICESIZEY,LATTICESIZEZ
      INTEGER  XPP,YPP,ZPP
      INTEGER  XP(CHAINN,CHAINL)
      INTEGER  YP(CHAINN,CHAINL)
      INTEGER  ZP(CHAINN,CHAINL)
      INTEGER  AX(CHAINN,CHAINL)
      INTEGER  AY(CHAINN,CHAINL)
      INTEGER  AZ(CHAINN,CHAINL)
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      INTEGER  bondvectorthirdshell(92,3)
      INTEGER  bondvectorfourthshell(162,3)
      INTEGER  BEADOCCUPANCY(0:50,0:50,0:150)
      INTEGER  FIRSTOCCUP(CHAINN,CHAINL)
      INTEGER  SECONDOCCUP(CHAINN,CHAINL)
      INTEGER  THIRDOCCUP(CHAINN,CHAINL)
      INTEGER  FOURTHOCCUP(CHAINN,CHAINL)

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER,
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ


      DO II = 1, CHAINNUMBER
        DO I = 1, CHAINLENGTH

	AX(II,I) = XP(II,I)
	AY(II,I) = YP(II,I)
	AZ(II,I) = ZP(II,I)

  103   IF( AX(II,I) .GE. LATTICESIZEX ) THEN
          AX(II,I) = AX(II,I) - LATTICESIZEX
          GO TO 103
        END IF
  104   IF( AX(II,I) .LT. 0. ) THEN
          AX(II,I) = AX(II,I) + LATTICESIZEX
          GO TO 104
        END IF

  105   IF( AY(II,I) .GE. LATTICESIZEY ) THEN
          AY(II,I) = AY(II,I) - LATTICESIZEY
          GO TO 105
        END IF
  106   IF( AY(II,I) .LT. 0. ) THEN
          AY(II,I) = AY(II,I) + LATTICESIZEY
          GO TO 106
        END IF

  107   IF( AZ(II,I) .GE. LATTICESIZEZ ) THEN
          AZ(II,I) = AZ(II,I) - LATTICESIZEZ
          GO TO 107
        END IF
  108   IF( AZ(II,I) .LT. 0. ) THEN
          AZ(II,I) = AZ(II,I) + LATTICESIZEZ
          GO TO 108
        END IF

        ENDDO
      ENDDO

      DO I = 0, LATTICESIZEZ-1
        DO J = 0, LATTICESIZEY-1
          DO K = 0, LATTICESIZEX-1
            BEADOCCUPANCY(K,J,I) = 0
          ENDDO
        ENDDO
      ENDDO

       DO II = 1, CHAINNUMBER
         DO I = 1, CHAINLENGTH
           BEADOCCUPANCY(AX(II,I),AY(II,I),AZ(II,I)) = 1
	   FIRSTOCCUP(II,I) = 0
	   SECONDOCCUP(II,I) = 0
	   THIRDOCCUP(II,I) = 0
	   FOURTHOCCUP(II,I) = 0
         ENDDO
       ENDDO


       DO II = 1, CHAINNUMBER
         DO I = 1, CHAINLENGTH

	   DO K = 1, 12
		XPP = XP(II,I) + BONDVECTOR(K,1)
		YPP = YP(II,I) + BONDVECTOR(K,2)
		ZPP = ZP(II,I) + BONDVECTOR(K,3)

            IF( I .GE. 3 ) THEN
              IF( (XPP .EQ. XP(II,I-2)) .AND.
     >            (YPP .EQ. YP(II,I-2)) .AND.
     >            (ZPP .EQ. ZP(II,I-2)) ) THEN
               FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
              ENDIF
            ENDIF
            IF( I .LE. CHAINLENGTH-2 ) THEN
              IF( (XPP .EQ. XP(II,I+2)) .AND.
     >            (YPP .EQ. YP(II,I+2)) .AND.
     >            (ZPP .EQ. ZP(II,I+2)) ) THEN
                FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
              ENDIF
            ENDIF

	    IF(IFRING.EQ.1) THEN
	    IF((I.EQ.1).AND.
     >            (XPP .EQ. XP(II,CHAINLENGTH-1)) .AND.
     >            (YPP .EQ. YP(II,CHAINLENGTH-1)) .AND.
     >            (ZPP .EQ. ZP(II,CHAINLENGTH-1)) ) THEN
                FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
            ENDIF
 
	    IF((I.EQ.2).AND.
     >            (XPP .EQ. XP(II,CHAINLENGTH)) .AND.
     >            (YPP .EQ. YP(II,CHAINLENGTH)) .AND.
     >            (ZPP .EQ. ZP(II,CHAINLENGTH)) ) THEN
                FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
            ENDIF

	    IF((I.EQ.CHAINLENGTH).AND.
     >      	  (XPP .EQ. XP(II,2)) .AND.
     >            (YPP .EQ. YP(II,2)) .AND.
     >            (ZPP .EQ. ZP(II,2)) ) THEN
                FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
            ENDIF

	    IF((I.EQ.CHAINLENGTH-1).AND.
     >      	  (XPP .EQ. XP(II,1)) .AND.
     >            (YPP .EQ. YP(II,1)) .AND.
     >            (ZPP .EQ. ZP(II,1)) ) THEN
                FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
            ENDIF
	    ENDIF

  203   IF( XPP .GE. LATTICESIZEX ) THEN
          XPP = XPP - LATTICESIZEX
          GO TO 203
        END IF
  204   IF( XPP .LT. 0. ) THEN
          XPP = XPP + LATTICESIZEX
          GO TO 204
        END IF

  205   IF( YPP .GE. LATTICESIZEY ) THEN
          YPP = YPP - LATTICESIZEY
          GO TO 205
        END IF
  206   IF( YPP .LT. 0. ) THEN
          YPP = YPP + LATTICESIZEY
          GO TO 206
        END IF

  207   IF( ZPP .GE. LATTICESIZEZ ) THEN
          ZPP = ZPP - LATTICESIZEZ
          GO TO 207
        END IF
  208   IF( ZPP .LT. 0. ) THEN
          ZPP = ZPP + LATTICESIZEZ
          GO TO 208
        END IF

	IF(BEADOCCUPANCY(XPP,YPP,ZPP).EQ.1)
     >	  FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) + 1

	   ENDDO

	IF(IFRING.EQ.0)THEN
	IF(I.EQ.1.OR.I.EQ.CHAINLENGTH) THEN
	  FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 1
	ELSE
	  FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 2
	ENDIF
	ELSE
	  FIRSTOCCUP(II,I) = FIRSTOCCUP(II,I) - 2
	ENDIF

	   DO K = 1, 42
		XPP = XP(II,I) + BONDVECTORSECONDSHELL(K,1)
		YPP = YP(II,I) + BONDVECTORSECONDSHELL(K,2)
		ZPP = ZP(II,I) + BONDVECTORSECONDSHELL(K,3)

            IF( I .GE. 3 ) THEN
              IF( (XPP .EQ. XP(II,I-2)) .AND.
     >            (YPP .EQ. YP(II,I-2)) .AND.
     >            (ZPP .EQ. ZP(II,I-2)) ) THEN
                SECONDOCCUP(II,I) = SECONDOCCUP(II,I) - 1
              ENDIF
            ENDIF
            IF( I .LE. CHAINLENGTH-2 ) THEN
              IF( (XPP .EQ. XP(II,I+2)) .AND.
     >            (YPP .EQ. YP(II,I+2)) .AND.
     >            (ZPP .EQ. ZP(II,I+2)) ) THEN
                SECONDOCCUP(II,I) = SECONDOCCUP(II,I) - 1
              ENDIF
            ENDIF

	    IF(IFRING.EQ.1) THEN
	    IF((I.EQ.1).AND.
     >            (XPP .EQ. XP(II,CHAINLENGTH-1)) .AND.
     >            (YPP .EQ. YP(II,CHAINLENGTH-1)) .AND.
     >            (ZPP .EQ. ZP(II,CHAINLENGTH-1)) ) THEN
                SECONDOCCUP(II,I) = SECONDOCCUP(II,I) - 1
            ENDIF
 
	    IF((I.EQ.2).AND.
     >            (XPP .EQ. XP(II,CHAINLENGTH)) .AND.
     >            (YPP .EQ. YP(II,CHAINLENGTH)) .AND.
     >            (ZPP .EQ. ZP(II,CHAINLENGTH)) ) THEN
                SECONDOCCUP(II,I) = SECONDOCCUP(II,I) - 1
            ENDIF

	    IF((I.EQ.CHAINLENGTH).AND.
     >      	  (XPP .EQ. XP(II,2)) .AND.
     >            (YPP .EQ. YP(II,2)) .AND.
     >            (ZPP .EQ. ZP(II,2)) ) THEN
                SECONDOCCUP(II,I) = SECONDOCCUP(II,I) - 1
            ENDIF

	    IF((I.EQ.CHAINLENGTH-1).AND.
     >      	  (XPP .EQ. XP(II,1)) .AND.
     >            (YPP .EQ. YP(II,1)) .AND.
     >            (ZPP .EQ. ZP(II,1)) ) THEN
                SECONDOCCUP(II,I) = SECONDOCCUP(II,I) - 1
            ENDIF
	    ENDIF


  303   IF( XPP .GE. LATTICESIZEX ) THEN
          XPP = XPP - LATTICESIZEX
          GO TO 303
        END IF
  304   IF( XPP .LT. 0. ) THEN
          XPP = XPP + LATTICESIZEX
          GO TO 304
        END IF

  305   IF( YPP .GE. LATTICESIZEY ) THEN
          YPP = YPP - LATTICESIZEY
          GO TO 305
        END IF
  306   IF( YPP .LT. 0. ) THEN
          YPP = YPP + LATTICESIZEY
          GO TO 306
        END IF

  307   IF( ZPP .GE. LATTICESIZEZ ) THEN
          ZPP = ZPP - LATTICESIZEZ
          GO TO 307
        END IF
  308   IF( ZPP .LT. 0. ) THEN
          ZPP = ZPP + LATTICESIZEZ
          GO TO 308
        END IF

	IF(BEADOCCUPANCY(XPP,YPP,ZPP).EQ.1)
     >	  SECONDOCCUP(II,I) = SECONDOCCUP(II,I) + 1

	   ENDDO


	   DO K = 1, 92
		XPP = AX(II,I) + BONDVECTORTHIRDSHELL(K,1)
		YPP = AY(II,I) + BONDVECTORTHIRDSHELL(K,2)
		ZPP = AZ(II,I) + BONDVECTORTHIRDSHELL(K,3)

  403   IF( XPP .GE. LATTICESIZEX ) THEN
          XPP = XPP - LATTICESIZEX
          GO TO 403
        END IF
  404   IF( XPP .LT. 0. ) THEN
          XPP = XPP + LATTICESIZEX
          GO TO 404
        END IF

  405   IF( YPP .GE. LATTICESIZEY ) THEN
          YPP = YPP - LATTICESIZEY
          GO TO 405
        END IF
  406   IF( YPP .LT. 0. ) THEN
          YPP = YPP + LATTICESIZEY
          GO TO 406
        END IF

  407   IF( ZPP .GE. LATTICESIZEZ ) THEN
          ZPP = ZPP - LATTICESIZEZ
          GO TO 407
        END IF
  408   IF( ZPP .LT. 0. ) THEN
          ZPP = ZPP + LATTICESIZEZ
          GO TO 408
        END IF

	IF(BEADOCCUPANCY(XPP,YPP,ZPP).EQ.1)
     >	  THIRDOCCUP(II,I) = THIRDOCCUP(II,I) + 1

	   ENDDO

	IF(IFOCC4TH.EQ.1) THEN
	   DO K = 1, 162
		XPP = AX(II,I) + BONDVECTORFOURTHSHELL(K,1)
		YPP = AY(II,I) + BONDVECTORFOURTHSHELL(K,2)
		ZPP = AZ(II,I) + BONDVECTORFOURTHSHELL(K,3)

  503   IF( XPP .GE. LATTICESIZEX ) THEN
          XPP = XPP - LATTICESIZEX
          GO TO 503
        END IF
  504   IF( XPP .LT. 0. ) THEN
          XPP = XPP + LATTICESIZEX
          GO TO 504
        END IF

  505   IF( YPP .GE. LATTICESIZEY ) THEN
          YPP = YPP - LATTICESIZEY
          GO TO 505
        END IF
  506   IF( YPP .LT. 0. ) THEN
          YPP = YPP + LATTICESIZEY
          GO TO 506
        END IF

  507   IF( ZPP .GE. LATTICESIZEZ ) THEN
          ZPP = ZPP - LATTICESIZEZ
          GO TO 507
        END IF
  508   IF( ZPP .LT. 0. ) THEN
          ZPP = ZPP + LATTICESIZEZ
          GO TO 508
        END IF

	IF(BEADOCCUPANCY(XPP,YPP,ZPP).EQ.1)
     >	  FOURTHOCCUP(II,I) = FOURTHOCCUP(II,I) + 1

	   ENDDO
	ENDIF

      ENDDO 
      ENDDO

      RETURN
      END
