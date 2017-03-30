

! calculate the conduit node concentration by arithmetic equation
SUBROUTINE CALC_CONC_NODE()

    USE GLOBAL, ONLY: DENFR, DRHODC
    USE CON,    ONLY: NNODE, CONCCONN, CONCCONNLTP, CONCCONT, DENRHOCON, QCON, QEXCON, QSCON, QSCONCONC, CONLOC, RESHEADCON
    USE TRANS,  ONLY: CONC2D
    USE INIT,   ONLY: SCONC 
   
    INTEGER (KIND = 4)  :: CN
    REAL (KIND = 8)     :: CONT, CONQS, VCON

!    WRITE(*, "(10F10.3)") QCON
!    WRITE(*, "(10F10.3)") CONCCONT 

    CONCCONN = 0

    DO CN = 1, NNODE

        IF (CN .EQ. 1)   THEN

            CONCCONN(CN) = 35.0
!            CONCCONN(CN) = SCONC(CONLOC(2, CN), CONLOC(3, CN))
        ELSEIF(CN .EQ. NNODE)   THEN 
    
            CONCCONN(CN) = 0.0

        ELSE

            CONT = 0
            VCON = 0
            
            IF(QCON(CN-1) .GT. 0)   THEN
                CONT = CONT + QCON(CN-1) * CONCCONT(CN-1)
                VCON = VCON + QCON(CN-1)
            ENDIF

            IF (QCON(CN) .LT. 0)    THEN
                CONT = CONT - QCON(CN) * CONCCONT(CN+1)
                VCON = VCON - QCON(CN)
            ENDIF

            ! QEXCON > 0, flow from matrix to conduit
            IF(QEXCON(CN) .GT. 0)   THEN
                CONT = CONT + QEXCON(CN) * CONC2D(CONLOC(2, CN), CONLOC(3, CN))
                VCON = VCON + QEXCON(CN)
            ENDIF

            ! QSCON > 0, source term to conduit
            IF (QSCON(CN) .LT. 0)    THEN
                CONT = CONT + QSCON(CN) * QSCONCONC(CN)
                VCON = VCON + QSCON(CN)
            ENDIF

! test
!            IF(CN .EQ. 2)   THEN            
!                WRITE(*, "(10F10.3)")  QCON(CN-1), QCON(CN), QEXCON(CN)
!                WRITE(*, "(10F10.3)")  CONCCONT(CN-1), CONCCONT(CN)
!            ENDIF
! end test

            IF (VCON .GT. RESHEADCON)   THEN

!               IF(CN .EQ. 1)   THEN
!                CONT = - QCON(CN) * CONCCONT(CN)
!               ELSEIF(CN .EQ. NNODE)   THEN
!                CONT = QCON(CN-1) * CONCCONT(CN-1) 
!               ELSE
!               ENDIF

                CONCCONN(CN) = CONT / VCON 

            ELSE
                CONCCONN(CN) = CONCCONNLTP(CN)
            ENDIF

            DENRHOCON(CN) = DENFR + DRHODC * CONCCONN(CN)

! test
!            IF(CN .EQ. 49)  THEN
!                WRITE(*, "(F10.1)") QCON(CN-1), QCON(CN), CONCCONT(CN-1), CONCCONT(CN)
!                WRITE(*, "(F10.1)") QEXCON(CN), CONC2D(CONLOC(2, CN), CONLOC(3, CN))
!            ENDIF
! end test

        ENDIF

    ENDDO

    RETURN

END SUBROUTINE CALC_CONC_NODE


SUBROUTINE ADV_CON_CONC_IMP_TB()

    USE GLOBAL, ONLY: TPLEN, DENFR, DRHODC
    USE CON,    ONLY: TNODE, TUBELEN, QCONF, CONCCONN, CONCCONT, CONCCONNLTP, CONCCONTLTP, CONCCONNLIMP, &
                      & CONCCONTLIMP, DENCONT

    INTEGER (KIND = 4) :: TUBSEG, TNODESEG
    INTEGER (KIND = 4) :: CT, SG, N, NRHS, LDAS, LDBS, INFO
    INTEGER (KIND = 4), ALLOCATABLE, DIMENSION(:) :: IPIV
    REAL (KIND = 8)     :: CR, SUMCON
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: ADVCONMAT, A
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: B
    
    TUBSEG = 10

    TNODESEG = TUBSEG * TNODE

    ALLOCATE(ADVCONMAT(TNODESEG, TNODESEG))
    ALLOCATE(A(TNODESEG, TNODESEG))
    ALLOCATE(B(TNODESEG))
    ALLOCATE(IPIV(TNODESEG))

    N = TNODESEG
    NRHS = 1
    LDAS = TNODESEG
    LDBS = TNODESEG

    ADVCONMAT = 0

!    B = CONCCONTLIMP

    DO CT = 1, TNODESEG
        
        B(CT) = CONCCONTLTP(INT((CT+9)/10.0))

    ENDDO

!    WRITE(*, *) "QCONF"
!    WRITE(*, "(10F15.5)") QCONF

    DO CT = 1, TNODESEG

        SG = INT((CT+9)/TUBSEG)
        CR = 0.5 * QCONF(SG) * TPLEN / (TUBELEN(SG)/TUBSEG) 

        IF (CT .EQ. 1)  THEN

            ADVCONMAT(CT, CT) = 1
            ADVCONMAT(CT, CT+1) = CR
            B(CT) = B(CT) + CR * CONCCONNLTP(1)

        ELSEIF (CT .EQ. TNODESEG)  THEN

            ADVCONMAT(CT, CT-1) = - CR
            ADVCONMAT(CT, CT) = 1
            B(CT) = B(CT) - CR * CONCCONNLTP(TNODE+1) 

        ELSE

            ADVCONMAT(CT, CT-1) = - CR
            ADVCONMAT(CT, CT) = 1
            ADVCONMAT(CT, CT+1) = CR
    
        ENDIF

    ENDDO
    
    A = ADVCONMAT

    CALL DGESV(N, NRHS, A, LDAS, IPIV, B, LDBS, INFO)

    IF(INFO .NE. 0) THEN
        STOP
    ENDIF

    DO CT = 1, TNODE
        SUMCON = 0        
        DO SG = 1, TUBSEG
            SUMCON = SUMCON + B((CT-1)*TUBSEG+SG)
        ENDDO
        CONCCONT(CT) = SUMCON/TUBSEG
    ENDDO

    DO CT = 1, TNODE
        DENCONT(CT) = DENFR + DRHODC * CONCCONT(CT)
    ENDDO

    DEALLOCATE(ADVCONMAT)
    DEALLOCATE(A)
    DEALLOCATE(B)
    DEALLOCATE(IPIV)

    RETURN

END SUBROUTINE ADV_CON_CONC_IMP_TB


