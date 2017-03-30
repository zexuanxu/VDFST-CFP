! calculate mass budget
! do this subroutine at the end of each time step
SUBROUTINE MASSBUDGET()

    USE GLOBAL, ONLY: NCOL, NLAY, IBOUND, DELR, DELC, DZ, VBLOCK, TPLEN, DENFR
    USE GWF,    ONLY: SF, RCH, DRCH, CCHL, CCHR, QFLV, QFLH, HEAD2D, HEAD2DLTP
    USE CON,    ONLY: GCON, TNODE, NNODE, QSCON, QCON, QEXCON, DENRHOCON, CONLOC, ICON, DENCONT, DENCONTLTP &
                      &, AREACON, TUBELEN
    USE TRANS,  ONLY: DENRHO, DENRHOLTP, CONC2D, CONC2DLTP, PRSITI
    USE BUDGET, ONLY: DIFF, PERC, TOTIN, TOTOUT, RCHBUD, CONSTBUD, CONSTBUDIN, CONSTBUDOUT, STOBUDGWF, DCDT &
                      & ,QEXCONBUD, QSCONBUD, STOBUDTRANS, CONSTOBUD, DIFFCON, PERCCON, TOTINCON, TOTOUTCON

    INTEGER (KIND = 4) :: M, N, CN, CT

    CALL ASSIVALBUD()

    DO M = 1, NLAY
        DO N = 1, NCOL

            IF(N .EQ. 1)    THEN
            IF((IBOUND(N, M) .EQ. -1) .AND. (IBOUND(N+1, M) .GT. 0))    THEN
                IF (QFLH(N, M) .GT. 0)  THEN

                    CONSTBUDIN = CONSTBUDIN + DENRHO(N, M) * QFLH(N, M) * DELR(N, M) * DZ(N, M)
!                    CONSTBUDIN = CONSTBUDIN + DENRHO(N, M) * QFLH(N, M) * DELR(N, M) * DELC(N, M)

                ELSE

                    CONSTBUDOUT = CONSTBUDOUT - DENRHO(N+1, M) * QFLH(N, M) * DELR(N, M) * DZ(N, M)
!                    CONSTBUDOUT = CONSTBUDOUT - DENRHO(N+1, M) * QFLH(N, M) * DELR(N, M) * DELC(N, M)

                ENDIF
            ENDIF
            ENDIF

            IF(N .EQ. NCOL) THEN
            IF((IBOUND(N, M) .EQ. -1) .AND. (IBOUND(N-1, M) .GT. 0))    THEN
                IF (QFLH(N-1, M) .GT. 0)  THEN

                    CONSTBUDOUT = CONSTBUDOUT + DENRHO(N-1, M) * QFLH(N-1, M) * DELR(N, M) * DZ(N, M)
!                    CONSTBUDOUT = CONSTBUDOUT + DENRHO(N-1, M) * QFLH(N-1, M) * DELR(N, M) * DELC(N, M)

                ELSE

                    CONSTBUDIN = CONSTBUDIN - DENRHO(N, M) * QFLH(N-1, M) * DELR(N, M) * DZ(N, M)
!                    CONSTBUDIN = CONSTBUDIN - DENRHO(N, M) * QFLH(N-1, M) * DELR(N, M) * DELC(N, M)

                ENDIF
            ENDIF
            ENDIF

            RCHBUD = RCHBUD + DRCH(N, M) * RCH(N, M)

            IF (IBOUND(N, M) .GT. 0)    THEN

                STOBUDGWF = STOBUDGWF + DENRHO(N, M) * SF(N, M) * VBLOCK(N, M) * (HEAD2D(N, M) - HEAD2DLTP(N, M)) / TPLEN

!                STOBUDTRANS = STOBUDTRANS + (CONC2D(N, M) - CONC2DLTP(N, M)) * PRSITI(N, M) * VBLOCK(N, M) / TPLEN
!                STOBUDTRANS = STOBUDTRANS + (CONC2D(N, M) - CONC2DLTP(N, M)) * PRSITI(N, M) * VBLOCK(N, M)

                DCDT = DCDT + PRSITI(N, M) * (DENRHO(N, M) - DENRHOLTP(N, M)) * VBLOCK(N, M) / TPLEN

            ENDIF

        ENDDO
    ENDDO

    TOTIN = TOTIN + CONSTBUDIN
    TOTOUT = TOTOUT + CONSTBUDOUT

!    IF (STOBUD .GT. 0)  THEN
        TOTOUT = TOTOUT + STOBUDGWF
!    ELSE
!        TOTIN = TOTIN - STOBUD
!    ENDIF

    TOTIN = TOTIN + RCHBUD
    TOTOUT = TOTOUT + DCDT

    ! QSCON > 0, conduit inflow (source)
    QSCON(1) = GCON(1) + QCON(1) - QEXCON(1)
    QSCON(NNODE) = GCON(NNODE) - QCON(NNODE-1) - QEXCON(NNODE)

    ! the density at QSCON(1) needs to be reset
    QSCONBUD = QSCON(1) * DENRHOCON(1) + QSCON(NNODE) * DENFR

    ! QEXCONBUD > 0: groundwater flow from matrix system into conduit system
    DO CN = 1, NNODE
        IF (QEXCON(CN) .GT. 0)  THEN
            QEXCONBUD = QEXCONBUD + QEXCON(CN) * DENRHO(CONLOC(2, CN), CONLOC(3, CN))
        ELSE
            QEXCONBUD = QEXCONBUD + QEXCON(CN) * DENRHOCON(CN)
        ENDIF
    ENDDO

    DO CT = 1, TNODE
        CONSTOBUD = CONSTOBUD + (DENCONT(CT) - DENCONTLTP(CT)) * AREACON(CT) * TUBELEN(CT) / TPLEN
    ENDDO

    IF (QEXCONBUD .GT. 0)  THEN
        TOTOUT = TOTOUT + QEXCONBUD
    ELSE
        TOTIN = TOTIN - QEXCONBUD
    ENDIF

    TOTINCON = QSCONBUD
    TOTOUTCON = - QEXCONBUD + CONSTOBUD

    DIFFCON = TOTINCON - TOTOUTCON
    PERCCON = 2*100*DIFFCON/(TOTINCON+TOTOUTCON)

    DIFF = TOTIN - TOTOUT
    PERC = 2*100*DIFF/(TOTIN+TOTOUT)

!   IF (ABS(PERC) .GT. 0.01) THEN
!        STOP
!    ENDIF

! test
!    WRITE(*, *) "CONSTBUDIN"
!    WRITE(*, "(F15.2)") CONSTBUDIN
!    WRITE(*, *) "CONSTBUDOUT"
!    WRITE(*, "(F15.2)") CONSTBUDOUT
!    WRITE(*, *) "STOBUDGWF"
!    WRITE(*, "(F15.2)") STOBUDGWF
!    WRITE(*, *) "STOBUDTRANS"
!    WRITE(*, "(F15.2)") STOBUDTRANS
!    WRITE(*, *) "RCHBUD"
!    WRITE(*, "(F15.2)") RCHBUD
!    WRITE(*, *) "DCDT"
!    WRITE(*, "(F15.2)") DCDT
!    WRITE(*, *) "QEXCONBUD"
!    WRITE(*, "(F15.2)") QEXCONBUD
!    WRITE(*, *) "QSCONBUD"
!    WRITE(*, "(F15.2)") QSCONBUD
!    WRITE(*, *) "TOTIN"
!    WRITE(*, "(F15.2)") TOTIN
!    WRITE(*, *) "TOTOUT"
!    WRITE(*, "(F15.2)") TOTOUT
!    WRITE(*, *) "DIFF"
!    WRITE(*, "(F15.2)") DIFF
!    WRITE(*, "(F10.5)") PERC
! end test

    RETURN

END SUBROUTINE MASSBUDGET

