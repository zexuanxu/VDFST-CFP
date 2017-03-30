! VDFST-CFP: Variable Density Flow and Solute Transport - Conduit Flow Process
! Author: Zexuan Xu, Department of Earth, Ocean and Atmosphere Science, Florida State University
! Email: xuzexuan@gmail.com
!
! All Right Resersed


PROGRAM VDFST_CFP

    USE GLOBAL, ONLY: NPER, NSTP, TPLEN, IMPCON, IMPAG, IMPNO, KONV
! test used only
    USE GWF,    ONLY: QFLH, HEAD2D
    USE CON,    ONLY: HEADCON, HEADCONLTP, QEXCON, CONCCONN, CONCCONNLTP, CONCCONT, CONCCONTLTP, DENRHOCON, DENRHOCONLTP
    USE TRANS,  ONLY: CONC, CONC2D, CONCLTP, CONC2DLTP

    INTEGER (KIND = 4) :: KPER, KSTP, OCPI
    REAL (KIND = 8)    :: total_time, start_time, end_time

    CALL CPU_TIME(start_time)

    ! assign all parameters
    CALL ASSIVAL()

    ! assign initial conditions for head, concentration and density
    CALL ASSIVALINIT()
    
    OPEN(30, FILE='vdfst_cfp.lst', STATUS='REPLACE')
    OPEN(31, FILE='mediahead.txt', STATUS='REPLACE')
    OPEN(32, FILE='conduithead.txt', STATUS='REPLACE')
    OPEN(33, FILE='mediaconc.txt', STATUS='REPLACE')
    OPEN(34, FILE='conduitconc.txt', STATUS='REPLACE')

!    WRITE(*,*) "TEST"

    KPER = 0
    KSTP = 0

! OCPI: determine initial condition, print or not print
! OCPI = 1: print head/conc in both list file and head/conc file
! OCPI = -1: print head/conc in list file only

    OCPI = -1

    CALL PRINTOUT(KSTP, KPER, OCPI)

    OCPI = 1

    KONV = 0

    DO KPER = 1, NPER
        DO KSTP = 1, NSTP
   
            WRITE(30, "('')")
            WRITE(30, "('Time step 'I6', Stress period 'I3', Elapsed time: 'F10.6'')") KSTP, KPER, TPLEN*KSTP
 
            WRITE(*, "('Time step 'I6', Stress period 'I3', Elapsed time: 'F10.6'')") KSTP, KPER, TPLEN*KSTP

            IMPAG = 0
            IMPCON = 0

            IMPNO = 0

            DO WHILE((IMPAG .EQ. 0) .OR. (IMPCON .EQ. 0))

                IMPNO = IMPNO + 1

                WRITE(30, "('Implicit loop 'I6'')")  IMPNO

                IF (IMPNO .EQ. 1)   THEN    
  
                    CALL INITIMESTEP()

                ENDIF
    
                CALL INITIMPLOOP(IMPNO)
    
! test (only for no GWF_CON case)
!                IMPCON = 1
!                QEXCON = 0
! end test

                CALL GWF_CON(KSTP, KPER)
                CALL GWF_MEDIA(KSTP, KPER)
!                CALL GWF_CON(KSTP, KPER)

!                WRITE(*, "(10F15.5)") QEXCON

                CALL TRANS_CON(KSTP, KPER)
                CALL TRANS_MEDIA(KSTP, KPER)        
!                CALL TRANS_CON(KSTP, KPER)

                CALL DETERIMPDEN()

                IF (IMPNO .GT. 200) THEN
                    STOP
                ENDIF
            ENDDO


!            WRITE(*, "(10F15.5)") QEXCON

!            WRITE(*, "(10F10.5)")  HEADCON
!            WRITE(*, "(10F10.5)")  HEAD2D(:,11)
!            WRITE(*, "(10F10.5)")  HEADCON - HEAD2D(:,11)

!            CALL GWF_CON(KSTP, KPER)
!            CALL TRANS_CON(KSTP, KPER)

            CALL MASSBUDGET()

            IF(NSTP .LE. 50)    THEN
                CALL PRINTOUT(KSTP, KPER, OCPI)
            ELSE
                IF(MOD(KSTP, NSTP/200) .EQ. 0)    THEN
                    CALL PRINTOUT(KSTP, KPER, OCPI)
                ENDIF
            ENDIF
       
 
            WRITE(30, "('End of Time step 'I6', Stress period 'I3'')") KSTP, KPER
            
        ENDDO
    ENDDO

    ! Deallocate memory
    CALL FREEMEMORY()

    CLOSE(31)
    CLOSE(32)
    CLOSE(33)
    CLOSE(34)

    CALL CPU_TIME(end_time)

    total_time = end_time - start_time

    WRITE(30, "('Total CPU time: 'I5' minutes, 'F7.3' second')") INT(total_time/60.0), MOD(total_time, 60.0)
    CLOSE(30)

    WRITE(*, "('VDFST_CFP PROGRAM DONE!')")

!       STOP

END PROGRAM VDFST_CFP


! formulate FD matrix and rhs for groundwater flow in porous media
SUBROUTINE GWF_MEDIA(KSTP, KPER)

    USE GLOBAL, ONLY: TPLEN    
    USE GWF,    ONLY: HEAD, HEAD2D, HEADLTP, HEAD2DLTP
    USE CON,    ONLY: QEXCON
    USE INIT,   ONLY: STRT
    USE BUDGET, ONLY: STOBUDGWF

    INTEGER (KIND = 4) :: IMPCON
    INTEGER (KIND = 4)  :: KSTP, KPER
    REAL (KIND = 8) :: HEADSUM

    ! calculate FD matrix for groundwater flow in porous medium
    ! CALL GFW_MEDIA_MATRIX(GMAT, DENRHO)
    CALL GWF_MEDIA_MATRIX()       

    ! calculate RHS of head for groundwater flow in porous medium
    ! CALL GWF_MEDIA_RHS(GRHS, GMAT, DENH, DENV, CONC, CONLTP)
    CALL GWF_MEDIA_RHS()        

    ! calculate head from the FD matrix and RHS
    ! CALL SOLVE_MEDIA_HEAD(HEAD, HEADLTP, QFLH, QFLV, GMAT, GRHS)
    CALL SOLVE_MEDIA_HEAD()

    CALL VEC2MAT(HEAD2D, HEAD, STRT)    

    CALL FLOWMEDIA()

!    WRITE(*,*) "TEST1"

    RETURN

END SUBROUTINE GWF_MEDIA


! solve conduit flow 
SUBROUTINE GWF_CON(KSTP, KPER)

    INTEGER (KIND = 4)  :: KSTP, KPER, IMPCON

    ! calculate the jacobian for Newton-Rasphon method for groundwater flow in conduits.
    ! calculate head of conduit
    ! all calculation has been inclued in the function below 

!    CALL SOLVE_CON_HEAD()

    CALL SOLVE_CON_HEAD_DGESV()

    RETURN

END SUBROUTINE GWF_CON


! formulate FD matrix and rhs for transport in porous media
SUBROUTINE TRANS_MEDIA(KSTP, KPER)

    USE TRANS,  ONLY: CONC, CONC2D, DENRHO, DENRHOLTP
    USE INIT,   ONLY: SCONC

    INTEGER (KIND = 4)  :: KSTP, KPER

    ! calculate FD matrix for solute transport in porous medium
    CALL TRANS_MEDIA_MATRIX()        

    ! calculate RHS of solute concentration for transport in porous medium
    CALL TRANS_MEDIA_RHS()

    ! calculate solute concentration from the FD matrix and RHS in porous medium
    CALL SOLVE_MEDIA_CONC()

    CALL VEC2MAT(CONC2D, CONC, SCONC)

    CALL CONC2DEN() 

    RETURN

END SUBROUTINE TRANS_MEDIA


! formulate FD matrix and rhs for transport in conduit 
SUBROUTINE TRANS_CON(KSTP, KPER)

    INTEGER (KIND = 4)  :: KSTP, KPER

    CALL ADV_CON_CONC_IMP_TB()

    CALL CALC_CONC_NODE()

    RETURN

END SUBROUTINE TRANS_CON




! print out results in output files (both flow and transport)
SUBROUTINE PRINTOUT(KSTP, KPER, OCPI)

    USE GLOBAL, ONLY: NCOL, NLAY, DENFR
    USE GWF,    ONLY: HEAD, HEAD2D
    USE CON,    ONLY: NNODE, HEADCON, CONCCONN, DENRHOCON, ZCON
    USE TRANS,  ONLY: CONC, CONC2D
    USE INIT,   ONLY: SCONC, STRT
    USE BUDGET, ONLY: CONSTBUDIN, CONSTBUDOUT, RCHBUD, STOBUDGWF, DCDT, DIFF, PERC, TOTIN, TOTOUT, & 
                      QEXCONBUD, QSCONBUD, STOBUDTRANS, CONSTOBUD, DIFFCON, PERCCON, TOTINCON, TOTOUTCON
    
    INTEGER (KIND = 4) :: K, CN, KSTP, KPER, OCPI
    REAL (KIND = 8), DIMENSION(NCOL, NLAY) :: MEA
    REAL (KIND = 8), DIMENSION(NNODE)   :: MEACON

    CALL HF2MEA(HEAD2D, MEA)

    WRITE(30, "('')")
    WRITE(30, "('Measured head in porous media for Time step 'I6', Stress period 'I3'')")  KSTP, KPER

    DO K = 1, NLAY

        WRITE(30, "('Measured head in Layer 'I5'')")  K
        WRITE(30, '(10F10.3)') MEA(:, K)
!        WRITE(30, '(10F10.3)') HEAD2D(:, K)   
        
        IF (OCPI .GT. 0)    THEN    
            WRITE(31, '(10F10.3)') MEA(:, K)
!            WRITE(31, '(10F10.3)') HEAD2D(:, K)  
        ENDIF

! test
!        WRITE(31, '(10F10.3)') HEAD2D(:, K)   
! end test

    ENDDO

    WRITE(30, "('')")
    WRITE(30, "('Measured head in conduit for Time step 'I6', Stress period 'I3'')")  KSTP, KPER

    DO CN = 1, NNODE

        MEACON(CN) = DENFR * HEADCON(CN) / DENRHOCON(CN) + (DENRHOCON(CN) - DENFR) * ZCON(CN) / DENRHOCON(CN)
    
    ENDDO

    WRITE(30, '(10F10.3)') MEACON 

    IF (OCPI .GT. 0)    THEN 
        WRITE(32, '(10F10.3)') MEACON
    ENDIF

    WRITE(30, "('')")
    WRITE(30, "('Concentration in porous media for Time step 'I6', Stress period 'I3'')")  KSTP, KPER

    DO K = 1, NLAY

        WRITE(30, "('Concentration in Layer 'I5'')")  K
        WRITE(30, '(10F10.3)') CONC2D(:, K)

        IF (OCPI .GT. 0)    THEN
            WRITE(33, '(10F10.3)') CONC2D(:, K) 
        ENDIF
    ENDDO

    WRITE(30, "('')")
    WRITE(30, "('Concentration in conduit for Time step 'I6', Stress period 'I3'')")  KSTP, KPER
    WRITE(30, '(10F10.3)') CONCCONN

    IF (OCPI .GT. 0)    THEN
        WRITE(34, '(10F10.3)') CONCCONN 
    ENDIF    
    
    WRITE(30, "('')")
    WRITE(30, "('Conduit mass budget for Time step 'I6', Stress period 'I3'')")  KSTP, KPER
    WRITE(30, "('Rates for this time step M**3/T')")
    WRITE(30, "('IN                            :')")
 
    IF (QSCONBUD .GT. 0)    THEN
        WRITE(30, "('Conduit sink/source           : 'F15.3' ')") QSCONBUD
    ENDIF
    IF (QEXCONBUD .GT. 0)    THEN
        WRITE(30, "('Conduit/Porous media exchange : 'F15.3' ')") QEXCONBUD
    ENDIF
    WRITE(30, "('Total IN                      : 'F15.3' ')")  TOTINCON

    WRITE(30, "('OUT                           :')")
    
    IF(CONSTOBUD .GT. 0)    THEN
        WRITE(30, "('Conduit storage (density)     : 'F15.3' ')") CONSTOBUD
    ENDIF
    IF (QSCONBUD .LT. 0)    THEN
        WRITE(30, "('Conduit sink/source           : 'F15.3' ')") -QSCONBUD
    ENDIF
    IF(QEXCONBUD .LT. 0)    THEN
        WRITE(30, "('Conduit/Porous media exchange : 'F15.3' ')") -QEXCONBUD
    ENDIF
    WRITE(30, "('Total OUT                     : 'F15.3' ')")  TOTOUTCON

    WRITE(30, "('')")
    WRITE(30, "('IN - OUT                      : 'F15.3' ')")  DIFFCON
    WRITE(30, "('PRECENT                       : 'F15.6' ')")  PERCCON

    WRITE(30, "('')")
    WRITE(30, "('Porous medium mass budget for Time step 'I6', Stress period 'I3'')")  KSTP, KPER
    WRITE(30, "('Rates for this time step M**3/T')")
    WRITE(30, "('IN                            :')")
    WRITE(30, "('Constant                      : 'F15.3' ')")  CONSTBUDIN
    WRITE(30, "('Recharge                      : 'F15.3' ')")  RCHBUD

    IF (QEXCONBUD .LE. 0)  THEN
        WRITE(30, "('Conduit/Porous media exchange : 'F15.3' ')") -QEXCONBUD
    ENDIF

    WRITE(30, "('Total IN                      : 'F15.3' ')")  TOTIN

    WRITE(30, "('')")
    WRITE(30, "('OUT                           :')")
    WRITE(30, "('Constant                      : 'F15.3' ')")  CONSTBUDOUT
    WRITE(30, "('Storage (Flow)                : 'F15.3' ')")  STOBUDGWF
    WRITE(30, "('Storage (Transport)           : 'F15.3' ')")  STOBUDTRANS
    WRITE(30, "('DCDT                          : 'F15.3' ')")  DCDT
    
    IF (QEXCONBUD .GT. 0)  THEN
        WRITE(30, "('Conduit/Porous media exchange : 'F15.3' ')") QEXCONBUD
    ENDIF

    WRITE(30, "('Total OUT                     : 'F15.3' ')")  TOTOUT

    WRITE(30, "('')")
    WRITE(30, "('IN - OUT                      : 'F15.3' ')")  DIFF
    WRITE(30, "('PRECENT                       : 'F15.6' ')")  PERC

    RETURN

END SUBROUTINE PRINTOUT



SUBROUTINE DETERIMPDEN()

    USE GLOBAL, ONLY: NCOL, NLAY, IMPAG
    USE GWF,    ONLY: HEAD2D, HEAD2DLIMP
    USE CON,    ONLY: NNODE, HEADCON, HEADCONLIMP, DENRHOCON, DENRHOCONLIMP, CONCCONN
    USE TRANS,  ONLY: DENRHO, DENRHOLIMP

    REAL (KIND = 8) :: HEADMAX, HEADCONMAX, HEADCRE, HEADCONCRE
    REAL (KIND = 8) :: DENMAX, DENCONMAX, DENCRE, DENCONCRE
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)    :: DENDIFF, HEAD2DIFF
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)      :: DENCONDIFF, HEADCONDIFF

    ALLOCATE(HEAD2DIFF(NCOL, NLAY))
    ALLOCATE(HEADCONDIFF(NNODE))
    ALLOCATE(DENDIFF(NCOL, NLAY))
    ALLOCATE(DENCONDIFF(NNODE))

    HEAD2DIFF = HEAD2D - HEAD2DLIMP
    HEADCONDIFF = HEADCON - HEADCONLIMP
    DENDIFF = DENRHO - DENRHOLIMP
    DENCONDIFF = DENRHOCON - DENRHOCONLIMP

    HEADMAX = MAXVAL(ABS(HEAD2DIFF))
    HEADCONMAX = MAXVAL(ABS(HEADCONDIFF))
    DENMAX = MAXVAL(ABS(DENDIFF))
    DENCONMAX = MAXVAL(ABS(DENCONDIFF))

!    HEADCRE = 0.0005
!    HEADCONCRE = 0.0005
!    DENCRE = 0.0001
!    DENCONCRE = 0.0001

    HEADCRE = 0.0001
    HEADCONCRE = 0.0001
    DENCRE = 0.00002
    DENCONCRE = 0.00002


    IF ((DENMAX .LT. DENCRE) .AND. (DENCONMAX .LT. DENCONCRE) .AND. (HEADMAX .LT. HEADCRE) & 
        & .AND. (HEADCONMAX .LT. HEADCONCRE))  THEN
        
        ! largest density is less than criterion, IMPAG = 1 to stop the do-while loop
        IMPAG = 1
    ELSE
    ENDIF    

    WRITE(30, "('HEADMAX    : 'F10.6' ')") HEADMAX
    WRITE(30, "('HEADCONMAX : 'F10.6' ')") HEADCONMAX
    WRITE(30, "('DENMAX     : 'F10.6' ')") DENMAX
    WRITE(30, "('DENCONMAX  : 'F10.6' ')") DENCONMAX

!    WRITE(*, "('HEADMAX    : 'F10.6' ')") HEADMAX
!    WRITE(*, "('HEADCONMAX : 'F10.6' ')") HEADCONMAX
!    WRITE(*, "('DENMAX     : 'F10.6' ')") DENMAX
!    WRITE(*, "('DENCONMAX  : 'F10.6' ')") DENCONMAX

    DEALLOCATE(HEAD2DIFF)
    DEALLOCATE(HEADCONDIFF)
    DEALLOCATE(DENCONDIFF)
    DEALLOCATE(DENDIFF)
    
    RETURN

END SUBROUTINE DETERIMPDEN


SUBROUTINE INITIMESTEP()
               
    USE GWF,    ONLY: HEAD, HEAD2D, HEADLTP, HEAD2DLTP
    USE CON,    ONLY: CONCCONN, CONCCONNLTP, CONCCONT, CONCCONTLTP, DENRHOCON, DENRHOCONLTP, CONCCONNLIMP, CONCCONTLIMP, &
                      & DENRHOCONLIMP, HEADCON, HEADCONLTP, HEADCONLIMP, DENCONT, DENCONTLTP
    USE TRANS,  ONLY: CONCLTPL, CONCLTP, CONC, CONC2D, CONC2DLTP, CONC2DLIMP, DENRHO, DENRHOLTP, CONCLIMP, CONCLIMPL, DENRHOLIMP

    HEADLTP = HEAD
    HEAD2DLTP = HEAD2D  
    HEADCONLTP = HEADCON                

    CONCLTPL = CONCLTP
    CONCLIMPL = CONCLTP

    CONCLTP = CONC
    CONCLIMP = CONC

    CONCLTP = CONC
    CONC2DLTP = CONC2D
    DENRHOLTP = DENRHO

    CONCCONNLTP = CONCCONN
    CONCCONTLTP = CONCCONT
    DENRHOCONLTP = DENRHOCON

    DENCONTLTP = DENCONT

    RETURN

END SUBROUTINE INITIMESTEP


SUBROUTINE INITIMPLOOP(IMPNO)

    USE GWF,    ONLY: HEAD, HEAD2D, HEADLTP, HEAD2DLTP, HEADLIMP, HEAD2DLIMP
    USE CON,    ONLY: HEADCON, HEADCONLIMP, CONCCONN, CONCCONT, CONCCONNLIMP, CONCCONTLIMP, DENRHOCON, DENRHOCONLIMP 
    USE TRANS,  ONLY: CONC, CONC2D, CONCLTP, CONC2DLTP, CONC2DLIMP, CONCLIMP, CONCLIMPL, DENRHO, DENRHOLIMP

    INTEGER :: IMPNO

    HEADLIMP = HEAD
    HEAD2DLIMP = HEAD2D
    HEADCONLIMP = HEADCON

    IF(IMPNO .GE. 2)    THEN
        CONCLIMPL = CONCLIMP
        CONCLIMP = CONC
    ENDIF

    CONC2DLIMP = CONC2D
    DENRHOLIMP = DENRHO

    CONCCONNLIMP = CONCCONN
    CONCCONTLIMP = CONCCONT
    DENRHOCONLIMP = DENRHOCON

    RETURN

END SUBROUTINE




