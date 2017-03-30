! calculate jacobian matrix for conduit flow (Newton-Raphson method)
SUBROUTINE GWF_CON_JACOBIAN()

    USE GLOBAL, ONLY: DENFR
    USE CON,    ONLY: TVHNODE, NNODE, CONJAC, HAD, KCCOND, HEADCON, DENRHOCON, QEXCON, ZCON, QCON, QFLCON, &
                      & TUBELEN, RESHEADCON, DENCONT, FRAT 

    INTEGER (KIND = 4) :: CN, CT 
    REAL (KIND = 8) :: DENTC    

    CONJAC = 0

    DO CN = 1, NNODE 

        IF((CN .EQ. 1) .OR. (CN .EQ. NNODE))    THEN

            CONJAC(CN, CN) = 1
        
        ELSE
    
! could be simplified to the function of QCON
!                    DENRHOI = (DENRHOCON(CN-1) + DENRHOCON(CN))/2
!                    CONJAC(CN, CN-1) = HAD(CN-1) / (2 * SQRT((ABS(HEADCON(CN)-HEADCON(CN-1)) + &
!                                     & (DENRHOI - DENFR) * (ZCON(CN-1) - ZCON(CN)) / DENFR) / TUBELEN(CN)))
!                    CONJAC(CN, CN) = CONJAC(CN, CN) - HAD(CN-1) / (2 * SQRT((ABS(HEADCON(CN)-HEADCON(CN-1)) + &
!                                     & (DENRHOI - DENFR) * (ZCON(CN-1) - ZCON(CN)) / DENFR) / TUBELEN(CN)))
! end simplification

            IF (ABS(QFLCON(CN-1)) .GT. RESHEADCON)   THEN

!                    CONJAC(CN, CN-1) = HAD(CN-1) / (2 * SQRT(ABS(QFLCON(CN-1))))
!                    CONJAC(CN, CN) = CONJAC(CN, CN) - HAD(CN-1) / (2 * SQRT(ABS(QFLCON(CN-1))))

                    DENTC = DENFR/DENCONT(CN-1)
!                    DENTC = 1

                    CONJAC(CN, CN-1) = DENTC * HAD(CN-1) / (2 * SQRT(ABS(QFLCON(CN-1)) * TUBELEN(CN-1))) / SQRT(FRAT(CN-1))
                    CONJAC(CN, CN) = CONJAC(CN, CN) - DENTC * HAD(CN-1) / (2 * SQRT(ABS(QFLCON(CN-1)) * TUBELEN(CN-1))) &
                                     & / SQRT(FRAT(CN-1))
 
            ENDIF

            IF (ABS(QFLCON(CN)) .GT. RESHEADCON)   THEN

! simplification
!                DENRHOI = (DENRHOCON(CN) + DENRHOCON(CN+1))/2
!                CONJAC(CN, CN) = CONJAC(CN, CN) - HAD(CN) / (2 * SQRT((ABS(HEADCON(CN)-HEADCON(CN+1)) + &
!                                & (DENRHOI - DENFR) * (ZCON(CN) - ZCON(CN+1)) / DENFR) / TUBELEN(CN)))
!                CONJAC(CN, CN+1) =  HAD(CN) / (2 * SQRT((ABS(HEADCON(CN)-HEADCON(CN+1)) + &
!                                & (DENRHOI - DENFR) * (ZCON(CN) - ZCON(CN+1)) / DENFR) / TUBELEN(CN))) 
! end simplification 

!                CONJAC(CN, CN) = CONJAC(CN, CN) - HAD(CN) / (2 * SQRT(ABS(QFLCON(CN))))
!                CONJAC(CN, CN+1) = HAD(CN) / (2 * SQRT(ABS(QFLCON(CN))))

                DENTC = DENFR/DENCONT(CN)
!                DENTC = 1

                CONJAC(CN, CN) = CONJAC(CN, CN) - DENTC * HAD(CN) / (2 * SQRT(ABS(QFLCON(CN)) * TUBELEN(CN))) / SQRT(FRAT(CN))
                CONJAC(CN, CN+1) = DENTC * HAD(CN) / (2 * SQRT(ABS(QFLCON(CN)) * TUBELEN(CN))) / SQRT(FRAT(CN))

           ENDIF
       ENDIF
                 
        CONJAC(CN, CN) = CONJAC(CN, CN) - KCCOND(CN) 
!        CONJAC(CN, CN) = CONJAC(CN, CN) + KCCOND(CN) 
 
    ENDDO
    
! test
!    WRITE(*, "(I5)") 1 
!    WRITE(*, "(F15.1)")  CONJAC(1,1), CONJAC(1,2)
!    WRITE(*, "(I5)") 54
!    WRITE(*, "(F15.1)")  CONJAC(54,53), CONJAC(54,54), CONJAC(54,55)
!    WRITE(*, "(I5)") 55 
!    WRITE(*, "(F15.1)")  CONJAC(55,54), CONJAC(55,55), CONJAC(55,56)
!    WRITE(*, "(I5)") 56 
!    WRITE(*, "(F15.1)")  CONJAC(56,55), CONJAC(56,56), CONJAC(56,57)
!    WRITE(*, "(I5)") NNODE-1 
!    WRITE(*, "(F15.1)")  CONJAC(NNODE-1,NNODE-2), CONJAC(NNODE-1,NNODE-1), CONJAC(NNODE-1,NNODE)
!    WRITE(*, "(I5)") NNODE 
!    WRITE(*, "(F15.1)")  CONJAC(NNODE, NNODE-1), CONJAC(NNODE, NNODE)
! end test
    
    RETURN

END SUBROUTINE GWF_CON_JACOBIAN


! calculate net flux in each conduit node (Kirchhoff's Law)
SUBROUTINE GWF_CON_GCON
    
    USE GLOBAL, ONLY: DENFR
    USE GWF,    ONLY: HEAD2D
    USE CON,    ONLY: TVHNODE, NNODE, HAD, GCON, ZCON, HEADCON, QCON, HEADCONLTP, &
                      & DENRHOCON, CONLOC, QSCON, RESCON, TUBELEN, QEXCON

    INTEGER (KIND = 4) :: CN
    REAL (KIND = 8)    :: DENRHODI

    GCON = 0
    QSCON = 0

    DO CN = 1, NNODE
        IF ((CN .EQ. 1) .OR. (CN .EQ. NNODE))   THEN
!            GCON(CN) = HEADCON(CN) - HEADCONLTP(CN)
            GCON(CN) = 0 

        ELSE

            ! QEXCON > 0, flow from matrix to conduit
            ! QEXCON < 0, flow from conduit to matrix
            
            GCON(CN) = QCON(CN-1) - QCON(CN) + QEXCON(CN) + QSCON(CN)
    
        ENDIF

    ENDDO

    RETURN

END SUBROUTINE GWF_CON_GCON


! calculate flow and flux within conduit
! QCON: L**3/T
! QCONF: L/T
SUBROUTINE FLOWCON()

    USE GLOBAL, ONLY: DENFR
    USE CON,    ONLY: FRAT, NNODE, TNODE, AREACON, HEADCON, QCON, QCONF, QFLCON, HAD, ZCON, DENRHOCON, TUBELEN, DENCONT
    
    INTEGER (KIND = 4) :: CN, CT
    REAL (KIND = 8)    :: DENRHOI

    QCON = 0
    QCONF = 0
    QFLCON = 0

!    WRITE(*, "(10F10.4)") DENCONT
!    WRITE(*, "(10F10.4)") DENRHOCON    

    DO CT = 1, TNODE

!        DENRHOI = (DENRHOCON(CT) + DENRHOCON(CT+1)) / 2

!        QFLCON(CT) = HEADCON(CT) - HEADCON(CT+1) + (DENRHOI - DENFR) * (ZCON(CT) - ZCON(CT+1)) / DENFR
!        QFLCON(CT) = HEADCON(CT) - HEADCON(CT+1) + (DENCONT(CT) - DENFR) * (ZCON(CT) - ZCON(CT+1)) / DENFR

! test
        QFLCON(CT) = ((HEADCON(CT) - HEADCON(CT+1)) * DENFR + (DENCONT(CT) - DENFR) * (ZCON(CT) - ZCON(CT+1))) &
                     & / DENCONT(CT)
!        QFLCON(CT) = ((HEADCON(CT) - HEADCON(CT+1)) * DENFR + (DENRHOI - DENFR) * (ZCON(CT) - ZCON(CT+1))) &
!                     & / DENRHOI
! test    

        IF (QFLCON(CT) .GT. 0)  THEN

!           QCON(CT) = HAD(CT) * SQRT(ABS(QFLCON(CT)))
            QCON(CT) = HAD(CT) * SQRT(ABS(QFLCON(CT)) / (TUBELEN(CT))) / SQRT(FRAT(CT))
    
        ELSE

        
!           QCON(CT) = - HAD(CT) * SQRT(ABS(QFLCON(CT)))
            QCON(CT) = - HAD(CT) * SQRT(ABS(QFLCON(CT) / TUBELEN(CT))) / SQRT(FRAT(CT))

        ENDIF

        QCONF(CT) = QCON(CT) / AREACON(CT)

    ENDDO

    RETURN

END SUBROUTINE FLOWCON


! calculate flow exchange between conduit and matrix
! QEXCON (L**3/T)
! QEXCON > 0, flow direction is from matrix to conduit
SUBROUTINE QEX_MEDIA_CON()

    USE GWF,    ONLY: HEAD2D
    USE CON,    ONLY: NNODE, HEADCON, QEXCON, CONLOC, KCCOND

    INTEGER (KIND = 4)  :: CN

    QEXCON = 0

    DO CN = 1, NNODE

        QEXCON(CN) = (HEAD2D(CONLOC(2, CN), CONLOC(3, CN)) - HEADCON(CN)) * KCCOND(CN)         
!       QEXCON(CN) = (HEADCON(CN) - HEAD2D(CONLOC(2, CN), CONLOC(3, CN))) * KCCOND(CN)

    ENDDO

    RETURN

END SUBROUTINE QEX_MEDIA_CON


! calculate the exchange conductance between conduit and matrix
! KCCON (L**2/T)
SUBROUTINE CALC_KCCOND()

    USE CON,    ONLY: NNODE, KCC, KCCOND, TUBELEN 

    INTEGER (KIND = 4)  :: CN

    DO CN = 1, NNODE

        IF(CN .EQ. 1)   THEN        

            KCCOND(CN) = 3.1416 *  (KCC(CN) * TUBELEN(CN))/2

        ELSEIF(CN .EQ. NNODE)   THEN

            KCCOND(CN) = 3.1416 * (KCC(CN-1) * TUBELEN(CN-1))/2
        ELSE

            KCCOND(CN) = 3.1416 * (KCC(CN-1) * TUBELEN(CN-1) + KCC(CN) * TUBELEN(CN))/2
        ENDIF

    ENDDO

    RETURN

END SUBROUTINE CALC_KCCOND






! solve conduit head by Newton-Raphson method 
SUBROUTINE SOLVE_CON_HEAD_DGESV()

    USE GLOBAL, ONLY: DENFR, IMPCON, KONV
    USE GWF,    ONLY: HEAD2D, HEAD
    USE CON,    ONLY: NNODE, TVHNODE, HAD, CONJAC, ZCON, HEADCON, HEADCONLTP, GCON, QCON, QEXCON, &
                      & DENRHOCON, CONLOC, QSCON, RESCON, TUBELEN, HEADCONLIMP

    USE INIT,   ONLY: STRT
    USE BUDGET, ONLY: QEXCONBUD, QSCONBUD

    INTEGER (KIND=4) :: CN, TN, IMACON

    INTEGER (KIND=4)     ::  N,  LDA, LDB, LDC, NRHS, INFO

    INTEGER (KIND=4), DIMENSION(NNODE)   :: IPIV
    REAL (KIND=8), DIMENSION(NNODE, NNODE) :: A
    REAL (KIND=8), DIMENSION(NNODE) :: B

    REAL (KIND=8)    :: CRES

    ! HEADCONL: used in the iterative newton-rasphon method
    REAL (KIND=8), ALLOCATABLE, DIMENSION(:)  :: HEADCONL  

    ALLOCATE(HEADCONL(NNODE))

    N = NNODE
    NRHS = 1

    LDA = NNODE
    LDB = NNODE

    HEADCONL = HEADCON

    KONV = 0
    RESCON = 0.00000002

!    WRITE(*, *) "HEADCONT"
!    WRITE(*, "(10F10.3)") HEADCONT

! have been assigned in initialization subroutine 

!    CRES = ABS(NORM2(GCON))

    DO TN = 1, 10000
    
        CALL FLOWCON()

        CALL FRICTION()

        CALL QEX_MEDIA_CON()

        CALL GWF_CON_GCON()

        CALL GWF_CON_JACOBIAN()

        A = CONJAC
        
        B = -GCON

        CALL DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)

        IF(INFO .NE. 0) THEN
             STOP
        ENDIF

        HEADCON = HEADCONL + B
   
        CRES = MAXVAL(ABS(HEADCON - HEADCONL))

        ! iteration converged
        IF ((CRES .LT. RESCON) .AND. (TN .GT. 2))  THEN
            IMPCON = 1
            KONV = 1
            WRITE(*, "(I5)") TN
            EXIT
        ELSE
            IMPCON = 0
            KONV = 0
        ENDIF

        HEADCONL = HEADCON

    ENDDO

    IF (KONV .EQ. 0)    THEN
        WRITE(*, *) 'NEWTON-RAPHSON METHOD NOT CONVERGED'
        WRITE(30, *)    'NEWTON-RAPHSON METHOD NOT CONVERGED'

        WRITE(*, "(F10.3)") CRES

        ! test
!        WRITE(*, *) "GCON"
!        WRITE(*, "(10F10.3)") GCON
!        WRITE(*, *) "HEADCON"
!        WRITE(*, "(10F10.3)") HEADCON
!        WRITE(*, *) "QCON"
!        WRITE(*, "(10F15.5)") QCON
!        WRITE(*, *) "QEXCON"
!        WRITE(*, "(10F15.5)") QEXCON
        ! end test
    
!        STOP
        HEADCON = HEADCONLIMP
!        HEADCON = HEADCONLTP 

     ENDIF

!    WRITE(*, *) "GCON"
!    WRITE(*, "(10F10.3)") GCON
!    WRITE(*, *) "HEADCON"
!    WRITE(*, "(10F10.3)") HEADCON
!    WRITE(*, *) "QCON"
!    WRITE(*, "(10F15.5)") QCON

    DEALLOCATE(HEADCONL)

    RETURN

END SUBROUTINE SOLVE_CON_HEAD_DGESV

!SUBROUTINE FRICTION

!    USE CON,    ONLY: FRAT, TURBI, TNODE, HEADCON, TUBELEN, DIAM, GRAVAC 
!    INTEGER (KIND = 4)  :: CT
    
!    FRAT = 0
!    TURBI = 0

!    DO CT = 1: TNODE

!        FRAT(CT) = SQRT(0.5 * ABS(HEADCON(CT+1)-HEADCON(CT)) * GRAVAC * (DIAM(CN)**5) * (PI**2) / TUBELEN(CT))

!        VEL(CT) = QCON(CT) / (PI * (DIAM(CT)/2)**2)

!        TRUBI(CT) = FRAT(CT) * LOG10(2.51*)

!    ENDDO

!END SUBROUTINE FRICTION

SUBROUTINE FRICTION

    USE CON,    ONLY: FRAT, REYNOLDS, TNODE, QCON, QCONF, DENCONT, DIAM 
!    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) : VELO

!    ALLOCATE(VELO(TNODE))    
    INTEGER (KIND = 4)  :: CT
    REAL (KIND = 8) :: VISCWA

    FRAT = 0
    REYNOLDS = 0

    VISCWA = 1.004E-006 * (3.28084**2) * (60.0) * (60.0) * (24.0)

    DO CT = 1, TNODE
!        VELO(CT) = QCON(CT) / (PI * (DIAM(CT)/2)**2)
        REYNOLDS(CT) = ABS(QCONF(CT)) * DENCONT(CT) * DIAM(CT) / VISCWA

        IF (REYNOLDS(CT) < 2000)    THEN
            FRAT(CT) = 64 / REYNOLDS(CT)
        ELSE
            FRAT(CT) = 0.316 * (REYNOLDS(CT) ** (-0.25)) + 0.0075 
!            FRAT(CT) = 0.0178

        ENDIF
    ENDDO

!    DEALLOCATE(VELO)
    
!    WRITE(*, *) "FRAT"
!    WRITE(*, "(10F15.8)") FRAT

!    STOP

END SUBROUTINE FRICTION


