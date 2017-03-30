SUBROUTINE VEC2MAT(MAT, VEC, MATINI)

    ! may not use NOCL, NLAY in the subroutine

    USE GLOBAL, ONLY: NCOL, NLAY, TVH

    INTEGER (KIND = 4) :: K, I
    REAL (KIND = 8), DIMENSION(NCOL, NLAY), INTENT(OUT)    :: MAT
    REAL (KIND = 8), DIMENSION(TVH), INTENT(IN)            :: VEC 
    REAL (KIND = 8), DIMENSION(NCOL, NLAY), INTENT(IN)     :: MATINI
    
    DO K = 1, NLAY

        DO I = 1, NCOL 
        
            MAT(I, K) = VEC(I + (K - 1) * NCOL)

        ENDDO
    ENDDO

    RETURN

END SUBROUTINE VEC2MAT


SUBROUTINE MAT2VEC(MAT, VEC)

    USE GLOBAL, ONLY: NCOL, NLAY, TVH, IBOUND

    INTEGER (KIND = 4) :: K, I

    REAL (KIND = 8), DIMENSION(NCOL, NLAY), INTENT(IN)    :: MAT
    REAL (KIND = 8), DIMENSION(TVH), INTENT(OUT)      :: VEC

    DO K = 1, NLAY

        DO I = 1, NCOL 

                VEC(I + (K - 1) * NCOL) = MAT(I, K)
            
        ENDDO
    ENDDO
    
    RETURN

END SUBROUTINE MAT2VEC


SUBROUTINE CONC2DEN()

    USE GLOBAL, ONLY: NCOL, NLAY, DENFR, DRHODC, IBOUND
    USE TRANS,  ONLY: CONC2D, DENRHO
    USE INIT,   ONLY: SDRHO

    INTEGER (KIND = 4) :: K, I, CDC

    CDC = 0

    DO K = 1, NLAY

        DO I = 1, NCOL

            IF(IBOUND(I, K) .LE. 0) THEN

                CDC = CDC + 1
                DENRHO(I, K) = SDRHO(I, K)   
    
            ELSE

                DENRHO(I, K) = DENFR + DRHODC * CONC2D(I, K)
            ENDIF
             
        ENDDO

    ENDDO

END SUBROUTINE CONC2DEN

! not use anymore
SUBROUTINE VEC2MATDEN(MAT, VEC, MATINI)

    ! may not use NOCL, NLAY in the subroutine

    USE GLOBAL, ONLY: NCOL, NLAY, TVH, DENFR, DRHODC, IBOUND

    INTEGER (KIND = 4) :: K, I, CDC
    REAL (KIND = 8), DIMENSION(NCOL, NLAY), INTENT(OUT)    :: MAT
    REAL (KIND = 8), DIMENSION(TVH), INTENT(IN)      :: VEC 
    REAL (KIND = 8), DIMENSION(NCOL, NLAY), INTENT(IN)     :: MATINI

    CDC = 0

    DO K = 1, NLAY

        DO I = 1, NCOL 

            IF(IBOUND(I, K) .LE. 0) THEN
     
                CDC = CDC + 1
                MAT(I, K) = MATINI(I, K)

            ELSE

                MAT(I, K) = DENFR + DRHODC * (VEC(I + (K-1)*NCOL - CDC)) 
            ENDIF
        ENDDO
    ENDDO

    RETURN

END SUBROUTINE VEC2MATDEN


SUBROUTINE PRINTOUT2D(MAT, N, M)

   !  USE GLOBAL, ONLY: NLAY, NCOL

    INTEGER (KIND = 4) :: K, N, M
    REAL (KIND = 8), DIMENSION(N, M)    :: MAT

    DO K = 1, M

        WRITE(*, '(10F15.3, 1X)') MAT(:,K)

    ENDDO

END SUBROUTINE PRINTOUT2D


SUBROUTINE PRINTOUT1D(VEC, N)

    ! USE CON,    ONLY: NNODE

    INTEGER (KIND = 4) :: N
    REAL (KIND = 8), DIMENSION(N)  :: VEC

        WRITE(*, "(10F10.3, 1X)") VEC


END SUBROUTINE PRINTOUT1D


SUBROUTINE MEA2HF(MAT, HF)

    USE GLOBAL, ONLY: NCOL, NLAY, DENFR, ELAY
    USE TRANS,  ONLY: DENRHO
    REAL (KIND = 8), DIMENSION(NCOL, NLAY)    :: MAT, HF    
    
    INTEGER (KIND = 4) :: K, I

    DO K = 1, NLAY
        DO I = 1, NCOL
    
            HF(I, K) = DENRHO(I, K) * MAT(I, K) / DENFR -  (DENRHO(I, K) - DENFR) * ELAY(I, K) / DENFR
   
        ENDDO
    ENDDO


END SUBROUTINE MEA2HF


SUBROUTINE HF2MEA(HF, MAT)

    USE GLOBAL, ONLY: NCOL, NLAY, DENFR, ELAY
    USE TRANS,  ONLY: DENRHO

    INTEGER (KIND = 4):: K, I
    REAL (KIND = 8), DIMENSION(NCOL, NLAY) :: MAT, HF

    DO K = 1, NLAY
        DO I = 1, NCOL

            MAT(I, K) =  DENFR * HF(I, K) / DENRHO(I, K) + (DENRHO(I, K) - DENFR) * ELAY(I, K) / DENRHO(I, K)

        ENDDO
    ENDDO

END SUBROUTINE HF2MEA

SUBROUTINE FREEMEMORY()

    USE GLOBAL
    USE GWF
    USE CON
    USE TRANS
    USE INIT

! deallocate GLOBAL module

    DEALLOCATE(IBOUND)
    DEALLOCATE(VBLOCK)
    DEALLOCATE(HTOP)
    DEALLOCATE(DELC)
    DEALLOCATE(DELR)
    DEALLOCATE(DZ)
    DEALLOCATE(ELAY)

! deallocate GWF module
    DEALLOCATE(SF)
    DEALLOCATE(HY)
    DEALLOCATE(HV)
    DEALLOCATE(CC)
    DEALLOCATE(CV)
    DEALLOCATE(CCHL)
    DEALLOCATE(CCHR)
    DEALLOCATE(CRVT)
    DEALLOCATE(CRVB)
    DEALLOCATE(CVVT)
    DEALLOCATE(CVVB)
    DEALLOCATE(RCH)
    DEALLOCATE(DRCH)
    DEALLOCATE(HEAD)
    DEALLOCATE(HEADLTP)
    DEALLOCATE(QFLH)
    DEALLOCATE(QFLV)
    DEALLOCATE(GMAT)
    DEALLOCATE(DENH)
    DEALLOCATE(DENV)
    DEALLOCATE(DENCH)
    DEALLOCATE(DENCV)
    DEALLOCATE(GRHS)
    DEALLOCATE(HEAD2D)
    DEALLOCATE(HEAD2DLTP)

! deallocate CON module
    DEALLOCATE(DIAM)
    DEALLOCATE(AREACON)
    DEALLOCATE(TUBELEN)
    DEALLOCATE(KCC)
    DEALLOCATE(KCCOND)
    DEALLOCATE(HAD)
    DEALLOCATE(DSX)
    DEALLOCATE(ZCON)
    DEALLOCATE(CONLOC)
    DEALLOCATE(ICON)
    DEALLOCATE(HEADCON)
    DEALLOCATE(HEADCONLTP)
    DEALLOCATE(CONCCONN)
    DEALLOCATE(CONCCONNLTP)
    DEALLOCATE(CONCCONT)
    DEALLOCATE(CONCCONTLTP)    
    DEALLOCATE(DENRHOCON)
    DEALLOCATE(QCON)
    DEALLOCATE(QCONF)
    DEALLOCATE(QEXCON)
    DEALLOCATE(QSCON)
    DEALLOCATE(QSCONCONC)
    DEALLOCATE(CONJAC)
    DEALLOCATE(GCON)
    DEALLOCATE(TCMAT)
    DEALLOCATE(TCRHS)

! deallocate TRANS module
    DEALLOCATE(PRSITI)
    DEALLOCATE(PRDC)
    DEALLOCATE(PRDZ)
    DEALLOCATE(SQX)
    DEALLOCATE(SQZ)
    DEALLOCATE(SSM)
    DEALLOCATE(CONC)
    DEALLOCATE(CONCLTP)
    DEALLOCATE(CONC2D)
    DEALLOCATE(DENRHO)
    DEALLOCATE(DENRHOLTP)
    DEALLOCATE(TMAT)
    DEALLOCATE(TRHS)

! deallocate INIT module
    DEALLOCATE(STRT)
    DEALLOCATE(SCONC)
    DEALLOCATE(SDRHO)

    RETURN

END SUBROUTINE FREEMEMORY