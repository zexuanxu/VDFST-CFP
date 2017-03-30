! formuate the FD matrix for groundwater transport model
SUBROUTINE TRANS_MEDIA_MATRIX

    USE GLOBAL, ONLY: TPLEN, NLAY, NCOL, IBOUND, DELC, DELR, DZ, VBLOCK
    USE GWF,    ONLY: QFLH, QFLV, COEH, COEV, HEAD2D
    USE TRANS,  ONLY: DXX, DZZ, PRSITI, PRDC, PRDZ, SQX, SQZ, TMAT

    INTEGER (KIND = 4):: I, K, M, N, MNC

    REAL (KIND = 8)   :: TAMTBC, DISPXZ
    REAL (KIND = 8)   :: TMATL, TMATR, TMATU, TMATD
    REAL (KIND = 8)   :: DISPL, DISPR, DISPU, DISPD
    
    ! velocity field (used for dispersion coefficient calculation)
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)   :: VELDH, VELDV 

    ALLOCATE(VELDH(NCOL-1, NLAY))
    ALLOCATE(VELDV(NCOL, NLAY-1))

    VELDH = 0
    VELDV = 0

    TMAT = 0

    COEH = 0
    COEV = 0

!    WRITE(*, "(10F10.2)") QFLV

    ! calculate the upstream weighting scheme coefficient
    DO K = 1, NLAY
        DO I = 1, NCOL

            IF(I .NE. NCOL) THEN
                IF (QFLH(I, K) .GE. 0) THEN

                    COEH(I, K) = 1.0
                
                ELSE

                    COEH(I, K) = 0.0
                
                ENDIF
            ENDIF

            IF (K .NE. NLAY) THEN
                IF (QFLV(I, K) .GE. O)  THEN

                    COEV(I, K) = 1.0
                  
                ELSE

                    COEV(I, K) = 0.0
               
                ENDIF
            ENDIF

        ENDDO
    ENDDO

! test 
    ! center-in-space weighting scheme
!    COEV = 0.5
!    COEH = 0.5

    ! calculate the dispersion coefficient
    DO K = 1, NLAY 
        DO I = 1, NCOL

            IF (I .NE. NCOL)    THEN
        
!                VELDH(I, K) = QFLH(I, K) / (DELR * DZ(I, K) * PRSITI(I, K))
                VELDH(I, K) = QFLH(I, K) / PRSITI(I, K)
                SQX(I, K) =  DXX * ABS(VELDH(I, K)) / (DELC(I, K) * DELC(I, K))      

            ENDIF

            IF (K .NE. NLAY)    THEN

!                VELDV(I, K) =  QFLV(I, K) / (DELC(I, K) * DELR * PRSITI(I, K))
                VELDV(I, K) = QFLV(I, K) /  PRSITI(I, K)

! comment for horizontal 
!                SQZ(I, K) =  DZZ * ABS(VELDV(I, K)) / (DZ(I, K) * DZ(I, K))
                SQZ(I, K) =  DZZ * ABS(VELDV(I, K)) / (DELR(I, K) * DELR(I, K))

            ENDIF

        ENDDO
    ENDDO


    DO M = 1, NLAY
        DO N = 1, NCOL

            MNC = NCOL * (M-1) + N
        
            ! constant concentration boundary    
            IF(IBOUND(N, M) .EQ. -1) THEN
                
                TMAT(MNC, MNC) = 1

            ELSEIF (IBOUND(N, M) .EQ. 0)  THEN

            ELSE 

                TMATBC = 0
                TMATL = 0
                TMATR = 0
                TMATU = 0
                TMATD = 0

                DISPL = 0
                DISPR = 0
                DISPU = 0
                DISPD = 0
                DISPXZ = 0

                ! PART 1
                ! formulate the FD matrix for surrounding four cells
                IF (N .NE. 1)   THEN 
                    IF (IBOUND(N-1, M) .NE. 0)    THEN
          
                        TMAT(MNC, MNC-1) = COEH(N-1, M) * QFLH(N-1, M) / PRDC(N, M) + SQX(N-1, M)
                    
                        TMATL =  (1 - COEH(N-1, M)) * QFLH(N-1, M) / PRDC(N, M)     
                        DISPL = - SQX(N-1, M)    
 
                    ENDIF
                ENDIF

                IF (N .NE. NCOL)    THEN 
                    IF (IBOUND(N+1, M) .NE. 0)    THEN
            
                        TMAT(MNC, MNC+1) = -(1-COEH(N, M)) * QFLH(N, M) / PRDC(N,M) + SQX(N, M)
                    
                        TMATR = - COEH(N, M) * QFLH(N, M) / PRDC(N, M) 
                        DISPR = - SQX(N, M)

                    ENDIF
                ENDIF        
                 
                IF (M .NE. 1)   THEN
                    IF (IBOUND(N, M-1) .NE. 0)    THEN

                        TMAT(MNC, MNC - NCOL) = COEV(N, M-1) * QFLV(N, M-1) / PRDZ(N,M) + SQZ(N, M-1)
                    
                        TMATU = (1 - COEV(N, M-1)) * QFLV(N, M-1) / PRDZ(N, M) 
                        DISPU = - SQZ(N, M-1)

                    ENDIF
                ENDIF
                  
                IF (M .NE. NLAY)    THEN
                    IF (IBOUND(N, M+1) .NE. 0)    THEN

                        TMAT(MNC, MNC + NCOL) =  -(1-COEV(N, M)) * QFLV(N, M) / PRDZ(N,M) + SQZ(N, M)
                    
                        TMATD = - COEV(N, M) * QFLV(N, M) / PRDZ(N, M)                  
                        DISPD = - SQZ(N, M)

                    ENDIF
                ENDIF

! test
!                IF((MNC .EQ. 2) .OR. (MNC .EQ. 142))   THEN
!                    WRITE(*, *) "TMATBCpost"
!                    WRITE(*, "(I5)") MNC
!                    WRITE(*, "(F10.1)") TMATL, TMATR, TMATU, TMATD
!                    WRITE(*, "(F10.1)") DISPL, DISPR, DISPU, DISPD
!                ENDIF
! end test
            
                !TMAT(MNC, MNC) = ((-QFLH(N-1, M) + QFLH(N, M))/DELC(N, M) + (-QFLV(N-1, M) + QFLV(N, M))/DZ(N, M))  &
                ! &                /PRSITI(N, M) - 2*SQX(N, M) - 2*SQZ(N, M) - 1/TPLEN

                TMATBC = TMATL + TMATR + TMATU + TMATD
                DISPXZ = DISPL + DISPR + DISPU + DISPD
  
                TMAT(MNC, MNC) = TMATBC + DISPXZ - 1/TPLEN 
            
            ENDIF

        ENDDO
    ENDDO

    DEALLOCATE(VELDH)
    DEALLOCATE(VELDV)

    RETURN

END SUBROUTINE TRANS_MEDIA_MATRIX


! formulate the FD right hand side (rhs) of transport model
SUBROUTINE TRANS_MEDIA_RHS()

    USE GLOBAL, ONLY: TPLEN, NLAY, NCOL, IBOUND, DELC, VBLOCK, KONV
    USE GWF,    ONLY: QFLH, QFLV, RCH, COEH, COEV
    USE CON,    ONLY: NNODE, QEXCON, CONCCONN, ICON, KCCOND
    USE TRANS,  ONLY: PRSITI, SQX, SQZ, PRDC, PRDZ, SSM, CONC2D, CONCLTP, CONC, CONCLIMP, TMAT, TRHS
    USE INIT,   ONLY: SCONC

    INTEGER (KIND = 4)  :: M, N, MNC, T

    REAL (KIND = 8)     :: TRHSBC

    TRHS = 0

    DO M = 1, NLAY
        DO N = 1, NCOL

            MNC = NCOL * (M - 1) + N

            ! constant concentration boundary
            IF (IBOUND(N, M) .EQ. -1)    THEN

                TRHS(MNC) = SCONC(N, M)
        
            ELSEIF (IBOUND(N, M) .EQ. 0)    THEN

                TRHS(MNC) = 0

            ELSE

                ! source(rch) term only: there is no sink term in this model
                ! which should be added in the FD matrix
            
                TRHS(MNC) = - CONCLTP(MNC)/TPLEN - RCH(N, M) * SSM(N, M)/ (VBLOCK(N, M) * PRSITI(N, M)) 

!                TRHS(MNC) = - CONCLIMP(MNC)/TPLEN - RCH(N, M) * SSM(N, M)/ (VBLOCK(N, M) * PRSITI(N, M)) + TRHSBC
            
                IF ((ICON(N, M) .GT. 0) .AND. (KONV .EQ. 1))  THEN

                    ! from matrix to conduit            
                    IF (QEXCON(ICON(N, M)) .GT. 0)   THEN
            
                         TMAT(MNC, MNC) = TMAT(MNC, MNC) - QEXCON(ICON(N, M))/(VBLOCK(N, M) * PRSITI(N, M))

                        ! need to be checked!
!                         TRHS(MNC) = TRHS(MNC) + QEXCON(ICON(N, M)) * CONC2D(N, M)/ (VBLOCK(N, M) * PRSITI(N, M))
                   
                    ! from conduit to matrix
                    ELSE
        
                        TRHS(MNC) = TRHS(MNC) + QEXCON(ICON(N, M)) * CONCCONN(ICON(N, M)) / (VBLOCK(N, M) * PRSITI(N, M))

                    ENDIF
                ENDIF

! test
!                IF((MNC .EQ. 2) .OR. (MNC .EQ. 5042))  THEN
!                   ! WRITE(*, "(F15.5)")  COEH(N, M), COEH(N-1, M), COEV(N, M), COEH(N, M-1)
!                   ! WRITE(*, *) "TEST"
!                   ! WRITE(*, "(F15.5)")  QFLH(N-1, M), QFLH(N, M), QFLV(N, M-1), QFLV(N, M) 
!                   ! WRITE(*, *) "TEST"
!                    WRITE(*,*) "NOTRANS"
!                    WRITE(*, "(I5)") MNC
!                    DO T = 1, NLAY*NCOL
!                        IF(TMAT(MNC,T) .NE. 0)    THEN
!                           WRITE(*, "(I5)")    T
!                           WRITE(*, "(F15.3)") TMAT(MNC, T)
!                        ENDIF
!                    ENDDO   
!                    WRITE(*, *) "TRHS"
!                    WRITE(*, "(F15.2)")  TRHS(MNC)
!                 ENDIF
! end test
 
            ENDIF
        ENDDO
    ENDDO

!    WRITE(*, "(10F15.5)") CONCCONN

    
! test
!    MNC = 5041
!    WRITE(*, *) "NOTRANS"
!    WRITE(*, "(I5)") MNC
!    DO T = 1, NLAY*NCOL
!       IF(TMAT(MNC,T) .NE. 0)    THEN
!           WRITE(*, "(I5)")    T
!           WRITE(*, "(F15.3)") TMAT(MNC, T)
!       ENDIF
!    ENDDO   
!    WRITE(*, *) "TRHS"
!    WRITE(*, "(F15.2)")  TRHS(MNC)
! end test

    RETURN

END SUBROUTINE TRANS_MEDIA_RHS


SUBROUTINE SOLVE_MEDIA_CONC()
! SUBROUTINE SOLVE_MEDIA_CONC(CONC, CONCLTP, TMAT, TRHS, DENRHO, SCONC)

    USE GLOBAL, ONLY: TVH, NCOL, NLAY
    USE TRANS,  ONLY: CONCLTP, CONC, TMAT, TRHS, DENRHO, CONC2D
    USE INIT,   ONLY: SCONC, SDRHO

    INTEGER (KIND = 4) :: N, NRHS, LDA, LDB, INFO
    INTEGER (KIND = 4), DIMENSION(TVH)   :: IPIV  
 
    INTEGER (KIND = 4) :: K

    REAL (KIND = 8), DIMENSION(TVH, TVH) :: A
    REAL (KIND = 8), DIMENSION(TVH)      :: B

    ! DGESV: LAPACK linear algebra package function to solve A * X = B
    ! where A is an N-by-N matrix, X and B are N-by-NRHS matrices
    ! N: the number of linear equations, the order of the matrix A
    N = TVH        

    ! NRHS: the number of right hand sides, i.e., the number of columns of the matrix
    NRHS = 1 

    ! A: the N-by-N coefficient matrix A
    A = TMAT

    ! LDA: the leading dimension of the array A
    LDA = TVH 

    ! B: the N-by-NRHS matrix of right hand side
    B = TRHS
    
    ! LDB: the leading dimension of the array B
    LDB = TVH 

    ! INFO: = 0: successful exit

    CALL DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)

    IF(INFO .NE. 0) THEN
        STOP
    ENDIF

    CONC = B
 
    RETURN

END SUBROUTINE SOLVE_MEDIA_CONC



