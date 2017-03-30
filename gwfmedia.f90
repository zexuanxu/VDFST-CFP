! formulate FD matrix for groundwater flow model 
! SUBROUTINE GWF_MEDIA_MATRIX(GMAT, DENRHO)
SUBROUTINE GWF_MEDIA_MATRIX()

    USE GLOBAL, ONLY: TPLEN, NCOL, NLAY, TVH, IBOUND, VBLOCK
    USE GWF,    ONLY: SF, CC, CV, GMAT, DENH, DENV, DENCH, DENCV, CCHL, CCHR, CVVT, CVVB, HEAD2D
    USE TRANS,  ONLY: DENRHO, CONCLTP 

    INTEGER (KIND = 4):: K, I, M, N, MNC

    ! cumulative deleted cells
    REAL (KIND = 8) :: GMATBC
    REAL (KIND = 8) :: GMATL, GMATR, GMATU, GMATD

!    WRITE(*, "(10F10.2)") CV 

    ! calculate upstream density, and arithmetic mean of density
    DO K = 1, NLAY   
        DO I = 1, NCOL

            IF(K .NE. NLAY) THEN
        
                DENV(I, K) = (DENRHO(I, K) + DENRHO(I, K+1))/2
            
                IF(HEAD2D(I, K) .GT. HEAD2D(I, K+1))    THEN

                    DENCV(I, K) = DENRHO(I, K)
                
                ELSEIF(HEAD2D(I, K) .LT. HEAD2D(I, K+1))    THEN

                    DENCV(I, K) = DENRHO(I, K+1)

                ELSE

                    DENCV(I, K) = DENV(I, K)
                
                ENDIF
            ENDIF
       
            IF(I .NE. NCOL) THEN

                DENH(I, K) = (DENRHO(I, K) + DENRHO(I+1, K))/2 
  
                IF(HEAD2D(I, K) .GT. HEAD2D(I+1, K))    THEN

                    DENCH(I, K) = DENRHO(I, K)
                
                ELSEIF(HEAD2D(I, K) .LT. HEAD2D(I+1, K))    THEN

                    DENCH(I, K) = DENRHO(I+1, K)

                ELSE

                    DENCH(I, K) = DENH(I, K)
                ENDIF
            ENDIF

        ENDDO
    ENDDO
  
    GMAT = 0    

    ! the first row represents grid(1,1) -- (NCOL, NLAY)
    ! the second row represents grid(2,1) 

    ! h(i,k-1)
    ! GMAT(NLAY*(M-1)+N, NLAY*(M-2)+N)

    ! h(i-1, k)
    ! GMAT(NLAY*(M-1)+N, NLAY*(M-1)+N-1)

    ! h(i, k)
    ! GMAT(NLAY*(M-1)+N, NLAY*(M-1)+N)

    ! h(i+1, k)
    ! GMAT(NLAY*(M-1)+N, NLAY*(M-1)+N+1)

    ! h(i,k+1)
    ! GMAT(NLAY*(M-1)+N, NLAY*M+N)
 
!    WRITE(*, "(10F10.2)") DENCV 
!    WRITE(*, "(I5)") NLAY
 
    DO M = 1, NLAY
        DO N = 1, NCOL

            MNC = NCOL * (M-1) + N
    
            ! constant head boundary
            IF(IBOUND(N, M) .EQ. -1)  THEN
            
                GMAT(MNC, MNC) = 1

            ! no flow boundary: do nothing
            ELSEIF(IBOUND(N, M) .EQ. 0) THEN
            
            ELSE
            
                GMATBC = 0
                GMATL = 0
                GMATR = 0        
                GMATU = 0
                GMATD = 0
            
                IF (N .NE. 1)   THEN
                    ! flow exchange between no-flow cells is not considered
                    IF (IBOUND(N-1, M) .NE. 0)  THEN
                    
                        GMAT(MNC, MNC - 1) = DENCH(N-1, M) * CC(N-1, M) 

                        GMATL = - DENCH(N-1, M) * CC(N-1, M)

                    ENDIF
                ENDIF
                
                IF (N .NE. NCOL)   THEN
                    IF (IBOUND(N+1, M) .NE. 0)    THEN

                        GMAT(MNC, MNC + 1) = DENCH(N, M) * CC(N, M)

                        GMATR = - DENCH(N, M) * CC(N, M)

                    ENDIF
                ENDIF

                IF (M .NE. 1)   THEN
                    IF (IBOUND(N, M-1) .NE. 0)  THEN

                        GMAT(MNC, MNC - NCOL) = DENCV(N, M-1) * CV(N, M-1)
   
                        GMATU = - DENCV(N, M-1) * CV(N, M-1)
                 
                    ENDIF              
                ENDIF

                IF (M .NE. NLAY)     THEN 
                    IF (IBOUND(N, M+1) .NE. 0)   THEN

                        GMAT(MNC, MNC + NCOL) = DENCV(N, M) * CV(N, M)

                        GMATD = - DENCV(N, M) * CV(N, M)

                    ENDIF
                ENDIF

                GMATBC = GMATL + GMATR + GMATU + GMATD

                GMAT(MNC, MNC) = GMATBC - DENRHO(N, M) * SF(N, M) * VBLOCK(N, M) / TPLEN

!                WRITE(*, "(I5)") MNC

!                IF((MNC .EQ. 2) .OR. (MNC .EQ. 2402))   THEN
!                    WRITE(*, "(F15.0)") GMATBC 
!                ENDIF

            ENDIF

        ENDDO
    ENDDO

!    WRITE(*, "(F15.0)") GMAT(2402,2282), GMAT(2402,2402), GMAT(2402,2403)

    RETURN

END SUBROUTINE GWF_MEDIA_MATRIX


! formulate the FD right hand side (rhs) of groundwater flow model
! SUBROUTINE GWF_MEDIA_RHS (GRHS, GMAT, DENH, DENV, CONC, CONCLTP) 
SUBROUTINE GWF_MEDIA_RHS()

    USE GLOBAL, ONLY: KONV, TPLEN, NLAY, NCOL, TVH, IBOUND, DRHODC, VBLOCK, ELAY, DENFR
    USE GWF,    ONLY: SF, CC, CV, RCH, DRCH, GRHS, GMAT, DENH, DENV, DENCH, DENCV, HEAD, HEADLTP, HEAD2D
    USE CON,    ONLY: QEXCON, ICON, DENRHOCON, KCCOND, HEADCON
    USE TRANS,  ONLY: PRSITI, CONCLTP, CONCLTPL, DENRHO, CONCLIMP, CONCLIMPL
    USE INIT,   ONLY: STRT

    INTEGER(KIND = 4) :: IMPCON
    INTEGER (KIND = 4):: K, I, M, N, MNC, T
    REAL (KIND = 8)    :: GRHSD

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DRHS

    ! The Dik term on the rhs
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DRHSHM, DRHSHP, DRHSVM, DRHSVP

    REAL (KIND = 8) :: QEXTEST

    ALLOCATE(DRHS(NCOL, NLAY))
    ALLOCATE(DRHSHM(NCOL-1, NLAY))
    ALLOCATE(DRHSHP(NCOL-1, NLAY))
    ALLOCATE(DRHSVM(NCOL, NLAY-1))
    ALLOCATE(DRHSVP(NCOL, NLAY-1))

    DRHS = 0
    DRHSHM = 0
    DRHSHP = 0
    DRHSVM = 0
    DRHSVP = 0

    ! calculate the density division terms
    DO K = 1, NLAY
        DO I = 1, NCOL

            IF(I .NE. 1)    THEN

                DRHSHM(I-1, K) = DENCH(I-1, K) * CC(I-1, K) * (DENH(I-1, K) - DENFR) * (ELAY(I-1, K) - ELAY(I, K))/DENFR

            ENDIF

            IF(I .NE. NCOL)  THEN
               
                DRHSHP(I, K) = DENCH(I, K) * CC(I, K) * (DENH(I, K) - DENFR) * (ELAY(I+1, K) - ELAY(I, K))/DENFR
            
            ENDIF
        
            IF(K .NE. 1)    THEN

                DRHSVM(I, K-1) = DENCV(I, K-1) * CV(I, K-1) * (DENV(I, K-1) - DENFR) * (ELAY(I, K-1) - ELAY(I, K))/DENFR

            ENDIF

            IF(K .NE. NLAY) THEN

                DRHSVP(I, K) = DENCV(I, K) * CV(I, K) * (DENV(I, K) - DENFR) * (ELAY(I, K+1) - ELAY(I, K))/DENFR
            
            ENDIF
        ENDDO
    ENDDO
    
    ! calculate the relative desity-difference terms
    DO K = 1, NLAY
        DO I = 1, NCOL

            IF (I .NE. 1)  THEN

                DRHS(I, K) = DRHS(I, K) + DRHSHM(I-1, K)

            ENDIF

            IF (I .NE. NCOL) THEN

                DRHS(I, K) = DRHS(I, K) + DRHSHP(I, K)

            ENDIF

            IF (K .NE. 1)   THEN

                DRHS(I, K) = DRHS(I, K) + DRHSVM(I, K-1)

            ENDIF

            IF (K .NE. NLAY) THEN

                DRHS(I, K) = DRHS(I, K) + DRHSVP(I, K)
            
            ENDIF

        ENDDO
    ENDDO

    GRHS = 0

!    WRITE(*, *) "SUMQEXCON"
!    WRITE(*, "(F15.2)") SUM(QEXCON)
!    WRITE(*, "(10F15.5)") QEXCON

    DO M = 1, NLAY
        DO N = 1, NCOL

        MNC = NCOL * (M - 1) + N

            ! constant head boundary        
            IF(IBOUND(N, M) .EQ. -1) THEN

                GRHS(MNC) = STRT(N, M)
        
            ELSEIF(IBOUND(N, M) .EQ. 0) THEN

                GRHS(MNC) = 0

            ELSE

! test
!                GRHS(MNC) = (-DENRHO(N, M) * SF(N, M) * HEAD(MNC) / TPLEN + PRSITI(N, M) *  &
!                &           DRHODC * (CONCLIMP(MNC) - CONCLIMPL(MNC)) / TPLEN ) * VBLOCK(N, M) &
!                &            - RCH(N, M) * DRCH(N, M) - DRHS(N, M) 
! end test

! HEAD(MNC) or HEADLTP(MNC)?

                GRHS(MNC) = (-DENRHO(N, M) * SF(N, M) * HEADLTP(MNC) / TPLEN + PRSITI(N, M) *  &
                &           DRHODC * (CONCLIMP(MNC) - CONCLIMPL(MNC)) / TPLEN ) * VBLOCK(N, M) &
                &            - RCH(N, M) * DRCH(N, M) - DRHS(N, M) 


!                WRITE(*, "(F10.0)")  GRHS(MNC)

! test (freshwater)                 
!               GRHS(MNC) = (-DENRHO(N, M) * SF(N, M) * HEAD(MNC) / TPLEN ) * VBLOCK(N, M) &
!               &            - RCH(N, M) * DRCH(N, M) - DRHS(N, M) + GRHSBC
! end test

                IF ((ICON(N, M) .GT. 0) .AND. (KONV .EQ. 1))  THEN
!                IF (ICON(N, M) .GT. 0)  THEN 
                  ! QEXCON > 0, flow direction is from matrix to conduit node
                   ! Therefore, it's a sink term for matrix
                     IF(QEXCON(ICON(N, M)) .GT. 0) THEN

!                        GRHS(MNC) = GRHS(MNC) + QEXCON(ICON(N, M)) * DENRHO(N, M) 

! test (possibily right, however, conduit exchange should be the sink/source term in matrix simulation)                        
                        GMAT(MNC, MNC) = GMAT(MNC, MNC) - KCCOND(ICON(N, M)) * DENRHO(N, M)
                        GRHS(MNC) = GRHS(MNC) - HEADCON(ICON(N, M)) * KCCOND(ICON(N, M)) * DENRHO(N, M)
! end test                                
                    ELSE

!                        GRHS(MNC) = GRHS(MNC) + QEXCON(ICON(N, M)) * DENRHOCON(ICON(N, M))
! test                    
                        GMAT(MNC, MNC) = GMAT(MNC, MNC) - KCCOND(ICON(N, M)) * DENRHOCON(ICON(N, M))
                        GRHS(MNC) = GRHS(MNC) - HEADCON(ICON(N, M)) * KCCOND(ICON(N, M)) * DENRHOCON(ICON(N, M))
! end test

!                        WRITE(*, "(F15.3)") KCCOND(ICON(N, M)) * DENRHOCON(ICON(N, M))
!                        QEXTEST = QEXTEST + KCCOND(ICON(N, M)) * (HEAD(MNC) - HEADCON(ICON(N, M)))

                    ENDIF
                ENDIF

! test
!                IF((MNC .EQ. 122) .OR. (MNC .EQ. 1202) .OR. (MNC .EQ. 2282))  THEN
!                IF ((MNC .EQ. 162) .OR. (MNC .EQ. 22))    THEN
!                    WRITE(*, "(I5)") MNC
!                    WRITE(*, *) "HEAD"
!                    WRITE(*, "(F10.3)")  HEAD(MNC)
!                    WRITE(*, "(F10.3)")  HEADCON(ICON(N, M))
!                    WRITE(*, "(F10.5)") HEADLTP(MNC)
!                    WRITE(*, "(F10.5)") DENRHO(N, M)
!                    WRITE(*, "(F15.3)") DRHS(N, M)
!                    DO T = 1, NLAY*NCOL
!                        IF(GMAT(MNC,T) .NE. 0)    THEN
!                            WRITE(*, "(I5)")    T
!                            WRITE(*, "(F15.3)") GMAT(MNC, T)
!                        ENDIF
!                    ENDDO
!                    WRITE(*, *) "GRHS"
!                    WRITE(*, "(F15.2)")  GRHS(MNC)
!                ENDIF
! end test

            ENDIF
        ENDDO
    ENDDO

 !    WRITE(*, "(F15.2)")  GRHS

!    WRITE(*, *) "QEXTEST"
!    WRITE(*, "(F15.2)") QEXTEST

! test
!    MNC = 5041 
!    WRITE(*, *) "NOGWF"
!    WRITE(*, "(I5)") MNC
!    DO T = 1, NLAY*NCOL
!        IF(GMAT(MNC,T) .NE. 0)    THEN
!            WRITE(*, "(I5)")    T
!            WRITE(*, "(F15.3)") GMAT(MNC, T)
!        ENDIF
!    ENDDO
!    WRITE(*, *) "GRHS"
!    WRITE(*, "(F15.2)")  GRHS(MNC)
! end test

    DEALLOCATE(DRHSHM)
    DEALLOCATE(DRHSHP)
    DEALLOCATE(DRHSVM)
    DEALLOCATE(DRHSVP)
    DEALLOCATE(DRHS)

    RETURN

END SUBROUTINE GWF_MEDIA_RHS


! solve head using LAPACK numerical solver
SUBROUTINE SOLVE_MEDIA_HEAD()

    USE GLOBAL, ONLY: NCOL, NLAY, TVH, BANDSUB, BANDSUP
    USE GWF,    ONLY: HEAD, HEADLTP, HEAD2D, HEAD2DLTP, QFLH, QFLV, GMAT, GRHS 
    USE INIT,   ONLY: STRT  
 
    INTEGER(KIND=4) :: NS, NRHS, LDAS, LDBS, INFO
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: IPIV    

    INTEGER(KIND=4) :: K, LWORK, IS, JS
    
    INTEGER(KIND=4) :: MS, KL, KU, LDAB

    CHARACTER   :: TRANS

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A, AB 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: B
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: TAU, WORK

    REAL(KIND=8) R(TVH, TVH)
    REAL(KIND=8) Q(TVH, TVH)
    REAL(KIND=8) SO(TVH)

    ALLOCATE(IPIV(TVH))
    ALLOCATE(A(TVH, TVH))
    ALLOCATE(B(TVH))
    
    ALLOCATE(AB(2*BANDSUB+BANDSUP+1, TVH))
    ALLOCATE(TAU(TVH))
    ALLOCATE(WORK(TVH))

    ! DGESV: LAPACK linear algebra package function to solve A * X = B
    ! where A is an N-by-N matrix, X and B are N-by-NRHS matrices

    M = TVH

    ! for DGESV, N should be N = TVH
    ! N: the number of linear equations, the order of the matrix A
    ! N = TVH    
    ! for DGBTRF, N = 1

    N = TVH

    KL = BANDSUB
    KU = BANDSUP

    LDAB = 2 * KL + KU + 1 

    ! NRHS: the number of right hand sides, i.e., the number of columns of the matrix
    NRHS = 1

    ! A: the N-by-N coefficient matrix A
    A = GMAT

    TRANS = 'n'

    ! LDA: the leading dimension of the array A
    LDAS = TVH

    ! B: the N-by-NRHS matrix of right hand side
    B = GRHS

    ! LDB: the leading dimension of the array B
    LDBS = TVH

    ! INFO: = 0: successful exit

    LWORK = TVH

!  BEGIN QR decomposition method
!    CALL DGEQRF(NS, MS, A, LDAS, TAU, WORK, LWORK , INFO)
!   R = 0
!   KQR = TVH
!   DO I = 1, KQR 
!       R(I, I:KQR) = A(I, I:KQR)
!    ENDDO
!    CALL DORMQR('L', 'T', MS, NRHS, KQR, A, LDAS, TAU, B, MS, WORK, LWORK, INFO)
!    DO J = TVH, 1, -1
!        SO(J) = B(J)/R(J,J)
!        DO I = 1, TVH-1
!            B(I) = B(I) - R(I,J) * SO(J)
!        ENDDO
!    ENDDO
! END QR decomposition method

! BEGIN LU decomposition (directly solve)
    ! subroutine of LAPACK fortran to solve AX=B

    CALL DGESV(N, NRHS, A, LDAS, IPIV, B, LDBS, INFO)
 
! END LU decomposition (directly solve)

! BEGIN LU decomposition
!    CALL DGETRF(MS, NS, A, LDAS, IPIV, INFO)
!    CALL DGETRS(TRANS, NS, NRHS, A, LDAS, IPIV, B, LDBS, INFO)
! END LU decomposition

! BEGIN LU decomposition band solve
!    CALL DGBTRF(MS, NS, KL, KU, AB, LDAB, IPIV, INFO)  
!    CALL DGBTRS(TRANS, NS, KL, KU, NRHS, AB, LDAB, IPIV, B, LDBS, INFO)
! END LU decomposition band solve


    IF(INFO .NE. 0) THEN
        STOP
    ENDIF

    HEAD = B    

    DEALLOCATE(IPIV)
    DEALLOCATE(A)
    DEALLOCATE(B)
    
    DEALLOCATE(AB)
    DEALLOCATE(TAU)
    DEALLOCATE(WORK)

    RETURN

END SUBROUTINE SOLVE_MEDIA_HEAD


! calculate the specific flux(horizontal and vertical)
SUBROUTINE FLOWMEDIA()

    USE GLOBAL, ONLY: NCOL, NLAY, DELC, DELR, DZ, DENFR, ELAY, IBOUND
    USE GWF,    ONLY: HY, HV, HEAD2D, QFLH, QFLV
    USE TRANS,  ONLY: DENRHO

    INTEGER (KIND = 4):: K, I

    QFLH = 0
    QFLV = 0

!    WRITE(*, "(10F10.3)") HEAD2D


    DO K = 1, NLAY
        DO I = 1, NCOL
   
            ! QFLH(1, 1): horizontal flow between (1, 1) and (1, 2)
            ! flow from (1, 1) to (1, 2): positive
         
            IF ((I .NE. 1) .AND. (IBOUND(I, K) .NE. 0 ) .AND. (IBOUND(I-1, K) .NE. 0))     THEN
    
                QFLH(I-1, K) = ((HY(I-1, K) + HY(I, K))/2) * ((HEAD2D(I-1, K) - HEAD2D(I, K)) / ((DELC(I, K)+DELC(I-1, K))/2) + &
                & (((DENRHO(I-1, K) + DENRHO(I, K))/2) - DENFR) * (ELAY(I-1, K) - ELAY(I, K))/(DENFR * (DELC(I, K)+DELC(I-1, K))/2))

            ENDIF

            IF ((K .NE. 1) .AND. (IBOUND(I, K) .NE. 0) .AND. (IBOUND(I, K-1) .NE. 0))    THEN
            
!                QFLV(I, K-1) = ((HV(I, K-1) + HV(I, K))/2) * ((HEAD2D(I, K-1) - HEAD2D(I, K))/((DZ(I, K)+DZ(I, K-1))/2) +  &
!                & (((DENRHO(I, K-1) + DENRHO(I, K))/2) - DENFR) / DENFR)  

                QFLV(I, K-1) = ((HV(I, K-1) + HV(I, K))/2) * ((HEAD2D(I, K-1) - HEAD2D(I, K))/((DELR(I, K)+DELR(I, K-1))/2) +  &
               & (((DENRHO(I, K-1) + DENRHO(I, K))/2) - DENFR) * (ELAY(I, K-1) - ELAY(I, K))/(DENFR * (DELR(I, K)+DELR(I, K-1))/2))  

! test
!                IF( (K .EQ. 1) .AND. (I .EQ. 3))  THEN
!                    WRITE(*, *) "DENRHO"
!                    WRITE(*, "(F10.5)") HEAD2D(I, K-1), HEAD2D(I, K)
!                    WRITE(*, "(F10.5)") DENRHO(I, K-1), DENRHO(I, K) 
!                    WRITE(*, "(F10.5)") HEAD2D(I, K-1) - HEAD2D(I, K)
!                    WRITE(*, "(F10.5)") QFLV(I, K-1)
!                ENDIF
! end test

            ENDIF

        ENDDO
    ENDDO

! test
!    WRITE(*, "(10F10.5)") HEAD2D(:,1)
!    WRITE(*, *) "TEST"
!    WRITE(*, "(10F10.5)") HEAD2D(:,10)
!    WRITE(*, *) "TEST"
!    WRITE(*, "(10F10.5)") HEAD2D(:,30)
!    WRITE(*, *) "TEST"

!    WRITE(*, "(10F10.5)") DENRHO(:,1)
!    WRITE(*, "(10F10.5)") DENRHO(:,2)

!    WRITE(*, "(10F15.5)") QFLH(:,10)
!    WRITE(*, *) "TEST"
!    WRITE(*, "(10F15.5)") QFLH(1,:)
!    WRITE(*, "(10F15.5)") QFLV(:,13)
!    WRITE(*, "(10F15.5)") QFLV(:,12)
!    WRITE(*, "(10F15.5)") QFLV(:,11)
!    WRITE(*, "(10F15.5)") QFLV(:,36)
! end test

    RETURN  

END SUBROUTINE FLOWMEDIA






