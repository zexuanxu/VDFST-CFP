
! call four subroutines to assgin values
SUBROUTINE ASSIVAL()

    CALL ASSIVALGLO()
    CALL ASSIVALGWF()

    CALL ASSIVALCON()
    CALL ASSIVALTRANS()

    CALL ASSIVALBUD()

    RETURN

END SUBROUTINE ASSIVAL


! assgin values for global parameters
SUBROUTINE ASSIVALGLO()

    USE GLOBAL

    INTEGER (KIND = 4):: K, I
    
    NCOL = 120
    NROW = 21

    NLAY = NROW

    ALLOCATE(IBOUND(NCOL, NROW))
    ALLOCATE(VBLOCK(NCOL, NROW))
    ALLOCATE(HTOP(NCOL))
    ALLOCATE(DELC(NCOL, NROW))
    ALLOCATE(DELR(NCOL, NROW))
    ALLOCATE(DZ(NCOL, NROW))
    ALLOCATE(ELAY(NCOL, NROW))

    NPER = 1
!    NSTP = 10

    PERLEN = 1.0
!    PERLEN = 5.0
!    TPLEN = PERLEN/NSTP
! timestep length should be less than 0.0001, which is 8.64 second
    TPLEN =  0.001

    NSTP = PERLEN / TPLEN

    DZ = 50 
   
    ! value should be changed due to different unit in length
    DENFR = 64.00
    DRHODC = 0.0457
! test
!    DRHODC = 0.0

    TVH = NCOL * NROW

    ELELEFT =  -334.0834
    ELERIGHT = -319.0504

!    ELELEFT = -30.0
!    ELERIGHT = -30.0

    DO K = 1, NROW
        DO I = 1, NCOL

            ! constant BC (-1) 
            IF(I .EQ. 1)     THEN
                IBOUND(I, K) = -1
            ELSEIF(I .EQ. NCOL)   THEN
                IBOUND(I, K) = -1
            ELSE
                IBOUND(I, K) = 1
            ENDIF

            DELC(I, K) = 500 
            DELR(I, K) = 500
                  
            VBLOCK(I, K) = DELC(I, K) * DELR(I, K) * DZ(I, K) 

            IF(I .EQ. 1)    THEN
                HTOP(I) = ELELEFT
                ELAY(I, K) = HTOP(I)
            ELSE
                HTOP(I) = HTOP(I-1) + (ELERIGHT - ELELEFT) / (NCOL-1) 
                ELAY(I, K) = HTOP(I)
            ENDIF
 
        ENDDO
    ENDDO

    BANDSUB = 0
    BANDSUP = 0

!    WRITE(*, "(10F15.0)") VBLOCK

!        IF(IBOUND(1, K) .GT. 0) THEN
!            BANDSUB = NCOL
!        ENDIF

!        IF(IBOUND(NCOL, K) .GT. 0)  THEN
!            BANDSUP = NCOL
!        ENDIF

    RETURN

END SUBROUTINE ASSIVALGLO


! assign values for parameter in groundwater flow model 
SUBROUTINE ASSIVALGWF()


    USE GLOBAL,     ONLY: NCOL, NROW, TVH, DELR, DELC, DZ, BANDSUB, BANDSUP, DENFR
    USE GWF

    INTEGER (KIND = 4) :: K, I, T

    ALLOCATE(SF(NCOL, NROW))

    ALLOCATE(HY(NCOL, NROW))
    ALLOCATE(HV(NCOL, NROW))

    ALLOCATE(CC(NCOL-1, NROW))
    ALLOCATE(CR(NCOL, NROW-1))
    ALLOCATE(CV(NCOL, NROW-1))

    ALLOCATE(CCHL(NROW))
    ALLOCATE(CCHR(NROW))
    ALLOCATE(CRVT(NCOL))
    ALLOCATE(CRVB(NCOL))
    ALLOCATE(CVVT(NCOL))
    ALLOCATE(CVVB(NCOL))

    ALLOCATE(RCH(NCOL, NROW)) 
    ALLOCATE(DRCH(NCOL, NROW))

    ALLOCATE(HEAD(TVH))
    ALLOCATE(HEADLTP(TVH))

    ALLOCATE(QFLH(NCOL-1, NROW))
    ALLOCATE(QFLV(NCOL, NROW-1))
    
    ALLOCATE(COEH(NCOL-1, NROW))
    ALLOCATE(COEV(NCOL, NROW-1))

    ALLOCATE(GMAT(TVH, TVH))
    ALLOCATE(DENH(NCOL-1, NROW))
    ALLOCATE(DENV(NCOL, NROW-1))    

    ALLOCATE(DENCH(NCOL-1, NROW))
    ALLOCATE(DENCV(NCOL, NROW-1))

    ALLOCATE(GRHS(TVH))

    ALLOCATE(HEAD2D(NCOL, NROW))
    ALLOCATE(HEAD2DLTP(NCOL, NROW))

    SF = 5.0E-005
!    SF = 5.0E-006
!    SF(:,11) = 0.05

    HY = 7500
!    HY(:,11) = 2000000

    HV = HY

    DO K = 1, NROW
        DO I = 1, NCOL   

            IF(K .NE. NROW) THEN
                CR(I, K) = ((HY(I, K) + HY(I, K+1)) / 2) * DZ(I, K) * DELC(I, K) / DELR(I, K)
            ENDIF 
            
            IF(I .NE. NCOL) THEN
                CC(I, K) = ((HY(I, K) + HY(I+1, K)) / 2) * DZ(I, K) * DELR(I, K) / DELC(I, K)
            ENDIF
        ENDDO
    ENDDO

    CV = CR

    DO K = 1, NROW
        CCHL(K) = HY(1, K) * DELR(1, K) * DZ(1, K) / DELC(1, K)
        CCHR(K) = HY(NCOL, K) * DELR(NCOL, K) * DZ(NCOL, K) / DELC(NCOL, K)  
    ENDDO

    DO I = 1, NCOL
        CRVT(I) = HY(I, 1) * DELC(I, 1) * DZ(I, 1) / DELR(I, 1)
        CRVB(I) = HY(I, NROW) * DELC(I, NROW) * DZ(I, NROW) / DELR(I, NROW)
    ENDDO

    CVVT = CRVT
    CVVB = CRVB

    DO K = 1, NROW
        DO I = 1, NCOL
            RCH(I, K) = 0.03 * DELR(I, K) * DELC(I, K)
        ENDDO
    ENDDO

! for test only
    RCH = 0

    DRCH = DENFR

    QFLH = 0
    QFLV = 0

    COEH = 0
    COEV = 0

    GMAT = 0
    GRHS = 0

    GMATLU = 0

    DENH = 0
    DENV = 0

    RETURN

END SUBROUTINE ASSIVALGWF    


SUBROUTINE ASSIVALCON() 

    USE GLOBAL, ONLY: NCOL, NROW, DELC, DZ, ELAY
    USE CON

    INTEGER (KIND = 4) ::CN, CT
    REAL (KIND = 8) :: DZCON

    NNODE = NCOL
    TNODE = NNODE - 1   

    TVHNODE = NCOL - 2

    ALLOCATE(TVHEADCON(NNODE))

    DO CN = 1, NNODE
        IF (CN .EQ. 1)  THEN
            TVHEADCON(CN) = 1.0
        ELSEIF (CN .EQ. NNODE) THEN
            TVHEADCON(CN) = 5.0
        ELSE
            TVHEADCON(CN) = -1.0
        ENDIF    
    ENDDO

    ALLOCATE(DIAM(TNODE))
    ALLOCATE(AREACON(TNODE))
    ALLOCATE(TUBELEN(TNODE))

    ALLOCATE(KCC(NNODE))
    ALLOCATE(KCCOND(NNODE))
    ALLOCATE(HAD(TNODE))
    ALLOCATE(DSX(TNODE))

    ALLOCATE(ZCON(NNODE))

    ALLOCATE(CONLOC(3, NNODE))
    ALLOCATE(ICON(NCOL, NROW))

    ALLOCATE(HEADCON(NNODE))
    ALLOCATE(HEADCONLTP(NNODE))
    ALLOCATE(CONCCONN(NNODE))
    ALLOCATE(CONCCONNLTP(NNODE))
    ALLOCATE(CONCCONT(TNODE))
    ALLOCATE(CONCCONTLTP(TNODE))
    ALLOCATE(DENRHOCON(NNODE))
    ALLOCATE(DENRHOCONLTP(NNODE))
    ALLOCATE(DENRHOCONLIMP(NNODE))

    ALLOCATE(QCON(TNODE))
    ALLOCATE(QFLCON(TNODE))
    ALLOCATE(QCONF(TNODE))
    ALLOCATE(QEXCON(NNODE))
    ALLOCATE(QSCON(NNODE))
    ALLOCATE(QSCONCONC(NNODE))

    ALLOCATE(CONJAC(NNODE, NNODE))
    ALLOCATE(GCON(NNODE))

    ALLOCATE(TCMAT(NNODE, NNODE))
    ALLOCATE(TCRHS(NNODE))
   
    ALLOCATE(CONCCONNLIMP(NNODE))
    ALLOCATE(CONCCONTLIMP(TNODE))

    ALLOCATE(DENCONT(TNODE))
    ALLOCATE(DENCONTLTP(TNODE)) 

    ALLOCATE(FRAT(TNODE))
    ALLOCATE(REYNOLDS(TNODE))

    ! value should be changed due to different unit in length
    GRAVAC = 32.169 * (60.0**2) * (60.0**2) * (24.0**2) 

    ! roughness height = 0.1
    ! friction factor will be calculated by Colebrook-White formula
    ! FRIC = 0.0316
    ! FRIC = 0.0178
    ! FRIC = 0.200
    ! FRIC = 5.00

    ! UNSURE about the dispersivity in conduit
    DISCON = 1.00
 
! comment for friction update   
!    HCON = SQRT(2*GRAVAC/FRIC)
    HCON = SQRT(2*GRAVAC)    

    ! residence criterion for Newton-Rasphon method
    RESCON = 0.01

    DO CT = 1, TNODE    
        DIAM(CT) = 3.00
        AREACON(CT) = 3.1416 * (DIAM(CT)/2) * (DIAM(CT)/2)
        HAD(CT) = HCON * AREACON(CT) * SQRT(DIAM(CT))
    ENDDO

    HEADCON = 0
    HEADCONLTP = 0

    CONCCONN = 0
    CONCCONNLTP = 0

    CONCCONT = 0
    CONCCONTLTP = 0

    DENRHOCON = 0
    DENRHOCONLTP = 0
 
    ICON = -1

    DO CN = 1, NNODE
 
        CONLOC(1, CN) = CN
        CONLOC(2, CN) = CN
        CONLOC(3, CN) = 11

        ICON(CONLOC(2, CN), CONLOC(3, CN)) = CN  
        ZCON(CN) = ELAY(CONLOC(2, CN), CONLOC(3, CN))

    ENDDO

    DO CT = 1, TNODE

        DZCON = ELAY(CONLOC(2, CT), CONLOC(3, CT)) - ELAY(CONLOC(2, CT+1), CONLOC(3, CT+1)) 
        TUBELEN(CT) = SQRT(DELC(CONLOC(2, CT), CONLOC(3, CT))**2 + DZCON ** 2) 
        
        DSX(CT) = DISCON/(TUBELEN(CT) * TUBELEN(CT))

    ENDDO

!    KCC = 500
    KCC = 7500
    
    QCON = 0
    QFLCON = 0
    QCONF = 0
    QEXCON = 0
    QSCON = 0

    QSCONCONC = 0

    !QSCON > 0, directly recharge into conduit node
    !QSCON(NNODE) = - 5000000

    CONJAC = 0
    GCON = 0

    TCMAT = 0
    TCRHS = 0

    RESHEADCON = 1e-8 

    DENCONT = 0

    CALL CALC_KCCOND()

!   WRITE(*, "(10I5)") ICON

    FRAT = 0.01

    RETURN

END SUBROUTINE ASSIVALCON


! assign values for parameters in transport model
SUBROUTINE ASSIVALTRANS()

    USE GLOBAL, ONLY: NROW, NCOL, TVH, DELC, DELR, BANDSUB, BANDSUP
    USE TRANS

    INTEGER (KIND = 4) :: K, I

    ALLOCATE(PRSITI(NCOL, NROW))

    ALLOCATE(PRDC(NCOL, NROW))
    ALLOCATE(PRDZ(NCOL, NROW))

    ALLOCATE(SQX(NCOL, NROW))
    ALLOCATE(SQZ(NCOL, NROW))

    ALLOCATE(SSM(NCOL, NROW))

    ALLOCATE(CONC(TVH))
    ALLOCATE(CONCLTP(TVH))    
    ALLOCATE(CONCLTPL(TVH))   

    ALLOCATE(CONC2D(NCOL, NROW))
    ALLOCATE(CONC2DLTP(NCOL, NROW))
 
    ALLOCATE(DENRHO(NCOL, NROW))
    ALLOCATE(DENRHOLTP(NCOL, NROW))
    ALLOCATE(DENRHOLIMP(NCOL, NROW))

    ALLOCATE(CONCLIMP(TVH))
    ALLOCATE(CONCLIMPL(TVH))

    ALLOCATE(TMAT(TVH, TVH))

    ALLOCATE(TRHS(TVH))    

    PRSITI = 0.003
!   PRSITI(:,11) = 0.300

    DXX = 32.80
    DZZ = 32.80

! test
!    DXX = 500.0
!    DZZ = 500.0
! end test

    SSM = 0

    DO K = 1, NROW
        DO I = 1, NCOL

            PRDC (I, K) = PRSITI(I, K) * DELC(I, K)
            PRDZ (I, K) = PRSITI(I, K) * DELR(I, K)

        ENDDO
    ENDDO    

    SQX = 0
    SQZ = 0

    DENRHO = 0
    DENRHOLTP = 0

    TMAT = 0
    TRHS = 0

    TMATLU = 0

    RETURN

END SUBROUTINE ASSIVALTRANS

SUBROUTINE ASSIVALINIT()

    USE GLOBAL, ONLY: NCOL, NROW, DENFR, DRHODC, TVH
    USE GWF,    ONLY: HEAD, HEADLTP, HEAD2D, HEAD2DLTP 
    USE CON,    ONLY: NNODE, HEADCON, CONCCONN, CONCCONT, CONCCONTLTP, DENRHOCON, CONLOC
    USE TRANS,  ONLY: CONC, CONCLTP, CONCLTPL, DENRHO, DENRHOLTP, CONC2D
    USE INIT,   ONLY: STRT, SCONC, SDRHO

    INTEGER (KIND = 4) :: K, I, CN
    REAL (KIND = 8)    :: HLEFT, HRIGHT
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)    :: HF

    ALLOCATE(STRT(NCOL, NROW))
    ALLOCATE(SCONC(NCOL, NROW))
    ALLOCATE(SDRHO(NCOL, NROW))
    ALLOCATE(HF(NCOL, NROW))

    HLEFT =  0.0
!    HLEFT =  8.315
    HRIGHT = 5.0

    DO K = 1, NROW       
        DO I = 1, NCOL  
    
            IF(I .EQ. 1)   THEN
                SCONC(I, K) = 35.00 
!             ELSEIF((I .GT. 1) .AND. (I .LE. 35))   THEN
!                SCONC(I, K) = SCONC(I-1, K) - 1.0
            ELSE
                SCONC(I, K) = 0
            ENDIF
  
            IF(I .EQ. 1)    THEN
                STRT(I, K) = HLEFT
            ELSE
                STRT(I, K) = STRT(I-1, K) + (HRIGHT - HLEFT)/(NCOL-1)
            ENDIF

            SDRHO(I, K) = DENFR + DRHODC * SCONC(I, K) 
        
         ENDDO

! for density used only: to garantee hydrostatic condition
!         HLEFT = HLEFT - 0.25
    
    ENDDO

    DENRHO = SDRHO
    DENRHOLTP = DENRHO
    
    ! convert measured head to equivalent freshwater head
    CALL MEA2HF(STRT, HF)
    
    STRT = HF

!    WRITE(*, "(10F10.2)")  HF
!    STOP

    ! convert from 2D to 1D
    CALL MAT2VEC(STRT, HEAD)
    CALL MAT2VEC(STRT, HEADLTP)

    CALL MAT2VEC(SCONC, CONC)
    CALL MAT2VEC(SCONC, CONCLTP)
    CALL MAT2VEC(SCONC, CONCLTPL)

    HEAD2D = STRT
    CONC2D = SCONC

    HEADLTP = HEAD
    HEAD2DLTP = HEAD2D

    CALL ASSINITCON()

    DEALLOCATE(HF)

    RETURN

END SUBROUTINE ASSIVALINIT


SUBROUTINE ASSINITCON()

    USE GLOBAL, ONLY: DENFR, DRHODC
    USE GWF,    ONLY: HEAD2D
    USE CON,    ONLY: NNODE, TVHEADCON, HEADCON, CONCCONN, CONCCONT, DENRHOCON, QSCONCONC, CONCCONTLTP, CONLOC, ZCON, DENCONT
    USE TRANS,  ONLY: CONC2D

    INTEGER(KIND = 4) :: CN

    DO CN = 1, NNODE  
        
        IF (CN .EQ. 1)  THEN
            CONCCONN(CN) = 35.0 
            QSCONCONC(CN) = 35.0
        ELSE    
            CONCCONN(CN) = CONC2D(CONLOC(2, CN), CONLOC(3, CN))
        ENDIF

        DENRHOCON(CN) = DENFR + DRHODC * CONCCONN(CN)

        IF ((CN .EQ. 1) .OR. (CN .EQ. NNODE))  THEN
            HEADCON(CN) =  DENRHOCON(CN) * TVHEADCON(CN) / DENFR - (DENRHOCON(CN) - DENFR) * ZCON(CN) / DENFR
!            HEADCON(CN) = HEAD2D(CONLOC(2, CN), CONLOC(3, CN))
        ELSE
!            HEADCON(CN) = 0
!            HEADCON(CN) = DENRHOCON(CN) * TVHEADCON(CN) / DENFR - (DENRHOCON(CN) - DENFR) * ZCON(CN) / DENFR
            HEADCON(CN) = HEAD2D(CONLOC(2, CN), CONLOC(3, CN))
        ENDIF
    
        IF(CN .GE. 2)   THEN

            CONCCONT(CN-1) = (CONCCONN(CN-1) + CONCCONN(CN))/2
            CONCCONTLTP(CN-1) = CONCCONT(CN-1)        
    
            DENCONT(CN-1) = DENFR + DRHODC * CONCCONT(CN-1)
        ENDIF
        
    ENDDO

!    WRITE(*, "(10F10.2)") CONCCONN
!     HEADCON = HEAD2D(:,11)    
!    WRITE(*, "(10F10.5)") HEADCON 

!    WRITE(*,*) "TEST"

    RETURN

END SUBROUTINE ASSINITCON


SUBROUTINE ASSIVALBUD()

    USE BUDGET

    CONSTBUD = 0
    RCHBUD = 0
    STOBUDGWF = 0
    STOBUDTRANS = 0
    DCDT = 0

    CONSTBUDIN = 0
    CONSTBUDOUT = 0

    DIFF = 0 
    PERC = 0
    TOTIN = 0
    TOTOUT = 0

    DIFFCON = 0
    PERCCON = 0 
    TOTINCON = 0
    TOTOUTCON = 0

    QEXCONBUD = 0
    QSCONBUD = 0
    CONSTOBUD = 0

    RETURN

END SUBROUTINE ASSIVALBUD




