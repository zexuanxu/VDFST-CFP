MODULE GLOBAL

    IMPLICIT NONE

    ! NPER: number of stress period
    ! NSTP: number of timestep in each stress period
    ! NCOL: number of column
    ! NLAY: number of layer
    INTEGER (KIND = 4), SAVE :: NPER, NSTP, NCOL, NROW, NLAY
    INTEGER (KIND = 4), SAVE :: TVH

    INTEGER (KIND = 4), SAVE :: IMPCON, IMPAG, IMPNO, KONV

    ! PERLEN: length of each stress period
    ! TPLEN: length of each timestep
    REAL (KIND = 8), SAVE    :: PERLEN, TPLEN

    ! -1: constant head/constant concentration cells
    !  1: variable head/variable concentration cells
    ! BC should be the same for flow and transport in this model
    INTEGER (KIND = 4), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: IBOUND

    INTEGER (KIND = 4), SAVE    :: BANDSUB, BANDSUP
    ! the number of variable head/variable concentration cells in each layer
    ! REAL, SAVE, ALLOCATABLE, DIMENSION(:)     :: IBNLAY    

    ! cell width
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: DZ

    ! freshwater density
    REAL (KIND = 8), SAVE    :: DENFR    
    
    ! The slope of linear equation of state that relates fluid density to solute concentration
    REAL (KIND = 8), SAVE :: DRHODC

    ! cell volumn
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:) :: VBLOCK

    ! aquifer top elevation on the left and right hand side
    REAL (KIND = 8), SAVE  :: ELELEFT, ELERIGHT

    ! top of the model
	REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:) :: HTOP

    ! column width of each cell
	REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:) :: DELC

    ! thickness of each layer
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:) :: DELR

    ! elevation of each layer (center of cell)
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:) :: ELAY

END MODULE GLOBAL


MODULE GWF 

    IMPLICIT NONE

    ! specific storage (with density term)
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: SF

    ! HY: horizontal hydraulic conductivity
    ! HV: vertical hydraulic conductivity
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: HY, HV

    ! CC: horizontal conductance
    ! CV: vertical conductance
    ! NOTE: CC and CV are not hydraulic conductivity
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: CC, CR, CV

    ! used for groundwater flow matrix calculation
    ! the conductance of domain boundary
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CCHL, CCHR, CRVT, CRVB, CVVT, CVVB

    ! recharge and sink/source term in flow model 
    ! RCH here is the volume of recharge/sink/source term [L^3]
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: RCH

    ! density of recharge/sink/source term
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DRCH

    ! head and head at last timestep
    ! HEAD and HEADLTP are 1D array in this model
    ! which needs to convert for output
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: HEAD
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: HEADLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: HEADLIMP
 
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: HEAD2D
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: HEAD2DLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: HEAD2DLIMP

    ! horizontal and vertical flow
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: QFLH
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: QFLV

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: COEH, COEV
  
    ! FD matrix for groundwater flow model
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: GMAT

    ! density between two cells
    ! DENH: horizontal density (NCOL-1, NLAY) 
    ! DENV: vertical density (NCOL, NLAY-1)
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DENH, DENV
   
    ! Rho cap, which has different definition with DENV, DENH
    ! see SEAWAT manual, page 20
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DENCH, DENCV

    ! The calculated rhs vector
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: GRHS


END MODULE GWF


MODULE CON 

    IMPLICIT NONE

    ! total conduit node number
    INTEGER (KIND = 4)   :: NNODE

    ! total tube number
    INTEGER (KIND = 4)   :: TNODE
    
    ! total variable conduit head node
    INTEGER (KIND = 4)   :: TVHNODE

    ! gravitational acceralation 
    REAL (KIND = 8)    :: GRAVAC

    ! friction factor of conduit wall
    ! it's a constant in this model
    REAL (KIND = 8)    :: FRIC
 
    ! dispersivity within conduit
    REAL (KIND = 8)    :: DISCON

    ! The H constant in conduit, which equals to SQRT(2*GRAVCC/FRIC)    
    REAL (KIND = 8)   :: HCON

    ! residual of newton-rasphon method in conduit calculation
    REAL (KIND = 8)   :: RESCON

    REAL (KIND = 4), ALLOCATABLE, DIMENSION(:)    :: TVHEADCON

    ! conduit diameter
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: DIAM

    ! cross section of each conduit tube
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: AREACON
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: TUBELEN

    ! KCCON: the exchange coefficient between conduit and porous media
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: KCC

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: KCCOND

    ! HAD(CN) = HCON * AREACON(CN) * SQRT(DIAM(CN)))
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: HAD

    ! DSX(CN) = DISCON/(TUBELEN(CN) * TUBELEN(CN))
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: DSX

    ! conduit node height
    ! probably do not need it
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)  :: ZCON

    ! conduit node locations array, 3 columns
    ! column 1: nnode number
    ! column 2: column number of conduit in cells
    ! column 3: layer number of conduit in cells
    INTEGER (KIND = 4), ALLOCATABLE, DIMENSION(:,:)    :: CONLOC

    ! ICON(NCOL, NLAY)
    ! ICON = 1, cells with conduit
    ! ICON = -1, cells without conduit
    INTEGER (KIND = 4), ALLOCATABLE, DIMENSION(:,:)    :: ICON

    ! head and head at last timestep in conduit
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: HEADCON
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: HEADCONLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: HEADCONLIMP

    ! CONCCONN: concentration in conduit node 
    ! CONCCONT: concentration in conduit tube
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CONCCONN
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CONCCONNLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CONCCONNLIMP

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CONCCONT
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CONCCONTLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: CONCCONTLIMP

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: DENCONT
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: DENCONTLTP

    ! density and density at last timestep in conduit
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: DENRHOCON
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: DENRHOCONLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: DENRHOCONLIMP

    ! conduit flow rate
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: QCON
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: QFLCON
    ! conduit flow flux
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: QCONF

    ! flow exchange between conduit and porous media
    ! QEXCON > 0, flow direction is from matrix to conduit node
    ! QEXCON < 0, flow direction is from conduit node to matrix
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: QEXCON

    ! source term of conduit
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: QSCON
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: QSCONCONC

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: CONJAC   
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: GCON

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: TCMAT
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)   :: TCRHS

    REAL (KIND = 8) :: RESHEADCON

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: FRAT
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: REYNOLDS
 
END MODULE CON 


MODULE TRANS

    ! porosity
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)    :: PRSITI

    ! DXX: horizontal dispersivity
    ! DZZ: vertical dispersivity
    REAL (KIND = 8)   :: DXX, DZZ

    ! PRDX: PRSITI * DELC(N, M)
    ! PRDZ: PRSITI * DZ(N, M)
    ! which are going to be use in the FD matrix of transport simulation
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)    :: PRDC, PRDZ

    ! SQX(N, M) = DXX / (DELC(N, M) * DELC(N, M))
    ! SQZ(N, M) = DZZ / (DZ(N, M) * DZ(N, M))
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)    :: SQX, SQZ

    ! concentration of recharge and sink/source term
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)    :: SSM
   
    ! density and density at last timestep
    ! density terms are 2D array in this model
    ! which needs to convert from concentration in each timestep
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DENRHO
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DENRHOLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: DENRHOLIMP

    ! concentration and concentration at last timestep
    ! this is 1D array, which is the same as HEAD
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)    :: CONC
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)    :: CONCLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)    :: CONCLTPL
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)    :: CONCLIMP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)    :: CONCLIMPL

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)  :: CONC2D
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)  :: CONC2DLTP
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)  :: CONC2DLIMP

    ! concentration at two timesteps forward
    ! REAL, DIMENSION(:)    :: CONCLTPX  

    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:)  :: TMAT  
    REAL (KIND = 8), ALLOCATABLE, DIMENSION(:)    :: TRHS

END MODULE TRANS

MODULE INIT

    ! initial head
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: STRT

    ! initial concentration
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: SCONC

    ! initial density
    REAL (KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: SDRHO


END MODULE INIT


MODULE BUDGET

    REAL (KIND = 8) :: CONSTBUD, RCHBUD, STOBUDGWF, STOBUDTRANS, DCDT, CONSTBUDIN, CONSTBUDOUT

    REAL (KIND = 8) :: DIFF, PERC, TOTIN, TOTOUT
    REAL (KIND = 8) :: DIFFCON, PERCCON, TOTINCON, TOTOUTCON

    REAL (KIND = 8) :: QEXCONBUD, QSCONBUD, CONSTOBUD


END MODULE BUDGET

