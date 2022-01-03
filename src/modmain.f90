module modmain

integer Dimen ! # lattice dimension

!interger parmeters (input)
integer TauGridNum ! # of grid for k, energy space
integer iWnGridNum ! # of grid for matsubara frequency
integer K1GridNum   ! # of grid for k, energy space
integer K2GridNum   ! # of grid for k, energy space
integer K3GridNum   ! # of grid for k, energy space
integer R1GridNum   ! # of grid for k, energy space
integer R2GridNum   ! # of grid for k, energy space
integer R3GridNum   ! # of grid for k, energy space
integer SCFloop    ! # Maximum # of self-consistent loop
integer, parameter :: MaxOrbitN = 100 ! matrix dimension , dimension index spin
integer MatDim  ! matrix dimension , dimension index spin
integer OrbN  ! number of orbital
integer, parameter :: InitOccup = 1 ! inital Ocuupancy choice 1:para, 2:ferro, !3:anti-ferro

!for matrix eigenvalue and eigenvector solver
integer lwork

!for saving data
integer allgreenfSave, locselfESave, LogOutOn

!for doing Hartree-Fock
integer DoHF


!for doing CTQMC++
integer DoCTQMC, DoGW_CTQMC

!for check whether num of grid enough or not
integer iWnGridCheck, TauGridCheck1, TauGridCheck2


!real parameters (input)
real(8) a1(3), a2(3), a3(3) ! lattice vector
real(8) b1(3), b2(3), b3(3) ! reciprocal lattice distance

!real(8) t0(MaxOrbitN,MaxOrbitN) ! Onsite energy (orbital energy)
!real(8) t1(MaxOrbitN,MaxOrbitN) ! NN hopping parameter
!real(8) t2(MaxOrbitN,MaxOrbitN) ! 2NN hopping parameter
!real(8) t3(MaxOrbitN,MaxOrbitN) ! 3NN hopping parameter
!real(8) U(MaxOrbitN,MaxOrbitN)  ! Onsite U
!real(8) U1(MaxOrbitN,MaxOrbitN) ! NN intersite U
!real(8) U2(MaxOrbitN,MaxOrbitN) ! 2NN intersite U
!real(8) U3(MaxOrbitN,MaxOrbitN) ! 3NN intersite U
real(8) Tem ! Temperature
real(8) Ntot ! Total # of electrons ! 1.0 half filled
real(8) GWMixing_Pol, GWMixing_SelfE ! Mixing for GW Pol_k_iWn, SelfE_k_iWn
real(8) DMFTMixing_Pol, DMFTMixing_SelfE ! Mixing for DMFT Pol_iWn, SelfE_iWn
real(8) MuShift

real(8), parameter :: PI = 3.1415926535897932384626433d0
real(8), parameter :: boltzmanC=0.0000861733034d0 ! boltzman constant
real(8) Erange    ! Energy range 15.0 means range from -15 and 15
real(8) beta ! 
real(8) delta   ! small segment for green function

real(8) ChemTH ! threashold for convergency
real(8) ConvTH ! threashold for convergency

real(8) iWnGridCut! grid determination by taking |(2n+1)pi/beta|<iWnGridCut


!real(8), parameter :: iWnGridCut=50.0 ! grid determination by taking |(2n+1)pi/beta|<iWnGridCut
integer dumGNum ! # of grid for matsubara frequency
!integer, parameter :: iWnGridNum=dumGNum*2 ! # of grid for matsubara frequency



! allocatable cairable
real(8), allocatable :: SpinSplit(:) !for Spin spliting. This values are used for chemical potential shift of spin up(SpinSplit/2) an down(-SpinSplit/2) only at first SCF step. 
real(8), allocatable :: t0(:,:) ! Onsite energy (orbital energy)
real(8), allocatable :: t1(:,:) ! NN hopping parameter
real(8), allocatable :: t2(:,:) ! 2NN hopping parameter
real(8), allocatable :: t3(:,:) ! 3NN hopping parameter
real(8), allocatable :: U(:,:)  ! Onsite U
real(8), allocatable :: U1(:,:) ! NN intersite U
real(8), allocatable :: U2(:,:) ! 2NN intersite U
real(8), allocatable :: U3(:,:) ! 3NN intersite U
real(8), allocatable :: OrbitPosi(:,:) !fractional
real(8), allocatable :: OrbitPosi_car(:,:) !cartesian


integer, allocatable :: NeighborN(:,:) ! Numberof Neighbor
integer, allocatable :: distingAtom(:,:) ! Numberof Neighbor
real(8), allocatable :: NeighborDist(:,:) ! reciprocal lattice distance
real(8), allocatable :: distingPos(:,:) ! reciprocal lattice distance
real(8), allocatable :: NeighborInfo(:,:,:,:) 
real(8), allocatable :: K1list(:,:),K2list(:,:),K3list(:,:), Taulist(:), R1list(:,:),R2list(:,:),R3list(:,:),Taulist2(:)
complex(8), allocatable :: localE(:)
complex(8), allocatable :: Noccup(:)
complex(8), allocatable :: ChemOccup(:)

complex(8), allocatable :: tq(:,:,:,:,:), Vq(:,:,:,:,:)
complex(8), allocatable :: iWnlist(:), iWnlist2(:)
complex(8), allocatable :: Greenf_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: NonIntGreenf_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: NonIntGreenf0ChemP_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: Polarization_R_tau(:,:,:,:,:,:)
complex(8), allocatable :: SelfEcorr_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: Hartree_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: Fock_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: ScreenedWcorr_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: W_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: dummy_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: dummy_k_iWn2(:,:,:,:,:,:)
complex(8), allocatable :: dummy_k_tau(:,:,:,:,:,:)
complex(8), allocatable :: dummy_k_tauII(:,:,:,:,:,:)
complex(8), allocatable :: dummy_R_iWn(:,:,:,:,:,:)
complex(8), allocatable :: dummy_R_tau(:,:,:,:,:,:)
complex(8), allocatable :: dummy_R_tauInv(:,:,:,:,:,:)
complex(8), allocatable :: dummy_R_tau2(:,:,:,:,:,:)
complex(8), allocatable :: dummy_R_tauII(:,:,:,:,:,:)
complex(8), allocatable :: dummy_R_tauIII(:,:,:,:,:,:)
complex(8), allocatable :: locGreenf(:,:,:)
complex(8), allocatable :: locW(:,:,:)
complex(8), allocatable :: locPol(:,:,:)
complex(8), allocatable :: locSelfE(:,:,:)
complex(8), allocatable :: NonIntHeigVal(:,:,:,:,:)
complex(8), allocatable :: ImpuritySelfE(:,:), ImpurityPol(:,:)
complex(8), allocatable :: GWSelfE_k_iWn_Old(:,:,:,:,:,:)
complex(8), allocatable :: GWPol_k_iWn_Old(:,:,:,:,:,:)
complex(8), allocatable :: GWCTQMCSelfE_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: GWCTQMCPol_k_iWn(:,:,:,:,:,:)
complex(8), allocatable :: Hartree_loc(:,:,:)
complex(8), allocatable :: Fock_loc(:,:,:)


end module
