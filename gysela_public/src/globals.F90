!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!----------------------------------------
! file : globals.f90
! date : 16/09/2005
! global variable definition
!----------------------------------------
module globals
  use prec_const
  
  implicit none
      
  !**********************************************
  !   GLOBAL VARIABLES FOR MPI PARALLELISATION
  !********************************************** 
  integer, parameter :: cible        = 0    
  integer, parameter :: cible_mu     = 0
  integer, parameter :: cible_jstart = 0
  integer, parameter :: cible_istart = 0
  ! Local space domain size in r (i) and theta (j)
  integer :: istart, iend, inumber
  integer :: jstart, jend, jnumber
  ! Local space domain size in phi (k) and vpar (l)
  integer :: kstart, kend, knumber
  integer :: lstart, lend, lnumber
  ! Local space domain size in  r (i) and theta (j) 
  !  with an extra ghost zone : buf
  integer  :: jstart_buf, istart_buf, jend_buf, iend_buf
  ! Size of the local patch (domain decomposition) in 
  ! directions r, theta, phi and vpar
  integer  :: dom_r, dom_theta, dom_phi, dom_vpar
  ! Identity of the processor
  integer  :: pglobal_id
  ! Number of processors
  integer  :: Nbproc_phi, Nbproc_vpar 
  integer  :: Nbproc_tot
  ! Size of the ghost zone
  integer, parameter :: bufsize = 3
  ! Size of the stencil 
  integer, parameter :: stencil = 10
  ! Index of the local mu value
  integer  :: mu_id
  ! Communicator shared between processors 
  !   that owns the same mu value
  integer  :: mpi_comm_mu
  ! Communicator shared between processors 
  !   that have the same plocal_id
  integer  :: mpi_comm_intermu
  ! Communicator with all the processors having the same mu_id 
  !  and the same jstart
  integer  :: mpi_comm_column
  ! Communicator with all the processors having the same mu_id 
  !  and the same istart
  integer  :: mpi_comm_row
  ! Communicator with all the processors managing a phi value
  integer  :: mpi_comm_moments
  ! Identity of processor in the 'mpi_comm_mu' communicator
  integer  :: plocal_id
  ! Identity of processor in the 'mpi_comm_column' communicator
  integer  :: jstart_id
  ! Identity of processor in the 'mpi_comm_row' communicator
  integer  :: istart_id
  ! Number of procs in 'mpi_comm_mu' communicator
  integer  :: Nbproc_loc
  integer  :: Nbpoints_muid
  ! Indicate if diagnostics should be parallel or not
  logical  :: diag_para
  ! Variables that defines the processor targets used
  ! for diagnostic computation and savings
  integer,  PARAMETER :: nbrp  = 6
  integer             :: diag_targ(1:nbrp)
  ! Number of digits that ends the filenames of diagnostics
  integer             :: nbdigits = 5
  character(20)       :: numfmt     = "'_d',i5.5"
  character(20)       :: numfmt_rst = "'_r',i3.3"
  ! Identity of processor responsible for the output file "gysela_res.out"
  integer  :: outputproc
      
  integer     :: iter_run
  real(RKIND) :: dt_save
      
  !********************************************
  !   GLOBAL VARIABLES FOR INPUT DATAS
  !********************************************
  ! Parallel datas
  integer, save :: Nbproc_r        = 1
  integer, save :: Nbproc_theta    = 1
  integer, save :: Nbproc_mu       = 1
  integer, save :: bloc_phi        = 1
  integer, save :: Nbthread        = 1
  logical, save :: transpose4D     = .false.
  logical, save :: transp_velocity = .false.
      
  ! Geometry datas
  integer    , save :: Nr             = 50
  integer    , save :: Ntheta         = 50
  integer    , save :: Nphi           = 50
  integer    , save :: Nvpar          = 50
  integer    , save :: Nmu            = 0
  real(RKIND), save :: a              = 1._RKIND
  real(RKIND), save :: &
    Ltheta         = 6.283185307179586476925286766559005768394_RKIND
  real(RKIND), save :: &
    Lphi           = 6.283185307179586476925286766559005768394_RKIND
  real(RKIND), save :: aspect_ratio   = 3._RKIND
  real(RKIND), save :: nb_vth0        = 3._RKIND
  real(RKIND), save :: Lmu            = 0._RKIND
  !-> a*rhomin < r < a*rhomax
  real(RKIND), save :: rhomin         = 0._RKIND
  !-> a*rhomin < r < a*rhomax    
  real(RKIND), save :: rhomax         = 1._RKIND
  real(RKIND), save :: mumin          = 0._RKIND
  logical    , save :: integral_vperp = .false.
  !*** Temporary value
  integer    , save :: imu0           = 0
      
  ! Equilibrium datas
  logical    , save :: plasma_current = .false.
  logical    , save :: rbar_in_feq    = .false.
  logical    , save :: use_psibar     = .false.
  integer    , save :: Zi             = 1
  real(RKIND), save :: tau0           = 1._RKIND
  ! -> for temperature and density profiles
  !  = 1: 1/T(r) dT(r)/dr = -kappaT*func(cosh^-2)
  !  = 2: 1/T(r) dT(r)/dr = -kappaT*func(tanh)
  logical    , save :: read_n0        = .false.
  logical    , save :: read_Te        = .false.
  logical    , save :: read_Ti        = .false.
  logical    , save :: read_q         = .false.
  integer    , save :: profile_choice = 1       
  real(RKIND), save :: rpeak          = 0.5_RKIND
  real(RKIND), save :: kappan         = 0._RKIND
  real(RKIND), save :: kappaTi        = 1.0_RKIND/0.28_RKIND
  real(RKIND), save :: kappaTe        = 0._RKIND
  real(RKIND), save :: deltarn        = 1._RKIND/3._RKIND
  real(RKIND), save :: deltarTe       = 1._RKIND/3._RKIND
  real(RKIND), save :: deltarTi       = 1._RKIND/3._RKIND
  logical    , save :: magnetic_shear = .true.
  real(RKIND), save :: q0             = 1._RKIND
  real(RKIND), save :: deltaq         = 2.7858_RKIND
  real(RKIND), save :: alphaq         = 2.8_RKIND
  logical    , save :: reversed_shear = .false.
  real(RKIND), save :: qmin           = 1.5_RKIND
  real(RKIND), save :: qa             = 3._RKIND
  real(RKIND), save :: rho_qmin       = 0.5_RKIND
  !-> The case B_curvature=.false. corresponds to a 
  !     temporary test case where the effects of the
  !     curvature of the magnetic field are not 
  !     taken into account in the advection equations
  logical    , save :: B_curvature    = .true.  
  !  -> initial perturbation
  real(RKIND), save :: epsilon        = 0.001_RKIND
  integer    , save :: m              = 10
  integer    , save :: n              = 3
  ! -> Initialization of the perturbation of 
  !      the distribution function
  !    = 1 if initialization with the mode (m,n)
  !    = 2 if initialization with several 
  !        modes 1 < im < m and 1 < in < n 
  !    = 3 if initialization of several in m for one 
  !       mode in n (1 < m and m > n*q)         
  integer    , save :: init_perturb   = 2   
  logical    , save :: zonal_flow     = .true.
  integer    , save :: filter_choice  = 0
  real(RKIND), save :: filter_deltam  = 5._RKIND
  logical    , save :: ExB_shear      = .false.
  real(RKIND), save :: max_shear      = 1._RKIND
  real(RKIND), save :: rho_shear      = 0.5_RKIND
  real(RKIND), save :: deltar_shear   = 0.1_RKIND
      
  ! Algorithm datas
  !-> predictor-corrector or leap-frog
  logical    , save :: leapfrog       = .true.
  integer    , save :: LF_average     = 50
  integer    , save :: hffilterfreq   = 0  ! frequency
  integer    , save :: hffiltertype   = 0  ! (0=no,1=FFT,2=wavelet)
  integer    , save :: hffilterfref   = 1  ! ref function type
  real(RKIND), save :: advec2Ddt      = 0._RKIND
  real(RKIND), save :: reduceddt      = 0._RKIND
  real(RKIND), save :: reducedbegin   = -1e30_RKIND
  real(RKIND), save :: reducedend     = -1e30_RKIND
  !-> time algo specifications
  real(RKIND), save :: deltat         = 0.1_RKIND
  integer    , save :: nbiter         = 1
  integer    , save :: diag_nbstep    = 1
  integer    , save :: refresh_nbstep = 0
  ! Test datas
  logical    , save :: memory_test       = .false.
  logical    , save :: gyroaverage       = .true.
  logical    , save :: slab_geometry     = .false.
  logical    , save :: Rosenbluth_Hinton = .false. 
  logical    , save :: m0n0_eq0          = .false. 
  logical    , save :: single_m          = .false.
  logical    , save :: single_n          = .false.
  logical    , save :: solve_poiss       = .true.
  real(RKIND), save :: coeff_pol         = 1._RKIND
  logical    , save :: RHS_only          = .false.
  integer    , save :: feq_choice        = 1
  integer    , save :: vpar0_option      = 1
  logical    , save :: Phi00BCr_Neumann  = .false.
  logical    , save :: hvpar_in_fperturb = .true.
      
  ! Output datas
  logical   , save :: integration_CS    = .true.
  integer   , save :: diag_level        = 3
  logical   , save :: Phi3D_saving      = .false.
  logical   , save :: FMoments3D_saving = .false.
  logical   , save :: f5D_saving        = .false.
  logical   , save :: CalviExp_saving   = .false.
  integer   , save :: Phi3D_nbstep      = 1
  integer   , save :: FMoments3D_nbstep = 1
  integer   , save :: f5D_nbstep        = 1
  integer   , save :: CalviExp_nbstep   = 1
  logical   , save :: rst_saving        = .true.
      
  !********************************************
  !   GLOBAL VARIABLES 
  !********************************************
  logical            :: restart
  character(LEN=20)  :: file_name_res = 'gysela_res.out'
  integer, save      :: uout_res           = 100
  integer, save      :: nbiter_prev        = 0 
  integer, save      :: nb_restart         = 0
  integer, save      :: num_restart        = 0
  integer, save      :: nb_restart_diag    = 0
  integer, save      :: nb_diag            = 0
  integer, save      :: Phi3D_nb_diag      = 0
  integer, save      :: FMoments3D_nb_diag = 0
  integer, save      :: f5D_nb_diag        = 0
      
  integer, parameter :: nbmoments_GC   = 24
  integer, parameter :: nbmoments_part = 10
      
  !********************************************
  !  BOUNDARY CONDITIONS
  !   = 0 if Dirichlet condition (func=0)
  !   = 1 if Neumann (func'=0)
  !   = 2 if Hermite conditions 
  !       (func'=approximation of the derivate)
  !********************************************
  integer, parameter :: BC_Hermite    = 2
  integer, parameter :: Phi_BCr_left  = BC_Hermite
  integer, parameter :: Phi_BCr_right = BC_Hermite
      
  ! power for the integral computation
  !(ipow=0 if Nmu=0 and Lmu=0.,ipow=1 otherwise)
  integer           :: ipow
      
  ! physical datas
  real(RKIND) :: R0            !(major radius)
  ! Iteration following
  integer           :: iter_glob
  real(RKIND)       :: init_time
      
  ! Input arguments
  integer     , save :: num_args = 1
  character(LEN=100) :: string   
      
  ! used for equilbrium part 
  real(RKIND), dimension(:,:)  , pointer :: M2_equil
  real(RKIND), dimension(:,:)  , pointer :: nGieq_rtheta
  real(RKIND), dimension(:)    , pointer :: neq_r
  real(RKIND), dimension(:)    , pointer :: dneq_dr
  real(RKIND), dimension(:)    , pointer :: Tieq_r
  real(RKIND)                            :: nbions_eq
      
  ! used for number of electrons and ions 
  real(RKIND) :: nbelectrons_diag, nbions_diag
      
  ! min and max error (r-rbar(r,theta,vpar,mu))
  real(RKIND) :: min_err_rbar, max_err_rbar
      
  ! position for the output saving
  integer     :: ir_save, itheta_save, iphi_save
      
  ! position for the output f2D saving
  ! for trapped particles
  integer     :: ir_f2D_trapped, itheta_f2D_trapped, &
    iphi_f2D_trapped, ivpar_f2D_trapped, imu_f2D_trapped
      
  ! position for the output f2D saving
  ! for passing particles
  integer     :: ir_f2D_passing, itheta_f2D_passing, &
    iphi_f2D_passing, ivpar_f2D_passing, imu_f2D_passing
      
  ! position for the output Phi2D saving
  integer     :: ir_Phi, itheta_Phi, iphi_Phi
      
  ! For P, V and T computation
  real(RKIND), dimension(:,:,:), pointer :: Pmean_loc
  real(RKIND), dimension(:,:,:), pointer :: Vmean_loc
  real(RKIND), dimension(:,:,:), pointer :: Tmean_loc
      
  ! For precomputation of Bstar and contravariant 
  ! components of bstar
  real(RKIND), dimension(:,:,:), pointer :: Bstar_3D
  real(RKIND), dimension(:,:,:), pointer :: bstar_gradr_3D
  real(RKIND), dimension(:,:,:), pointer :: bstar_gradtheta_3D
  real(RKIND), dimension(:,:,:), pointer :: bstar_gradphi_3D
      
  ! Partial derivatives of J0.Phi(r,theta,phi)
  real(RKIND), dimension(:,:,:), pointer, public :: dJ0Phidr
  real(RKIND), dimension(:,:,:), pointer, public :: dJ0Phidtheta
  real(RKIND), dimension(:,:,:), pointer, public :: dJ0Phidphi
      
  ! Coefficients according to the case treated 
  real(RKIND) :: K_curv
      
  ! Time value for output saving
  real(RKIND) :: time_diag
      
  !********************************************
  !   GLOBAL VARIABLES FOR DISPLACEMENT STATISTICS
  !********************************************
  integer, dimension(:), pointer :: nbpart_rleft_iter
  integer, dimension(:), pointer :: nbpart_rright_iter
  integer, dimension(:), pointer :: nbpart_thleft_iter
  integer, dimension(:), pointer :: nbpart_thright_iter
  integer, dimension(:), pointer :: nbpart_rleft_muid
  integer, dimension(:), pointer :: nbpart_rright_muid
  integer, dimension(:), pointer :: nbpart_thleft_muid
  integer, dimension(:), pointer :: nbpart_thright_muid
      
  !********************************************
  !   GLOBAL VARIABLES FOR USEFUL ARRAYS
  !********************************************
  !*** -> for QN parallel solver (version 2 et 3) ***
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_NrNthetamuphi_nbM2
  real(RKIND), &
    dimension(:,:,:,:), pointer :: Rarray_PNrPNthetaNphi_nbM2
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_recv_nbM5
  real(RKIND), &
    dimension(:,:,:,:), pointer :: Rarray_send_nbM5
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_PNr_nbM5
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_Nr_nbM5
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_PNthetaNphi_nbM2
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_NthetaNphi_nbM2
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_PNr_nbM6
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_Nr_nbM6
  real(RKIND), &
    dimension(:,:)    , pointer :: Rarray_NrNthetamuphi_nbM8
  real(RKIND), &
    dimension(:,:,:,:), pointer :: Rarray_PNrPNthetaNphi_nbM8
      
  !*** real array allocations ***
  real(RKIND), dimension(:)    , pointer :: Rarray1_Nr
  real(RKIND), dimension(:)    , pointer :: Rarray2_Nr
  real(RKIND), dimension(:)    , pointer :: Rarray3_Nr
  real(RKIND), dimension(:)    , pointer :: Rarray4_Nr
  real(RKIND), dimension(:)    , pointer :: Rarray1_Ntheta
  real(RKIND), dimension(:)    , pointer :: Rarray2_Ntheta
  real(RKIND), dimension(:)    , pointer :: Rarray1_Nphi
  real(RKIND), dimension(:)    , pointer :: Rarray2_Nphi
  real(RKIND), dimension(:)    , pointer :: Rarray1_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray2_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray3_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray4_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray5_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray6_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray7_Nvpar
  real(RKIND), dimension(:)    , pointer :: Rarray_m1Nphip1
  real(RKIND), dimension(:)    , pointer :: Rarray_m1Nvparp1
  real(RKIND), dimension(:,:)  , pointer :: Rarray_1Ntheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray1_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray2_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray3_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray4_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray5_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray6_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray7_NrNtheta
  real(RKIND), dimension(:,:)  , pointer :: Rarray_NrNphi
  real(RKIND), dimension(:,:)  , pointer :: Rarray2_NrNphi
  real(RKIND), dimension(:,:)  , pointer :: Rarray_NrNvpar
  real(RKIND), dimension(:,:)  , pointer :: Rarray_NphiNvpar
  real(RKIND), dimension(:,:)  , pointer :: Rarray_NthetaNvpar
  real(RKIND), dimension(:,:)  , pointer :: Rarray_NvparNmu
  real(RKIND), dimension(:,:)  , pointer :: Rarray1_m1Nrp1m1Nthetap1
  real(RKIND), dimension(:,:)  , pointer :: Rarray2_m1Nrp1m1Nthetap1
  real(RKIND), dimension(:,:)  , pointer :: Rarray3_m1Nrp1m1Nthetap1
  real(RKIND), dimension(:,:)  , pointer :: Rarray_m1Nphip1m1Nvparp1
  real(RKIND), dimension(:,:,:), pointer :: Rarray1_NrNthetaDomphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray2_NrNthetaDomphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray3_NrNthetaDomphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray4_NrNthetaDomphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray5_NrNthetaDomphi
!baoter ????
  real(RKIND), dimension(:,:,:), pointer :: Rarray1_NrNthetaNphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray2_NrNthetaNphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray3_NrNthetaNphi
!eaoter
  real(RKIND), dimension(:,:,:), pointer :: Rarray1_PNrPNthetaNphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray2_PNrPNthetaNphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray3_PNrPNthetaNphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray4_PNrPNthetaNphi
  real(RKIND), dimension(:,:,:), pointer :: Rarray5_PNrPNthetaNphi
      
  complex(RKIND), &
    dimension(:,:,:), pointer :: Carray_1Nrp11Ntheta1Nphi
end module globals
