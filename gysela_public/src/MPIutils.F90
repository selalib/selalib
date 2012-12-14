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
      
!-------------------------------------------------------
! file : MPIutils.f90
! date : 20/03/2003
!  subroutines used for the MPI parallelization
! => Developped by G. LATU 
!---------------------------------------------
module MPIutils_module
  use prec_const
  use globals
  use mem_alloc_module
  implicit none
  include "mpiperso.h"
      
  private
  public :: ppdist1d, pp_mapproc
  public :: pptranspose_backward, pptranspose_forward
  public :: pptranspvel_backward, pptranspvel_forward
  public :: ppinit, ppinit_comm, ppdeallocate
  public :: ppexit, ppbarrier, ppbarrier_timer
  public :: comm_irecv, comm_isend, comm_ssend, comm_test
  public :: comm_wait_recv, comm_wait_send
  public :: comm_gather2D, comm_gather3D
      
  ! Index of the processors in the neiborhood : 
  !   North, South, West ...
  integer :: NN, SS, WW, EE, NW, NE, SW, SE
  ! MPI types to send and receive blocs of 
  !   the distribution function
  integer :: mpitype_WE, mpitype_NS, mpitype_corner
  ! MPI Receiving buffers
  real(RKIND), dimension (:,:), pointer, public :: rbufWW 
  real(RKIND), dimension (:,:), pointer, public :: rbufNN
  real(RKIND), dimension (:,:), pointer, public :: rbufSS
  real(RKIND), dimension (:,:), pointer, public :: rbufEE
  real(RKIND), dimension (:,:), pointer, public :: rbufNW 
  real(RKIND), dimension (:,:), pointer, public :: rbufNE
  real(RKIND), dimension (:,:), pointer, public :: rbufSW
  real(RKIND), dimension (:,:), pointer, public :: rbufSE
  ! MPI Sending buffers
  real(RKIND), dimension (:,:), pointer, public :: sbufWW 
  real(RKIND), dimension (:,:), pointer, public :: sbufNN
  real(RKIND), dimension (:,:), pointer, public :: sbufSS
  real(RKIND), dimension (:,:), pointer, public :: sbufEE
  real(RKIND), dimension (:,:), pointer, public :: sbufNW 
  real(RKIND), dimension (:,:), pointer, public :: sbufNE
  real(RKIND), dimension (:,:), pointer, public :: sbufSW
  real(RKIND), dimension (:,:), pointer, public :: sbufSE
  real(RKIND), dimension (:,:), pointer, public :: fillderiv
  ! request pending for isends
  integer, dimension(1:8) :: reqs
  ! request pending for irecvs
  integer, dimension(1:8) :: reqr
  ! status for irecvs
  integer, dimension(MPI_STATUS_SIZE,1:8) :: status
  ! for testing purpose
!baoter
  integer    , public :: mpitype_locbloc, mpitype_globbloc
  integer    , public :: mpitype_globbloc_rvpar
  real(RKIND), public, dimension(:,:,:)    , pointer :: tmprhs2
  real(RKIND), public, dimension(:,:,:,:)  , pointer :: f4D_transp
  real(RKIND), public, dimension(:,:,:,:)  , pointer :: f4D_send
  real(RKIND), public, dimension(:,:,:,:)  , pointer :: f4D_recv
  real(RKIND), public, &
    dimension(:,:,:,:,:), pointer :: f5D_transp
  real(RKIND), public, dimension(:)        , pointer :: f1D_send
  real(RKIND), public, dimension(:)        , pointer :: f1D_recv
!eaoter
  ! for moment parallelization
  integer                         , public :: dom_mapphi, iphistart
  integer                         , public :: dom_bisphi
  integer                         , public :: collstart, collend
  integer, dimension(:)  , pointer, public :: moments_mapphi
  integer, dimension(:)  , pointer, public :: moments_reqs
  integer, dimension(:)  , pointer, public :: moments_reqr
  integer, dimension(:)  , pointer, public :: bcast_reqs
  integer, dimension(:)  , pointer, public :: bcast_reqr
  integer, dimension(:,:), pointer, public :: moments_status
  integer                                  :: moments_status_size
      
  !******************************
  contains
  !******************************   
  
  !-----------------------------------------------------
  ! global parallel initialization 
  !-----------------------------------------------------
  subroutine ppinit
    use globals, only : Nbproc_tot, pglobal_id, diag_para, &
      nbrp, diag_targ, uout_res
      
    integer :: i, ierr, obt, resultlen
    character (MPI_MAX_PROCESSOR_NAME) :: name
      
#ifdef MPI2
    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, obt, ierr)
#else
    call MPI_INIT(ierr)
#endif
    call MPI_GET_PROCESSOR_NAME(name,resultlen,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,Nbproc_tot,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,pglobal_id,ierr)
      
    !*** define if the diagnostics must be done parallely or ***
    !***  sequentially                                       ***
    diag_para = (Nbproc_tot .ge. 8) 
      
    !*** Defines a set of processor targets used ***
    !*** for diagnostic computation and savings  ***
    do i = 1,nbrp
      diag_targ(i) = cible
    end do
    if (diag_para) then
      outputproc = cible
      if (pglobal_id.eq.outputproc) then
        print*, ' '
        print*, '---> Physical computation : PARALLEL '
        print*, ' '
      end if
      
      diag_targ(1) = 1
      diag_targ(2) = 2
      diag_targ(3) = 3
      diag_targ(4) = 4
      !***  HDF5_f2D_saving                              ***
      !***  performed by the proc. diag_targ(5) = cible  ***
      diag_targ(5) = cible
      diag_targ(6) = 5
    else
      outputproc = 0
      if (pglobal_id.eq.outputproc) then
        print*, ' '
        print*, '---> Physical computation : SEQUENTIAL '
        print*, ' '
      end if
    end if
  end subroutine ppinit
      
  !-----------------------------------------------------
  ! find the item of the global processor having the
  !  value istart, jstart and imu
  !-----------------------------------------------------
  subroutine pp_mapproc(pmu_id,pistart,pjstart,procid)
    integer, intent(in) :: pmu_id,pjstart,pistart
    integer, intent(out) :: procid
    integer  :: base, dep
    dep = pistart / dom_r
    base = pjstart / dom_theta
    procid = pmu_id * Nbproc_loc + Nbproc_r * base + dep
  end subroutine pp_mapproc
  
      
  !-----------------------------------------------------
  ! initialization of the communicators
  !-----------------------------------------------------
  subroutine ppinit_comm
    use globals, only : memory_test
    integer           :: base, dep, size, manage_phivalue
    integer           :: stride
    character(LEN=11) :: string
    character(LEN=80) :: message
    !-> variables for moment parallelisation
    integer           :: i, iphi, iproc, Nbt, ierr
      
    !*** partitioning in r and theta directions ***
    dom_r     = (Nr+1)/Nbproc_r
    dom_theta = Ntheta/Nbproc_theta
    
    !-> the size of the surrouding area must be lower or equal
    !   than the domain size (in r direction or theta direction)
    if ((dom_r < bufsize) .or. (dom_theta < bufsize)) then 
      print *, "(dom_theta or dom_r) < bufsize"
      call ppexit
      stop
    endif
      
    if ((dom_r < stencil+2) .or. (dom_theta < stencil+2)) then 
      print *, "(dom_theta or dom_r) < points needed", &
        " for the derivative"
      call ppexit
      stop
    endif
      
    !*** partitioning in phi and vpar
!baoter (temporary used for fluid moment computation even if 
!         transpose4D=.false.)
!VG!    if (transpose4D) then
      Nbproc_phi  = MIN(Nphi/bloc_phi, Nbproc_r*Nbproc_theta)
      Nbproc_vpar = (Nbproc_r*Nbproc_theta) / Nbproc_phi
      if ( (Nbproc_r*Nbproc_theta) .ne. &
        (Nbproc_phi*Nbproc_vpar)) then
        print *, '(Nbproc_r*Nbproc_theta) .ne. ', &
          '(Nbproc_phi*Nbproc_vpar)'
        print *, 'Nproc_phi = ', Nbproc_phi, &
          'Nproc_vpar = ', Nbproc_vpar
        call ppexit
        stop
      end if
      dom_phi  = Nphi/Nbproc_phi
      if ( mod(dom_phi,bloc_phi) .ne. 0 ) then
        print *, " mod(dom_phi,bloc_phi) .ne. 0 "
        print *, "dom_phi = ", dom_phi, "bloc_phi = ", bloc_phi
        call ppexit
        stop
      end if
      dom_vpar = ceiling((Nvpar+1)/float(Nbproc_vpar))
!VG!    end if
!eaoter
      
    !*** partitioning in mu direction ***
    stride  = Nbproc_tot / Nbproc_mu
    mu_id   = pglobal_id / stride
      
    if (pglobal_id .lt. Nphi) then 
      manage_phivalue = 1
    else
      manage_phivalue = 0
    endif
    !*** Split the global group of processors ***
    !-> communicator with all the processors managing a phi 
    !-> value in the moments computation
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, manage_phivalue, &
      pglobal_id,mpi_comm_moments, ierr)
      
    !*** Split the global group of processors into 'Nbproc_mu' ***
    !-> communicator with all the processors having the same mu
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mu_id, pglobal_id, &
      mpi_comm_mu, ierr)
    !-> get the processor plocal_identity into the new 
    !    mpi_comm_mu communicator
    call MPI_COMM_SIZE(mpi_comm_mu, Nbproc_loc,ierr)
    call MPI_COMM_RANK(mpi_comm_mu, plocal_id,ierr)
    !-> communicator with all the processors having 
    !    the same plocal_id
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, plocal_id, mu_id, &
      mpi_comm_intermu, ierr)
    Nbpoints_muid = (Nr+1)*(Ntheta+1)*(Nphi+1)*(Nvpar+1)
      
    !*** Definition of the blocks on each processors ***
    !-> local coordinates
    base       = (plocal_id/Nbproc_r)
    dep        = mod(plocal_id,Nbproc_r)              
    istart     = dep *  dom_r
    iend       = istart + dom_r - 1
    inumber    = iend - istart + 1
    jstart     = base * dom_theta
    jend       = jstart + dom_theta - 1
    jnumber    = jend - jstart + 1
    istart_buf = istart - bufsize
    iend_buf   = iend + bufsize
    jstart_buf = jstart - bufsize 
    jend_buf   = jend + bufsize
!baoter (temporary used for fluid moment computation
!        even if transpose4D=false)
!VG!    if (transpose4D) then
      base       = (plocal_id/Nbproc_phi)
      dep        = mod(plocal_id,Nbproc_phi)
      kstart     = dep *  dom_phi
      kend       = kstart + dom_phi - 1
      knumber    = kend - kstart + 1
      lstart     = (base * (Nvpar+1)) / Nbproc_vpar
      lend       = ((base+1) * (Nvpar+1)) / Nbproc_vpar - 1
      lnumber    = lend - lstart + 1
      call glob_allocate(f4D_transp,0,Nr,0,Ntheta, &
        kstart,kend,lstart,lend,'f4D_transp')
      call glob_allocate(f4D_send,0,dom_r-1,0,dom_theta-1, &
        0,dom_phi-1,0,dom_vpar-1,'f4D_send')
      call glob_allocate(f4D_recv,0,dom_r-1,0,dom_theta-1, &
        0,dom_phi-1,0,dom_vpar-1,'f4D_recv')
!VG!    end if
!eaoter
      
      if (transp_velocity) then
        dom_bisphi = (Nphi / (Nmu+1))
        if ((Nmu+1)*dom_bisphi .ne. Nphi) dom_bisphi = dom_bisphi+1
        collstart = min(mu_id * dom_bisphi, Nphi-1)
        collend   = min(collstart + dom_bisphi-1, Nphi-1)
        call glob_allocate(f5D_transp,istart,iend,jstart,jend, &
          collstart,collend,0,Nvpar,0,Nmu,'f5D_transp')
        call glob_allocate(f1D_send,0, &
          dom_r*dom_theta*(Nvpar+1)*dom_bisphi-1,'f1D_send')
        call glob_allocate(f1D_recv,0, &
          dom_r*dom_theta*(Nvpar+1)*dom_bisphi-1,'f1D_recv')
      end if
      
    !-> communicator with all the processors having the same 
    !   mu_id and the same jstart
    call MPI_COMM_SPLIT(mpi_comm_mu, jstart, plocal_id, &
      mpi_comm_column, ierr)
    call MPI_COMM_RANK(mpi_comm_column, jstart_id,ierr)
      
    !-> communicator with all the processors having the same
    !   mu_id and the same istart
    call MPI_COMM_SPLIT(mpi_comm_mu, istart, plocal_id, &
      mpi_comm_row, ierr)
    call MPI_COMM_RANK(mpi_comm_row, istart_id,ierr)
      
    !*** Define the local identity of the eight neighbours of ***
    !***  the local processor 'plocal_id'                     ***
    WW   = mod(plocal_id+Nbproc_loc-Nbproc_r, Nbproc_loc) ! West
    EE   = mod(plocal_id+Nbproc_r, Nbproc_loc)            ! East
    base = (plocal_id/Nbproc_r)*Nbproc_r           
    dep  = mod(plocal_id+1,Nbproc_r)              
    SS   = mod(base+dep, Nbproc_loc)                  ! South
    dep  = mod(Nbproc_r+plocal_id-1,Nbproc_r)             
    NN   = mod(base+dep, Nbproc_loc)                  ! North 
    NW   = mod(NN+Nbproc_loc-Nbproc_r, Nbproc_loc)    ! North-West
    NE   = mod(NN+Nbproc_r, Nbproc_loc)               ! North-East
    SW   = mod(SS+Nbproc_loc-Nbproc_r, Nbproc_loc)    ! South-West
    SE   = mod(SS+Nbproc_r, Nbproc_loc)               ! South-East
      
    !*** Define three MPI types needed to transmit blocks ***
    !***   of data to processors in the neighborhood      ***
    !***     mpitype_WE     : block to send to processors ***
    !***                      on the West or East         ***
    !***     mpitype_NS     : block to send to processors ***
    !***                      on the North or South       ***
    !***     mpitype_corner : block to send to processors ***
    !***                      on NE, NW, SW, SE           ***
    !-> mpitype_WE
    if (transpose4D) then
      stride = (3 * (Nr+1) + 6) 
    else
      stride = (3 * dom_r + 6) 
    end if
    size = stride * 3*bloc_phi
    call mpi_type_vector(size, 1, 1, MPI_REAL8, mpitype_WE, ierr)
    call mpi_type_commit(mpitype_WE, ierr)
    call glob_allocate(rbufWW,0,stride-1,0,3*bloc_phi-1,'rbufWW')
    call glob_allocate(rbufEE,0,stride-1,0,3*bloc_phi-1,'rbufEE')
    call glob_allocate(sbufWW,0,stride-1,0,3*bloc_phi-1,'sbufWW')
    call glob_allocate(sbufEE,0,stride-1,0,3*bloc_phi-1,'sbufEE')
    !-> mpitype_NS
    if (transpose4D) then
      stride = (3 * Ntheta + 6) 
    else
      stride = (3 * dom_theta + 6) 
    end if
    size   = stride * 3*bloc_phi
    call mpi_type_vector(size, 1, 1, MPI_REAL8, mpitype_NS, ierr)
    call mpi_type_commit(mpitype_NS, ierr)
    call glob_allocate(rbufNN,0,stride-1,0,3*bloc_phi-1,'rbufNN')
    call glob_allocate(rbufSS,0,stride-1,0,3*bloc_phi-1,'rbufSS')
    call glob_allocate(sbufNN,0,stride-1,0,3*bloc_phi-1,'sbufNN')
    call glob_allocate(sbufSS,0,stride-1,0,3*bloc_phi-1,'sbufSS')
    !-> mpitype_corner
    stride = 9
    size   = stride * 3*bloc_phi
    call mpi_type_vector(size, 1, 1, MPI_REAL8, &
      mpitype_corner, ierr)
    call mpi_type_commit(mpitype_corner, ierr)
    call glob_allocate(rbufNW,0,stride-1,0,3*bloc_phi-1,'rbufNW')
    call glob_allocate(rbufSW,0,stride-1,0,3*bloc_phi-1,'rbufSW')
    call glob_allocate(rbufNE,0,stride-1,0,3*bloc_phi-1,'rbufNE')
    call glob_allocate(rbufSE,0,stride-1,0,3*bloc_phi-1,'rbufSE')
    call glob_allocate(sbufNW,0,stride-1,0,3*bloc_phi-1,'sbufNW')
    call glob_allocate(sbufSW,0,stride-1,0,3*bloc_phi-1,'sbufSW')
    call glob_allocate(sbufNE,0,stride-1,0,3*bloc_phi-1,'sbufNE')
    call glob_allocate(sbufSE,0,stride-1,0,3*bloc_phi-1,'sbufSE')
    call glob_allocate(fillderiv,1,12,0,3*bloc_phi-1,'fillderiv')
      
!baoter
    if (transpose4D) then
      call glob_allocate(tmprhs2,-bufsize,Nr+bufsize, &
        -bufsize,Ntheta-1+bufsize,0,3*bloc_phi-1,'tmprhs2')
    else
      call glob_allocate(tmprhs2,istart_buf,iend_buf, &
        jstart_buf,jend_buf,0,6*bloc_phi-1,'tmprhs2')
    end if
      
    stride = Nr+1 ! Because it starts at 0 and ends at Nr
    call MPI_TYPE_VECTOR( dom_theta, dom_r, stride,&
      MPI_REAL8, mpitype_locbloc, ierr )
    call MPI_TYPE_COMMIT( mpitype_locbloc, ierr )
    
    stride = dom_r 
    call MPI_TYPE_VECTOR( dom_theta, dom_r, stride,&
      MPI_REAL8, mpitype_globbloc, ierr )
    call MPI_TYPE_COMMIT( mpitype_globbloc, ierr )
      
    stride = dom_r 
    call MPI_TYPE_VECTOR( (Nvpar+1), dom_r, stride,&
      MPI_REAL8, mpitype_globbloc_rvpar, ierr )
    call MPI_TYPE_COMMIT( mpitype_globbloc_rvpar, ierr )
!eaoter
      
    !*****************************************************
    ! FOR MOMENTS
    !*****************************************************
    call glob_allocate(moments_mapphi,0,Nphi-1,'moments_mapphi')
    dom_mapphi = ceiling(1.0*Nphi/Nbproc_tot)
    if ((dom_mapphi.gt.1).and.(dom_mapphi*Nbproc_tot.ne.Nphi)) &
      then
      print*, 'dom_mapphi = ',dom_mapphi
      print*, 'dom_mapphi*Nbproc_tot = ',dom_mapphi*Nbproc_tot
      STOP "condition not checked : dom_mapphi*Nbproc_tot .ne. Nphi"
    endif
    if (.not.memory_test) then
      iphi      = 0
      iphistart = -1
      do iproc = 0, Nbproc_tot
        do i = 1, dom_mapphi
          if (iphi .le. Nphi-1) then
            moments_mapphi(iphi) = iproc
            if ((iphistart.eq.-1).and.(iproc.eq.pglobal_id)) then
              iphistart = iphi
            endif
            iphi = iphi+1
          endif
        enddo
      enddo
    end if
    call glob_allocate(moments_reqs,0,Nphi-1,'moments_reqs')
    call glob_allocate(moments_reqr,0, &
      Nbproc_loc*(Nmu+1)*dom_mapphi,'moments_reqr')
    call glob_allocate(bcast_reqs,0,dom_mapphi*Nbproc_tot-1, &
      'bcast_reqs')
    call glob_allocate(bcast_reqr,0,Nphi-1,'bcast_reqr')
    moments_status_size = max(Nphi-1,Nbproc_loc*(Nmu+1)*dom_mapphi)
    call glob_allocate(moments_status,1,MPI_STATUS_SIZE,0, &
      moments_status_size,'moments_status')    
   call ppbarrier()
  end subroutine ppinit_comm
      
  !---------------------------------------------------------------
  ! Transposition of a 4D structure (r,theta,phi,vpar) 
  !  distributed in r and theta into a 4D structure distributed
  !  in phi and vpar, i.e:
  !   fval(istart:iend,jstart:jend,0:Nphi,0:Nvpar)
  !   -> fval_transp(0:Nr,0:Ntheta,kstart:kend,lstart:lend)
  !---------------------------------------------------------------
  subroutine pptranspose_forward(fval)
    real*8, dimension (istart_buf:iend_buf,&
      jstart_buf:jend_buf,0:Nphi,0:Nvpar), intent(in) :: fval
      
    integer :: ir, itheta, iphi, ivpar, zz, k, q, tag, ierr
    integer :: base, dep
    integer :: dististart, distiend, distjstart, distjend
    integer :: distkstart, distkend, distlstart, distlend
    integer :: distknumber, distlnumber
    integer :: offset, sendsize, send_id, recv_id
!R3 #include "r3_info.h" !R3
      
!R3 call r3_info_begin (r3_info_index_1, 'Comm_transp_forw') !R3
    sendsize = dom_vpar * dom_phi * dom_theta * dom_r
    do offset = 0,Nbproc_loc-1
      send_id     = mod(plocal_id + offset, Nbproc_loc)
      base        = (send_id/Nbproc_phi)
      dep         = mod(send_id,Nbproc_phi)
      distkstart  = dep *  dom_phi
      distkend    = distkstart + dom_phi - 1
      distknumber = distkend - distkstart + 1
      distlstart  = (base * (Nvpar+1)) / Nbproc_vpar
      distlend    = ((base+1) * (Nvpar+1)) / Nbproc_vpar - 1
      distlnumber = distlend - distlstart + 1
      do ivpar = 0,distlnumber-1
        do iphi = 0,distknumber-1
          do itheta = 0,dom_theta-1
            do ir = 0,dom_r-1
              f4D_send(ir,itheta,iphi,ivpar) = &
                fval(istart+ir,jstart+itheta, &
                distkstart+iphi,distlstart+ivpar)
            end do
          end do
        end do
      end do
      recv_id = mod(plocal_id - offset + Nbproc_loc, Nbproc_loc)
      tag     = offset + 10
      call MPI_IRECV(f4D_recv,sendsize,MPI_REAL8,recv_id, &
        tag,mpi_comm_mu,reqr(1),ierr)
      call MPI_SEND(f4D_send,sendsize,MPI_REAL8,send_id, &
        tag,mpi_comm_mu,ierr)
      call MPI_WAIT(reqr(1),status,ierr)
      
      base       = (recv_id/Nbproc_r)
      dep        = mod(recv_id,Nbproc_r)              
      dististart = dep *  dom_r
      distiend   = dististart + dom_r - 1
      distjstart = base * dom_theta
      distjend   = distjstart + dom_theta - 1
      do ivpar = 0,lnumber-1
        do iphi = 0,knumber-1
          do itheta = 0,dom_theta-1
            do ir = 0,dom_r-1
              f4D_transp(dististart+ir,distjstart+itheta, &
                kstart+iphi,lstart+ivpar) = &
                f4D_recv(ir,itheta,iphi,ivpar) 
            end do
          end do
        end do
      end do
    end do
!R3 call r3_info_end(r3_info_index_1) !R3
  end subroutine pptranspose_forward
      
  !---------------------------------------------------------------
  ! Transposition of a 4D structure (r,theta,phi,vpar) 
  !  distributed in phi and vpar into a 4D structure distributed
  !  in r and theta, i.e:
  !   fval_transp(0:Nr,0:Ntheta,kstart:kend,lstart:lend)
  !   -> fval(istart:iend,jstart:jend,0:Nphi,0:Nvpar)
  !
  ! ( Rk : Opposite transposition of the one performed by 
  !   pptranspose_backward )
  !---------------------------------------------------------------
  subroutine pptranspose_backward(fval)
    real*8, dimension (istart_buf:iend_buf, &
      jstart_buf:jend_buf,0:Nphi,0:Nvpar), intent(inout) :: fval
      
    integer :: ir, itheta, iphi, ivpar, zz, k, q, tag, ierr
    integer :: base, dep
    integer :: dististart, distiend, distjstart, distjend
    integer :: distkstart, distkend, distlstart, distlend
    integer :: distknumber, distlnumber
    integer :: offset, sendsize, send_id, recv_id
!R3 #include "r3_info.h" !R3
      
!R3 call r3_info_begin (r3_info_index_1, 'Comm_transp_back') !R3
    sendsize = dom_vpar * dom_phi * dom_theta * dom_r
    do offset = 0,Nbproc_loc-1
      send_id    = mod(plocal_id + offset, Nbproc_loc)
      base       = (send_id/Nbproc_r)
      dep        = mod(send_id,Nbproc_r)              
      dististart = dep *  dom_r
      distiend   = dististart + dom_r - 1
      distjstart = base * dom_theta
      distjend   = distjstart + dom_theta - 1
      do ivpar = 0,lnumber-1
        do iphi = 0,knumber-1
          do itheta = 0,dom_theta-1
            do ir = 0,dom_r-1
              f4D_send(ir,itheta,iphi,ivpar)  = & 
                f4D_transp(dististart+ir,distjstart+itheta,&
                kstart+iphi,lstart+ivpar) 
      
            end do
          end do
        end do
      end do
      
      recv_id = mod(plocal_id - offset + Nbproc_loc, Nbproc_loc)
      tag     = offset + 10
      call MPI_IRECV(f4D_recv,sendsize,MPI_REAL8,recv_id, &
        tag,mpi_comm_mu,reqr(1),ierr)
      call MPI_SEND(f4D_send,sendsize,MPI_REAL8,send_id, &
        tag,mpi_comm_mu,ierr)
      call MPI_WAIT(reqr(1),status,ierr)
      
      base        = (recv_id/Nbproc_phi)
      dep         = mod(recv_id,Nbproc_phi)
      distkstart  = dep * dom_phi
      distkend    = distkstart + dom_phi - 1
      distknumber = distkend - distkstart + 1
      distlstart  = (base * (Nvpar+1)) / Nbproc_vpar
      distlend    = ((base+1) * (Nvpar+1)) / Nbproc_vpar - 1
      distlnumber = distlend - distlstart + 1
      do ivpar = 0,distlnumber-1
        do iphi = 0,distknumber-1
          do itheta = 0,dom_theta-1
            do ir = 0,dom_r-1
              fval(istart+ir,jstart+itheta, &
                distkstart+iphi,distlstart+ivpar) = &
                f4D_recv(ir,itheta,iphi,ivpar)
            end do
          end do
        end do
      end do
    end do
!R3 call r3_info_end(r3_info_index_1) !R3
  end subroutine pptranspose_backward
      
  !---------------------------------------------------------------
  ! Transposition of a 5D structure (r,theta,phi,vpar) 
  !  distributed in (r, theta, mu) into a 5D structure distributed
  !  in (r, theta, phi), i.e:
  !   fval(istart:iend,jstart:jend,0:Nphi,0:Nvpar,mu_id)
  !   -> fval_transp(istart:iend,jstart:jend,
  !                 cstart:cend,0:Nvpar,0:Nmu)
  !---------------------------------------------------------------
  subroutine pptranspvel_forward(fval)
    real*8, dimension (istart_buf:iend_buf,&
      jstart_buf:jend_buf,0:Nphi,0:Nvpar), intent(in) :: fval
      
    integer :: ir, itheta, iphi, ivpar, zz, k, q, tag, ierr
    integer :: indx
    integer :: base, dep
    integer :: dististart, distiend, distjstart, distjend
    integer :: offset, sendsize, send_id, recv_id
    integer :: iproc, mudist
!R3 #include "r3_info.h" !R3
      
!R3 call r3_info_begin (r3_info_index_1, 'Comm_transpvel_forw') !R3
    sendsize = dom_bisphi * dom_theta * dom_r * (Nvpar+1)
      
    do offset = 0, Nmu
      mudist      = mod(mu_id + offset, Nmu+1)
      send_id     = plocal_id + mudist * Nbproc_loc
      base        = min(Nphi-1, mudist * dom_bisphi)
      if (base .eq. mudist * dom_bisphi) then
        do ivpar = 0, Nvpar
          do iphi = 0, min(Nphi-1-base,dom_bisphi-1)
            do itheta = 0, dom_theta-1
              do ir = 0, dom_r-1
                indx = ir + dom_r * &
                  (itheta+dom_theta*(iphi+dom_bisphi*ivpar))
                f1D_send(indx) = &
                  fval(istart+ir,jstart+itheta, &
                  base+iphi,ivpar)
              end do
            end do
          end do
        end do
      end if
      mudist      = mod(Nmu+1 + mu_id - offset, Nmu+1)
      recv_id     = plocal_id + mudist * Nbproc_loc
      tag     = offset + 10
      call MPI_IRECV(f1D_recv,sendsize,MPI_REAL8,recv_id, &
        tag,mpi_comm_world,reqr(1),ierr)
      call MPI_SEND(f1D_send,sendsize,MPI_REAL8,send_id, &
        tag,mpi_comm_world,ierr)
      call MPI_WAIT(reqr(1),status,ierr)
      
      if (collstart .eq. (mu_id * dom_bisphi)) then
        do ivpar = 0, Nvpar
          do iphi = 0, collend-collstart
            do itheta = 0, dom_theta-1
              do ir = 0, dom_r-1
                indx = ir + dom_r * &
                  (itheta+dom_theta*(iphi+dom_bisphi*ivpar))
                f5D_transp(istart+ir,jstart+itheta, &
                  collstart+iphi,ivpar,mudist) = &
                  f1D_recv(indx) 
              end do
            end do
          end do
        end do
      end if
    end do
!R3 call r3_info_end(r3_info_index_1) !R3
  end subroutine pptranspvel_forward
      
  !---------------------------------------------------------------
  ! Transposition of a 5D structure (r,theta,phi,vpar) 
  !  distributed in (r, theta, phi) into a 5D structure distributed
  !  in (r, theta, mu), i.e:
  !  fval_transp(istart:iend,jstart:jend,cstart:cend,0:Nvpar,0:Nmu)
  !    -> fval(istart:iend,jstart:jend,0:Nphi,0:Nvpar,mu_id)
  !---------------------------------------------------------------
  subroutine pptranspvel_backward(fval)
    real*8, dimension (istart_buf:iend_buf,&
      jstart_buf:jend_buf,0:Nphi,0:Nvpar), intent(out) :: fval
      
    integer :: ir, itheta, iphi, ivpar, zz, k, q, tag, ierr
    integer :: indx
    integer :: base, dep
    integer :: dististart, distiend, distjstart, distjend
    integer :: offset, sendsize, send_id, recv_id
    integer :: iproc, mudist
!R3 #include "r3_info.h" !R3
      
!R3 call r3_info_begin (r3_info_index_1, 'Comm_transpvel_back') !R3
    sendsize = dom_bisphi * dom_theta * dom_r * (Nvpar+1)
      
    do offset = 0, Nmu
      mudist      = mod(mu_id + offset, Nmu+1)
      send_id     = plocal_id + mudist * Nbproc_loc
      
      if (collstart .eq. mu_id * dom_bisphi) then
        do ivpar = 0, Nvpar
          do iphi = 0, collend-collstart
            do itheta = 0, dom_theta-1
              do ir = 0, dom_r-1
                indx = ir + dom_r * &
                  (itheta+dom_theta*(iphi+dom_bisphi*ivpar))
                f1D_send(indx) = &
                  f5D_transp(istart+ir,jstart+itheta, &
                  collstart+iphi,ivpar,mudist) 
              end do
            end do
          end do
        end do
      end if
      mudist  = mod(Nmu+1 + mu_id - offset, Nmu+1)
      recv_id = plocal_id + mudist * Nbproc_loc
      tag     = offset + 10
      call MPI_IRECV(f1D_recv,sendsize,MPI_REAL8,recv_id, &
        tag,mpi_comm_world,reqr(1),ierr)
      call MPI_SEND(f1D_send,sendsize,MPI_REAL8,send_id, &
        tag,mpi_comm_world,ierr)
      call MPI_WAIT(reqr(1),status,ierr)
      
      base        = min(Nphi-1, mudist * dom_bisphi)
      if (base .eq. mudist * dom_bisphi) then
        do ivpar = 0, Nvpar
          do iphi = 0, min(Nphi-1-base,dom_bisphi-1)
            do itheta = 0, dom_theta-1
              do ir = 0, dom_r-1
                indx = ir + dom_r * &
                  (itheta+dom_theta*(iphi+dom_bisphi*ivpar))
                fval(istart+ir,jstart+itheta, &
                  base+iphi,ivpar) = f1D_recv(indx)
              end do
            end do
          end do
        end do
      end if
    end do
!R3 call r3_info_end (r3_info_index_1) !R3
  end subroutine pptranspvel_backward
      
  !-----------------------------------------------------
  ! Delete all the arrays associated to the parallel
  !  buffers
  !-----------------------------------------------------
  subroutine ppdeallocate
    call glob_deallocate(rbufWW)
    call glob_deallocate(rbufEE)
    call glob_deallocate(sbufWW)
    call glob_deallocate(sbufEE)
    call glob_deallocate(rbufNN)
    call glob_deallocate(rbufSS)
    call glob_deallocate(sbufNN)
    call glob_deallocate(sbufSS)
    call glob_deallocate(rbufNW)
    call glob_deallocate(rbufSW)
    call glob_deallocate(rbufNE)
    call glob_deallocate(rbufSE)
    call glob_deallocate(sbufNW)
    call glob_deallocate(sbufSW)
    call glob_deallocate(sbufNE)
    call glob_deallocate(sbufSE)
!baoter
    call glob_deallocate(tmprhs2)
    if (transpose4D) then
      call glob_deallocate(f4D_transp)
      call glob_deallocate(f4D_send)
      call glob_deallocate(f4D_recv)
    end if
    if (transp_velocity) then
      call glob_deallocate(f5D_transp)
      call glob_deallocate(f1D_send)
      call glob_deallocate(f1D_recv)
    end if
!eaoter
  end subroutine ppdeallocate
      
  !-----------------------------------------------------
  ! send all datas on the local processors 
  !  (using MPI_ISEND directive)
  !-----------------------------------------------------
  subroutine comm_isend 
    integer :: ierr
      
    call MPI_ISEND(sbufWW,1,mpitype_WE,WW,10,mpi_comm_mu, &
      reqs(1),ierr)
    call MPI_ISEND(sbufEE,1,mpitype_WE,EE,11,mpi_comm_mu, &
      reqs(2),ierr)
    call MPI_ISEND(sbufNW,1,mpitype_corner,NW,20,mpi_comm_mu, &
      reqs(3),ierr)
    call MPI_ISEND(sbufNE,1,mpitype_corner,NE,21,mpi_comm_mu, &
      reqs(4),ierr)
    call MPI_ISEND(sbufSE,1,mpitype_corner,SE,23,mpi_comm_mu, &
      reqs(5),ierr)
    call MPI_ISEND(sbufSW,1,mpitype_corner,SW,22,mpi_comm_mu, &
      reqs(6),ierr)
    call MPI_ISEND(sbufSS,1,mpitype_NS,SS,13,mpi_comm_mu, &
      reqs(7),ierr)
    call MPI_ISEND(sbufNN,1,mpitype_NS,NN,12,mpi_comm_mu, &
      reqs(8),ierr)
  end subroutine comm_isend
      
  !-----------------------------------------------------
  ! send all datas on the local processors 
  !  (using MPI_SSEND directive)
  !-----------------------------------------------------
  subroutine comm_ssend 
    integer :: ierr
      
    call MPI_SSEND(sbufWW,1,mpitype_WE,WW,10,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufEE,1,mpitype_WE,EE,11,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufNW,1,mpitype_corner,NW,20,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufNE,1,mpitype_corner,NE,21,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufSE,1,mpitype_corner,SE,23,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufSW,1,mpitype_corner,SW,22,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufSS,1,mpitype_NS,SS,13,mpi_comm_mu,ierr)
    call MPI_SSEND(sbufNN,1,mpitype_NS,NN,12,mpi_comm_mu,ierr)
  end subroutine comm_ssend
      
  !-----------------------------------------------------------------
  ! receive all datas on the local processor sended by other procs.
  !-----------------------------------------------------------------
  subroutine comm_irecv 
    integer :: ierr
      
   call MPI_IRECV(rbufEE,1,mpitype_WE,EE,10,mpi_comm_mu, &
     reqr(1),ierr)
   call MPI_IRECV(rbufWW,1,mpitype_WE,WW,11,mpi_comm_mu, &
     reqr(2),ierr)
   call MPI_IRECV(rbufSE,1,mpitype_corner,SE,20,mpi_comm_mu, &
     reqr(3),ierr)
   call MPI_IRECV(rbufSW,1,mpitype_corner,SW,21,mpi_comm_mu, &
     reqr(4),ierr)
   call MPI_IRECV(rbufNW,1,mpitype_corner,NW,23,mpi_comm_mu, &
     reqr(5),ierr)
   call MPI_IRECV(rbufNE,1,mpitype_corner,NE,22,mpi_comm_mu, &
     reqr(6),ierr)
   call MPI_IRECV(rbufNN,1,mpitype_NS,NN,13,mpi_comm_mu, &
     reqr(7),ierr)
   call MPI_IRECV(rbufSS,1,mpitype_NS,SS,12,mpi_comm_mu, &
     reqr(8),ierr)
 end subroutine comm_irecv
      
  !-----------------------------------------------------
  ! used for test
  !-----------------------------------------------------
  subroutine comm_wait_recv
    integer :: ierr
      
    call MPI_WAITALL(8,reqr,status, ierr)
  end subroutine comm_wait_recv
      
  !-----------------------------------------------------
  ! used for test
  !-----------------------------------------------------
  subroutine comm_wait_send
   integer :: ierr
      
   call MPI_WAITALL(8,reqs,status, ierr)
  end subroutine comm_wait_send
      
  !-----------------------------------------------------
  ! used for test 
  !-----------------------------------------------------
  subroutine comm_test
    integer :: ierr
    logical flag
      
    call mpi_testall(8,reqs,flag,status, ierr)
    call mpi_testall(8,reqr,flag,status, ierr)
  end subroutine comm_test
      
  !-----------------------------------------------------
  ! exit 
  !-----------------------------------------------------
  subroutine ppexit
    integer :: ierr
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
  end subroutine ppexit
  
  
  !-----------------------------------------------------
  ! MPI_BARRIER
  !-----------------------------------------------------  
  subroutine ppbarrier
    integer :: ierr
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine ppbarrier
      
  !-----------------------------------------------------
  ! MPI_BARRIER
  !-----------------------------------------------------  
  subroutine ppbarrier_timer
    integer :: ierr
    
#ifdef TIMER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
  end subroutine ppbarrier_timer
      
  !-----------------------------------------------------
  ! writing of the character string 'string'  
  !-----------------------------------------------------  
  subroutine ppwrite(string,proc_num)
    character(len=*), intent(in) :: string
    integer         , intent(in) :: proc_num
      
    integer :: ierr
      
    if (proc_num.eq.0) then
      write(6,*) string(1:len(string))
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)    
  end subroutine ppwrite
  
  
  !-----------------------------------------------------
  ! distribution of a direction according to the 
  !  processor number
  !-----------------------------------------------------  
  subroutine ppdist1d(begin_tot, nb_tot, begin_loc, end_loc, nb_loc)
    use globals, only : pglobal_id, Nbproc_tot
    integer, intent(in)  :: begin_tot, nb_tot
    integer, intent(out) :: begin_loc, end_loc, nb_loc
  
    integer :: naver, rem
      
    naver     = nb_tot/Nbproc_tot
    rem       = modulo(nb_tot,Nbproc_tot)
    begin_loc = begin_tot + min(rem,pglobal_id) + pglobal_id*naver
    nb_loc    = naver
    if ( pglobal_id.lt.rem ) nb_loc = nb_loc+1
    end_loc   = begin_loc + nb_loc - 1
  end subroutine ppdist1d
      
  !-----------------------------------------------------
  ! definition of the generic MPI directive MPI_GATHER
  !  directly with the simple MPI directive MPI_ISEND
  !   and MPI_RECV for a 2D structure
  !-----------------------------------------------------  
  subroutine comm_gather2D(buf_source,nbelts,start,stride, &
    cible,buf_cible,comm) 
    real(RKIND), dimension(0:,0:), intent(in)  :: buf_source
    integer                      , intent(in)  :: nbelts
    integer                      , intent(in)  :: start
    integer                      , intent(in)  :: stride
    integer                      , intent(in)  :: cible
    real(RKIND), dimension(0:,0:), intent(out) :: buf_cible
    integer                      , intent(in)  :: comm
      
    integer :: gath_reqs, tag, nbreqr, ierr, pid, curs, nbp, mypid
    integer, &
      dimension(0:Nbproc_tot-1)                   :: gath_reqr
    integer, &
      dimension(1:MPI_STATUS_SIZE,0:Nbproc_tot-1) :: gath_status
      
    call MPI_COMM_SIZE(comm,nbp,ierr)
    call MPI_COMM_RANK(comm,mypid,ierr)
    tag = 4001 
    if (mypid .eq. cible) then
      nbreqr = 0
      do pid = 0,nbp-1
        curs = pid*stride
        call MPI_IRECV(buf_cible(0,curs),nbelts,MPI_REAL8,&
          pid,tag,comm,gath_reqr(pid),ierr)
      enddo
    endif
    call MPI_ISEND(buf_source(0,start),nbelts,MPI_REAL8, &
      cible,tag,comm, gath_reqs, ierr)
      
    if (mypid .eq. cible) &
      call MPI_WAITALL(nbp,gath_reqr(0),gath_status(1,0),ierr)
    call MPI_WAIT(gath_reqs,gath_status(1,0),ierr)
  end subroutine comm_gather2D
      
  !-----------------------------------------------------
  ! definition of the generic MPI directive MPI_GATHER
  !  directly with the simple MPI directive MPI_ISEND
  !   and MPI_RECV for a 3D structure
  !-----------------------------------------------------  
  subroutine comm_gather3D(buf_source,nbelts,start,stride, &
    cible,buf_cible,comm) 
    real(RKIND), dimension(:,:,:)   , pointer     :: buf_source
    integer                         , intent(in)  :: nbelts
    integer                         , intent(in)  :: start
    integer                         , intent(in)  :: stride
    integer                         , intent(in)  :: cible
    real(RKIND), dimension(0:,0:,0:), intent(out) :: buf_cible
    integer                         , intent(in)  :: comm
      
    integer :: gath_reqs, tag, nbreqr, ierr, pid, curs, nbp, mypid
    integer, &
      dimension(0:Nbproc_tot-1)                   :: gath_reqr
    integer, &
      dimension(1:MPI_STATUS_SIZE,0:Nbproc_tot-1) :: gath_status
      
    call MPI_COMM_SIZE(comm,nbp,ierr)
    call MPI_COMM_RANK(comm,mypid,ierr)
    tag = 2001 
    if (mypid .eq. cible) then
      nbreqr = 0
      do pid = 0,nbp-1
        curs = pid*stride
        call MPI_IRECV(buf_cible(0,0,curs),nbelts,MPI_REAL8,&
          pid,tag,comm,gath_reqr(pid),ierr)
      enddo
    endif
    call MPI_ISEND(buf_source(0,0,start),nbelts,MPI_REAL8, &
      cible,tag,comm,gath_reqs,ierr)
      
    if (mypid .eq. cible) &
      call MPI_WAITALL(nbp,gath_reqr(0),gath_status(1,0),ierr)
    call MPI_WAIT(gath_reqs,gath_status(1,0),ierr)
  end subroutine comm_gather3D
end module MPIutils_module
