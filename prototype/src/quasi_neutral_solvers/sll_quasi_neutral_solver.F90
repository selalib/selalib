  !---------------------------------------------------------------------------
  !  Module for solving the quasi-neutral equation in polar coordinates using 
  !  arbitrary order spline finite elements
  !  Eric Sonnendrucker 2011-05-05
  !---------------------------------------------------------------------------
  
  ! Some 'fake' input information that eventually must come from somewhere 
  ! else: These are the dimensions of the 'processor' mesh, i.e., the number 
  ! of processors that contain the decomposed domain in each direction.
  ! r, z, t = r, z, theta. This is a 'logical' organization of the processors
  ! and
  
#define NUMP_X1    1
#define NUMP_X2    2
#define NUMP_X3    1
  
  module sll_quasi_neutral_solver
    use sll_collective
#include "sll_memory.h"
#include "sll_working_precision.h"
!#include "sll_mesh_types.h"
#include "sll_remap.h"
    use sll_constants
    implicit none
    
    type quasi_neutral_plan
       integer   :: spline_degree  ! degree of splines used by F.E. solver
       integer   :: num_pts_r
       integer   :: num_pts_t   
       integer   :: num_pts_z
       real(f64) :: dr, dtheta     ! cell size in r and theta
       real(f64) :: rmin
       sll_real64, dimension(:,:), allocatable :: F ! RHS in matrix form
       sll_real64, dimension(:,:), allocatable :: U ! solution
       real(f64), dimension(:), pointer :: knotsr, knotsth  ! knot vectors for 
                                                            ! splines in r and 
                                                            ! theta direction
       real(f64), dimension(:,:), pointer :: Kar, Mar, Mcr  ! symmetric banded 
                                                            ! matrices in r 
                                                            ! direction
       real(f64), dimension(:), pointer ::  Kth, Mth        ! coefficients of 
                                                            ! symmetric,
                                                            ! circulant 
                                                            ! matrices in theta 
                                                            ! direction
       real(f64), dimension(:), pointer :: xgauss, wgauss   ! Gauss points and 
                                                            ! weights
       ! The data is distributed in such a way that sequential 
       ! algorithms can be used. The redistribution is made by the 
       ! remapper utility. The QN plan includes the required information,
       ! in terms of remap plans, to transform the data from one layout
       ! to another.
       type(remap_plan_3D_t), pointer :: seq_theta_to_seq_r
       type(remap_plan_3D_t), pointer :: seq_r_to_seq_theta
    end type quasi_neutral_plan
    
    
  contains
    
    
    ! FIXME: The name of this function should express the fact that this is the 
    ! plan for a 3D mesh...
    ! Here we explore the idea of having all the information about the data
    ! layouts and remap plans internally, without any exposure in the interface.
    function new_qn_plan( &
      spline_degree, &
      num_pts_r, &
      num_pts_theta, &
      r_min, &
      r_max, &
      theta_min, &
      theta_max )
 
      type(quasi_neutral_plan), pointer  :: new_qn_plan
      sll_int32, intent(in) :: spline_degree   ! degree of spline basis funcs
      sll_int32, intent(in) :: num_pts_r, num_pts_theta
      sll_real64, intent(in) :: r_min, r_max, theta_min, theta_max
      ! point dimensions in r, theta and z
      integer   :: num_pts_r
      integer   :: num_pts_theta   
      real(f64) :: dr, dtheta      ! cell size in r and theta
      real(f64) :: rprof           ! r profile for Kar and Mar matrices
      integer   :: i,j,k,ig,j1,j2  ! loop indices
      integer   :: ierr            ! for error codes
      real(f64) :: xg              ! gauss point
      real(f64) :: ri,ti           ! grid point, knot
      ! work array for de Boor routines:
      real(f64), dimension(spline_degree+1,spline_degree+1) :: work  
      ! values of splines and derivatives at point x
      real(f64), dimension(spline_degree+1,2) :: dbiatx 
      ! Variables related to the parallel data layout.
      type(layout_3D_t), pointer :: sequential_r
      type(layout_3D_t), pointer :: sequential_theta
      ! The following determines the data layout.
      sll_int32                  :: nproc_t   ! number of processors in theta
      sll_int32                  :: nproc_r   ! number of processors in r
      sll_int32                  :: nproc_z   ! number of processors in z
      sll_int32                  :: node
      sll_int32                  :: local_r
      sll_int32                  :: local_t
      sll_int32                  :: local_z
      sll_int32                  :: i_min, i_max
      sll_int32                  :: j_min, j_max
      sll_int32                  :: k_min, k_max
      sequential_theta => new_layout_3D( sll_world_collective )
      sequential_r     => new_layout_3D( sll_world_collective )
      
      SLL_ALLOCATE( new_qn_plan, ierr )
      ! set scalars in quasi_neutral_plan object
!!$      rmin                      = GET_MESH_RMIN(rtz_mesh)
!!$      npr                       = GET_MESH_NCR(rtz_mesh)+1   ! num of points
!!$      npt                       = GET_MESH_NCTHETA(rtz_mesh) ! periodic
!!$      npz                       = GET_MESH_NCZ(rtz_mesh)+1   ! num of points
!!$      dr                        = GET_MESH_DELTA_R(rtz_mesh)
!!$      dtheta                    = GET_MESH_DELTA_THETA(rtz_mesh)
      new_qn_plan%spline_degree = spline_degree
      new_qn_plan%num_pts_r     = num_pts_r
      new_qn_plan%num_pts_t     = num_pts_theta
      new_qn_plan%dr            = (r_max-r_min)/(num_pts_r-1)
      new_qn_plan%dtheta        = (theta_max-theta_min)/(num_pts_theta-1)
      new_qn_plan%r_min         = r_min
      ! Start with a layout for sequential operations in theta
      nproc_t = NUMP_X1
      nproc_r = NUMP_X2
      nproc_z = NUMP_X3
      ! Initialize the layout
      local_t = npt/nproc_t
      local_r = npr/nproc_r
      local_z = npz/nproc_z
      write (*,'(a, i8)') 'local_t = ', local_t
      write (*,'(a, i8)') 'local_r = ', local_r
      write (*,'(a, i8)') 'local_z = ', local_z
      write (*,'(a, i8)') 'npt = ', npt

      ! Loop over the "processor mesh" to assign the min/max limits in the 
      ! layout. This tends to get complicated if the number of points is not 
      ! equally partitioned amond the available processes. Hence the 
      ! disgusting "if's". This should be improved.
      do k=0, nproc_z-1
         do j=0, nproc_r-1
            do i=0, nproc_t-1
               node = i+nproc_t*(j+nproc_r*k)  ! global node index
               ! theta values are contiguous, thus they are the first coord.
               i_min = i*local_t + 1
               if(i .eq. nproc_t-1) then
                  i_max = npt
               else
                  i_max = i*local_t + local_t
               end if
               ! r values follow
               j_min = j*local_r + 1
               if(j .eq. nproc_r-1) then
                  j_max = npr
               else
                  j_max = j*local_r + local_r
               end if
               ! z values finally
               k_min = k*local_z + 1
               if(k .eq. nproc_z-1) then
                  k_max = npz
               else
                  k_max = k*local_z + local_z
               end if
               call set_layout_i_min(sequential_theta, node, i_min )
               call set_layout_i_max(sequential_theta, node, i_max )
               call set_layout_j_min(sequential_theta, node, j_min )
               call set_layout_j_max(sequential_theta, node, j_max )
               call set_layout_k_min(sequential_theta, node, k_min )
               call set_layout_k_max(sequential_theta, node, k_max )
            end do
         end do
      end do
      ! Repeat the same process but with a layout that allows 
      ! sequential operations in 'r' (i.e.: the number of processors in
      ! 'r' is equal to 1. The rest stays the same, we can think of 
      ! abstracting this into a subroutine it seems...
      nproc_t = NUMP_X2
      nproc_r = NUMP_X1
      nproc_z = NUMP_X3
      ! Initialize the layout
      local_t = npt/nproc_t
      local_r = npr/nproc_r
      local_z = npz/nproc_z

      ! Loop over the "processor mesh" to assign the min/max limits in the 
      ! layout. This tends to get complicated if the number of points is not 
      ! equally partitioned amond the available processes. Hence the 
      ! disgusting "if's". This should be improved.
      do k=0, nproc_z-1
         do j=0, nproc_r-1
            do i=0, nproc_t-1
               node = i+nproc_t*(j+nproc_r*k)  ! global node index
               ! theta values are contiguous, thus they are the first coord.
               i_min = i*local_t + 1
               if(i .eq. nproc_t-1) then
                  i_max = npt
               else
                  i_max = i*local_t + local_t
               end if
               ! r values follow
               j_min = j*local_r + 1
               if(j .eq. nproc_r-1) then
                  j_max = npr
               else
                  j_max = j*local_r + local_r
               end if
               ! z values finally
               k_min = k*local_z + 1
               if(k .eq. nproc_z-1) then
                  k_max = npz
               else
                  k_max = k*local_z + local_z
               end if
               call set_layout_i_min(sequential_r, node, i_min )
               call set_layout_i_max(sequential_r, node, i_max )
               call set_layout_j_min(sequential_r, node, j_min )
               call set_layout_j_max(sequential_r, node, j_max )
               call set_layout_k_min(sequential_r, node, k_min )
               call set_layout_k_max(sequential_r, node, k_max )
            end do
         end do
      end do
#if 0
      print *, 'View limits. Sequential THETA'
      call sll_view_lims_3D(sequential_theta)
      print *, 'View limits. Sequential R'
      call sll_view_lims_3D(sequential_r)
#endif
      new_qn_plan%seq_theta_to_seq_r => &
           NEW_REMAPPER_PLAN_3D(sequential_theta, sequential_r, rtz_mesh%data)
      new_qn_plan%seq_r_to_seq_theta => &
           NEW_REMAPPER_PLAN_3D(sequential_r, sequential_theta, rtz_mesh%data)
      ! POR AQUI!!!!
      ! allocate arrays
      SLL_ALLOCATE(new_qn_plan%Kar(spline_degree+1,npr+spline_degree-1),ierr)
      SLL_ALLOCATE(new_qn_plan%Mar(spline_degree+1,npr+spline_degree-1),ierr)
      SLL_ALLOCATE(new_qn_plan%Mcr(spline_degree+1,npr+spline_degree-1),ierr)
      SLL_ALLOCATE(new_qn_plan%Kth(spline_degree+1),ierr)
      SLL_ALLOCATE(new_qn_plan%Mth(spline_degree+1),ierr)
      SLL_ALLOCATE(new_qn_plan%xgauss(spline_degree+1),ierr)
      SLL_ALLOCATE(new_qn_plan%wgauss(spline_degree+1),ierr)
      SLL_ALLOCATE(new_qn_plan%knotsr(npr+2*spline_degree),ierr) 
      SLL_ALLOCATE(new_qn_plan%knotsth(2*spline_degree+2),ierr)
      
      ! set Gauss points and weights
      select case(new_qn_plan%spline_degree)
      case(1) 
         new_qn_plan%xgauss(1) = -1.0_f64/sqrt(3._f64)
         new_qn_plan%xgauss(2) =  1.0_f64/sqrt(3.0_f64)
         new_qn_plan%wgauss(1) =  1.0_f64 
         new_qn_plan%wgauss(2) =  1.0_f64
      case(2)
         new_qn_plan%xgauss(1) = -sqrt(3.0_f64/5.0_f64)
         new_qn_plan%xgauss(2) = 0. 
         new_qn_plan%xgauss(3) = sqrt(3.0_f64/5.0_f64)
         new_qn_plan%wgauss(1) = 5.0_f64/9.0_f64
         new_qn_plan%wgauss(2) = 8.0_f64/9.0_f64
         new_qn_plan%wgauss(3) = new_qn_plan%wgauss(1)
      case(3)
         new_qn_plan%xgauss(4) = &
              sqrt((3.0_f64+2.0_f64*sqrt(6.0_f64/5.0_f64))/7.0_f64)
         new_qn_plan%xgauss(3) = &
              sqrt((3.0_f64-2.0_f64*sqrt(6.0_f64/5.0_f64))/7.0_f64)
         new_qn_plan%xgauss(2) = -new_qn_plan%xgauss(3) 
         new_qn_plan%xgauss(1) = -new_qn_plan%xgauss(4) 
         new_qn_plan%wgauss(1) = (18.0_f64-sqrt(30.0_f64))/36.0_f64
         new_qn_plan%wgauss(2) = (18.0_f64+sqrt(30.0_f64))/36.0_f64  
         new_qn_plan%wgauss(3) = new_qn_plan%wgauss(2)  
         new_qn_plan%wgauss(4) = new_qn_plan%wgauss(1)
      case(4)
         new_qn_plan%xgauss(1) = -0.90617984593866374_f64
         new_qn_plan%xgauss(2) = -0.53846931010568311_f64
         new_qn_plan%xgauss(3) =  0._f64
         new_qn_plan%xgauss(4) =  0.53846931010568311_f64
         new_qn_plan%xgauss(5) =  0.90617984593866374_f64
         new_qn_plan%wgauss(1) =  0.23692688505618875_f64  
         new_qn_plan%wgauss(2) =  0.47862867049936653_f64
         new_qn_plan%wgauss(3) =  0.568888888888888888_f64
         new_qn_plan%wgauss(4) =  0.47862867049936653_f64
         new_qn_plan%wgauss(5) =  0.23692688505618875_f64
      case(5)
         new_qn_plan%xgauss(1) = -0.93246951420315249_f64 
         new_qn_plan%xgauss(2) = -0.66120938646626437_f64
         new_qn_plan%xgauss(3) = -0.23861918608319668_f64
         new_qn_plan%xgauss(4) =  0.23861918608319668_f64  
         new_qn_plan%xgauss(5) =  0.66120938646626437_f64  
         new_qn_plan%xgauss(6) =  0.93246951420315249_f64
         new_qn_plan%wgauss(1) =  0.17132449237917155_f64 
         new_qn_plan%wgauss(2) =  0.36076157304813761_f64
         new_qn_plan%wgauss(3) =  0.46791393457269115_f64
         new_qn_plan%wgauss(4) =  0.46791393457269115_f64 
         new_qn_plan%wgauss(5) =  0.36076157304813761_f64  
         new_qn_plan%wgauss(6) =  0.17132449237917155_f64
      case(6)
         new_qn_plan%xgauss(1) = -0.94910791234275871_f64
         new_qn_plan%xgauss(2) = -0.74153118559939302_f64 
         new_qn_plan%xgauss(3) = -0.4058451513773973_f64
         new_qn_plan%xgauss(4) =  0._f64
         new_qn_plan%xgauss(5) =  0.4058451513773973_f64
         new_qn_plan%xgauss(6) =  0.74153118559939302_f64
         new_qn_plan%xgauss(7) =  0.94910791234275871_f64
         new_qn_plan%wgauss(1) =  0.12948496616886976_f64  
         new_qn_plan%wgauss(2) =  0.27970539148927637_f64 
         new_qn_plan%wgauss(3) =  0.38183005050511931_f64
         new_qn_plan%wgauss(4) =  0.41795918367346824_f64  
         new_qn_plan%wgauss(5) =  0.38183005050511931_f64 
         new_qn_plan%wgauss(6) =  0.27970539148927637_f64
         new_qn_plan%wgauss(7) =  0.12948496616886976_f64
      case(7)
         new_qn_plan%xgauss(1) = -0.96028985649753706_f64 
         new_qn_plan%xgauss(2) = -0.79666647741362673_f64
         new_qn_plan%xgauss(3) = -0.52553240991632888_f64
         new_qn_plan%xgauss(4) = -0.18343464249565036_f64
         new_qn_plan%xgauss(5) =  0.18343464249565036_f64 
         new_qn_plan%xgauss(6) =  0.52553240991632888_f64
         new_qn_plan%xgauss(7) =  0.79666647741362673_f64 
         new_qn_plan%xgauss(8) =  0.96028985649753706_f64
         new_qn_plan%wgauss(1) =  0.10122853629037747_f64 
         new_qn_plan%wgauss(2) =  0.22238103445337337_f64 
         new_qn_plan%wgauss(3) =  0.31370664587788677_f64
         new_qn_plan%wgauss(4) =  0.36268378337836155_f64  
         new_qn_plan%wgauss(5) =  0.36268378337836155_f64 
         new_qn_plan%wgauss(6) =  0.31370664587788677_f64
         new_qn_plan%wgauss(7) =  0.22238103445337337_f64 
         new_qn_plan%wgauss(8) =  0.10122853629037747_f64
      case default
         print*, 'spline degree ', spline_degree, ' not implemented'
         stop
      end select
      
      ! Compute circulant matrices in theta direction
      !----------------------------------------------
      ! set knot vector
      do k=-spline_degree,spline_degree+1
         new_qn_plan%knotsth(k+spline_degree+1)=dtheta*k 
      enddo
      ! initialize matrices
      new_qn_plan%Mth(:) = 0.0
      new_qn_plan%Kth(:) = 0.0
      do ig=1,spline_degree+1   ! sum over Gauss points
         ! rescale Gauss points to be in interval [0,h]
         xg = 0.5_f64*dtheta*(new_qn_plan%xgauss(ig)+1.)  
         !print*, xg, knots(spline_degree+1),knots(spline_degree+2)
         call bsplvd( new_qn_plan%knotsth,    &
              spline_degree+1, &
              xg,              &
              spline_degree+1, &
              work,            &
              dbiatx,          &
              2 )
         ! sum over bsplines in interval [x_i, x_i+1]
         do j = 0,spline_degree
            do i = 0,spline_degree 
               if (i+j < spline_degree+1) then
                  new_qn_plan%Mth(j+1) = new_qn_plan%Mth(j+1) + &
                       0.5_f64*dtheta*new_qn_plan%wgauss(ig)*dbiatx(i+1,1)*dbiatx(i+j+1,1)
                  new_qn_plan%Kth(j+1) = new_qn_plan%Kth(j+1) + &
                       0.5_f64*dtheta*new_qn_plan%wgauss(ig)*dbiatx(i+1,2)*dbiatx(i+j+1,2)
               end if
            end do
         end do
      end do
      
      ! Compute banded matrices in r direction (Dirichlet BC)
      !------------------------------------------------------
      ! set knot vector
      do i = 1, spline_degree + 1
         new_qn_plan%knotsr(i) = rmin
      enddo
      ti = rmin
      do i = spline_degree + 2, Npr + spline_degree
         ti = ti+dr
         new_qn_plan%knotsr(i) = ti
      enddo
      do i = Npr + spline_degree + 1, Npr + 2 * spline_degree
         new_qn_plan%knotsr(i) = ti
      enddo
      
      ! intialize arrays
      new_qn_plan%Mar(:,:) = 0.0
      new_qn_plan%Mcr(:,:) = 0.0
      new_qn_plan%Kar(:,:) = 0.0
      ! loop over cells
      do i = 0, Npr-2
         ri = rmin + i*dr
         ! sum over Gauss points
         do ig = 1, spline_degree+1
            ! rescale Gauss points to be in interval [ri,ri+dr]
            xg = ri + 0.5_f64*dr*(new_qn_plan%xgauss(ig)+1.0_f64)  
            rprof = 1.0_f64   ! r-profile at Gauss point xg 
            call bsplvd( new_qn_plan%knotsr,      &
                 spline_degree+1,  &
                 xg,               &
                 spline_degree+1+i,&
                 work,             &
                 dbiatx,           &
                 2)            
            do j1 = 0, spline_degree
               do j2 = j1, spline_degree
                  new_qn_plan%Mar(j2-j1+1,i+j1+1) = new_qn_plan%Mar(j2-j1+1,i+j1+1) &
                       + 0.5_f64*rprof/xg*dr*new_qn_plan%wgauss(ig)*dbiatx(j1+1,1)*dbiatx(j2+1,1)
                  new_qn_plan%Kar(j2-j1+1,i+j1+1) = new_qn_plan%Kar(j2-j1+1,i+j1+1) &
                       + 0.5_f64*rprof*xg*dr*new_qn_plan%wgauss(ig)*dbiatx(j1+1,2)*dbiatx(j2+1,2)
               enddo
            enddo
         enddo
      enddo
    end function new_qn_plan
    
    subroutine apply_quasi_neutral_solver_plan( plan, F, U )
      ! Solve 2D QN equation on structured grid in polar coordinates with 
      ! Dirichlet BC in r and periodic in theta
      ! Tensor product technique used to bring the solution back to 1D problems
      ! Discrete problem is of the form Kar U Mtheta +Mar U K th +Mcr U Mth = F
      ! where the unknow U and the RHS F are in matrix form and U is multiplied
      ! on left or right by 1D matrices
      ! See Crouseilles, Ratnani, Sonnendrucker 2011
      !------------------------------------------------------------------------
      ! Input :
      type(quasi_neutral_plan), pointer        :: plan
      real(f64), dimension(:,:), intent(inout) :: F !  RHS in matrix form
      ! Output : 
      real(f64), dimension(:,:), intent(out)   :: U ! solution in matrix form
      ! local variables 
      integer :: info, i,k, nband
      ! banded matrix for problem in r :
      real(f64), dimension(plan%spline_degree+1, &
           plan%num_pts_r + plan%spline_degree - 3) :: A 
      ! eigenvalues of circulant mass and stiffness matrix in theta direction 
      real(f64), dimension(plan%num_pts_t/2+1) :: valpm, valpk 
      real(f64), dimension(2*plan%num_pts_t+15) :: wsave  ! help array for FFTPACK
      
      ! intialize dffft
      call dffti(plan%num_pts_t, wsave) 
      
      ! forward FFT of lines of F with FFTPACK
      do i=2,plan%num_pts_r+plan%spline_degree-2
         call dfftf(plan%num_pts_t, F(i,:), wsave) 
      end do
      
      ! compute eigenvalues of circulant matrices
      do k=1,plan%num_pts_t/2+1
         valpm(k) = plan%Mth(1)
         valpk(k) = plan%Kth(1)
         do i=2,plan%spline_degree+1
            valpm(k) = valpm(k) + &
                 2*plan%Mth(i)*cos(2*sll_pi*(i-1)*(k-1)/plan%num_pts_t)
            valpk(k) = valpk(k) + &
                 2*plan%Kth(i)*cos(2*sll_pi*(i-1)*(k-1)/plan%num_pts_t)
         end do
      end do
      
      ! Banded solves in x direction
      nband = plan%spline_degree+1
      U(:,:) = F(:,:)  ! copy rhs into solution
      A(:,:) = valpm(1)*(plan%Kar(:,2:plan%num_pts_r+plan%spline_degree-2) + &
           plan%Mcr(:,2:plan%num_pts_r+plan%spline_degree-2)) + &
           valpk(1)*plan%Mar(:,2:plan%num_pts_r+plan%spline_degree-2)
      ! Cholesky factorisation of A
      call DPBTRF( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, A, nband, info )
      call DPBTRS( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband, &
           U(2:plan%num_pts_r+plan%spline_degree-2,1), plan%num_pts_r+plan%spline_degree-3, &
           info ) ! Solution
      do k = 1, plan%num_pts_t/2-1
         A(:,:) = valpm(k+1)*(plan%Kar(:,2:plan%num_pts_r+plan%spline_degree-2) + &
              plan%Mcr(:,2:plan%num_pts_r+plan%spline_degree-2)) + &
              valpk(k+1)*plan%Mar(:,2:plan%num_pts_r+plan%spline_degree-2)
         ! Cholesky factorisation of A
         call DPBTRF( 'L',plan%num_pts_r+plan%spline_degree-3, nband-1, A, nband, info ) 
         call DPBTRS( 'L',plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband, &
              U(2:plan%num_pts_r+plan%spline_degree-2,2*k), &
              plan%num_pts_r+plan%spline_degree-3, info )
         call DPBTRS( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband,&
              U(2:plan%num_pts_r+plan%spline_degree-2,2*k+1), &
              plan%num_pts_r+plan%spline_degree-3, info )
      end do
      A(:,:) = &
           valpm(plan%num_pts_t/2+1)*(plan%Kar(:,2:plan%num_pts_r+plan%spline_degree-2) + &
           plan%Mcr(:,2:plan%num_pts_r+plan%spline_degree-2)) + &
           valpk(plan%num_pts_t/2+1)*plan%Mar(:,2:plan%num_pts_r+plan%spline_degree-2) 
      ! Cholesky factorisation of A
      call DPBTRF( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, A, nband, info ) 
      call DPBTRS( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband, &
           U(2:plan%num_pts_r+plan%spline_degree-2,plan%num_pts_t), &
           plan%num_pts_r+plan%spline_degree-3, info ) ! Solution
      
      ! backward FFT of lines of U with FFTPACK
      U(1,:) = 0.0 
      U(plan%num_pts_r+plan%spline_degree-1,:) = 0.0
      do i=1,plan%num_pts_r+plan%spline_degree-1
         call dfftb(plan%num_pts_t, U(i,:), wsave) 
      end do
      U(:,:) = U(:,:)/plan%num_pts_t    ! normalization
    end subroutine apply_quasi_neutral_solver_plan

    subroutine solve_quasi_neutral_equation( plan, rho, phi ) !F, U )
      ! Solve 2D QN equation on structured grid in polar coordinates with 
      ! Dirichlet BC in r and periodic in theta
      ! Tensor product technique used to bring the solution back to 1D problems
      ! Discrete problem is of the form Kar U Mtheta +Mar U K th +Mcr U Mth = F
      ! where the unknow U and the RHS F are in matrix form and U is multiplied
      ! on left or right by 1D matrices
      ! See Crouseilles, Ratnani, Sonnendrucker 2011
      !------------------------------------------------------------------------
      ! Input :
      type(quasi_neutral_plan), pointer        :: plan
      sll_real64, dimension(:,:), intent(in)   :: rho
      sll_real64, dimension(:,:), intent(out)  :: phi
      real(f64), dimension(:,:) :: F !  RHS in matrix form
      ! Output : 
      real(f64), dimension(:,:) :: U ! solution in matrix form
      ! local variables 
      integer :: info, i,k, nband
      ! banded matrix for problem in r :
      real(f64), dimension(plan%spline_degree+1, &
           plan%num_pts_r + plan%spline_degree - 3) :: A 
      ! eigenvalues of circulant mass and stiffness matrix in theta direction 
      real(f64), dimension(plan%num_pts_t/2+1) :: valpm, valpk 
      real(f64), dimension(2*plan%num_pts_t+15) :: wsave  ! help array for FFTPACK
      
      ! intialize dffft
      call dffti(plan%num_pts_t, wsave) 
      
      ! forward FFT of lines of F with FFTPACK
      do i=2,plan%num_pts_r+plan%spline_degree-2
         call dfftf(plan%num_pts_t, F(i,:), wsave) 
      end do
      
      ! compute eigenvalues of circulant matrices
      do k=1,plan%num_pts_t/2+1
         valpm(k) = plan%Mth(1)
         valpk(k) = plan%Kth(1)
         do i=2,plan%spline_degree+1
            valpm(k) = valpm(k) + &
                 2*plan%Mth(i)*cos(2*sll_pi*(i-1)*(k-1)/plan%num_pts_t)
            valpk(k) = valpk(k) + &
                 2*plan%Kth(i)*cos(2*sll_pi*(i-1)*(k-1)/plan%num_pts_t)
         end do
      end do
      
      ! Banded solves in x direction
      nband = plan%spline_degree+1
      U(:,:) = F(:,:)  ! copy rhs into solution
      A(:,:) = valpm(1)*(plan%Kar(:,2:plan%num_pts_r+plan%spline_degree-2) + &
           plan%Mcr(:,2:plan%num_pts_r+plan%spline_degree-2)) + &
           valpk(1)*plan%Mar(:,2:plan%num_pts_r+plan%spline_degree-2)
      ! Cholesky factorisation of A
      call DPBTRF( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, A, nband, info )
      call DPBTRS( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband, &
           U(2:plan%num_pts_r+plan%spline_degree-2,1), plan%num_pts_r+plan%spline_degree-3, &
           info ) ! Solution
      do k = 1, plan%num_pts_t/2-1
         A(:,:) = valpm(k+1)*(plan%Kar(:,2:plan%num_pts_r+plan%spline_degree-2) + &
              plan%Mcr(:,2:plan%num_pts_r+plan%spline_degree-2)) + &
              valpk(k+1)*plan%Mar(:,2:plan%num_pts_r+plan%spline_degree-2)
         ! Cholesky factorisation of A
         call DPBTRF( 'L',plan%num_pts_r+plan%spline_degree-3, nband-1, A, nband, info ) 
         call DPBTRS( 'L',plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband, &
              U(2:plan%num_pts_r+plan%spline_degree-2,2*k), &
              plan%num_pts_r+plan%spline_degree-3, info )
         call DPBTRS( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband,&
              U(2:plan%num_pts_r+plan%spline_degree-2,2*k+1), &
              plan%num_pts_r+plan%spline_degree-3, info )
      end do
      A(:,:) = &
           valpm(plan%num_pts_t/2+1)*(plan%Kar(:,2:plan%num_pts_r+plan%spline_degree-2) + &
           plan%Mcr(:,2:plan%num_pts_r+plan%spline_degree-2)) + &
           valpk(plan%num_pts_t/2+1)*plan%Mar(:,2:plan%num_pts_r+plan%spline_degree-2) 
      ! Cholesky factorisation of A
      call DPBTRF( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, A, nband, info ) 
      call DPBTRS( 'L', plan%num_pts_r+plan%spline_degree-3, nband-1, 1, A, nband, &
           U(2:plan%num_pts_r+plan%spline_degree-2,plan%num_pts_t), &
           plan%num_pts_r+plan%spline_degree-3, info ) ! Solution
      
      ! backward FFT of lines of U with FFTPACK
      U(1,:) = 0.0 
      U(plan%num_pts_r+plan%spline_degree-1,:) = 0.0
      do i=1,plan%num_pts_r+plan%spline_degree-1
         call dfftb(plan%num_pts_t, U(i,:), wsave) 
      end do
      U(:,:) = U(:,:)/plan%num_pts_t    ! normalization
    end subroutine solve_quasi_neutral_equation

    
    subroutine evalsplgrid(plan,C,U)
      ! evaluate spline defined by coefficient array C at all grid points
      !------------------------------------------------------------------
      type (quasi_neutral_plan), pointer :: plan
      real(f64), dimension(:,:) :: C  ! spline coefficients
      real(f64), dimension(:,:) :: U  ! Value of function at grid points
      ! local variables
      real(f64), dimension(plan%spline_degree+1) :: biatr,biatth
      real(f64) :: r
      integer   :: i,j,ii,jj
      
      U = 0.0
      ! Values of periodic spline at grid points (need only one call as it is 
      ! uniform)
      call bsplvb(plan%knotsth,plan%spline_degree+1,1,0.0,plan%spline_degree+1,&
           biatth)
      do i = 0, plan%num_pts_r-2
         do ii=0,plan%spline_degree
            r= plan%rmin + i*plan%dr
            call bsplvb(plan%knotsr,plan%spline_degree+1,1,r,&
                 plan%spline_degree+1+i,biatr)
            do j = 0, plan%num_pts_t-1
               do jj=0, plan%spline_degree              
                  U(i+1,j+1) = U(i+1,j+1) + &
                       biatr(ii+1)*biatth(jj+1)*C(i+ii+1,mod(j+jj,plan%num_pts_t)+1)
               enddo
            enddo
         enddo
      enddo
    end subroutine evalsplgrid
    
  end module sll_quasi_neutral_solver
  

