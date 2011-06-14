!------------------------------------------------------------------------------
!  Module for solving the quasi-neutral equation in polar coordinates using 
!  arbitrary order spline finite elements
!  Eric Sonnendrucker 2011-05-05
!------------------------------------------------------------------------------
module sll_quasi_neutral_solver
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

  type quasi_neutral_plan
     integer   :: spline_degree  ! degree of splines used by finite element 
                                 ! solver
     integer   :: nr, ntheta     ! dimensions in r and theta
     real(f64) :: dr, dtheta     ! cell size in r and theta
     real(f64) :: rmin
     real(f64), dimension(:), pointer :: knotsr, knotsth  ! knot vectors for 
                                                          ! splines in r and 
                                                          ! theta direction
     real(f64), dimension(:,:), pointer :: Kar, Mar, Mcr  ! symmetric banded 
                                                          ! matrices in r 
                                                          ! direction
     real(f64), dimension(:), pointer ::  Kth, Mth        ! coefficients of 
                                                          ! symmetric circulant 
                                                          ! matrices in theta 
                                                          ! direction
     real(f64), dimension(:), pointer :: xgauss, wgauss   ! Gauss points and 
                                                          ! weights
  end type quasi_neutral_plan


contains



  function new_qn_plan(spline_degree, rmin, nr, ntheta, dr, dtheta)
    type(quasi_neutral_plan), pointer :: new_qn_plan
    integer   :: spline_degree   ! degree of spline basis functions
    real(f64) :: rmin            ! Minimum value of r
    integer   :: nr, ntheta      ! dimensions in r and theta
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

    SLL_ALLOCATE( new_qn_plan, ierr )
    ! set scalars in quasi_neutral_plan object
    new_qn_plan%spline_degree = spline_degree
    new_qn_plan%nr            = nr
    new_qn_plan%ntheta        = ntheta
    new_qn_plan%dr            = dr
    new_qn_plan%dtheta        = dtheta
    new_qn_plan%rmin          = rmin

    ! allocate arrays
    SLL_ALLOCATE(new_qn_plan%Kar(spline_degree+1,nr+spline_degree-1),ierr)
    SLL_ALLOCATE(new_qn_plan%Mar(spline_degree+1,nr+spline_degree-1),ierr)
    SLL_ALLOCATE(new_qn_plan%Mcr(spline_degree+1,nr+spline_degree-1),ierr)
    SLL_ALLOCATE(new_qn_plan%Kth(spline_degree+1),ierr)
    SLL_ALLOCATE(new_qn_plan%Mth(spline_degree+1),ierr)
    SLL_ALLOCATE(new_qn_plan%xgauss(spline_degree+1),ierr)
    SLL_ALLOCATE(new_qn_plan%wgauss(spline_degree+1),ierr)
    SLL_ALLOCATE(new_qn_plan%knotsr(nr+2*spline_degree),ierr) 
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
    do i = spline_degree + 2, Nr + spline_degree
       ti = ti+dr
       new_qn_plan%knotsr(i) = ti
     enddo
    do i = Nr + spline_degree + 1, Nr + 2 * spline_degree
       new_qn_plan%knotsr(i) = ti
    enddo
    
    ! intialize arrays
    new_qn_plan%Mar(:,:) = 0.0
    new_qn_plan%Mcr(:,:) = 0.0
    new_qn_plan%Kar(:,:) = 0.0
    ! loop over cells
    do i = 0, Nr-2
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
                         plan%nr + plan%spline_degree - 3) :: A 
    ! eigenvalues of circulant mass and stiffness matrix in theta direction 
    real(f64), dimension(plan%ntheta/2+1) :: valpm, valpk 
    real(f64), parameter :: pi = 3.1415926535897931_f64
    real(f64), dimension(2*plan%ntheta+15) :: wsave  ! help array for FFTPACK

    ! intialize dffft
    call dffti(plan%ntheta, wsave) 

    ! forward FFT of lines of F with FFTPACK
    do i=2,plan%nr+plan%spline_degree-2
       call dfftf(plan%ntheta, F(i,:), wsave) 
    end do

    ! compute eigenvalues of circulant matrices
    do k=1,plan%ntheta/2+1
       valpm(k) = plan%Mth(1)
       valpk(k) = plan%Kth(1)
       do i=2,plan%spline_degree+1
          valpm(k) = valpm(k) + 2*plan%Mth(i)*cos(2*pi*(i-1)*(k-1)/plan%ntheta)
          valpk(k) = valpk(k) + 2*plan%Kth(i)*cos(2*pi*(i-1)*(k-1)/plan%ntheta)
       end do
    end do

    ! Banded solves in x direction
    nband = plan%spline_degree+1
    U(:,:) = F(:,:)  ! copy rhs into solution
    A(:,:) = valpm(1)*(plan%Kar(:,2:plan%nr+plan%spline_degree-2) + &
             plan%Mcr(:,2:plan%nr+plan%spline_degree-2)) + &
             valpk(1)*plan%Mar(:,2:plan%nr+plan%spline_degree-2)
   ! Cholesky factorisation of A
    call DPBTRF( 'L', plan%nr+plan%spline_degree-3, nband-1, A, nband, info )
    call DPBTRS( 'L', plan%nr+plan%spline_degree-3, nband-1, 1, A, nband, &
         U(2:plan%nr+plan%spline_degree-2,1), plan%nr+plan%spline_degree-3, &
         info ) ! Solution
    do k = 1, plan%ntheta/2-1
       A(:,:) = valpm(k+1)*(plan%Kar(:,2:plan%nr+plan%spline_degree-2) + &
            plan%Mcr(:,2:plan%nr+plan%spline_degree-2)) + &
            valpk(k+1)*plan%Mar(:,2:plan%nr+plan%spline_degree-2)
       ! Cholesky factorisation of A
       call DPBTRF( 'L',plan%nr+plan%spline_degree-3, nband-1, A, nband, info ) 
       call DPBTRS( 'L',plan%nr+plan%spline_degree-3, nband-1, 1, A, nband, &
            U(2:plan%nr+plan%spline_degree-2,2*k), &
            plan%nr+plan%spline_degree-3, info )
       call DPBTRS( 'L', plan%nr+plan%spline_degree-3, nband-1, 1, A, nband,&
            U(2:plan%nr+plan%spline_degree-2,2*k+1), &
            plan%nr+plan%spline_degree-3, info )
    end do
    A(:,:) = &
         valpm(plan%ntheta/2+1)*(plan%Kar(:,2:plan%nr+plan%spline_degree-2) + &
         plan%Mcr(:,2:plan%nr+plan%spline_degree-2)) + &
         valpk(plan%ntheta/2+1)*plan%Mar(:,2:plan%nr+plan%spline_degree-2) 
    ! Cholesky factorisation of A
    call DPBTRF( 'L', plan%nr+plan%spline_degree-3, nband-1, A, nband, info ) 
    call DPBTRS( 'L', plan%nr+plan%spline_degree-3, nband-1, 1, A, nband, &
         U(2:plan%nr+plan%spline_degree-2,plan%ntheta), &
         plan%nr+plan%spline_degree-3, info ) ! Solution

    ! backward FFT of lines of U with FFTPACK
    U(1,:) = 0.0 
    U(plan%nr+plan%spline_degree-1,:) = 0.0
    do i=1,plan%nr+plan%spline_degree-1
       call dfftb(plan%ntheta, U(i,:), wsave) 
    end do
    U(:,:) = U(:,:)/plan%ntheta    ! normalization
  end subroutine apply_quasi_neutral_solver_plan

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
    do i = 0, plan%nr-2
       do ii=0,plan%spline_degree
          r= plan%rmin + i*plan%dr
          call bsplvb(plan%knotsr,plan%spline_degree+1,1,r,&
               plan%spline_degree+1+i,biatr)
          do j = 0, plan%ntheta-1
             do jj=0, plan%spline_degree              
                U(i+1,j+1) = U(i+1,j+1) + &
                     biatr(ii+1)*biatth(jj+1)*C(i+ii+1,mod(j+jj,plan%ntheta)+1)
             enddo
          enddo
       enddo
    enddo
  end subroutine evalsplgrid

end module sll_quasi_neutral_solver


