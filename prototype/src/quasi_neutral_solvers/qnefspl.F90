!-----------------------------------------------------------------------------------------------------------------
!  Module for solving the quasi-neutral equation in polar coordinates using arbitrary order spline finite elements
!  Eric Sonnendrucker 2011-05-05
!-----------------------------------------------------------------------------------------------------------------
module qnefspl_module
  implicit none
  type :: qndata
     integer :: spline_degree  ! degree of splines used by finite element solver
     integer :: nr, ntheta  ! dimensions in r and theta
     real(8) :: dr, dtheta  ! cell size in r and theta
     real(8) :: rmin
     real(8), dimension(:), pointer :: knotsr, knotsth  ! knot vectors for splines in r and theta direction
     real(8), dimension(:,:), pointer :: Kar, Mar, Mcr ! symmetric banded matrices in r direction
     real(8), dimension(:), pointer ::  Kth, Mth  ! coefficients of symmetric circulant matrices in theta direction
     real(8), dimension(:), pointer :: xgauss, wgauss ! Gauss points and weights
  end type qndata
contains
  subroutine new_qn(this, spline_degree, rmin, nr, ntheta, dr, dtheta)
    type (qndata) :: this
    integer :: spline_degree ! degree of spline basis functions
    real(8) :: rmin   ! Minimum value of r
    integer :: nr, ntheta  ! dimensions in r and theta
    real(8) :: dr, dtheta  ! cell size in r and theta

    ! local variables
    real(8) :: rprof   ! r profile for Kar and Mar matrices
    integer :: i,j,k,ig,j1,j2  ! loop indices
    real(8) :: xg  ! gauss point
    real(8) :: ri,ti  ! grid point, knot
    real(8), dimension(spline_degree+1,spline_degree+1) :: work  ! work array for de Boor routines
    real(8), dimension(spline_degree+1,2) :: dbiatx  ! values of splines and derivatives at point x

    !print*, 'new_qn', spline_degree, nr, ntheta, dr, dtheta

    ! set scalars in qndata object
    this%spline_degree = spline_degree
    this%nr = nr
    this%ntheta = ntheta
    this%dr = dr
    this%dtheta = dtheta
    this%rmin = rmin

    ! allocate arrays
    allocate(this%Kar(spline_degree+1,nr+spline_degree-1))
    allocate(this%Mar(spline_degree+1,nr+spline_degree-1))
    allocate(this%Mcr(spline_degree+1,nr+spline_degree-1))
    allocate(this%Kth(spline_degree+1))
    allocate(this%Mth(spline_degree+1))
    allocate(this%xgauss(spline_degree+1))
    allocate(this%wgauss(spline_degree+1))
    allocate(this%knotsr(nr+2*spline_degree)) 
    allocate(this%knotsth(2*spline_degree+2))

    ! set Gauss points and weights
    select case(this%spline_degree)
    case(1) 
       this%xgauss(1) = -1.0_8/sqrt(3._8); this%xgauss(2) = 1.0_8/sqrt(3.0_8)
       this%wgauss(1) = 1.0_8 ; this%wgauss(2) = 1.0_8
    case(2)
       this%xgauss(1) = -sqrt(3.0_8/5.0_8); this%xgauss(2) = 0. ;  this%xgauss(3) = sqrt(3.0_8/5.0_8)
       this%wgauss(1) = 5.0_8/9.0_8; this%wgauss(2) = 8.0_8/9.0_8
       this%wgauss(3) = this%wgauss(1)
    case(3)
       this%xgauss(4) = sqrt((3.0_8+2.0_8*sqrt(6.0_8/5.0_8))/7.0_8); 
       this%xgauss(3) = sqrt((3.0_8-2.0_8*sqrt(6.0_8/5.0_8))/7.0_8);
       this%xgauss(2) = -this%xgauss(3) ;  this%xgauss(1) = -this%xgauss(4) ;
       this%wgauss(1) = (18.0_8-sqrt(30.0_8))/36.0_8;  this%wgauss(2) = (18.0_8+sqrt(30.0_8))/36.0_8  
       this%wgauss(3) = this%wgauss(2);  this%wgauss(4) = this%wgauss(1);
    case(4)
       this%xgauss(1) = -0.90617984593866374_8; this%xgauss(2) = -0.53846931010568311_8;this%xgauss(3) = 0._8;
       this%xgauss(4) = 0.53846931010568311_8;  this%xgauss(5) = 0.90617984593866374_8;
       this%wgauss(1) = 0.23692688505618875_8;  this%wgauss(2) = 0.47862867049936653_8; this%wgauss(3) = 0.568888888888888888_8;
       this%wgauss(4) = 0.47862867049936653_8;  this%wgauss(5) = 0.23692688505618875_8;
    case(5)
       this%xgauss(1) = -0.93246951420315249_8; this%xgauss(2) = -0.66120938646626437_8;this%xgauss(3) = -0.23861918608319668_8;
       this%xgauss(4) = 0.23861918608319668_8;  this%xgauss(5) = 0.66120938646626437_8;  this%xgauss(6) = 0.93246951420315249_8;
       this%wgauss(1) = 0.17132449237917155_8;  this%wgauss(2) = 0.36076157304813761_8;  this%wgauss(3) = 0.46791393457269115_8;
       this%wgauss(4) = 0.46791393457269115_8;  this%wgauss(5) = 0.36076157304813761_8;  this%wgauss(6) = 0.17132449237917155_8;
    case(6)
       this%xgauss(1) = -0.94910791234275871_8; this%xgauss(2) = -0.74153118559939302_8; this%xgauss(3) = -0.4058451513773973_8;
       this%xgauss(4) = 0._8;                   this%xgauss(5) = 0.4058451513773973_8;   this%xgauss(6) = 0.74153118559939302_8;
       this%xgauss(7) = 0.94910791234275871_8;
       this%wgauss(1) = 0.12948496616886976_8;  this%wgauss(2) = 0.27970539148927637_8; this%wgauss(3) = 0.38183005050511931_8;
       this%wgauss(4) = 0.41795918367346824_8;  this%wgauss(5) = 0.38183005050511931_8; this%wgauss(6) = 0.27970539148927637_8;
       this%wgauss(7) = 0.12948496616886976_8;
    case(7)
       this%xgauss(1) = -0.96028985649753706_8; this%xgauss(2) = -0.79666647741362673_8;this%xgauss(3) = -0.52553240991632888_8;
       this%xgauss(4) = -0.18343464249565036_8; this%xgauss(5) = 0.18343464249565036_8; this%xgauss(6) = 0.52553240991632888_8;
       this%xgauss(7) = 0.79666647741362673_8; this%xgauss(8) = 0.96028985649753706_8;
       this%wgauss(1) = 0.10122853629037747_8;  this%wgauss(2) = 0.22238103445337337_8; this%wgauss(3) = 0.31370664587788677_8;
       this%wgauss(4) = 0.36268378337836155_8;  this%wgauss(5) = 0.36268378337836155_8; this%wgauss(6) = 0.31370664587788677_8;
       this%wgauss(7) = 0.22238103445337337_8; this%wgauss(8) = 0.10122853629037747_8;
    case default
       print*, 'spline degree ', spline_degree, ' not implemented'
       stop
    end select

    ! Compute circulant matrices in theta direction
    !----------------------------------------------
    ! set knot vector
    do k=-spline_degree,spline_degree+1
       this%knotsth(k+spline_degree+1)=dtheta*k 
    enddo
    ! initialize matrices
    this%Mth(:) = 0.0
    this%Kth(:) = 0.0
    do ig=1,spline_degree+1   ! sum over Gauss points
       xg = 0.5_8*dtheta*(this%xgauss(ig)+1.)  ! rescale Gauss points to be in interval [0,h]
       !print*, xg, knots(spline_degree+1),knots(spline_degree+2)
       call bsplvd(this%knotsth,spline_degree+1,xg,spline_degree+1,work,dbiatx,2)
       ! sum over bsplines in interval [x_i, x_i+1]
       do j = 0,spline_degree
          do i = 0,spline_degree 
             if (i+j < spline_degree+1) then
                this%Mth(j+1) = this%Mth(j+1) + 0.5_8*dtheta*this%wgauss(ig)*dbiatx(i+1,1)*dbiatx(i+j+1,1)
                this%Kth(j+1) = this%Kth(j+1) + 0.5_8*dtheta*this%wgauss(ig)*dbiatx(i+1,2)*dbiatx(i+j+1,2)
             end if
          end do
       end do
    end do
   
    ! Compute banded matrices in r direction (Dirichlet BC)
    !------------------------------------------------------
    ! set knot vector
    do i = 1, spline_degree + 1
       this%knotsr(i) = rmin
     enddo
    ti = rmin
    do i = spline_degree + 2, Nr + spline_degree
       ti = ti+dr
       this%knotsr(i) = ti
     enddo
    do i = Nr + spline_degree + 1, Nr + 2 * spline_degree
       this%knotsr(i) = ti
    enddo
    
    ! intialize arrays
    this%Mar(:,:) = 0.0
    this%Mcr(:,:) = 0.0
    this%Kar(:,:) = 0.0
    ! loop over cells
    do i = 0, Nr-2
       ri = rmin + i*dr
       ! sum over Gauss points
       do ig = 1, spline_degree+1
          xg = ri + 0.5_8*dr*(this%xgauss(ig)+1.0_8)  ! rescale Gauss points to be in interval [ri,ri+dr]
          rprof = 1.0_8   ! r-profile at Gauss point xg 
          call bsplvd(this%knotsr,spline_degree+1,xg,spline_degree+1+i,work,dbiatx,2)            
          do j1 = 0, spline_degree
             do j2 = j1, spline_degree
                this%Mar(j2-j1+1,i+j1+1) = this%Mar(j2-j1+1,i+j1+1) &
                     + 0.5_8*rprof/xg*dr*this%wgauss(ig)*dbiatx(j1+1,1)*dbiatx(j2+1,1)
                this%Kar(j2-j1+1,i+j1+1) = this%Kar(j2-j1+1,i+j1+1) &
                     + 0.5_8*rprof*xg*dr*this%wgauss(ig)*dbiatx(j1+1,2)*dbiatx(j2+1,2)
             enddo
          enddo
       enddo
    enddo

  end subroutine new_qn

  subroutine solve2Dtensor(this,U,F)
    ! Solve 2D QN equation on structured grid in polar coordinates with Dirichlet BC in r and periodic in theta
    ! Tensor product technique used to bring the solution back to 1D problems
    ! Discrete problem is of the form Kar U Mtheta +Mar U K th +Mcr U Mth = F
    ! where the unknow U and the RHS F are in matrix form and U is multiplied on left or right by 1D matrices
    ! See Crouseilles, Ratnani, Sonnendrucker 2011
    !--------------------------------------------------------------------------------------------------------
    ! Input :
    type(qndata), intent(inout) :: this
    real(8), dimension(:,:), intent(inout) :: F !  RHS in matrix form
    ! Output : 
    real(8), dimension(:,:), intent(out) :: U ! solution in matrix form
    ! local variables 
    integer :: info, i,k, nband
    real(8), dimension(this%spline_degree+1,this%nr+this%spline_degree-3) :: A ! banded matrix for problem in r 
    real(8), dimension(2,this%nr+this%spline_degree-3) :: V ! auxiliary matrix
    real(8), dimension(this%ntheta/2+1) :: valpm, valpk ! eigenvalues of circulant mass and stiffness matrix in theta direction 
    real(8), parameter :: pi = 3.1415926535897931_8
    real(8), dimension(2*this%ntheta+15) :: wsave  ! auxiliary array for FFTPACK

    ! intialize dffft
    call dffti(this%ntheta, wsave) 

    ! forward FFT of lines of F with FFTPACK
    do i=2,this%nr+this%spline_degree-2
       call dfftf(this%ntheta, F(i,:), wsave) 
    end do

    ! compute eigenvalues of circulant matrices
    do k=1,this%ntheta/2+1
       valpm(k) = this%Mth(1)
       valpk(k) = this%Kth(1)
       do i=2,this%spline_degree+1
          valpm(k) = valpm(k) + 2*this%Mth(i)*cos(2*pi*(i-1)*(k-1)/this%ntheta)
          valpk(k) = valpk(k) + 2*this%Kth(i)*cos(2*pi*(i-1)*(k-1)/this%ntheta)
       end do
    end do


    ! Banded solves in x direction
    nband = this%spline_degree+1
    U(:,:) = F(:,:)  ! copy rhs into solution
    A(:,:) = valpm(1)*(this%Kar(:,2:this%nr+this%spline_degree-2) + this%Mcr(:,2:this%nr+this%spline_degree-2))&
         + valpk(1)*this%Mar(:,2:this%nr+this%spline_degree-2)  
    call DPBTRF( 'L', this%nr+this%spline_degree-3, nband-1, A, nband, info ) ! Cholesky factorisation of A
    call DPBTRS( 'L', this%nr+this%spline_degree-3, nband-1, 1, A, nband, U(2:this%nr+this%spline_degree-2,1), &
         this%nr+this%spline_degree-3, info ) ! Solution
    do k = 1, this%ntheta/2-1
       A(:,:) = valpm(k+1)*(this%Kar(:,2:this%nr+this%spline_degree-2) + this%Mcr(:,2:this%nr+this%spline_degree-2)) &
               + valpk(k+1)*this%Mar(:,2:this%nr+this%spline_degree-2)
       call DPBTRF( 'L', this%nr+this%spline_degree-3, nband-1, A, nband, info ) ! Cholesky factorisation of A
       call DPBTRS( 'L', this%nr+this%spline_degree-3, nband-1, 1, A, nband, U(2:this%nr+this%spline_degree-2,2*k), &
         this%nr+this%spline_degree-3, info )
       call DPBTRS( 'L', this%nr+this%spline_degree-3, nband-1, 1, A, nband, U(2:this%nr+this%spline_degree-2,2*k+1), &
            this%nr+this%spline_degree-3, info )
    end do
    A(:,:) = valpm(this%ntheta/2+1)*(this%Kar(:,2:this%nr+this%spline_degree-2) + this%Mcr(:,2:this%nr+this%spline_degree-2)) &
           + valpk(this%ntheta/2+1)*this%Mar(:,2:this%nr+this%spline_degree-2) 
    call DPBTRF( 'L', this%nr+this%spline_degree-3, nband-1, A, nband, info ) ! Cholesky factorisation of A
    call DPBTRS( 'L', this%nr+this%spline_degree-3, nband-1, 1, A, nband, U(2:this%nr+this%spline_degree-2,this%ntheta), &
         this%nr+this%spline_degree-3, info ) ! Solution

    ! backward FFT of lines of U with FFTPACK
    U(1,:) = 0.0 
    U(this%nr+this%spline_degree-1,:) = 0.0
    do i=1,this%nr+this%spline_degree-1
       call dfftb(this%ntheta, U(i,:), wsave) 
    end do
    U(:,:) = U(:,:)/this%ntheta    ! normalization
 
  end subroutine solve2Dtensor

  subroutine evalsplgrid(this,C,U)
    ! evaluate spline defined by coefficient array C at all grid points
    !------------------------------------------------------------------
    type (qndata) :: this
    real(8), dimension(:,:) :: C  ! spline coefficients
    real(8), dimension(:,:) :: U  ! Value of function at grid points
    ! local variables
    real(8), dimension(this%spline_degree+1) :: biatr,biatth
    real(8) :: r
    integer   :: i,j,ii,jj

    U = 0.0
    ! Values of periodic spline at grid points (need only one call as it is uniform)
    call bsplvb(this%knotsth,this%spline_degree+1,1,0.0,this%spline_degree+1,biatth)
    do i = 0, this%nr-2
       do ii=0,this%spline_degree
          r= this%rmin + i*this%dr
          call bsplvb(this%knotsr,this%spline_degree+1,1,r,this%spline_degree+1+i,biatr)
          do j = 0, this%ntheta-1
             do jj=0, this%spline_degree              
                U(i+1,j+1) = U(i+1,j+1) + biatr(ii+1)*biatth(jj+1)*C(i+ii+1,mod(j+jj,this%ntheta)+1)
             enddo
          enddo
       enddo
    enddo
  end subroutine evalsplgrid
end module qnefspl_module


