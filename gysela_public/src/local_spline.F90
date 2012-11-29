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
      
!---------------------------------------------
! file : local_spline
! date : 28/02/2006
! cubic spline interpolation by using local
!  cubic spline bases
! => Developped by G. LATU and N. Crouseilles
!---------------------------------------------
module local_spline_module
  use prec_const
  use mem_alloc_module
  implicit none
      
  public :: new_splinehh, del_splinehh 
      
  type, public :: splinehh
     ! number of points in x and y directions
     integer                              :: nx, ny
     ! buffer size
     integer                              :: buf
     ! information on mesh in x
     real(RKIND)                          :: x0, dx, invhx
     ! information on mesh in y 
     real(RKIND)                          :: y0, dy, invhy
     ! grid in x and y directions
     real(RKIND), dimension(:)  , pointer :: xgrid, ygrid
     ! diagonal terms of L in LU decomposition
     real(RKIND), dimension(:)  , pointer :: lmatx, lmaty
     ! diagonal terms of U in LU decomposition
     real(RKIND), dimension(:)  , pointer :: umatx, umaty
     ! coefficients des splines
     real(RKIND), dimension(:,:), pointer :: coef, bcoef
     real(RKIND), dimension(:,:), pointer :: aux
     real(RKIND), dimension(:,:), pointer :: auy
  end type splinehh
      
  !*** public variables ***
  real(RKIND) rightconst, leftconst
      
  !******************************
  contains
  !******************************
      
  !-------------------------------------------------------- 
  ! Computes the derivative of 10th order
  !--------------------------------------------------------
  function derive10(vect,dh)
    real(RKIND), dimension (-10:10), intent(in) :: vect
    real(RKIND)                    , intent(in) :: dh
    real(RKIND)                                 :: derive10
      
    real(RKIND) :: tpr
      
    tpr = &
      .2214309755e-5*vect(-10) - &
      .1771447804e-4*vect(-9) + &
      .7971515119e-4*vect(-8) - &
      .3011461267e-3*vect(-7) + &
      .1113797807e-2*vect(-6) - &
      .4145187862e-2*vect(-5) + &
      .1546473933e-1*vect(-4) - &
      .5771376946e-1*vect(-3) + &
      .2153903385*vect(-2) - &
      .8038475846*vect(-1)
    tpr = tpr - &
      .2214309755e-5*vect(10) + &
      .1771447804e-4*vect(9) - &
      .7971515119e-4*vect(8) + &
      .3011461267e-3*vect(7) - &
      .1113797807e-2*vect(6) + &
      .4145187862e-2*vect(5) - &
      .1546473933e-1*vect(4) + &
      .5771376946e-1*vect(3) - &
      .2153903385*vect(2) + &
      .8038475846*vect(1)
    derive10 = tpr/dh
  end function derive10
  
      
  !-------------------------------------------------------- 
  ! Computes the right contribution of the derivative of
  !  10th order
  !--------------------------------------------------------
  function rightderive10(vect,dh)
    real(RKIND), dimension(1:10), intent(in) :: vect
    real(RKIND)                 , intent(in) :: dh
    real(RKIND)                              :: rightderive10
      
    rightderive10 = &
      -.2214309755e-5*vect(10) + &
      .1771447804e-4*vect(9) - &
      .7971515119e-4*vect(8) + &
      .3011461267e-3*vect(7) - &
      .1113797807e-2*vect(6) + &
      .4145187862e-2*vect(5) - &
      .1546473933e-1*vect(4) + &
      .5771376946e-1*vect(3) - &
      .2153903385*vect(2) + &
      .8038475846*vect(1)
    rightderive10 = rightderive10/dh
  end function rightderive10
      
  !-------------------------------------------------------- 
  ! Computes the left contribution of the derivative of
  !  10th order
  !--------------------------------------------------------
  function leftderive10(vect,dh)
    real(RKIND), dimension(-10:-1), intent(in) :: vect
    real(RKIND)                   , intent(in) :: dh
    real(RKIND)                                :: leftderive10
      
    leftderive10 = &
      .2214309755e-5*vect(-10) - &
      .1771447804e-4*vect(-9) + &
      .7971515119e-4*vect(-8) - &
      .3011461267e-3*vect(-7) + &
      .1113797807e-2*vect(-6) - &
      .4145187862e-2*vect(-5) + &
      .1546473933e-1*vect(-4) - &
      .5771376946e-1*vect(-3) + &
      .2153903385*vect(-2) - &
      .8038475846*vect(-1)
    leftderive10 = leftderive10/dh
  end function leftderive10
      
  !------------------------------------------------------ 
  ! Constructor
  !------------------------------------------------------
  subroutine new_splinehh(hsthis,nx,ny,x0,y0,dx,dy)
    use globals
    type(splinehh), intent(inout) :: hsthis
    ! number of points in x and y directions
    integer       , intent(in)    :: nx, ny 
    ! coordinates of mesh origin in each direction
    real(RKIND)   , intent(in)    :: x0, y0
    ! cell size in each direction 
    real(RKIND)   , intent(in)    :: dx, dy
      
    !*** local variables ***
    real(RKIND) :: vect(-stencil:stencil) 
    integer     :: i, j, info
      
    ! initialization of the mesh
    hsthis%nx = nx
    hsthis%ny = ny
    hsthis%x0 = x0
    hsthis%dx = dx
    hsthis%y0 = y0
    hsthis%dy = dy
      
      call glob_allocate(hsthis%xgrid,1,nx+3,'hsthis%xgrid')
      call glob_allocate(hsthis%ygrid,1,ny+3,'hsthis%ygrid')
      if (.not.memory_test) then
        do i = 1,nx+3
          hsthis%xgrid(i)= x0+(i-1)*dx
        end do
        do i = 1,ny+3
          hsthis%ygrid(i)= y0+(i-1)*dy
        end do
      end if
      
      ! memory allocation of matrices
      call glob_allocate(hsthis%coef,-2,nx+1, &
        -2,ny+1,'hsthis%coef')
      call glob_allocate(hsthis%bcoef,-2,nx+1, &
        -2,ny+1,'hsthis%bcoef')
      call glob_allocate(hsthis%aux,-2,nx+1,-2,ny+1,'hsthis%aux')
      call glob_allocate(hsthis%auy,-2,nx+1,-2,ny+1,'hsthis%auy')
      call glob_allocate(hsthis%lmatx,1,nx+2,'hsthis%lmatx')
      call glob_allocate(hsthis%lmaty,1,ny+2,'hsthis%lmaty')
      call glob_allocate(hsthis%umatx,1,3,'hsthis%umatx')
      call glob_allocate(hsthis%umaty,1,3,'hsthis%umaty')
      
      if (.not.memory_test) then
        ! Compute LU decomposition in x and y directions
        call LUHermiteSpline(hsthis%lmatx,hsthis%umatx,nx)
        call LUHermiteSpline(hsthis%lmaty,hsthis%umaty,ny)
        
        hsthis%invhx = 1._RKIND/dx
        hsthis%invhy = 1._RKIND/dy
        vect(-stencil:stencil) = 1._RKIND
        
        rightconst = rightderive10(vect(1:stencil),1.0_RKIND)
        leftconst  = leftderive10(vect(-stencil:-1),1.0_RKIND)
      end if
  end subroutine new_splinehh
      
  !------------------------------------------------------ 
  ! Destructor
  !------------------------------------------------------
  subroutine del_splinehh(hsthis)
    type(splinehh), intent(inout) :: hsthis                
      
    call glob_deallocate(hsthis%aux)
    call glob_deallocate(hsthis%auy)
    call glob_deallocate(hsthis%xgrid)
    call glob_deallocate(hsthis%ygrid)
    call glob_deallocate(hsthis%lmatx)
    call glob_deallocate(hsthis%lmaty)
    call glob_deallocate(hsthis%umatx)
    call glob_deallocate(hsthis%umaty)
    call glob_deallocate(hsthis%bcoef)
    call glob_deallocate(hsthis%coef)
  end subroutine del_splinehh
      
  !*****************************************************
  ! Computation of the spline coefficients
  !*****************************************************
  !----------------------------------------------------------
  ! Compute LU decomposition for splines with 
  !  Hermite boundary conditions :
  !   -> initialization of li and di
  !      which are the inverse of diagonal terms of U matrix 
  !      in the LU decomposition
  !----------------------------------------------------------  
  subroutine LUHermiteSpline(li,di,n)
    real(RKIND), dimension(1:), intent(out) :: li,di 
    integer                   , intent(in)  :: n
      
    !*** local variables ***
    integer                     :: i
    real(RKIND)                 :: d
    real(RKIND), dimension(1:n) :: Dinv
   
    di    = 0._RKIND
    li    = 0._RKIND
    d     = 7._RKIND/2._RKIND
    li(1) = 0.25_RKIND
    li(2) = 2._RKIND/7._RKIND
    do i = 2,n-2
      li(i) = 1._RKIND/d
      d     = 4._RKIND-li(i)
    enddo
    di(1)   = d
    li(n-1) = 1._RKIND/di(1)
    di(2)   = 4._RKIND-li(n-1)
    li(n)   = 1._RKIND/(di(1)*di(2)) 
    di(3)   = 1._RKIND-li(n)
    li(n+1) = 1._RKIND/di(3)    
  end subroutine LUHermiteSpline
      
  !----------------------------------------------------------
  ! Compute LU decomposition for splines with 
  !  Hermite second derivative boundary conditions :
  !   -> initialization of li and di
  !      which are the inverse of diagonal terms of U matrix 
  !      in the LU decomposition
  !----------------------------------------------------------  
  subroutine LUHermiteSpline2(li,di,h,is,n)
    real(RKIND), dimension(is:n), intent(out) :: li
    real(RKIND), dimension(is:n), intent(out) :: di 
    integer                     , intent(in)  :: n,is
      
    !*** local variables ***
    integer                     :: i
    real(RKIND)                 :: h,a,b,c,d
    real(RKIND), dimension(1:n) :: Dinv
      
    !*** Second derivative is imposed ***   
    a  =  6._RKIND/h**2
    b  = -2._RKIND*a
    c  =  a
      
    di =  0._RKIND
    li =  0._RKIND
      
    di(is)   = a
    li(is)   = 1._RKIND/di(is)
    di(is+1) = 4._RKIND-b*li(is)
    li(is+1) = 1._RKIND/di(is+1)
    di(is+2) = 4._RKIND-(1._RKIND-c/a)*li(is+1)
      
    do i = is+2,n-2
      li(i)   = 1._RKIND/di(i)
      di(i+1) = 4._RKIND-li(i)
    enddo
      
    li(n)   = a/di(n-2) 
    li(n-1) = (b-li(n))/di(n-1)
    di(n)   = c-li(n-1)   
  end subroutine LUHermiteSpline2
      
  !------------------------------------------------------------- 
  ! Compute Hermite spline coefficients
  !   rhs(:,1:ntx) contains values of function to be 
  !   interpolated on the full patch and rhs(:0) and 
  !   rhs(:,nx+1) contains the values of the
  !   derivatives of the function at points x_0 and x_N 
  !   respectively.
  !   The transposed function is being used for better 
  !   data locality
  !------------------------------------------------------------- 
  subroutine hermite(hsthis,rhs)
    use globals, only : bufsize
    type(splinehh)                  , intent(inout) :: hsthis
    real(RKIND), &
      dimension(-bufsize:,-bufsize:), intent(inout) :: rhs
    
    !*** local variables ***
    integer     :: i, j, k 
    integer     :: nx, ny
    real(RKIND) :: threeOverh, hOverThree, twoh
    
    !*** initializations for the x direction ***
    nx         = hsthis%nx
    ny         = hsthis%ny
    threeOverh = 3._RKIND/hsthis%dx
    hOverThree = hsthis%dx/3._RKIND
    twoh       = 2._RKIND*hsthis%dx
      
    !-> Sweep down 
    do j = -2,ny+1
      hsthis%aux(-1,j) = rhs(-1,j)
      hsthis%aux(0,j)  = rhs(0,j) + hOverThree * hsthis%aux(-1,j)
      do i = 1,nx-1
        hsthis%aux(i,j) = rhs(i,j) - hsthis%lmatx(i) * &
          hsthis%aux(i-1,j)
      end do
      hsthis%aux(nx,j) = rhs(nx,j)  &
        + threeOverh * hsthis%lmatx(nx-1) * hsthis%aux(nx-2,j) &
        - threeOverh * hsthis%lmatx(nx) * hsthis%aux(nx-1,j)
    enddo
      
    !-> Sweep up 
    do j = -2,ny+1 
      hsthis%bcoef(nx,j)   = twoh*hsthis%aux(nx,j)/hsthis%umatx(3)
      hsthis%bcoef(nx-1,j) = (6._RKIND*hsthis%aux(nx-1,j) - &
        hsthis%bcoef(nx,j)) / hsthis%umatx(2)
      do i = nx-2,1,-1
        hsthis%bcoef(i,j) = (6._RKIND*hsthis%aux(i,j) - &
          hsthis%bcoef(i+1,j)) * hsthis%lmatx(i+1)
      end do
      hsthis%bcoef(0,j)    = (6._RKIND*hsthis%aux(0,j) - &
        2._RKIND*hsthis%bcoef(1,j)) * hsthis%lmatx(1)
      hsthis%bcoef(-1,j)   = hsthis%bcoef(1,j)-twoh*hsthis%aux(-1,j)
    enddo
    do j=-2,ny+1
      hsthis%bcoef(-2,j)   = 6._RKIND*rhs(-2,j) - &
        4._RKIND*hsthis%bcoef(-1,j) - hsthis%bcoef(0,j)
      hsthis%bcoef(nx+1,j) = 6._RKIND*rhs(nx+1,j) - &
        4._RKIND*hsthis%bcoef(nx,j) - hsthis%bcoef(nx-1,j)
    enddo
      
    !*** initializations for the y direction ***
    threeOverh = 3._RKIND/hsthis%dy
    hOverThree = hsthis%dy/3._RKIND
    twoh       = 2._RKIND*hsthis%dy
      
    !-> Sweep down
    do i = -2,nx+1
      hsthis%auy(i,-1) = hsthis%bcoef(i,-1)
      hsthis%auy(i,0)  = hsthis%bcoef(i,0) + &
        hOverThree * hsthis%auy(i,-1)
      do j = 1, ny-1
        hsthis%auy(i,j) = hsthis%bcoef(i,j) - &
          hsthis%lmaty(j) * hsthis%auy(i,j-1)
      enddo
      hsthis%auy(i,ny) = hsthis%bcoef(i,ny) & 
        + threeOverh * hsthis%lmaty(ny-1) * hsthis%auy(i,ny-2) & 
        - threeOverh * hsthis%lmaty(ny) * hsthis%auy(i,ny-1)
    enddo
      
    !-> Sweep up
    do i = -2,nx+1
      hsthis%coef(i,ny)   = twoh*hsthis%auy(i,ny)/hsthis%umaty(3)
      hsthis%coef(i,ny-1) = (6._RKIND*hsthis%auy(i,ny-1) - &
        hsthis%coef(i,ny)) / hsthis%umaty(2)
      do j = ny-2,1,-1
        hsthis%coef(i,j) = (6._RKIND*hsthis%auy(i,j) - &
          hsthis%coef(i,j+1)) * hsthis%lmaty(j+1)
      enddo
      hsthis%coef(i,0)    = (6._RKIND*hsthis%auy(i,0) - &
        2._RKIND*hsthis%coef(i,1)) * hsthis%lmaty(1)
      hsthis%coef(i,-1)   = hsthis%coef(i,1)-twoh*hsthis%auy(i,-1)
    enddo
    do i=-2,nx+1
      hsthis%coef(i,-2)   = 6._RKIND*hsthis%bcoef(i,-2) - &
        4._RKIND*hsthis%coef(i,-1) - hsthis%coef(i,0)
      hsthis%coef(i,ny+1) = 6._RKIND*hsthis%bcoef(i,ny+1) - &
        4._RKIND*hsthis%coef(i,ny) - hsthis%coef(i,ny-1)
    enddo
  end subroutine hermite
      
  !------------------------------------------------------ 
  ! Compute spline basis in x direction
  !------------------------------------------------------  
  subroutine splinex_basis(hsthis,xprev,xstar,xnext,sbase)
    type(splinehh)              , intent(in)    :: hsthis          
    real(RKIND)                 , intent(in)    :: xprev, xstar
    real(RKIND)                 , intent(in)    :: xnext
    real(RKIND), dimension(-1:2), intent(inout) :: sbase
      
    real(RKIND) :: coeff, d_prev, d_next
      
    coeff     = hsthis%invhx*hsthis%invhx*hsthis%invhx
    d_prev    = xstar-xprev
    d_next    = xnext-xstar
      
    sbase(2)  = coeff*d_prev*d_prev*d_prev*0.166666666666667_RKIND
    sbase(1)  = (1._RKIND + coeff * &
      (3._RKIND*d_prev * &
      (hsthis%dx*hsthis%dx+hsthis%dx*d_prev-d_prev*d_prev))) * &
      0.166666666666667_RKIND
    sbase(0)  = (1._RKIND + coeff * &
      (3._RKIND*d_next * &
      (hsthis%dx*hsthis%dx+hsthis%dx*d_next-d_next*d_next))) * &
      0.166666666666667_RKIND
    sbase(-1) = (coeff*d_next*d_next*d_next) * &
      0.166666666666667_RKIND  
  end subroutine splinex_basis
      
  !------------------------------------------------------ 
  ! Compute spline basis in y direction
  !------------------------------------------------------  
  subroutine spliney_basis(hsthis,xprev,xstar,xnext,sbase)
    type(splinehh)              , intent(in)    :: hsthis
    real(RKIND)                 , intent(in)    :: xprev, xstar
    real(RKIND)                 , intent(in)    :: xnext
    real(RKIND), dimension(-1:2), intent(inout) :: sbase
      
    real(RKIND) :: coeff, d_prev, d_next
      
    coeff     = hsthis%invhy*hsthis%invhy*hsthis%invhy
    d_prev    = xstar-xprev
    d_next    = xnext-xstar
      
    sbase(2)  = coeff*d_prev*d_prev*d_prev*0.166666666666667_RKIND
    sbase(1)  = (1._RKIND + coeff * &
      (3._RKIND*d_prev*(hsthis%dy*hsthis%dy + &
      hsthis%dy*d_prev-d_prev*d_prev))) * &
      0.166666666666667_RKIND
    sbase(0)  = (1._RKIND + coeff * &
      (3._RKIND*d_next*(hsthis%dy*hsthis%dy + &
      hsthis%dy*d_next-d_next*d_next))) * &
      0.166666666666667_RKIND
    sbase(-1) = (coeff*d_next*d_next*d_next) * &
      0.166666666666667_RKIND
  end subroutine spliney_basis
  
      
  !------------------------------------------------------ 
  ! 2D interpolation
  !------------------------------------------------------  
  subroutine interpol2d(hsthis,xstar1,xstar2,finterpol)
    type(splinehh), intent(in)    :: hsthis                
    real(RKIND)   , intent(in)    :: xstar1, xstar2
    real(RKIND)   , intent(inout) :: finterpol
      
    integer                      :: ipos1, ipos2
    integer                      :: i, j
    real(RKIND), dimension(-1:2) :: sbase1, sbase2
      
    !*** array position location in r ***
    ! -> if the particle is out of the domain 
    !    is put on the boundaries
    ipos1 = 1 + hsthis%invhx*(xstar1-hsthis%xgrid(1))
    ipos1 = max(min(ipos1,hsthis%nx+1),1)
      
    !*** array position location in theta ***
    ! -> if the particle is out of the domain 
    !    is put on the boundaries
    ipos2 = 1 + hsthis%invhy*(xstar2-hsthis%ygrid(1))
    ipos2 = max(min(ipos2,hsthis%ny+1),1)
      
    !*** calculation of cubic spline basis ***
    call splinex_basis(hsthis,hsthis%xgrid(ipos1), &
      xstar1,hsthis%xgrid(ipos1+1),sbase1)
    call spliney_basis(hsthis,hsthis%ygrid(ipos2), &
      xstar2,hsthis%ygrid(ipos2+1),sbase2)
      
    !*** computation of f(x*,y*) ***
    finterpol = 0._RKIND
    ipos1     = ipos1-2
    ipos2     = ipos2-2
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2-1)*sbase1(-1)*sbase2(-1) + &
      hsthis%coef(ipos1,ipos2-1)*sbase1(0)*sbase2(-1) + &
      hsthis%coef(ipos1+1,ipos2-1)*sbase1(1)*sbase2(-1) + &
      hsthis%coef(ipos1+2,ipos2-1)*sbase1(2)*sbase2(-1)
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2)*sbase1(-1)*sbase2(0) + &
      hsthis%coef(ipos1,ipos2)*sbase1(0)*sbase2(0) + &
      hsthis%coef(ipos1+1,ipos2)*sbase1(1)*sbase2(0) + &
      hsthis%coef(ipos1+2,ipos2)*sbase1(2)*sbase2(0)
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2+1)*sbase1(-1)*sbase2(1) + &
      hsthis%coef(ipos1,ipos2+1)*sbase1(0)*sbase2(1) + &
      hsthis%coef(ipos1+1,ipos2+1)*sbase1(1)*sbase2(1) + &
      hsthis%coef(ipos1+2,ipos2+1)*sbase1(2)*sbase2(1)
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2+2)*sbase1(-1)*sbase2(2) + &
      hsthis%coef(ipos1,ipos2+2)*sbase1(0)*sbase2(2) + &
      hsthis%coef(ipos1+1,ipos2+2)*sbase1(1)*sbase2(2) + &
      hsthis%coef(ipos1+2,ipos2+2)*sbase1(2)*sbase2(2)
  end subroutine interpol2d
      
  !----------------------------------------------------------
  ! Verify r coordinate, i.e.
  !  . if r<rmin set r=rmin 
  !  . if r>rmax set r=rmax 
  !----------------------------------------------------------
  subroutine local_r_verif(hsthis,r)
    use globals, only : mu_id, nbpart_rleft_iter, nbpart_rright_iter
    type(splinehh), intent(in)    :: hsthis
    real(RKIND)   , intent(inout) :: r
    
    real(RKIND) :: r_min_loc, r_max_loc 
    
    r_min_loc = hsthis%xgrid(1)
    r_max_loc = hsthis%xgrid(hsthis%nx+2)
    if (r.lt.r_min_loc) then
      r = r_min_loc
      nbpart_rleft_iter(mu_id) = nbpart_rleft_iter(mu_id) + 1
    else 
      if (r.gt.r_max_loc) then
        r = r_max_loc
        nbpart_rright_iter(mu_id) = nbpart_rright_iter(mu_id) + 1
      end if
    end if
  end subroutine local_r_verif
      
  !----------------------------------------------------------
  ! Verify theta coordinate, i.e.
  !  . if theta<theta_min or theta>theta_max use theta 
  !     periodicity
  !----------------------------------------------------------
   subroutine local_theta_verif(hsthis,geom,theta)
    use globals, only : mu_id, Nbproc_theta, &
      nbpart_thleft_iter, nbpart_thright_iter
    use geometry_class
    include "mpiperso.h"
    type(splinehh), intent(in)    :: hsthis
    type(geometry), intent(in)    :: geom
    real(RKIND)   , intent(inout) :: theta
    
    integer     :: ierr
    real(RKIND) :: theta_min_loc, theta_max_loc
    real(RKIND) :: theta_mid_glob
      
    theta_min_loc  = hsthis%ygrid(1)
    theta_max_loc  = hsthis%ygrid(hsthis%ny+1)
    theta_mid_glob = geom%thetag(0)+ geom%Ltheta/2._RKIND
    if (theta.lt.theta_min_loc) then
      if ((theta_max_loc.gt.geom%thetag(geom%Ntheta-1)).and. &
        (theta.lt.theta_mid_glob)) then
        !-> processor with the last patch and theta<theta_mid_glob
        theta                      = min(theta+TWOPI,theta_max_loc)
        nbpart_thright_iter(mu_id) = nbpart_thright_iter(mu_id) + 1
      else
        theta                     = theta_min_loc
        nbpart_thleft_iter(mu_id) = nbpart_thleft_iter(mu_id) + 1
      end if
    elseif (theta.gt.theta_max_loc) then
      if ((theta_min_loc.lt.geom%thetag(0)).and. &
        (theta.gt.theta_mid_glob)) then
        !-> processor with the first patch and theta>theta_mid_glob
        theta                     = max(theta-TWOPI,theta_min_loc)
        nbpart_thleft_iter(mu_id) = nbpart_thleft_iter(mu_id) + 1
      else
        theta                      = theta_max_loc
        nbpart_thright_iter(mu_id) = nbpart_thright_iter(mu_id) + 1
      end if
    end if
  end subroutine local_theta_verif
end module local_spline_module
