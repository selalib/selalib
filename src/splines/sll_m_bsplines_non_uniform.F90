module sll_m_bsplines_non_uniform
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64
  use sll_m_bsplines_base    , only: sll_c_bsplines

  implicit none

  public :: &
    sll_t_bsplines_non_uniform

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type, B-splines
  type, extends(sll_c_bsplines) :: sll_t_bsplines_non_uniform

  contains
    procedure :: init                    => s_bsplines_non_uniform__init
    procedure :: free                    => s_bsplines_non_uniform__free
    procedure :: find_cell               => f_bsplines_non_uniform__find_cell
    procedure :: eval_basis              => s_bsplines_non_uniform__eval_basis
    procedure :: eval_deriv              => s_bsplines_non_uniform__eval_deriv
    procedure :: eval_basis_and_n_derivs => s_bsplines_non_uniform__eval_basis_and_n_derivs

  end type sll_t_bsplines_non_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine s_bsplines_non_uniform__init( self, degree, periodic, breaks ) 
    class(sll_t_bsplines_non_uniform), intent(  out) :: self
    integer                          , intent(in   ) :: degree
    logical                          , intent(in   ) :: periodic
    real(wp),             allocatable, intent(in   ) :: breaks(:)

    integer  :: i
    real(wp) :: period ! length of period for periodic B-splines

    ! Public attributes
    self % degree   = degree
    self % periodic = periodic
    self % uniform  = .false.
    self % ncells   = size(breaks) - 1
    self % nbasis   = merge( self%ncells, self%ncells+degree, periodic )
    self % offset   = merge( degree/2   , 0                 , periodic )
    self % xmin     = breaks(1)
    self % xmax     = breaks(self%ncells+1)

    ! Create the knots array from the grid points. Here take the grid points
    ! as knots and simply add to the left and the right the
    ! amount of knots that depends on the degree of the requested 
    ! spline. We aim at setting up the indexing in such a way that the 
    ! original indexing of 'grid' is preserved, i.e.: grid(i) = knot(i), at
    ! least whenever the scope of the indices defined here is active.
    associate( npoints => self%ncells+1 )

      allocate( self % knots(1-degree:npoints+degree) )

      do i = 1, npoints
        self%knots(i) = breaks(i)
      end do

      ! Fill out the extra nodes 
      if ( periodic ) then
        period = breaks(npoints) - breaks(1)
        do i = 1, degree
          self%knots(1-i)       = breaks(npoints-i) - period
          self%knots(npoints+i) = breaks(i+1)       + period
        end do
      else ! open
        do i = 1, degree
          self%knots(1-i)       = breaks(1)
          self%knots(npoints+i) = breaks(npoints)
        end do
      end if

    end associate

  end subroutine s_bsplines_non_uniform__init

  !-----------------------------------------------------------------------------
  !> Free storage
  subroutine s_bsplines_non_uniform__free( self )
    class(sll_t_bsplines_non_uniform), intent(inout) :: self

    deallocate( self % knots )

  end subroutine s_bsplines_non_uniform__free

  !----------------------------------------------------------------------------
  function f_bsplines_non_uniform__find_cell( self, x ) result( icell )
    class(sll_t_bsplines_non_uniform), intent(in) :: self
    real(wp)                         , intent(in) :: x
    integer :: icell

    integer :: low, high

    associate( npoints => self%ncells+1 )

      ! check if point is outside of grid
      if (x > self%knots(npoints)) then; icell = -1; return; end if
      if (x < self%knots(1)      ) then; icell = -1; return; end if

      ! check if point is exactly on left/right boundary
      if (x == self%knots(1)      ) then; icell = 1          ; return; end if
      if (x == self%knots(npoints)) then; icell = self%ncells; return; end if

      low   = 1
      high  = npoints
      icell = (low + high) / 2
      do while (x <  self%knots(icell) .or. x >= self%knots(icell+1))
         if (x < self%knots(icell)) then
           high = icell
         else
           low  = icell
         end if
         icell = (low + high) / 2
      end do

    end associate

  end function

  !-----------------------------------------------------------------------------
  !> Evaluate value at x of all basis functions with support in local cell
  !> (jmin identifies index of basis functions)
  SLL_PURE subroutine s_bsplines_non_uniform__eval_basis( self, x, values, jmin )
    class(sll_t_bsplines_non_uniform), intent(in   ) :: self
    real(wp)                         , intent(in   ) :: x
    real(wp)                         , intent(  out) :: values(0:)
    integer                          , intent(  out) :: jmin

    integer  :: icell

    integer  :: j, r
    real(wp) :: temp
    real(wp) :: saved

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:self%degree)
    real(wp) :: right(1:self%degree)

    ! Check on inputs
    SLL_ASSERT( x > self%xmin - 1.0d-14 )
    SLL_ASSERT( x < self%xmax + 1.0d-14 )
    SLL_ASSERT( size(values) == 1+self%degree )

    ! 1. Compute cell index 'icell'
!    icell = f_bsplines_non_uniform__find_cell( self, x )
    icell = self % find_cell( x )
    SLL_ASSERT( icell >= 1               )
    SLL_ASSERT( icell <= self%ncells     )
    SLL_ASSERT( self%knots(icell) <= x   )
    SLL_ASSERT( x <= self%knots(icell+1) ) 

    ! 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell - self%offset

    ! 3. Compute values of aforementioned B-splines
    values(0) = 1.0_wp
    do j = 1, self%degree
       left (j) = x - self%knots(icell+1-j)
       right(j) = self%knots(icell+j) - x
       saved    = 0.0_wp
       do r = 0, j-1
          temp      = values(r) / (right(r+1) + left(j-r))
          values(r) = saved + right(r+1) * temp
          saved     = left(j-r) * temp
       end do
      values(j) = saved
    end do

  end subroutine s_bsplines_non_uniform__eval_basis

  !-----------------------------------------------------------------------------
  !> Evaluate derivative at x of all basis functions with support in local cell
  SLL_PURE subroutine s_bsplines_non_uniform__eval_deriv( self, x, derivs, jmin )
    class(sll_t_bsplines_non_uniform), intent(in   ) :: self
    real(wp)                         , intent(in   ) :: x
    real(wp)                         , intent(  out) :: derivs(0:)
    integer                          , intent(  out) :: jmin

    integer  :: icell

    integer  :: j, r
    real(wp) :: temp
    real(wp) :: saved

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:self%degree)
    real(wp) :: right(1:self%degree)

    ! Check on inputs
    SLL_ASSERT( x > self%xmin - 1.0d-14 )
    SLL_ASSERT( x < self%xmax + 1.0d-14 )
    SLL_ASSERT( size(derivs) == 1+self%degree )

    ! 1. Compute cell index 'icell' and x_offset
!    icell = f_bsplines_non_uniform__find_cell( self, x )
    icell = self % find_cell( x )
    SLL_ASSERT( icell >= 1               )
    SLL_ASSERT( icell <= self%ncells     )
    SLL_ASSERT( self%knots(icell) <= x   )
    SLL_ASSERT( x <= self%knots(icell+1) )

    ! 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell - self%offset

    ! 3. Compute derivatives of aforementioned B-splines
    associate( degree => self%degree, degree_real => real(self%degree,wp) )

      ! Compute nonzero basis functions and knot differences
      ! for splines up to degree deg-1 which are needed to compute derivative
      ! First part of Algorithm  A3.2 of NURBS book 
      derivs(0) = 1.0_wp
      do j = 1, degree-1
         left (j) = x - self%knots(icell+1-j)
         right(j) = self%knots(icell+j) - x
         saved    = 0.0_wp
         do r = 0, j-1
            ! compute and save bspline values
            temp      = derivs(r)/(right(r+1) + left(j-r))
            derivs(r) = saved + right(r+1) * temp
            saved     = left(j-r) * temp
         end do
         derivs(j) = saved
      end do

      ! Compute derivatives at x using values stored in bsdx and formula
      ! formula for spline derivative based on difference of splines of 
      ! degree deg-1
      ! -------
      ! j = 0
      saved = degree_real * derivs(0) / (self%knots(icell+1)-self%knots(icell+1-degree))
      derivs(0) = -saved
      do j = 1, degree-1
         temp      = saved 
         saved     = degree_real * derivs(j) / (self%knots(icell+j+1)-self%knots(icell+j+1-degree))
         derivs(j) = temp - saved
      end do
      ! j = degree
      derivs(degree) =  saved

    end associate

  end subroutine s_bsplines_non_uniform__eval_deriv

  !-----------------------------------------------------------------------------
  !> Evaluate value and n derivatives at x of all basis functions with support in local cell
  SLL_PURE subroutine s_bsplines_non_uniform__eval_basis_and_n_derivs( self, x, n, derivs, jmin )
    class(sll_t_bsplines_non_uniform), intent(in   ) :: self
    real(wp)                         , intent(in   ) :: x
    integer                          , intent(in   ) :: n
    real(wp)                         , intent(  out) :: derivs(0:,0:)
    integer                          , intent(  out) :: jmin

    integer  :: icell

    integer  :: j, r
    real(wp) :: temp
    real(wp) :: saved

    integer  :: k, s1, s2, rk, pk, j1, j2
    real(wp) :: d

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:self%degree)
    real(wp) :: right(1:self%degree)
    real(wp) :: ndu  (0:self%degree,0:self%degree)
    real(wp) :: a    (0:1          ,0:self%degree)

    ! Check on inputs
    SLL_ASSERT( x > self%xmin - 1.0d-14 )
    SLL_ASSERT( x < self%xmax + 1.0d-14 )
    SLL_ASSERT( n >= 0 )
    SLL_ASSERT( n <= self%degree )
    SLL_ASSERT( size(derivs,1) == 1+n           )
    SLL_ASSERT( size(derivs,2) == 1+self%degree )

    ! 1. Compute cell index 'icell' and x_offset
!    icell = f_bsplines_non_uniform__find_cell( self, x )
    icell = self % find_cell( x )
    SLL_ASSERT( icell >= 1               )
    SLL_ASSERT( icell <= self%ncells     )
    SLL_ASSERT( self%knots(icell) <= x   )
    SLL_ASSERT( x <= self%knots(icell+1) )

    ! 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell - self%offset

    ! 3. Compute nonzero basis functions and knot differences for splines
    !    up to degree (degree-1) which are needed to compute derivative
    !    Algorithm  A2.3 of NURBS book 
    !
    !    21.08.2017: save inverse of knot differences to avoid unnecessary divisions
    !                [Yaman Güçlü, Edoardo Zoni]
    associate( degree => self%degree )

      ndu(0,0) = 1.0_wp
      do j = 1, degree 
         left(j)  = x - self%knots(icell+1-j)
         right(j) = self%knots(icell+j) - x
         saved    = 0.0_wp
         do r = 0, j-1
            ! compute inverse of knot differences and save them into lower triangular part of ndu
            ndu(j,r) = 1.0_wp / (right(r+1) + left(j-r))
            ! compute basis functions and save them into upper triangular part of ndu
            temp     = ndu(r,j-1) * ndu(j,r)
            ndu(r,j) = saved + right(r+1) * temp
            saved    = left(j-r) * temp
         end do
         ndu(j,j) = saved
      end do
      derivs(0,:) = ndu(:,degree)

      do r = 0, degree
         s1 = 0
         s2 = 1
         a(0,0) = 1.0_wp
         do k = 1, n
            d  = 0.0_wp
            rk = r-k
            pk = degree-k
            if (r >= k) then
               a(s2,0) = a(s1,0) * ndu(pk+1,rk)
               d = a(s2,0) * ndu(rk,pk)
            end if
            if (rk > -1) then
               j1 = 1
            else
               j1 = -rk
            end if
            if (r-1 <= pk) then
               j2 = k-1
            else
               j2 = degree-r
            end if
            do j = j1, j2
               a(s2,j) = (a(s1,j) - a(s1,j-1)) * ndu(pk+1,rk+j)
               d = d + a(s2,j) * ndu(rk+j,pk)
            end do
            if (r <= pk) then
               a(s2,k) = - a(s1,k-1) * ndu(pk+1,r)
               d = d + a(s2,k) * ndu(r,pk)
            end if
            derivs(k,r) = d
            j  = s1
            s1 = s2
            s2 = j
         end do
      end do
      r = degree
      do k = 1, n
         derivs(k,:) = derivs(k,:) * r
         r = r * (degree-k)
      end do

    end associate

  end subroutine s_bsplines_non_uniform__eval_basis_and_n_derivs

end module sll_m_bsplines_non_uniform
