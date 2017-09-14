!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @ingroup splines
!> @brief 
!> Low level arbitrary degree splines
!> @details
!> This module defines low level algorithms of the arbitrary
!> degree spline.
!> For non periodic boundary conditions, knots defined from a grid
!> are either duplicated (open knot sequence) or mirror outside boundary.
!> It is a selalib implemnation of the classical algorithms
!> found in the de Boor book "Practical guide to splines"
!> or the NURBS book by Piegl and Tiller.  

module sll_m_low_level_bsplines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: &
    f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_open,     &
    sll_p_mirror

  implicit none

  public :: &
    sll_t_bsplines, &
    sll_s_bsplines_init_from_grid, &
    sll_s_bsplines_init_from_knots, &
    sll_s_bsplines_free, &
    sll_s_bsplines_eval_basis, &
    sll_s_bsplines_eval_deriv, &
    sll_s_bsplines_eval_basis_and_deriv, &
    sll_s_bsplines_eval_basis_and_n_derivs, &
    ! Binary search algorithm (should be in 'meshes' or 'low_level_utilities'):
    sll_f_find_cell, &
    ! Michel Mehrenberger's algorithms (slightly faster but not fully robust):
    sll_s_bsplines_eval_basis_mm, &
    sll_s_bsplines_eval_basis_and_deriv_mm, &
    ! Evaluation of uniform spline basis (faster, should go to separate module):
    sll_s_uniform_bsplines_eval_basis             , &
    sll_s_uniform_bsplines_eval_deriv             , &
    sll_s_uniform_bsplines_eval_basis_and_deriv   , &
    sll_s_uniform_bsplines_eval_basis_and_n_derivs, &
    ! Evaluation of uniform periodic spline (use 'sll_m_spline_1d' instead):
    sll_s_eval_uniform_periodic_spline_curve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> List of allowed boundary conditions
  integer, parameter :: allowed_bcs(*) = [sll_p_periodic, sll_p_open, sll_p_mirror]

  !> Information for evaluation of B-splines on non-uniform grid
  type :: sll_t_bsplines
     integer               :: num_pts
     integer               :: deg
     real(wp)              :: xmin
     real(wp)              :: xmax
     integer               :: n        ! dimension of spline space
     real(wp), allocatable :: knots(:) ! knots array
  end type sll_t_bsplines

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief
  !> Build new sll_t_bsplines object 
  !> based on a on a grid of strictly increasing points including the last point also
  !> for periodic domains.
  !> @param[inout] basis arbitrary_degree_spline_1d object
  !> @param[in] degree  spline_degree
  !> @param[in] grid    x coordinates of grid points
  !> @param[in] num_pts number of grid points
  !> @param[in] bc_xmin boundary condition at xmin [periodic/open/mirror]
  !> @param[in] bc_xmax boundary condition at xmax [periodic/open/mirror]
  !>
  !> @details
  !>
  !> The logic behind the periodic boundary condition is the following.
  !> The given grid array has minimum (grid(1)) and maximum (grid(n))
  !> values at either end. This defines a length 'L'. If interpreted
  !> as a periodic space, this is also the period. Thus, as we extend
  !> the number of knots at both ends of the given array, we use the
  !> periodicity condition to fill out the new values:
  !>
  !>                    .
  !>                    .
  !>                    .
  !>           knots(-1) = knots(n-2) - L
  !>           knots( 0) = knots(n-1) - L
  !>                    .
  !>                    .
  !>                    .
  !>           knots(n+1) = knots(1) + L
  !>           knots(n+2) = knots(2) + L
  !>                    .
  !>                    .
  !>                    .
  !>
  !> The 'open' boundary condition simply extends the new values
  !> of the local array at both ends with repeated endpoint values.
  !> That is
  !>
  !>     ... = knots(-2) = knots(-1) = knots(0) = knots(1)
  !>
  !> and
  !>
  !>    knots(n+1) = knots(n+2) = knots(n+3) = ... =  knots(n)
  !>
  !>
  !> The mirror boundary condition mirrors the knot values on
  !> each side of the grid
  !> of the local array at both ends with repeated endpoint values.
  !> That is
  !>
  !>     ... =  = knots(-1) = knots(0) = knots(1)
  !>
  !> and
  !>
  !>    knots(n+1) = knots(n+2) = knots(n+3) = ... =  knots(n)
  !-----------------------------------------------------------------------------
  subroutine sll_s_bsplines_init_from_grid( &
      basis  , &
      degree , &
      grid   , &
      bc_xmin, &
      bc_xmax )

    type(sll_t_bsplines), intent(  out) :: basis
    integer             , intent(in   ) :: degree
    real(wp)            , intent(in   ) :: grid(:)
    integer             , intent(in   ) :: bc_xmin
    integer             , intent(in   ) :: bc_xmax

    character(len=*), parameter :: this_sub_name = "sll_s_bsplines_init_from_grid"
    character(len=256) :: err_msg

    integer  :: i
    integer  :: num_pts
    real(wp) :: period ! length of period
    real(wp) :: min_cell_width

    ! Check that polynomial degree is at least 1
    if( degree < 1 ) then
      write(err_msg,"('Minimum degree = 1, given ',i0,' instead.')") degree
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    ! Check that grid contains at least two points
    num_pts = size( grid )
    if( num_pts < 2 ) then
      write(err_msg,"('Minimum size(grid) = 2, given ',i0,' instead.')") num_pts
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    ! Check whether grid points are in strictly increasing order
    min_cell_width = grid(2)-grid(1)
    do i = 3, num_pts
      min_cell_width = min( min_cell_width, grid(i)-grid(i-1) )
    end do
    if (min_cell_width <= 0.0_wp) then
      err_msg = "Grid points must be in strictly increasing order."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    else if (min_cell_width < 1.e-10_wp) then
      err_msg = "Length of grid cells must be >= 1e-10."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    ! Check that boundary conditions are OK
    if (.not. any( bc_xmin == allowed_bcs )) then
      err_msg = "Unrecognized boundary condition at xmin: " // &
                "possible values are [sll_p_periodic, sll_p_open, sll_p_mirror]."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if
    if (.not. any( bc_xmax == allowed_bcs )) then
      err_msg = "Unrecognized boundary condition at xmax: " // &
                "possible values are [sll_p_periodic, sll_p_open, sll_p_mirror]."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if
    if (any( [bc_xmin,bc_xmax]== sll_p_periodic ) .and. bc_xmin /= bc_xmax) then
      err_msg = "Periodic boundary conditions mismatch: bc_xmin /= bc_xmax."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    ! Periodic case: check that there are enough grid points for a given degree
    if (any( [bc_xmin,bc_xmax]== sll_p_periodic ) .and. num_pts < 2+degree) then
      err_msg = "Insufficient number of grid points for periodic spline: " // &
                "condition num_pts >= 2+degree not satisfied."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    basis%num_pts = num_pts
    basis%deg     = degree
    basis%xmin    = grid(1)
    basis%xmax    = grid(num_pts)

    ! 'period' is useful in the case of periodic boundary conditions.
    period = grid(num_pts) - grid(1)

    ! Create the knots array from the grid points. Here take the grid points
    ! as knots and simply add to the left and the right the
    ! amount of knots that depends on the degree of the requested 
    ! spline. We aim at setting up the indexing in such a way that the 
    ! original indexing of 'grid' is preserved, i.e.: grid(i) = knot(i), at
    ! least whenever the scope of the indices defined here is active.
    allocate( basis%knots (1-degree:num_pts+degree) )
    do i=1,num_pts
       basis%knots(i) = grid(i)
    end do

    ! Allocate array for spline coefficients
    if (bc_xmin == sll_p_periodic) then
       basis%n = num_pts - 1 ! dimension of periodic spline space
    else 
       basis%n = num_pts + degree - 1 ! dimension of non periodic spline space
    end if

    ! Fill out the extra points at both ends of the local knot array with
    ! values proper to the boundary condition requested.
    !
    ! Fill out the extra nodes on the left
    select case (bc_xmin)
    case (sll_p_periodic)
      do i = 1, degree
        basis%knots(1-i) = grid(num_pts-i) - period
      end do
    case (sll_p_open)
      do i = 1, degree
        basis%knots(1-i) = grid(1)
      end do
    case (sll_p_mirror)
      do i = 1, degree
        basis%knots(1-i) = 2*grid(1) - grid(1+i)
      end do
    end select
    !
    ! Fill out the extra nodes on the right
    select case (bc_xmax)
    case (sll_p_periodic)
      do i = 1, degree
        basis%knots(num_pts+i) = grid(i+1) + period
      end do
    case (sll_p_open)
      do i = 1, degree
        basis%knots(num_pts+i) = grid(num_pts)
      end do
    case (sll_p_mirror)
      do i = 1, degree
        basis%knots(num_pts+i) = 2*grid(num_pts) - grid(num_pts-i)
      end do
    end select

  end subroutine sll_s_bsplines_init_from_grid

  !-----------------------------------------------------------------------------
  !> @brief
  !> Build new sll_t_bsplines object 
  !> @details
  !> based on a given array of knots
  !> which is a non decreasing array
  !> @param[inout] basis spline object
  !> @param[in] degree spline degree
  !> @param[in] n dimension of spline space
  !> @param[in] knots array of give knots
  !-----------------------------------------------------------------------------
  subroutine sll_s_bsplines_init_from_knots( basis, degree, n, knots )

    type(sll_t_bsplines), intent(  out) :: basis
    integer             , intent(in   ) :: degree
    integer             , intent(in   ) :: n
    real(wp)            , intent(in   ) :: knots(:)
  
    integer  :: i
    real(wp) :: knotmin
    
    ! Error checking
    SLL_ASSERT( size(knots) == n+2*degree )
    if( degree < 1 ) then
       print *, 'ERROR. sll_s_bsplines_init_from_knots: ', &
            'only strictly positive integer values for degree are allowed, ', &
            'given: ', degree
    end if
    ! Check whether grid points are in strictly increasing order
    knotmin = knots(2)-knots(1)
    do i=2, n+2*degree
       knotmin = min(knotmin, knots(i)-knots(i-1)) 
    end do
    if (knotmin < 0.0_wp) then
       print*, 'ERROR. sll_s_bsplines_init_from_knots: ', &
            ' we require that knots  are in non decreasing.'
    end if

    ! Initialise variables
    basis%n       = n
    basis%deg     = degree
    basis%num_pts = n - degree + 1
    basis%xmin    = knots(degree+1)
    basis%xmax    = knots(n+1)

    ! We aim at setting up the indexing in such a way that the 
    ! the computation grid is [knot(1),knots(num_pts)],  at
    ! least whenever the scope of the indices defined here is active.
    allocate( basis%knots (1-degree:n+1) )
    basis%knots(1-degree:n+1) = knots
    
  end subroutine sll_s_bsplines_init_from_knots
  
  !-----------------------------------------------------------------------------
  !> @brief return index i of grid cell such that:
  !> basis%knots(i) <= x <= basis%knots(i+1).
  !>
  !> @detail
  !> If x is not between basis%knots(1) and basis%knots(basis%num_pts),
  !> then the value -1 is returned.
  !-----------------------------------------------------------------------------
  pure function sll_f_find_cell( basis, x ) result( icell )

    type(sll_t_bsplines), intent(in) :: basis
    real(wp)            , intent(in) :: x
    integer :: icell

    integer :: low
    integer :: high
    integer :: n

    n = basis%num_pts

    ! check if point is outside of grid
    if (x > basis%knots(n)) then
       icell = -1
       return
    end if
    if (x < basis%knots(1)) then
       icell = -1
       return
    end if

    ! check if point is exactly on left boundary
    if (x == basis%knots(1)) then
       icell = 1
       return
    end if

    ! check if point is exactly on right boundary
    if (x == basis%knots(n)) then
       icell = n-1
       return
    end if

    low  = 1
    high = n
    icell = (low + high) / 2
    do while (x <  basis%knots(icell) &
         .or. x >= basis%knots(icell+1))
       if (x < basis%knots(icell)) then
          high = icell
       else
          low  = icell
       end if
       icell = (low + high) / 2
    end do

  end function sll_f_find_cell

  !-----------------------------------------------------------------------------
  ! sll_f_splines_at_x returns the values of all the B-splines of a given 
  ! degree that have support in cell 'cell' and evaluated at the point 'x'. 
  ! In other words, if B[j,i](x) is the spline of degree 'j' whose leftmost 
  ! support is at cell 'i' and evaluated at 'x', then sll_f_splines_at_x returns 
  ! the sequence (in the form of an array):
  ! 
  ! B[j,i-degree](x), B[j,i-degree+1](x), B[j,i-degree+2](x), ..., B[j,i](x)
  !
  ! Implementation notes: 
  !
  ! The algorithm is based on Algorithm A2.2 from 
  ! (L. Piegl, W. Tiller, The NURBS book p. 70)
  ! A variant of the same Algorithm is implemented in the de Boor splines.
  ! It is known as the Cox - de Boor algorithm

  !> @brief
  !> Evaluates B-spline values at a point x in a given cell.
  !> @details
  !> sll_s_bsplines_eval_basis( basis, cell, x, splines_at_x ) computes the values of all the
  !> splines which have support in 'cell' and evaluates them at point 'x',
  !> which is supposed to be in cell. The spline object should have already
  !> been initialized and will contain information on the spline degree
  !> to use and the type of boundary condition desired.
  !> The algorithm implemented is numerically stable and known as The
  !> Cox - de Boor algorithm, which is a generalisation to splines of the
  !> de Casteljau algorithm for Bezier curves.
  !> @return b_spline_at_x B-spline values
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bsplines_eval_basis( basis, icell, x, splines_at_x )

    type(sll_t_bsplines), intent(in   ) :: basis
    integer             , intent(in   ) :: icell
    real(wp)            , intent(in   ) :: x
    real(wp)            , intent(  out) :: splines_at_x(0:basis%deg)

    real(wp) :: saved
    real(wp) :: temp
    integer  :: j
    integer  :: r

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:basis%deg)
    real(wp) :: right(1:basis%deg)

    ! Run some checks on the arguments.
    SLL_ASSERT( x > basis%xmin - 1.0d-14 )
    SLL_ASSERT( x < basis%xmax + 1.0d-14 )
    SLL_ASSERT( icell >= 1 )
    SLL_ASSERT( icell <= basis%num_pts - 1 )
    SLL_ASSERT( basis%knots(icell) <= x .and. x <= basis%knots(icell+1) )

    splines_at_x(0) = 1.0_wp
    do j = 1, basis%deg
       left (j) = x - basis%knots(icell+1-j)
       right(j) = basis%knots(icell+j) - x
       saved    = 0.0_wp
       do r = 0, j-1
          temp = splines_at_x(r) / (right(r+1) + left(j-r))
          splines_at_x(r) = saved + right(r+1) * temp
          saved = left(j-r) * temp
       end do
      splines_at_x(j) = saved
    end do

  end subroutine sll_s_bsplines_eval_basis

  !-----------------------------------------------------------------------------
  !> @brief
  !> returns first derivative values at x of all b-splines with support in cell
  !> @details
  !> sll_s_bsplines_eval_deriv returns an array with the derivative values of 
  !> the B-splines of a requested order that are supported in 'cell' and 
  !> evaluated at 'x'. 
  !> Algorithm derived from algorithm A3.2 of NURBS book 
  !> The return value has the format:
  !> \f[
  !> B'[deg,i-deg](x), B'[deg,i-deg+1](x), ..., B'[deg,i](x)
  !> \f]
  !> where 'deg' is the degree of the spline.
  !> @param[out] bsdx  B-spline derivatives
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bsplines_eval_deriv( basis, icell, x, bsdx )

    type(sll_t_bsplines), intent(in   ) :: basis
    integer                             , intent(in   ) :: icell
    real(wp)                            , intent(in   ) :: x
    real(wp)                            , intent(  out) :: bsdx(0:basis%deg)

    integer  :: deg
    integer  :: num_pts
    integer  :: r
    integer  :: j
    real(wp) :: saved
    real(wp) :: temp
    real(wp) :: rdeg

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:basis%deg)
    real(wp) :: right(1:basis%deg)

    ! Run some checks on the arguments.
    SLL_ASSERT( allocated( basis%knots ) )
    SLL_ASSERT( x > basis%xmin - 1.0d-14 )
    SLL_ASSERT( x < basis%xmax + 1.0d-14 )
    SLL_ASSERT( icell >= 1 )
    SLL_ASSERT( icell <= basis%num_pts - 1 )
    SLL_ASSERT( basis%knots(icell) <= x .and. x <= basis%knots(icell+1) )

    deg     = basis%deg
    rdeg    = real(deg,wp)
    num_pts = basis%num_pts

    ! compute nonzero basis functions and knot differences
    ! for splines up to degree deg-1 which are needed to compute derivative
    ! First part of Algorithm  A3.2 of NURBS book 
    bsdx(0) = 1.0_wp
    do j = 1, deg-1
       left (j) = x - basis%knots(icell+1-j)
       right(j) = basis%knots(icell+j) - x
       saved    = 0.0_wp
       do r = 0, j-1
          ! compute and save bspline values
          temp    = bsdx(r)/(right(r+1) + left(j-r))
          bsdx(r) = saved + right(r+1) * temp
          saved   = left(j-r) * temp
       end do
       bsdx(j) = saved
    end do

    ! Compute derivatives at x using values stored in bsdx and formula
    ! formula for spline derivative based on difference of splines of 
    ! degree deg-1
    ! -------
    ! j = 0
    saved = rdeg*bsdx(0) / (basis%knots(icell+1) - basis%knots(icell+1-deg))
    bsdx(0) = -saved
    do j = 1, deg-1
       temp    = saved 
       saved   = rdeg*bsdx(j) / (basis%knots(icell+j+1)-basis%knots(icell+j+1-deg))
       bsdx(j) = temp - saved
    end do
    ! j = deg
    bsdx(deg) =  saved

  end subroutine sll_s_bsplines_eval_deriv

  !-----------------------------------------------------------------------------
  !> @brief 
  !> returns splines and first derivatives
  !> @details
  !> See sll_s_bsplines_eval_deriv and  sll_f_splines_at_x
  !> @return b_spline_and_derivs_at_x B-spline values and derivatives
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bsplines_eval_basis_and_deriv( basis, icell, x, bsdx )

    type(sll_t_bsplines), intent(in   ) :: basis
    integer             , intent(in   ) :: icell
    real(wp)            , intent(in   ) :: x
    real(wp)            , intent(  out) :: bsdx(2,0:basis%deg)

    integer  :: deg
    integer  :: num_pts
    integer  :: r
    integer  :: j
    real(wp) :: saved
    real(wp) :: temp
    real(wp) :: rdeg

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:basis%deg)
    real(wp) :: right(1:basis%deg)

    ! Run some checks on the arguments.
    SLL_ASSERT( allocated( basis%knots ) )
    SLL_ASSERT( x > basis%xmin - 1.0d-14 )
    SLL_ASSERT( x < basis%xmax + 1.0d-14 )
    SLL_ASSERT( icell >= 1 )
    SLL_ASSERT( icell <= basis%num_pts - 1 )
    SLL_ASSERT( basis%knots(icell) <= x .and. x <= basis%knots(icell+1) )

    deg = basis%deg
    num_pts = basis%num_pts
    rdeg = real(deg,wp) 

    ! compute nonzero basis functions and knot differences
    ! for splines up to degree deg-1 which are needed to compute derivative
    ! First part of Algorithm  A2.3 of NURBS book 
    bsdx(1,0) = 1.0_wp
    do j = 1, deg-1
       left (j) = x - basis%knots(icell+1-j)
       right(j) = basis%knots(icell+j) - x
       saved = 0.0_wp
       do r = 0, j-1
          ! compute and save knot differences
          temp = bsdx(1,r)/(right(r+1) + left(j-r))
          bsdx(1,r) = saved + right(r+1) * temp
          saved = left(j-r) * temp
       end do
       bsdx(1,j) = saved
    end do

    ! Compute derivatives at x using values stored in bsdx and formula
    ! formula for spline derivative based on difference of splines of 
    ! degree deg-1
    ! -------
    ! j = 0
    saved = rdeg * bsdx(1,0) / &
         (basis%knots(icell+1) - basis%knots(icell+1-deg))
    bsdx(2,0) = -saved
    do j = 1, deg-1
       temp = saved 
       saved =  rdeg*bsdx(1,j) / &
            (basis%knots(icell+j+1) - basis%knots(icell+j+1-deg))
       bsdx(2,j) = temp - saved
    end do
    ! j = deg
    bsdx(2,deg) =  saved  

    ! Compute values of splines of degree deg
    !----------------------------------------
    j = deg
    left(j)  = x - basis%knots(icell+1-j)
    right(j) = basis%knots(icell+j) - x
    saved    = 0.0_wp
    do r = 0, j-1
       ! compute and save knot differences
       temp = bsdx(1,r)/(right(r+1) + left(j-r))
       bsdx(1,r) = saved + right(r+1) * temp
       saved = left(j-r) * temp
    end do
    bsdx(1,j) = saved 

  end subroutine sll_s_bsplines_eval_basis_and_deriv

  !-----------------------------------------------------------------------------
  !> @brief 
  !> returns splines and first derivatives
  !> @details
  !> See sll_s_bsplines_eval_deriv and  sll_f_splines_at_x
  !> @param[in] basis  bspline object
  !> @param[in] icell cell where bsplines are to be computed
  !> @param[in] x value of point where bsplines are to be computed
  !> @param[in] n number of derivatives to be computed
  !> @param[out] bsdx B-spline values and the first n derivatives
  !
  ! TODO: transpose output array 'bsdx'
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bsplines_eval_basis_and_n_derivs( basis, icell, x, n, bsdx )

    type(sll_t_bsplines), intent(in   ) :: basis
    integer             , intent(in   ) :: icell
    real(wp)            , intent(in   ) :: x
    integer             , intent(in   ) :: n
    real(wp)            , intent(  out) :: bsdx(0:n,0:basis%deg)

    integer  :: deg
    integer  :: num_pts
    integer  :: r
    integer  :: k
    integer  :: j
    integer  :: j1
    integer  :: j2
    integer  :: s1
    integer  :: s2
    integer  :: rk
    integer  :: pk
    real(wp) :: saved
    real(wp) :: temp
    real(wp) :: rdeg
    real(wp) :: d

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: left (1:basis%deg)
    real(wp) :: right(1:basis%deg)
    real(wp) :: ndu  (0:basis%deg,0:basis%deg)
    real(wp) :: a    (0:1        ,0:basis%deg)

    ! Run some checks on the arguments.
    SLL_ASSERT( allocated( basis%knots ) )
    SLL_ASSERT( x > basis%xmin - 1.0d-14 )
    SLL_ASSERT( x < basis%xmax + 1.0d-14 )
    SLL_ASSERT( icell >= 1 )
    SLL_ASSERT( icell <= basis%num_pts - 1 )
    SLL_ASSERT( n >= 0 )
    SLL_ASSERT( n <= basis%deg )
    SLL_ASSERT( basis%knots(icell) <= x .and. x <= basis%knots(icell+1) )

    deg     = basis%deg
    num_pts = basis%num_pts
    rdeg    = real(deg,wp) 

    ! compute nonzero basis functions and knot differences
    ! for splines up to degree deg-1 which are needed to compute derivative
    ! Algorithm  A2.3 of NURBS book 
    !
    ! 21.08.2017: save inverse of knot differences to avoid unnecessary divisions
    !             [Yaman Güçlü, Edoardo Zoni]

    ndu(0,0) = 1.0_wp
    do j = 1, deg 
       left(j)  = x - basis%knots(icell+1-j)
       right(j) = basis%knots(icell+j) - x
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
    bsdx(0,:) = ndu(:,deg)

    do r = 0, deg
       s1 = 0
       s2 = 1
       a(0,0) = 1.0_wp
       do k = 1, n
          d  = 0.0_wp
          rk = r-k
          pk = deg-k
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
             j2 = deg-r
          end if
          do j = j1, j2
             a(s2,j) = (a(s1,j) - a(s1,j-1)) * ndu(pk+1,rk+j)
             d = d + a(s2,j) * ndu(rk+j,pk)
          end do
          if (r <= pk) then
             a(s2,k) = - a(s1,k-1) * ndu(pk+1,r)
             d = d + a(s2,k) * ndu(r,pk)
          end if
          bsdx(k,r) = d
          j  = s1
          s1 = s2
          s2 = j
       end do
    end do
    r = deg
    do k = 1, n
       bsdx(k,:) = bsdx(k,:) * r
       r = r * (deg-k)
    end do

  end subroutine sll_s_bsplines_eval_basis_and_n_derivs

  !-----------------------------------------------------------------------------
  !> @brief Alternative direct implentation of recursion formula. 
  !>
  !> @detail
  !> This provides an evaluation of B-splines directly based on the recurrence
  !> formula. It is 10% faster than the classical Cox - de Boor formula 
  !> that is implented in sll_f_splines_at_x, but can have numerical stability issues.
  !> For this reason the Cox - de Boor formula should be the default implementation 
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bsplines_eval_basis_mm( &
       knots, &
       cell, &
       x, &
       degree, &
       out )

    integer , intent(in   ) :: degree
    real(wp), intent(in   ) :: knots(1-degree:)
    integer , intent(in   ) :: cell
    real(wp), intent(in   ) :: x
    real(wp), intent(  out) :: out(:)

    real(wp) :: tmp1
    real(wp) :: tmp2
    integer  :: ell
    integer  :: k

    out(1) = 1._wp
    do ell = 1, degree
       tmp1 = (x-knots(cell+1-ell))/(knots(cell+1)-knots(cell+1-ell))*out(1)
       out(1) = out(1) -tmp1
       do k = 2, ell
          tmp2 = (x-knots(cell+k-ell))/(knots(cell+k)-knots(cell+k-ell))*out(k)
          out(k) = out(k)+tmp1-tmp2
          tmp1 = tmp2
       enddo
       out(ell+1) = tmp1
    enddo

  end subroutine sll_s_bsplines_eval_basis_mm

  !-----------------------------------------------------------------------------
  !> @brief Alternative direct implentation of recursion formula. 
  !>
  !> @detail
  !> This provides an evaluation of B-splines directly based on the recurrence
  !> formula. It is about 80% faster than the classical Cox - de Boor formula 
  !> that is implented in sll_f_splines_at_x, but can have numerical stability issues.
  !> For this reason the Cox - de Boor formula should be the default implementation
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bsplines_eval_basis_and_deriv_mm( &
       knots, &
       cell, &
       x, &
       degree, &
       out )

    integer , intent(in   ) :: degree
    real(wp), intent(in   ) :: knots(1-degree:)
    integer , intent(in   ) :: cell
    real(wp), intent(in   ) :: x
    real(wp), intent(  out) :: out(:,:)

    real(wp) :: tmp1
    real(wp) :: tmp2
    integer  :: ell
    integer  :: k

    out(1,1) = 1._wp
    do ell = 1, degree
       tmp1 = (x-knots(cell+1-ell))/(knots(cell+1)-knots(cell+1-ell))*out(1,1)
       out(1,1) = out(1,1) -tmp1
       do k = 2, ell
          tmp2 = (x-knots(cell+k-ell))/(knots(cell+k)-knots(cell+k-ell))*out(1,k)
          out(1,k) = out(1,k)+tmp1-tmp2
          tmp1 = tmp2
       enddo
       out(2,ell+1) = tmp1
       if(ell==degree-1)then
          !compute the derivatives
          tmp1 = real(degree,wp)/(knots(cell+1)-knots(cell+1-degree))*out(1,1)
          out(2,1) = -tmp1
          do k = 2, degree
             out(2,k) = tmp1
             tmp1 = real(degree,wp)/(knots(cell+k)-knots(cell+k-degree))*out(2,k)
             out(2,k) = out(2,k)-tmp1
          enddo
          out(2,degree+1) = tmp1
       endif
    enddo

  end subroutine sll_s_bsplines_eval_basis_and_deriv_mm

  !-----------------------------------------------------------------------------
  subroutine sll_s_bsplines_free( spline )

    type(sll_t_bsplines), intent(inout) :: spline

    character(len=*), parameter :: this_sub_name = "sll_s_bsplines_free()"

    if( .not. allocated(spline%knots) ) then
       SLL_ERROR( this_sub_name, 'knots array is not allocated.' )
    end if
    deallocate( spline%knots )

  end subroutine sll_s_bsplines_free

  ! *************************************************************************
  !
  !                    UNIFORM B-SPLINE FUNCTIONS
  !
  ! *************************************************************************

  !-----------------------------------------------------------------------------
  !> @brief Evaluate all non vanishing uniform B-Splines in unit cell. 
  !>
  !> @detail Returns an array with the values of the b-splines of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, thus the offset given must be a number
  !> between 0 and 1.
  !> Output: bspl(1:d+1)= B_d(-(d+1)/2+d+x),...,B_d(-(d+1)/2+x) 
  !> with d=spline_degree and x=normalized_offset
  !> where B_d=B_{d-1}*B_0 and B_0=1_[-1/2,1/2] and * is convolution
  !> the following code can be used for comparison with deboor
  !> do i=-d,d+1
  !> t(i+d+1)=real(i,8)
  !> enddo
  !> call bsplvb(t,d+1,1,normalized_offset,d+1,out)
  !> We also have the property (from the symmetry of the B-spline)
  !> out(1:d+1)= B_d(-(d+1)/2+xx),...,B_d(-(d+1)/2+d+xx),..., 
  !> where xx=1-normalized_offset
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_uniform_bsplines_eval_basis( &
    spline_degree,     &
    normalized_offset, &
    bspl               )

    integer , intent(in   ) :: spline_degree
    real(wp), intent(in   ) :: normalized_offset
    real(wp), intent(  out) :: bspl(0:spline_degree)

    integer  :: j, r
    real(wp) :: inv_j
    real(wp) :: j_real
    real(wp) :: xx
    real(wp) :: temp
    real(wp) :: saved

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_wp )
    SLL_ASSERT( normalized_offset <= 1.0_wp )

    bspl(0) = 1.0_wp
    do j = 1, spline_degree
       xx     = -normalized_offset
       j_real = real(j,wp)
       inv_j  = 1.0_wp / j_real
       saved  = 0.0_wp
       do r = 0, j-1
          xx      = xx + 1.0_wp
          temp    = bspl(r) * inv_j
          bspl(r) = saved + xx * temp
          saved   = (j_real - xx) * temp
       end do
       bspl(j) = saved
    end do

  end subroutine sll_s_uniform_bsplines_eval_basis

  !-----------------------------------------------------------------------------
  !> @brief Evaluate all derivatives of non vanishing uniform B-Splines 
  !> in unit cell. 
  !>
  !> @detail
  !> Returns an array with the values of the b-spline derivatives of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, hence the results must be divided by the
  !> real cell size to scale back the results.
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_uniform_bsplines_eval_deriv( &
    spline_degree,     &
    normalized_offset, & 
    bspl               )

    integer , intent(in   ) :: spline_degree
    real(wp), intent(in   ) :: normalized_offset
    real(wp), intent(  out) :: bspl(0:spline_degree)

    real(wp) :: inv_j
    real(wp) :: x, xx
    real(wp) :: j_real
    integer  :: j, r
    real(wp) :: temp
    real(wp) :: saved
    real(wp) :: bj, bjm1 

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_wp )
    SLL_ASSERT( normalized_offset <= 1.0_wp )

    x = normalized_offset
    bspl(0) = 1.0_wp

    ! only need splines of lower degree to compute derivatives
    do j = 1, spline_degree - 1
       saved = 0.0_wp
       j_real = real(j,wp)
       inv_j = 1.0_wp / j_real
       xx = - x
       do r = 0, j-1
          xx = xx + 1.0_wp
          temp = bspl(r) * inv_j
          bspl(r) = saved + xx * temp
          saved = (j_real - xx) * temp
       end do
       bspl(j) = saved
    end do

    ! compute derivatives
    bjm1 = bspl(0)
    bj = bjm1
    bspl(0) = -bjm1
    do j=1, spline_degree - 1
       bj = bspl(j)
       bspl(j) = bjm1 - bj
       bjm1 = bj
    end do
    bspl(spline_degree) = bj

  end subroutine sll_s_uniform_bsplines_eval_deriv

  !-----------------------------------------------------------------------------
  !> @brief Evaluate all values and derivatives of non vanishing uniform B-Splines 
  !> in unit cell. 
  !>
  !> @detail
  !> returns an array with the values of the b-spline derivatives of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, hence the results must be divided by the
  !> real cell size to scale back the results.
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_uniform_bsplines_eval_basis_and_deriv( &
    degree,            &
    normalized_offset, &
    bspl               )

    integer , intent(in   ) :: degree
    real(wp), intent(in   ) :: normalized_offset
    real(wp), intent(  out) :: bspl(2,0:degree)

    real(wp) :: inv_j
    real(wp) :: xx
    real(wp) :: j_real
    integer  :: j, r
    real(wp) :: temp
    real(wp) :: saved

    SLL_ASSERT( degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_wp )
    SLL_ASSERT( normalized_offset <= 1.0_wp )

    ! compute splines up to degree spline_degree -1
    bspl(1,0) = 1.0_wp
    do j = 1, degree - 1
       saved = 0.0_wp
       j_real = real(j,wp)
       inv_j = 1.0_wp / j_real
       xx = - normalized_offset
       do r = 0, j-1
          xx = xx + 1.0_wp
          temp = bspl(1,r) * inv_j
          bspl(1,r) = saved + xx * temp
          saved = (j_real - xx) * temp
       end do
       bspl(1,j) = saved
    end do

    ! compute derivatives
    bspl(2,0) = -bspl(1,0)
    do j = 1, degree - 1
       bspl(2,j) = bspl(1,j-1) - bspl(1,j)
    end do
    bspl(2,degree) = bspl(1,degree-1)

    ! continue the Cox-De Boor algorithm to evaluate splines
    j = degree
    saved = 0.0_wp
    j_real = real(j,wp)
    inv_j = 1.0_wp / j_real
    xx = - normalized_offset
    do r = 0, j-1
       xx = xx + 1.0_wp
       temp = bspl(1,r) * inv_j
       bspl(1,r) = saved + xx * temp
       saved = (j_real - xx) * temp
    end do
    bspl(1,j) = saved

  end subroutine sll_s_uniform_bsplines_eval_basis_and_deriv

  !-----------------------------------------------------------------------------
  ! TODO: transpose output array 'bspl'
  SLL_PURE subroutine sll_s_uniform_bsplines_eval_basis_and_n_derivs( &
      spline_degree,     &
      normalized_offset, &
      n,                 &
      bspl               )

    integer , intent(in   ) :: spline_degree
    real(wp), intent(in   ) :: normalized_offset
    integer , intent(in   ) :: n
    real(wp), intent(  out) :: bspl(0:n,0:spline_degree)

    integer  :: j, r
    real(wp) :: j_real
    real(wp) :: xx
    real(wp) :: temp
    real(wp) :: saved

    integer  :: k, s1, s2, rk, pk, j1, j2
    real(wp) :: d

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: ndu    (0:spline_degree,0:spline_degree)
    real(wp) :: a      (0:1            ,0:spline_degree)

    ! Inverse of integers for later use (max spline degree = 32)
    real(wp), parameter :: inv_idx(*) = [(1.0_wp/real(j,wp), j=1,32)]

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( n >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_wp )
    SLL_ASSERT( normalized_offset <= 1.0_wp )

    ! Evaluate all basis splines (see "sll_s_uniform_bsplines_eval_basis")
    ndu(0,0) = 1.0_wp
    do j = 1, spline_degree
       xx     = -normalized_offset
       j_real = real(j,wp)
       saved  = 0.0_wp
       do r = 0, j-1
          xx       = xx + 1.0_wp
          temp     = ndu(r,j-1) * inv_idx(j)
          ndu(r,j) = saved + xx * temp
          saved    = (j_real - xx) * temp
       end do
       ndu(j,j) = saved
    end do
    bspl(0,:) = ndu(:,spline_degree)

    ! Use equation 2.10 in "The NURBS Book" to compute n derivatives
    associate( deg => spline_degree, bsdx => bspl )

    do r = 0, deg
       s1 = 0
       s2 = 1
       a(0,0) = 1.0_wp
       do k = 1, n
          d  = 0.0_wp
          rk = r-k
          pk = deg-k
          if (r >= k) then
             a(s2,0) = a(s1,0) * inv_idx(pk+1)
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
             j2 = deg-r
          end if
          do j = j1, j2
             a(s2,j) = (a(s1,j) - a(s1,j-1)) * inv_idx(pk+1)
             d = d + a(s2,j) * ndu(rk+j,pk)
          end do
          if (r <= pk) then
             a(s2,k) = - a(s1,k-1) * inv_idx(pk+1)
             d = d + a(s2,k) * ndu(r,pk)
          end if
          bsdx(k,r) = d
          j  = s1
          s1 = s2
          s2 = j
       end do
    end do

    ! Multiply result by correct factors:
    ! deg!/(deg-n)! = deg*(deg-1)*...*(deg-n+1)
    r = deg
    do k = 1, n
       bsdx(k,:) = bsdx(k,:) * r
       r = r * (deg-k)
    end do

    end associate

  end subroutine sll_s_uniform_bsplines_eval_basis_and_n_derivs

  !-----------------------------------------------------------------------------
  !> @brief
  !> Evaluate uniform periodic spline curve defined by coefficients scoef at 
  !> knots (which are the grid points) 
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_eval_uniform_periodic_spline_curve( degree, scoef, sval )

    integer , intent(in   ) :: degree
    real(wp), intent(in   ) :: scoef(:) 
    real(wp), intent(  out) :: sval (:) 
    
    real(wp) :: bspl(degree+1)
    real(wp) :: val 
    integer  :: i, j, imj, n

    ! get bspline values at knots
    call sll_s_uniform_bsplines_eval_basis(degree, 0.0_wp, bspl)
    n = size(scoef)
    do i= 1, n
       val = 0.0_wp
       do j=1, degree
          imj = mod(i-1-j+n,n) + 1 
          val = val + bspl(j)*scoef(imj)
          !print*, i,j, imj, bspl(j),val
       enddo
       sval(i) = val 
    end do

  end subroutine sll_s_eval_uniform_periodic_spline_curve

end module sll_m_low_level_bsplines
