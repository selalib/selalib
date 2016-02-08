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
!> It is a selalib implemnation of the classical algorithms
!> found in the de Boor book "Practical guide to splines"
!> or the NURBS book by Piegl and Tiller.  
module sll_m_arbitrary_degree_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

implicit none

  public :: &
    sll_t_arbitrary_degree_spline_1d, &
    sll_f_spline_derivatives_at_x, &
    sll_f_splines_and_derivs_at_x, &
    sll_f_splines_at_x, &
    sll_s_compute_b_spline_and_deriv_at_x_mm, &
    sll_s_compute_b_spline_at_x_mm, &
    sll_f_eval_uniform_periodic_spline_curve, &
    sll_f_find_cell, &
    sll_f_new_arbitrary_degree_spline_1d, &
    sll_p_open_arbitrary_deg_spline, &
    sll_p_periodic_arbitrary_deg_spline, &
    sll_o_delete, &
    sll_s_uniform_b_spline_derivatives_at_x, &
    sll_s_uniform_b_splines_and_derivs_at_x, &
    sll_s_uniform_b_splines_at_x

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! With these read-only parameters, we mimic the behavior of an enumerator
  ! or object macros. These are the alternatives with a different 
  ! implementation.
  ! This should be in sll_m_boundary_condition_descriptors ?
  !> Spline constructed from a grid with periodic boundary conditions
  sll_int32, parameter :: sll_p_periodic_arbitrary_deg_spline = 0
  !> Spline constructed from a grid with open boundary conditions
  !> This means that the end points are duplicated so as to have multiplicity
  !> degree + 1 in the knot sequence, which are otherwise the grid points
  sll_int32, parameter :: sll_p_open_arbitrary_deg_spline     = 1 

  type :: sll_t_arbitrary_degree_spline_1d
     sll_int32  :: num_pts
     sll_int32  :: bc_type
     sll_int32  :: degree
     sll_real64 :: xmin
     sll_real64 :: xmax
     sll_real64, dimension(:), pointer :: knots   ! knots array
     sll_real64, dimension(:), pointer :: left    ! needed for evaluation
     sll_real64, dimension(:), pointer :: right   ! needed for evaluation
     sll_real64, dimension(:,:), pointer :: ndu   ! needed for derivative
  end type sll_t_arbitrary_degree_spline_1d

  interface sll_o_delete
     module procedure delete_arbitrary_order_spline_1d
  end interface sll_o_delete


contains
  !> @brief
  !> Build new sll_t_arbitrary_degree_spline_1d object 
  !> @details
  !> based on a on a grid of 
  !> strictly increasing points including the last point also for periodic
  !> domains.
  function sll_f_new_arbitrary_degree_spline_1d( degree, grid, num_pts, bc_type )
    type(sll_t_arbitrary_degree_spline_1d), pointer :: sll_f_new_arbitrary_degree_spline_1d
    sll_int32, intent(in)                    :: degree !< spline degree
    sll_real64, dimension(:), intent(in)     :: grid   !< grid points 
    sll_int32,  intent(in)                   :: num_pts !< number of grid points
    sll_int32,  intent(in)                   :: bc_type !< boundary condition
    ! local variables
    sll_int32                                :: i
    sll_int32                                :: ierr
    sll_real64                               :: period ! length of period
    sll_real64                               :: gridmin

    if( size(grid) < num_pts ) then
       print *, 'ERROR. sll_f_new_arbitrary_degree_spline_1d(): ', &
            'size of given grid array is smaller than the stated ', &
            'number of points.'
       print *, 'size(grid) = ', size(grid), 'num_pts = ', num_pts
       STOP
    end if
    if( degree < 1 ) then
       print *, 'ERROR. sll_f_new_arbitrary_degree_spline_1d(): ', &
            'only strictly positive integer values for degree are allowed, ', &
            'given: ', degree
    end if
    ! Check whether grid points are in strictly increasing order
    gridmin = grid(2)-grid(1)
    do i=2, num_pts
       gridmin = min(gridmin, grid(i)-grid(i-1)) 
    end do
    if (gridmin < 1.d-10) then
       print*, 'ERROR. sll_f_new_arbitrary_degree_spline_1d(): ', &
            ' we require that grid points are in strictly increasing', &
            ' order and that each cell length is at least 1.e-10.'
    end if
    select case( bc_type )
    case(:-1)
       print *, 'ERROR. sll_f_new_arbitrary_degree_spline_1d(): ', &
            'invalid boundary condition type supplied.'
       STOP
    case(2:)
       print *, 'ERROR. sll_f_new_arbitrary_degree_spline_1d(): ', &
            'invalid boundary condition type supplied.'
       STOP
    end select

    SLL_ALLOCATE( sll_f_new_arbitrary_degree_spline_1d, ierr )
    sll_f_new_arbitrary_degree_spline_1d%num_pts  = num_pts
    sll_f_new_arbitrary_degree_spline_1d%bc_type  = bc_type
    sll_f_new_arbitrary_degree_spline_1d%degree   = degree
    sll_f_new_arbitrary_degree_spline_1d%xmin     = grid(1)
    sll_f_new_arbitrary_degree_spline_1d%xmax     = grid(num_pts)

    ! 'period' is useful in the case of periodic boundary conditions.
    period = grid(num_pts) - grid(1)

    ! Create the knots array from the grid points. Here take the grid points
    ! as konts and simply add to the left and the right the
    ! amount of knots that depends on the degree of the requested 
    ! spline. We aim at setting up the indexing in such a way that the 
    ! original indexing of 'grid' is preserved, i.e.: grid(i) = knot(i), at
    ! least whenever the scope of the indices defined here is active.
    SLL_ALLOCATE(sll_f_new_arbitrary_degree_spline_1d%knots(1-degree:num_pts+degree),ierr)
    do i=1,num_pts
       sll_f_new_arbitrary_degree_spline_1d%knots(i) = grid(i)
    end do

    ! Allocate work arrays for evaluation
    SLL_ALLOCATE(sll_f_new_arbitrary_degree_spline_1d%left(degree),ierr)
    SLL_ALLOCATE(sll_f_new_arbitrary_degree_spline_1d%right(degree),ierr)
    SLL_ALLOCATE(sll_f_new_arbitrary_degree_spline_1d%ndu(0:degree,0:degree),ierr)

    ! Fill out the extra points at both ends of the local knot array with
    ! values proper to the boundary condition requested.
    select case( bc_type )
    case(sll_p_periodic_arbitrary_deg_spline)
       ! The logic behind the periodic boundary condition is the following.
       ! The given grid array has minimum (grid(1)) and maximum (grid(n))
       ! values at either end. This defines a length 'L'. If interpreted 
       ! as a periodic space, this is also the period. Thus, as we extend
       ! the number of knots at both ends of the given array, we use the
       ! periodicity condition to fill out the new values:
       !
       !                    .
       !                    .
       !                    .
       !           knots(-1) = knots(n-2) - L
       !           knots( 0) = knots(n-1) - L
       !                    .
       !                    .
       !                    .
       !           knots(n+1) = knots(1) + L
       !           knots(n+2) = knots(2) + L
       !                    .
       !                    .
       !                    .
       !

       ! Check that there are enough grid points for a given degree
       if (num_pts < degree) then
          print *, 'ERROR. sll_f_new_arbitrary_degree_spline_1d(): ', &
               'new at leas ', degree, 'grid points for a periodic spline', &
               'of degree', degree
          stop 
       end if
       do i=1,degree
          ! Fill out the extra nodes on the left
          sll_f_new_arbitrary_degree_spline_1d%knots(1-i) = &
               grid(num_pts-i) - period
          ! Fill out the extra nodes on the right
          sll_f_new_arbitrary_degree_spline_1d%knots(num_pts+i) = &
               grid(i+1) + period
       end do

    case(sll_p_open_arbitrary_deg_spline)
       ! The 'open' boundary condition simply extends the new values
       ! of the local array at both ends with repeated endpoint values.
       ! That is
       !
       !     ... = knots(-2) = knots(-1) = knots(0) = knots(1)
       !
       ! and
       !
       !    knots(n+1) = knots(n+2) = knots(n+3) = ... =  knots(n)
       do i=1-degree,0
          sll_f_new_arbitrary_degree_spline_1d%knots(i) = grid(1)
       end do
       do i=num_pts+1,num_pts+degree
          sll_f_new_arbitrary_degree_spline_1d%knots(i) = grid(num_pts)
       end do
    end select
  end function sll_f_new_arbitrary_degree_spline_1d

  !>@brief find cell returns the index i of the grid cell such that:
  !> spline_obj%knots(i) <= x <= spline_obj%knots(i+1).
  !>
  !>@detail
  !> If x is not between spline_obj%knots(1) and 
  !> spline_obj%knots(spline_obj%num_pts),  then the value -1 is returned.
  function sll_f_find_cell( spline_obj, x)
    type(sll_t_arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_real64, intent(in) :: x
    sll_int32 :: sll_f_find_cell
    sll_int32 :: low
    sll_int32 :: high
    sll_int32 :: n

    n = spline_obj%num_pts

    ! check if point is outside of grid
    if (x > spline_obj%knots(n)) then
       sll_f_find_cell = -1
       return
    end if
    if (x < spline_obj%knots(1)) then
       sll_f_find_cell = -1
       return
    end if
    ! check is point is exactly on right boundary
    if (x == spline_obj%knots(n)) then
       sll_f_find_cell = n-1
       !print*, 'sll_f_find_cell=', sll_f_find_cell
       return
    end if

    low  = 1
    high = n
    sll_f_find_cell = (low + high) / 2
    do while (x < spline_obj%knots(sll_f_find_cell) &
         .or. x >= spline_obj%knots(sll_f_find_cell+1))
       if (x < spline_obj%knots(sll_f_find_cell)) then
          high = sll_f_find_cell
       else
          low  = sll_f_find_cell
       end if
       sll_f_find_cell = (low + high) / 2
    end do
  end function sll_f_find_cell

  ! sll_f_splines_at_x() returns the values of all the B-splines of a given 
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
  !> sll_f_splines_at_x( spline_obj, cell, x ) computes the values of all the
  !> splines which have support in 'cell' and evaluates them at point 'x',
  !> which is supposed to be in cell. The spline object should have already
  !> been initialized and will contain information on the spline degree
  !> to use and the type of boundary condition desired.
  !> The algorithm implemented is numerically stable and known as The
  !> Cox - de Boor algorithm, which is a generalisation to splines of the
  !> de Casteljau algorithm for Bezier curves.
  !> @return b_spline_at_x B-spline values
  function sll_f_splines_at_x( spline_obj, icell, x )
    type(sll_t_arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_int32, intent(in)                          :: icell
    sll_real64, intent(in)                         :: x
    sll_real64, dimension(0:spline_obj%degree)     :: sll_f_splines_at_x
    ! local variables
    sll_int32                                      :: deg
    sll_real64                                     :: saved
    sll_real64                                     :: temp
    sll_int32                                      :: j
    sll_int32                                      :: r

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(spline_obj))
    SLL_ASSERT(x >= spline_obj%xmin)
    SLL_ASSERT(x <= spline_obj%xmax)
    SLL_ASSERT(icell >= 1)
    SLL_ASSERT(icell <= spline_obj%num_pts - 1)
    ! This is checked always. If it appear to penalise too much the performances
    ! one could consider replacing this by ASSERT
    if( .not. ((x >= spline_obj%knots(icell)) &
         .and.(x <= spline_obj%knots(icell+1)))) then
       print *, 'ERROR. sll_f_splines_at_x(): the given value of x is not ', &
            'inside the specified cell.'
       STOP
    end if

    deg = spline_obj%degree

    sll_f_splines_at_x(0) = 1.0_f64
    do j = 1, deg
       spline_obj%left(j)  = x - spline_obj%knots(icell+1-j)
       spline_obj%right(j) = spline_obj%knots(icell+j) - x
       saved = 0.0_f64
       do r = 0, j-1
          temp = sll_f_splines_at_x(r) / (spline_obj%right(r+1) &
               + spline_obj%left(j-r))
          sll_f_splines_at_x(r) = saved + spline_obj%right(r+1) * temp
          saved = spline_obj%left(j-r) * temp
       end do
       sll_f_splines_at_x(j) = saved
    end do

  end function sll_f_splines_at_x

  !> @brief
  !> returns first derivative values at x of all b-splines with support in cell
  !> @details
  !> sll_f_spline_derivatives_at_x() returns an array with the derivative values of 
  !> the B-splines of a requested order that are supported in 'cell' and 
  !> evaluated at 'x'. 
  !> Algorithm derived from algorithm A3.2 of NURBS book 
  !> The return value has the format:
  !> \f[
  !> B'[deg,i-deg](x), B'[deg,i-deg+1](x), ..., B'[deg,i](x)
  !> \f]
  !> where 'deg' is the degree of the spline.
  !> @return sll_f_spline_derivatives_at_x B-spline derivatives
  function sll_f_spline_derivatives_at_x( spline_obj, icell, x ) result(bsdx)
    type(sll_t_arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_int32, intent(in)                          :: icell
    sll_real64, intent(in)                         :: x
    sll_real64, dimension(0:spline_obj%degree)     :: bsdx !sll_f_spline_derivatives_at_x
    ! local variables
    sll_int32                                      :: deg
    sll_int32                                      :: num_pts
    sll_int32                                      :: r
    sll_int32                                      :: j
    sll_real64                                     :: saved
    sll_real64                                     :: temp
    sll_real64                                     :: rdeg

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(spline_obj))
    SLL_ASSERT(x >= spline_obj%xmin)
    SLL_ASSERT(x <= spline_obj%xmax)
    SLL_ASSERT(icell >= 1)
    SLL_ASSERT(icell <= spline_obj%num_pts - 1)
    ! This is checked always.
    if( .not. (x >= spline_obj%knots(icell)) &
         .and. (x <= spline_obj%knots(icell+1))) then
       print *, 'ERROR. sll_f_splines_at_x(): the given value of x is not ', &
            'inside the specified cell.'
       STOP
    end if
    deg = spline_obj%degree
    rdeg = real(deg,f64)
    num_pts = spline_obj%num_pts

    ! compute nonzero basis functions and knot differences
    ! for splines up to degree deg-1 which are needed to compute derivative
    ! First part of Algorithm  A3.2 of NURBS book 
    bsdx(0) = 1.0_f64
    do j = 1, deg - 1
       spline_obj%left(j)  = x - spline_obj%knots(icell+1-j)
       spline_obj%right(j) = spline_obj%knots(icell+j) - x
       saved = 0.0_f64
       do r = 0, j-1
          ! compute and save bspline values
          temp = bsdx(r)/(spline_obj%right(r+1) + spline_obj%left(j-r))
          bsdx(r) = saved + spline_obj%right(r+1) * temp
          saved = spline_obj%left(j-r) * temp
       end do
       bsdx(j) = saved
    end do
    ! Compute derivatives at x using values stored in ndu and formula
    ! formula for spline derivative based on difference of splines of 
    ! degree deg-1
    ! -------
    ! j = 0
    saved = rdeg *bsdx(0) / &
         (spline_obj%knots(icell+1) - spline_obj%knots(icell+1-deg)) 
    bsdx(0) = -saved
    do j = 1, deg-1
       temp = saved 
       saved =  rdeg*bsdx(j) / &
            (spline_obj%knots(icell+j+1) - spline_obj%knots(icell+j+1-deg))
       bsdx(j) = temp - saved
    end do
    ! j = deg
    bsdx(deg) =  saved

  end function sll_f_spline_derivatives_at_x


  !> @brief 
  !> returns splines and first derivatives
  !> @details
  !> See sll_f_spline_derivatives_at_x and  sll_f_splines_at_x
  !> @return b_spline_and_derivs_at_x B-spline values and derivatives
  function sll_f_splines_and_derivs_at_x( spline_obj, icell, x ) result(bsdx)
    type(sll_t_arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_int32, intent(in)                          :: icell
    sll_real64, intent(in)                         :: x
    sll_real64, dimension(2,0:spline_obj%degree)   :: bsdx
    ! local variables
    sll_int32                                      :: deg
    sll_int32                                      :: num_pts
    sll_int32                                      :: r
    sll_int32                                      :: j
    sll_real64                                     :: saved
    sll_real64                                     :: temp
    sll_real64                                     :: rdeg

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(spline_obj))
    SLL_ASSERT(x >= spline_obj%xmin)
    SLL_ASSERT(x <= spline_obj%xmax)
    SLL_ASSERT(icell >= 1)
    SLL_ASSERT(icell <= spline_obj%num_pts - 1)
    ! This is checked always.
    if( .not. (x >= spline_obj%knots(icell)) &
         .and. (x <= spline_obj%knots(icell+1))) then
       print *, 'ERROR. sll_f_splines_at_x(): the given value of x is not ', &
            'inside the specified cell.'
       STOP
    end if
    deg = spline_obj%degree
    num_pts = spline_obj%num_pts
    rdeg = real(deg,f64) 

    ! compute nonzero basis functions and knot differences
    ! for splines up to degree deg-1 which are needed to compute derivative
    ! First part of Algorithm  A3.2 of NURBS book 
    bsdx(1,0) = 1.0_f64
    do j = 1, deg - 1
       spline_obj%left(j)  = x - spline_obj%knots(icell+1-j)
       spline_obj%right(j) = spline_obj%knots(icell+j) - x
       saved = 0.0_f64
       do r = 0, j-1
          ! compute and save knot differences
          temp = bsdx(1,r)/(spline_obj%right(r+1) + spline_obj%left(j-r))
          bsdx(1,r) = saved + spline_obj%right(r+1) * temp
          saved = spline_obj%left(j-r) * temp
       end do
       bsdx(1,j) = saved
    end do

    ! Compute derivatives at x using values stored in ndu and formula
    ! formula for spline derivative based on difference of splines of 
    ! degree deg-1
    ! -------
    ! j = 0
    saved = rdeg * bsdx(1,0) / &
         (spline_obj%knots(icell+1) - spline_obj%knots(icell+1-deg)) 
    bsdx(2,0) = -saved
    do j = 1, deg-1
       temp = saved 
       saved =  rdeg*bsdx(1,j) / &
            (spline_obj%knots(icell+j+1) - spline_obj%knots(icell+j+1-deg))
       bsdx(2,j) = temp - saved
    end do
    ! j = deg
    bsdx(2,deg) =  saved  
    ! Compute values of splines of degree deg
    !----------------------------------------
    j = deg
    spline_obj%left(j)  = x - spline_obj%knots(icell+1-j)
    spline_obj%right(j) = spline_obj%knots(icell+j) - x
    saved = 0.0_f64
    do r = 0, j-1
       ! compute and save knot differences
       temp = bsdx(1,r)/(spline_obj%right(r+1) + spline_obj%left(j-r))
       bsdx(1,r) = saved + spline_obj%right(r+1) * temp
       saved = spline_obj%left(j-r) * temp
    end do
    bsdx(1,j) = saved 
  end function sll_f_splines_and_derivs_at_x

  !> @brief Alternative direct implentation of recursion formula. 
  !>
  !> @detail
  !> This provides an evaluation of B-splines directly based on the recurrence
  !> formula. It is 10% faster than the classical Cox - de Boor formula 
  !> that is implented in sll_f_splines_at_x, but can have numerical stability issues.
  !>For this reason the Cox - de Boor formula should be the default implementation 
  subroutine sll_s_compute_b_spline_at_x_mm( &
       knots, &
       cell, &
       x, &
       degree, &
       out)
    sll_int32, intent(in) :: degree
    sll_real64, dimension(1-degree:), intent(in) :: knots
    sll_int32, intent(in) :: cell
    sll_real64, intent(in) :: x
    sll_real64, dimension(:), intent(out) :: out

    sll_real64 :: tmp1
    sll_real64 :: tmp2
    sll_int32 :: ell
    sll_int32 :: k

    out(1) = 1._f64
    do ell=1,degree
       tmp1 = (x-knots(cell+1-ell))/(knots(cell+1)-knots(cell+1-ell))*out(1)
       out(1) = out(1) -tmp1
       do k=2,ell
          tmp2 = (x-knots(cell+k-ell))/(knots(cell+k)-knots(cell+k-ell))*out(k)
          out(k) = out(k)+tmp1-tmp2
          tmp1 = tmp2
       enddo
       out(ell+1) = tmp1
    enddo

  end subroutine sll_s_compute_b_spline_at_x_mm

  !> @brief Alternative direct implentation of recursion formula. 
  !>
  !> @detail
  !> This provides an evaluation of B-splines directly based on the recurrence
  !> formula. It is about 80% faster than the classical Cox - de Boor formula 
  !> that is implented in sll_f_splines_at_x, but can have numerical stability issues.
  !> For this reason the Cox - de Boor formula should be the default implementation
  subroutine sll_s_compute_b_spline_and_deriv_at_x_mm( &
       knots, &
       cell, &
       x, &
       degree, &
       out)
    sll_int32, intent(in) :: degree
    sll_real64, dimension(1-degree:), intent(in) :: knots
    sll_int32, intent(in) :: cell
    sll_real64, intent(in) :: x
    sll_real64, dimension(:,:), intent(out) :: out

    sll_real64 :: tmp1
    sll_real64 :: tmp2
    sll_int32 :: ell
    sll_int32 :: k

    out(1,1) = 1._f64
    do ell=1,degree
       tmp1 = (x-knots(cell+1-ell))/(knots(cell+1)-knots(cell+1-ell))*out(1,1)
       out(1,1) = out(1,1) -tmp1
       do k=2,ell
          tmp2 = (x-knots(cell+k-ell))/(knots(cell+k)-knots(cell+k-ell))*out(1,k)
          out(1,k) = out(1,k)+tmp1-tmp2
          tmp1 = tmp2
       enddo
       out(2,ell+1) = tmp1
       if(ell==degree-1)then
          !compute the derivatives
          tmp1 = real(degree,f64)/(knots(cell+1)-knots(cell+1-degree))*out(1,1)
          out(2,1) = -tmp1
          do k=2,degree
             out(2,k) = tmp1
             tmp1 = real(degree,f64)/(knots(cell+k)-knots(cell+k-degree))*out(2,k)
             out(2,k) = out(2,k)-tmp1
          enddo
          out(2,degree+1) = tmp1
       endif
    enddo
  end subroutine sll_s_compute_b_spline_and_deriv_at_x_mm

  subroutine delete_arbitrary_order_spline_1d( spline )
    type(sll_t_arbitrary_degree_spline_1d), pointer :: spline
    sll_int32                    :: ierr
    if( .not. associated(spline) ) then
       print *, 'ERROR. delete_arbitrary_order_spline_1d(): given spline ', &
            'pointer is not associated.'
       STOP
    end if
    SLL_DEALLOCATE( spline%knots, ierr )
    SLL_DEALLOCATE( spline, ierr )
  end subroutine delete_arbitrary_order_spline_1d

  ! *************************************************************************
  !
  !                    UNIFORM B-SPLINE FUNCTIONS
  !
  ! *************************************************************************

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
  subroutine sll_s_uniform_b_splines_at_x( spline_degree,     &
                                           normalized_offset, &
                                           bspl               )

    implicit none
    sll_int32, intent(in)                      :: spline_degree
    sll_real64, intent(in)                     :: normalized_offset
    sll_real64, dimension(0:spline_degree)     :: bspl
    ! local variables
    sll_real64                                 :: inv_j
    sll_real64                                 :: x,xx
    sll_real64                                 :: j_real
    sll_int32                                  :: j, r
    sll_real64                                 :: temp
    sll_real64                                 :: saved

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_f64 )
    SLL_ASSERT( normalized_offset <= 1.0_f64 )

    x = normalized_offset
    bspl(0) = 1.0_f64
    do j = 1, spline_degree
       saved = 0.0_f64
       j_real = real(j,f64)
       inv_j = 1.0_f64 / j_real
       xx = - x
       do r = 0, j-1
          xx = xx + 1.0_f64
          temp = bspl(r) * inv_j
          bspl(r) = saved + xx * temp
          saved = (j_real - xx) * temp
       end do
       bspl(j) = saved
    end do
  end subroutine sll_s_uniform_b_splines_at_x

  !> @brief Evaluate all derivatives of non vanishing uniform B-Splines 
  !> in unit cell. 
  !>
  !> @detail
  !> Returns an array with the values of the b-spline derivatives of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, hence the results must be divided by the
  !> real cell size to scale back the results.
  subroutine sll_s_uniform_b_spline_derivatives_at_x( spline_degree,     &
                                                      normalized_offset, & 
                                                      bspl               )

    sll_int32, intent(in)                      :: spline_degree
    sll_real64, dimension(0:spline_degree)     :: bspl
    sll_real64, intent(in)                     :: normalized_offset
    ! local variables   
    sll_real64                                 :: inv_j
    sll_real64                                 :: x,xx
    sll_real64                                 :: j_real
    sll_int32                                  :: j, r
    sll_real64                                 :: temp
    sll_real64                                 :: saved
    sll_real64                                 :: bj,bjm1 

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_f64 )
    SLL_ASSERT( normalized_offset <= 1.0_f64 )

    x = normalized_offset
    bspl(0) = 1.0_f64
    ! only need splines of lower degree to compute derivatives
    do j = 1, spline_degree - 1
       saved = 0.0_f64
       j_real = real(j,f64)
       inv_j = 1.0_f64 / j_real
       xx = - x
       do r = 0, j-1
          xx = xx + 1.0_f64
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
  end subroutine sll_s_uniform_b_spline_derivatives_at_x

  !> @brief Evaluate all values and derivatives of non vanishing uniform B-Splines 
  !> in unit cell. 
  !>
  !> @detail
  !> returns an array with the values of the b-spline derivatives of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, hence the results must be divided by the
  !> real cell size to scale back the results.
  subroutine sll_s_uniform_b_splines_and_derivs_at_x( degree,            &
                                                      normalized_offset, &
                                                      bspl               )

    sll_int32, intent(in)                 :: degree
    sll_real64, dimension(2,0:degree)     :: bspl
    sll_real64, intent(in)                :: normalized_offset
    ! local variables
    sll_real64                            :: inv_j
    sll_real64                            :: xx
    sll_real64                            :: j_real
    sll_int32                             :: j, r
    sll_real64                            :: temp
    sll_real64                            :: saved

    SLL_ASSERT( degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_f64 )
    SLL_ASSERT( normalized_offset <= 1.0_f64 )

    bspl(1,0) = 1.0_f64
    ! compute splines up to degree spline_degree -1
    do j = 1, degree - 1
       saved = 0.0_f64
       j_real = real(j,f64)
       inv_j = 1.0_f64 / j_real
       xx = - normalized_offset
       do r = 0, j-1
          xx = xx + 1.0_f64
          temp = bspl(1,r) * inv_j
          bspl(1,r) = saved + xx * temp
          saved = (j_real - xx) * temp
       end do
       bspl(1,j) = saved
    end do
    ! compute derivatives
    bspl(2,0) = -bspl(1,0)
    do j=1, degree - 1
       bspl(2,j) = bspl(1,j-1) - bspl(1,j)
    end do
    bspl(2,degree) = bspl(1,degree-1)
    ! continue the Cox-De Boor algorithm to evaluate splines
    j= degree
    saved = 0.0_f64
    j_real = real(j,f64)
    inv_j = 1.0_f64 / j_real
    xx = - normalized_offset
    do r = 0, j-1
       xx = xx + 1.0_f64
       temp = bspl(1,r) * inv_j
       bspl(1,r) = saved + xx * temp
       saved = (j_real - xx) * temp
    end do
    bspl(1,j) = saved

  end subroutine sll_s_uniform_b_splines_and_derivs_at_x

  !> @brief
  !> Evaluate uniform periodic spline curve defined by coefficients scoef at 
  !> knots (which are the grid points) 
  function sll_f_eval_uniform_periodic_spline_curve(degree, scoef) result(sval)
    sll_int32  :: degree   ! spline degree
    !sll_int32  :: npoints  ! number of points where spline is evaluated
    sll_real64 :: scoef(:) 
    sll_real64, allocatable :: sval(:) 
    ! local variables
    sll_real64 :: bspl(degree+1)
    sll_real64 :: val 
    sll_int32 :: i,j, imj, ierr, n

    ! get bspline values at knots
    call sll_s_uniform_b_splines_at_x(degree, 0.0_f64, bspl)
    n = size(scoef)
    SLL_ALLOCATE(sval(n), ierr)
    do i= 1, n
       val = 0.0_f64
       do j=1, degree
          imj = mod(i-1-j+n,n) + 1 
          val = val + bspl(j)*scoef(imj)
          !print*, i,j, imj, bspl(j),val
       enddo
       sval(i) = val 
    end do
  end function sll_f_eval_uniform_periodic_spline_curve

end module sll_m_arbitrary_degree_splines
