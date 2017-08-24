!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_bsplines
!>
!> @author Eric Sonnendrücker - IPP Garching
!> @author Yaman Güçlü        - IPP Garching
!> @author Edoardo Zoni       - IPP Garching

module sll_m_spline_2d_non_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
  sll_p_periodic, &
  sll_p_hermite, &
  sll_p_greville

use sll_m_bsplines, only: &
  sll_f_find_cell, &
  sll_s_bsplines_eval_basis, &
  sll_s_bsplines_eval_deriv

use schur_complement, only: &
  schur_complement_solver, &
  schur_complement_fac   , &
  schur_complement_slv   , &
  schur_complement_free

use sll_m_spline_1d_non_uniform, only: &
  sll_t_spline_1d_non_uniform

implicit none

public :: &
  sll_t_spline_2d_boundary_data,                   &
  sll_t_spline_2d_non_uniform,                     &
  sll_s_spline_2d_non_uniform_init,                &
  sll_s_spline_2d_non_uniform_free,                &
  sll_s_spline_2d_non_uniform_compute_interpolant, &
  sll_f_spline_2d_non_uniform_eval,                & ! scalar functions for evaluation
  sll_f_spline_2d_non_uniform_eval_deriv_x1,       &
  sll_f_spline_2d_non_uniform_eval_deriv_x2,       &
  sll_s_spline_2d_non_uniform_eval_array,          & ! vector subroutines for evaluation
  sll_s_spline_2d_non_uniform_eval_array_deriv_x1, &
  sll_s_spline_2d_non_uniform_eval_array_deriv_x2

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief
!> basic type for two-dimensional B-spline data.
!> @details
!> treated as an opaque type. No access to its internals is directly allowed.
type :: sll_t_spline_2d_non_uniform

  integer ::     deg(2)
  integer ::       n(2)
  integer :: bc_xmin(2)
  integer :: bc_xmax(2)

  type(sll_t_spline_1d_non_uniform) :: bs1
  type(sll_t_spline_1d_non_uniform) :: bs2
  sll_real64, allocatable           :: bcoef(:,:)

  sll_real64, allocatable, private  :: bwork(:,:)
  integer                , private  :: nbc_xmin(2)
  integer                , private  :: nbc_xmax(2)

end type sll_t_spline_2d_non_uniform


!> Container for boundary condition data
!>
!>  x2_max  ____________
!>         |            |
!>         | c        d |
!>         |            |
!>         |            |
!>         |            |
!>         | a        b |
!>  x2_min |____________|
!>       x1_min       x1_max
!>
type :: sll_t_spline_2d_boundary_data
  sll_real64, allocatable :: derivs_x1_min (:,:)
  sll_real64, allocatable :: derivs_x1_max (:,:)
  sll_real64, allocatable :: derivs_x2_min (:,:)
  sll_real64, allocatable :: derivs_x2_max (:,:)
  sll_real64, allocatable :: mixed_derivs_a(:,:)
  sll_real64, allocatable :: mixed_derivs_b(:,:)
  sll_real64, allocatable :: mixed_derivs_c(:,:)
  sll_real64, allocatable :: mixed_derivs_d(:,:)
end type sll_t_spline_2d_boundary_data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief Initialize a 2D spline interpolation object
  !-----------------------------------------------------------------------------
  subroutine sll_s_spline_2d_non_uniform_init( &
    self   , &
    degree , &
    breaks1, &
    breaks2, &
    bc_xmin, &
    bc_xmax )

    type(sll_t_spline_2d_non_uniform), intent(  out) :: self
    sll_int32                        , intent(in   ) :: degree (2)
    sll_real64                       , intent(in   ) :: breaks1(:)
    sll_real64                       , intent(in   ) :: breaks2(:)
    sll_int32                        , intent(in   ) :: bc_xmin(2)
    sll_int32                        , intent(in   ) :: bc_xmax(2)

    call self%bs1%init( degree(1), breaks1, bc_xmin(1), bc_xmax(1) )
    call self%bs2%init( degree(2), breaks2, bc_xmin(2), bc_xmax(2) )

    associate( n1 => self%bs1%n, &
               n2 => self%bs2%n, &
               g1 => merge( 1+degree(1)/2, 0, bc_xmin(1) == sll_p_periodic ), &
               g2 => merge( 1+degree(2)/2, 0, bc_xmin(2) == sll_p_periodic ) )

      ! Save data into type
      self%deg     = degree
      self%n       = [n1,n2]
      self%bc_xmin = bc_xmin
      self%bc_xmax = bc_xmax

      ! Allocate array of spline coefficients, and work array for interpolation.
      ! in case of periodic BCs, a larger array of coefficients is used in order
      ! to avoid a loop with calls to the "mod( , )" function at evaluation.
      allocate( self%bcoef(1-g1:n1+g1,1-g2:n2+g2) );  self%bcoef = 0.0_f64
      allocate( self%bwork(1-g2:n2+g2,1-g1:n1+g1) );  self%bwork = 0.0_f64

      ! Calculate number of additional boundary data on each side of domain
      ! (i.e. derivatives for Hermite BC)
      self%nbc_xmin = merge( degree/2, 0, bc_xmin == sll_p_hermite )
      self%nbc_xmax = merge( degree/2, 0, bc_xmax == sll_p_hermite )

    end associate

  end subroutine sll_s_spline_2d_non_uniform_init

  !-----------------------------------------------------------------------------
  subroutine sll_s_spline_2d_non_uniform_compute_interpolant( &
    self, &
    gtau, &
    boundary_data )

    type(sll_t_spline_2d_non_uniform)  , intent(inout)           :: self
    sll_real64                         , intent(in   )           :: gtau(:,:)
    type(sll_t_spline_2d_boundary_data), intent(in   ), optional :: boundary_data

    character(len=*), parameter :: &
      this_sub_name = "sll_s_spline_2d_non_uniform_compute_interpolant"

    sll_int32 :: i1, i2

    ! User must provide derivatives at boundary in case of Hermite BC
    if (.not. present(boundary_data)) then
      if (self%bc_xmin(1) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x1_min requires derivatives" )
      else if (self%bc_xmax(1) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x1_max requires derivatives" )
      else if (self%bc_xmin(2) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x2_min requires derivatives" )
      else if (self%bc_xmax(2) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x2_max requires derivatives" )
      end if
    end if

    associate( w  => self%bcoef       , &
               wt => self%bwork       , &
               n1 => self%n(1)        , &
               a1 => self%nbc_xmin(1) , &
               b1 => self%nbc_xmax(1) , &
               n2 => self%n(2)        , &
               a2 => self%nbc_xmin(2) , &
               b2 => self%nbc_xmax(2) )

    SLL_ASSERT( all( shape(gtau) == [n1-a1-b1,n2-a2-b2] ) )

    ! DEBUG mode, Hermite boundary conditions:
    ! Verify that required arrays are passed and allocated with correct shape
    if (a1 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x1_min )  )
      SLL_ASSERT(      size( boundary_data % derivs_x1_min, 1 ) == a1       )
      SLL_ASSERT(      size( boundary_data % derivs_x1_min, 2 ) == n2-a2-b2 )
    end if

    if (b1 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x1_max )  )
      SLL_ASSERT(      size( boundary_data % derivs_x1_max, 1 ) == b1       )
      SLL_ASSERT(      size( boundary_data % derivs_x1_max, 2 ) == n2-a2-b2 )
    end if

    if (a2 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x2_min )  )
      SLL_ASSERT(      size( boundary_data % derivs_x2_min, 1 ) == a2       )
      SLL_ASSERT(      size( boundary_data % derivs_x2_min, 2 ) == n1-a1-b1 )
    end if

    if (b2 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x2_max )  )
      SLL_ASSERT(      size( boundary_data % derivs_x2_max, 1 ) == b2       )
      SLL_ASSERT(      size( boundary_data % derivs_x2_max, 2 ) == n1-a1-b1 )
    end if

    if (a1 > 0 .and. a2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_a )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_a, 1 ) == a1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_a, 2 ) == a2 )
    end if

    if (b1 > 0 .and. a2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_b )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_b, 1 ) == b1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_b, 2 ) == a2 )
    end if

    if (a1 > 0 .and. b2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_c )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_c, 1 ) == a1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_c, 2 ) == b2 )
    end if

    if (b1 > 0 .and. b2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_d )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_d, 1 ) == b1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_d, 2 ) == b2 )
    end if

    ! Copy interpolation data onto w array
    w(1+a1:n1-b1,1+a2:n2-b2) = gtau(:,:)

    ! Hermite BCs: store boundary data in appropriate chunk of w array
    if (present( boundary_data )) then

      if (a1 > 0) w(      1:a1,1+a2:n2-b2) = boundary_data % derivs_x1_min(1:a1,:)
      if (b1 > 0) w(n1-b1+1:n1,1+a2:n2-b2) = boundary_data % derivs_x1_max(1:b1,:)

      if (a2 > 0) w(1+a1:n1-b1,      1:a2) = transpose( boundary_data % derivs_x2_min(1:a2,:) )
      if (b2 > 0) w(1+a1:n1-b1,n2-b2+1:n2) = transpose( boundary_data % derivs_x2_max(1:b2,:) )

      if (a1 > 0 .and. a2 > 0) w(      1:a1,      1:a2) = boundary_data % mixed_derivs_a(:,:)
      if (b1 > 0 .and. a2 > 0) w(n1-b1+1:n1,      1:a2) = boundary_data % mixed_derivs_b(:,:)
      if (a1 > 0 .and. b2 > 0) w(      1:a1,n2-b2+1:n2) = boundary_data % mixed_derivs_c(:,:)
      if (b1 > 0 .and. b2 > 0) w(n1-b1+1:n1,n2-b2+1:n2) = boundary_data % mixed_derivs_d(:,:)

    end if

    ! Cycle over x2 position (or order of x2-derivative at boundary)
    ! and interpolate f along x1 direction. Store coefficients in bwork
    do i2 = 1, n2

      call self%bs1%compute_interpolant( &
        gtau        = w(   1+a1:n1-b1,i2) , &
        derivs_xmin = w(      1:a1   ,i2) , &
        derivs_xmax = w(n1-b1+1:n1   ,i2) )

      w(1:n1,i2) = self%bs1%bcoef(1:n1)

    end do

    wt = transpose( w )

    ! Cycle over x1 position (or order of x1-derivative at boundary)
    ! and interpolate w along x2 direction. Store coefficients in bcoef
    do i1 = 1, n1

      call self%bs2%compute_interpolant( &
        gtau        = wt(   1+a2:n2-b2,i1) , &
        derivs_xmin = wt(      1:a2   ,i1) , &
        derivs_xmax = wt(n2-b2+1:n2   ,i1) )

      wt(1:n2,i1) = self%bs2%bcoef(1:n2)

    end do

    ! x1-periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin(1) == sll_p_periodic) then
      associate( g1 => 1+self%deg(1)/2 )
        wt(:,1-g1:0    ) = wt(:,n1-g1+1:n1)
        wt(:,n1+1:n1+g1) = wt(:,      1:g1)
      end associate
    end if

    w = transpose( wt )

    ! x2-periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin(2) == sll_p_periodic) then
      associate( g2 => 1+self%deg(2)/2 )
        w(:,1-g2:0    ) = w(:,n2-g2+1:n2)
        w(:,n2+1:n2+g2) = w(:,      1:g2)
      end associate
    end if

    end associate

  end subroutine sll_s_spline_2d_non_uniform_compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function sll_f_spline_2d_non_uniform_eval( self, x1, x2 ) result( y )

    type(sll_t_spline_2d_non_uniform), intent(in) :: self
    sll_real64                       , intent(in) :: x1
    sll_real64                       , intent(in) :: x2
    sll_real64 :: y

    sll_int32 :: d1, d2
    sll_int32 :: a1, a2
    sll_int32 :: b1, b2
    sll_int32 :: k1, k2
    sll_int32 :: ic1, ic2

    ! Automatic arrays
    sll_real64 :: values1(1+self%bs1%deg)
    sll_real64 :: values2(1+self%bs2%deg)

    d1 = self%bs1%deg+1
    d2 = self%bs2%deg+1

    ! Find 2D cell (ic1,ic2)
    ic1 = sll_f_find_cell( self%bs1%bsp, x1 )
    ic2 = sll_f_find_cell( self%bs2%bsp, x2 )

    ! Compute arrays v1 and v2 of B-spline values
    call sll_s_bsplines_eval_basis( self%bs1%bsp, ic1, x1, values1(1:d1) )
    call sll_s_bsplines_eval_basis( self%bs2%bsp, ic2, x2, values2(1:d2) )

    ! Calculate (matrix) block C of B-spline coefficients to be used
    a1 = ic1 - self%bs1%offset;  b1 = ic1 - self%bs1%offset + self%bs1%deg
    a2 = ic2 - self%bs2%offset;  b2 = ic2 - self%bs2%offset + self%bs2%deg

    ! Compute scalar product <v1,v2> = (v1^T)*C*v2
    associate( bcoef => self%bcoef(a1:b1,a2:b2) )
      y = 0.0_f64
      do k2 = 1, d2
        do k1 = 1, d1
          y  = y + bcoef(k1,k2) * values1(k1) * values2(k2)
        end do
      end do
    end associate

  end function sll_f_spline_2d_non_uniform_eval


  !-----------------------------------------------------------------------------
  SLL_PURE function sll_f_spline_2d_non_uniform_eval_deriv_x1( self, x1, x2 ) result( y )

    type(sll_t_spline_2d_non_uniform), intent(in) :: self
    sll_real64                       , intent(in) :: x1
    sll_real64                       , intent(in) :: x2
    sll_real64 :: y

    sll_int32 :: d1, d2
    sll_int32 :: a1, a2
    sll_int32 :: b1, b2
    sll_int32 :: k1, k2
    sll_int32 :: ic1, ic2

    ! Automatic arrays
    sll_real64 :: derivs1(1+self%bs1%deg)
    sll_real64 :: values2(1+self%bs2%deg)

    d1 = self%bs1%deg+1
    d2 = self%bs2%deg+1

    ! Find 2D cell (ic1,ic2)
    ic1 = sll_f_find_cell( self%bs1%bsp, x1 )
    ic2 = sll_f_find_cell( self%bs2%bsp, x2 )

    ! Compute arrays d1 and v2 of B-spline derivatives/values
    call sll_s_bsplines_eval_deriv( self%bs1%bsp, ic1, x1, derivs1(1:d1) )
    call sll_s_bsplines_eval_basis( self%bs2%bsp, ic2, x2, values2(1:d2) )

    ! Calculate (matrix) block C of B-spline coefficients to be used
    a1 = ic1 - self%bs1%offset;  b1 = ic1 - self%bs1%offset + self%bs1%deg
    a2 = ic2 - self%bs2%offset;  b2 = ic2 - self%bs2%offset + self%bs2%deg

    ! Compute scalar product <d1,v2> = (d1^T)*C*v2
    associate( bcoef => self%bcoef(a1:b1,a2:b2) )
      y = 0.0_f64
      do k2 = 1, d2
        do k1 = 1, d1
          y  = y + bcoef(k1,k2) * derivs1(k1) * values2(k2)
        end do
      end do
    end associate

  end function sll_f_spline_2d_non_uniform_eval_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE function sll_f_spline_2d_non_uniform_eval_deriv_x2( self, x1, x2 ) result( y )

    type(sll_t_spline_2d_non_uniform), intent(in) :: self
    sll_real64                       , intent(in) :: x1
    sll_real64                       , intent(in) :: x2
    sll_real64 :: y

    sll_int32 :: d1, d2
    sll_int32 :: a1, a2
    sll_int32 :: b1, b2
    sll_int32 :: k1, k2
    sll_int32 :: ic1, ic2

    ! Automatic arrays
    sll_real64 :: values1(1+self%bs1%deg)
    sll_real64 :: derivs2(1+self%bs2%deg)

    d1 = self%bs1%deg+1
    d2 = self%bs2%deg+1

    ! Find 2D cell (ic1,ic2)
    ic1 = sll_f_find_cell( self%bs1%bsp, x1 )
    ic2 = sll_f_find_cell( self%bs2%bsp, x2 )

    ! Compute arrays v1 and d2 of B-spline values/derivatives
    call sll_s_bsplines_eval_basis( self%bs1%bsp, ic1, x1, values1(1:d1) )
    call sll_s_bsplines_eval_deriv( self%bs2%bsp, ic2, x2, derivs2(1:d2) )

    ! Calculate (matrix) block C of B-spline coefficients to be used
    a1 = ic1 - self%bs1%offset;  b1 = ic1 - self%bs1%offset + self%bs1%deg
    a2 = ic2 - self%bs2%offset;  b2 = ic2 - self%bs2%offset + self%bs2%deg

    ! Compute scalar product <v1,d2> = (v1^T)*C*d2
    associate( bcoef => self%bcoef(a1:b1,a2:b2) )
      y = 0.0_f64
      do k2 = 1, d2
        do k1 = 1, d1
          y  = y + bcoef(k1,k2) * values1(k1) * derivs2(k2)
        end do
      end do
    end associate

  end function sll_f_spline_2d_non_uniform_eval_deriv_x2

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_spline_2d_non_uniform_eval_array( self, x1, x2, y )

    type(sll_t_spline_2d_non_uniform), intent(in   ) :: self
    sll_real64                       , intent(in   ) :: x1(:,:)
    sll_real64                       , intent(in   ) :: x2(:,:)
    sll_real64                       , intent(  out) :: y (:,:)

    sll_int32 :: i1, i2
    sll_int32 :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_spline_2d_non_uniform_eval( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine sll_s_spline_2d_non_uniform_eval_array

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_spline_2d_non_uniform_eval_array_deriv_x1( self, x1, x2, y )

    type(sll_t_spline_2d_non_uniform), intent(in   ) :: self
    sll_real64                       , intent(in   ) :: x1(:,:)
    sll_real64                       , intent(in   ) :: x2(:,:)
    sll_real64                       , intent(  out) :: y (:,:)

    sll_int32 :: i1, i2
    sll_int32 :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_spline_2d_non_uniform_eval_deriv_x1( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine sll_s_spline_2d_non_uniform_eval_array_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_spline_2d_non_uniform_eval_array_deriv_x2( self, x1, x2, y )

    type(sll_t_spline_2d_non_uniform), intent(in   ) :: self
    sll_real64                       , intent(in   ) :: x1(:,:)
    sll_real64                       , intent(in   ) :: x2(:,:)
    sll_real64                       , intent(  out) :: y (:,:)

    sll_int32 :: i1, i2
    sll_int32 :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_spline_2d_non_uniform_eval_deriv_x2( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine sll_s_spline_2d_non_uniform_eval_array_deriv_x2

  !-----------------------------------------------------------------------------
  subroutine sll_s_spline_2d_non_uniform_free( self )

    type(sll_t_spline_2d_non_uniform), intent(inout) :: self

    ! Deallocate local 2D arrays
    deallocate( self%bwork )
    deallocate( self%bcoef )

    ! Free memory of 1D B-splines
    call self%bs1%free()
    call self%bs2%free()

  end subroutine sll_s_spline_2d_non_uniform_free

end module sll_m_spline_2d_non_uniform
