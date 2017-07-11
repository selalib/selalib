!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_arbitrary_degree_splines
module sll_m_bspline_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
     sll_p_periodic, &
     sll_p_hermite, &
     sll_p_greville, &
     sll_p_mirror

use sll_m_arbitrary_degree_splines, only: &
     sll_t_arbitrary_degree_spline_1d, &
     sll_s_arbitrary_degree_spline_1d_init, &
     sll_s_arbitrary_degree_spline_1d_free, &
     sll_f_find_cell, &
     sll_s_splines_at_x, &
     sll_s_spline_derivatives_at_x, &
     sll_s_splines_and_n_derivs_at_x, &
     sll_s_uniform_b_splines_at_x

use schur_complement, only: &
  schur_complement_solver, &
  schur_complement_fac   , &
  schur_complement_slv   , &
  schur_complement_free

use sll_m_bspline_1d, only: &
  sll_t_bspline_1d,         &
  sll_s_bspline_1d_init,    &
  sll_s_bspline_1d_free,    &
  sll_s_bspline_1d_compute_interpolant

implicit none

public :: &
  sll_t_bspline_2d,                     &
  sll_s_bspline_2d_init,                &
  sll_s_bspline_2d_free,                &
  sll_s_bspline_2d_compute_interpolant, &
  sll_f_bspline_2d_eval,                & ! scalar functions for evaluation
  sll_f_bspline_2d_eval_deriv_x1,       &
  sll_f_bspline_2d_eval_deriv_x2,       &
  sll_s_bspline_2d_eval_array,          & ! vector subroutines for evaluation
  sll_s_bspline_2d_eval_array_deriv_x1, &
  sll_s_bspline_2d_eval_array_deriv_x2

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief
!> basic type for two-dimensional B-spline data.
!> @details
!> treated as an opaque type. No access to its internals is directly allowed.
type :: sll_t_bspline_2d

  type(sll_t_bspline_1d)  :: bs1
  type(sll_t_bspline_1d)  :: bs2
  sll_real64, allocatable :: bcoef(:,:)
  sll_real64, allocatable :: bwork(:,:)

end type sll_t_bspline_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief Initialises a 2D spline interpolation object.
  !> @param[in] nx1 Number of points where the data to be interpolated are
  !>            represented.
  !> @param[in] degree1 Spline degree
  !> @param[in] x1_min Minimum value of the abscissae where the data are meant
  !> to be interpolated.
  !> @param[in] x1_max Maximum value of the abscissae where the data are meant
  !> to be interpolated.
  !> @param[in] bc1 A boundary condition specifier. Must be one of the
  !> symbols defined in the SLL_BOUNDARY_CONDITION_DESCRIPTORS module.
  !> @param[in] nx2 Number of points where the data to be interpolated are
  !>            represented.
  !> @param[in] degree2 Spline degree
  !> @param[in] x2_min Minimum value of the abscissae where the data are meant
  !> to be interpolated.
  !> @param[in] x2_max Maximum value of the abscissae where the data are meant
  !> to be interpolated.
  !> @param[in] bc2 A boundary condition specifier. Must be one of the
  !> symbols defined in the sll_m_boundary_condition_descriptors module.
  !> @param[in] spline_bc_type1 A boundary condition specifier (see sll_s_bspline_1d_init).
  !> @param[in] spline_bc_type2 A boundary condition specifier (see sll_s_bspline_1d_init).
  !> @param[in] bc_left1 value of function on the west boundary
  !> @param[in] bc_left2 value of function on the south boundary
  !> @param[in] bc_right1 value of function on the east boundary
  !> @param[in] bc_right2 value of the function on the north boundary
  !> @return a spline interpolation object.
  !-----------------------------------------------------------------------------
  subroutine sll_s_bspline_2d_init( &
    self,            &
    num_pts1,        &
    num_pts2,        &
    degree1,         &
    degree2,         &
    x1_min,          &
    x2_min,          &
    x1_max,          &
    x2_max,          &
    bc1_min,         &
    bc2_min,         &
    bc1_max,         &
    bc2_max,         &
    spline_bc_type1, &
    spline_bc_type2 )

    type(sll_t_bspline_2d), intent(  out) :: self
    sll_int32 ,             intent(in   ) :: num_pts1
    sll_int32 ,             intent(in   ) :: num_pts2
    sll_int32 ,             intent(in   ) :: degree1
    sll_int32 ,             intent(in   ) :: degree2
    sll_real64,             intent(in   ) :: x1_min
    sll_real64,             intent(in   ) :: x2_min
    sll_real64,             intent(in   ) :: x1_max
    sll_real64,             intent(in   ) :: x2_max
    sll_int32 ,             intent(in   ) :: bc1_min
    sll_int32 ,             intent(in   ) :: bc2_min
    sll_int32 ,             intent(in   ) :: bc1_max
    sll_int32 ,             intent(in   ) :: bc2_max
    sll_int32 , optional,   intent(in   ) :: spline_bc_type1
    sll_int32 , optional,   intent(in   ) :: spline_bc_type2

    sll_int32 :: n1
    sll_int32 :: n2

    ! NOTE: in the future different boundary conditions at xmin and xmax
    !       should be considered. For now we only check that bc_xmin==bc_xmax
    SLL_ASSERT( bc1_min == bc1_max )
    SLL_ASSERT( bc2_min == bc2_max )

    if (present(spline_bc_type1)) then
      call sll_s_bspline_1d_init( self%bs1, num_pts1, degree1, &
           x1_min, x1_max, bc1_min, bc1_max, spline_bc_type1)
    else
      call sll_s_bspline_1d_init( self%bs1, num_pts1, degree1, &
           x1_min, x1_max, bc1_min, bc1_max )
    end if
    if (present(spline_bc_type2)) then
      call sll_s_bspline_1d_init( self%bs2, num_pts2, degree2, &
           x2_min, x2_max, bc2_min, bc2_max, spline_bc_type2 )
    else
      call sll_s_bspline_1d_init( self%bs2, num_pts2, degree2, &
           x2_min, x2_max, bc2_min, bc2_max )
    end if

    n1 = self%bs1%n
    n2 = self%bs2%n
    allocate( self%bwork(1:n2,1:n1) ); self%bwork = 0.0_f64
    allocate( self%bcoef(1:n1,1:n2) ); self%bcoef = 0.0_f64

  end subroutine sll_s_bspline_2d_init

  !-----------------------------------------------------------------------------
  subroutine sll_s_bspline_2d_compute_interpolant( &
    self, &
    gtau, &
    derivs_x1_min, &
    derivs_x1_max, &
    derivs_x2_min, &
    derivs_x2_max, &
    derivs_corners )

    type(sll_t_bspline_2d), intent(inout)           :: self
    sll_real64            , intent(in   )           :: gtau(:,:)
    sll_real64            , intent(in   ), optional :: derivs_x1_min (:,:)
    sll_real64            , intent(in   ), optional :: derivs_x1_max (:,:)
    sll_real64            , intent(in   ), optional :: derivs_x2_min (:,:)
    sll_real64            , intent(in   ), optional :: derivs_x2_max (:,:)
    sll_real64            , intent(in   ), optional :: derivs_corners(:,:,:)

    sll_int32 :: i1, n1, ncond1
    sll_int32 :: i2, n2, ncond2, j2

    ! TODO: check size of gtau == tau
    ! TODO: check size of all input arrays

    ! Number of degrees of freedom
    n1 = self%bs1%n
    n2 = self%bs2%n

    ! Compute number of additional boundary values needed
    ncond1 = 0
    ncond2 = 0
    if (self%bs1%bc_type == sll_p_hermite)  ncond1 = self%bs1%deg/2
    if (self%bs2%bc_type == sll_p_hermite)  ncond2 = self%bs2%deg/2

    !--------------------------------------------
    ! Compute spline coefficients in x1 direction
    !--------------------------------------------

    ! Boundary points in x2 (only for Hermite):
    ! Additional x1-interpolation of x2-derivatives on boundaries is needed!
    if (self%bs2%bc_type == sll_p_hermite) then

      ! Hermite-Hermite case:
      if (self%bs1%bc_type == sll_p_hermite) then

        ! boundary conditions at x2_min
        if (present(derivs_x2_min)) then
          do j2 = 1, ncond2
            if( present(derivs_corners)) then
              call sll_s_bspline_1d_compute_interpolant( self%bs1, &
                derivs_x2_min (j2,:)  , &
                derivs_corners(:,j2,1), &
                derivs_corners(:,j2,2))
            else
              call sll_s_bspline_1d_compute_interpolant( self%bs1, derivs_x2_min(j2,:) )
            end if
            self%bwork(j2,:) = self%bs1%bcoef(:)
          end do
        else  ! set needed boundary values to 0
          self%bwork(1:ncond2,:) = 0.0_f64
        end if

        ! boundary conditions at x2_max
        if (present(derivs_x2_max)) then
          do j2 = 1, ncond2
            if (present(derivs_corners)) then
              call sll_s_bspline_1d_compute_interpolant( self%bs1, &
                derivs_x2_max (j2,:)  , &
                derivs_corners(:,j2,3), &
                derivs_corners(:,j2,4))
            else
              call sll_s_bspline_1d_compute_interpolant( self%bs1, derivs_x2_max(j2,:) )
            end if
            self%bwork(n2-ncond2+j2,:) = self%bs1%bcoef(:)
          end do
        else ! set needed boundary values to 0
          self%bwork(n2-ncond2+1:n2,:) = 0.0_f64
        end if

        ! Other-Hermite case:
      else

        ! boundary conditions at x2_min
        if (present(derivs_x2_min)) then
          do j2 = 1, ncond2
            call sll_s_bspline_1d_compute_interpolant( self%bs1, derivs_x2_min(j2,:) )
            self%bwork(j2,:) = self%bs1%bcoef(:)
          end do
        else  ! set needed boundary values to 0
          self%bwork(1:ncond2,:) = 0.0_f64
        end if

        ! boundary conditions at x2_max
        if (present(derivs_x2_max)) then
          do j2 = 1, ncond2
            call sll_s_bspline_1d_compute_interpolant( self%bs1, derivs_x2_max(j2,:) )
            self%bwork(n2-ncond2+j2,:) = self%bs1%bcoef(:)
          end do
        else ! set needed boundary values to 0
          self%bwork(n2-ncond2+1:n2,:) = 0.0_f64
        end if

      end if
    end if

    ! Interior points
    do i2 = 1, n2-2*ncond2
      if (present(derivs_x1_min) .and. present(derivs_x1_max)) then
        call sll_s_bspline_1d_compute_interpolant( self%bs1, &
          gtau         (:,i2), &
          derivs_x1_min(:,i2), &
          derivs_x1_max(:,i2) )
      else
        call sll_s_bspline_1d_compute_interpolant( self%bs1, gtau(:,i2) )
      end if
      self%bwork(i2+ncond2,:) = self%bs1%bcoef(:)
    end do

    !--------------------------------------------
    ! Compute spline coefficients in x2 direction
    !--------------------------------------------
    if (self%bs2%bc_type == sll_p_hermite) then

      if (present(derivs_x2_min) .and. present(derivs_x2_max)) then
        do i1 = 1, n1
          call sll_s_bspline_1d_compute_interpolant( self%bs2, &
               self%bwork(ncond2+1:n2-ncond2,i1), &
               self%bwork(1:ncond2,i1), &
               self%bwork(n2-ncond2+1:n2,i1) )
          self%bcoef(i1,:) = self%bs2%bcoef(:)
        end do
      else
        do i1 = 1, n1
          call sll_s_bspline_1d_compute_interpolant( self%bs2, self%bwork(ncond2+1:n2-ncond2,i1) )
          self%bcoef(i1,:) = self%bs2%bcoef(:)
        end do
      end if

    else

      do i1 = 1, n1
        call sll_s_bspline_1d_compute_interpolant( self%bs2, self%bwork(:,i1) )
        self%bcoef(i1,:) = self%bs2%bcoef(:)
      end do

    end if

  end subroutine sll_s_bspline_2d_compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function sll_f_bspline_2d_eval( self, x1, x2 ) result( y )

    type(sll_t_bspline_2d), intent(in) :: self
    sll_real64            , intent(in) :: x1
    sll_real64            , intent(in) :: x2
    sll_real64 :: y

    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: ib
    sll_int32 :: jb
    sll_int32 :: k1
    sll_int32 :: k2
    sll_int32 :: icell
    sll_int32 :: jcell

    ! Automatic arrays
    sll_real64 :: values(1+max(self%bs1%deg,self%bs2%deg))
    sll_real64 :: work  (1+self%bs2%deg)

    k1 = self%bs1%deg+1
    k2 = self%bs2%deg+1

    icell = sll_f_find_cell( self%bs1%bsp, x1 )
    jcell = sll_f_find_cell( self%bs2%bsp, x2 )

    work = 0.0_f64
    do j = 1, k2
      jb = mod( jcell+j-2-self%bs2%offset+self%bs2%n, self%bs2%n ) + 1
      call sll_s_splines_at_x( self%bs1%bsp, icell, x1, values(1:k1) )
      do k = 1, k1
        ib = mod( icell+k-2-self%bs1%offset+self%bs1%n, self%bs1%n ) + 1
        work(j) = work(j) + values(k)*self%bcoef(ib,jb)
      end do
    end do

    call sll_s_splines_at_x( self%bs2%bsp, jcell, x2, values(1:k2) )
    y = 0.0_f64
    do k = 1, k2
      y  = y + values(k)*work(k)
    end do

  end function sll_f_bspline_2d_eval


  !-----------------------------------------------------------------------------
  SLL_PURE function sll_f_bspline_2d_eval_deriv_x1( self, x1, x2 ) result( y )

    type(sll_t_bspline_2d), intent(in) :: self
    sll_real64            , intent(in) :: x1
    sll_real64            , intent(in) :: x2
    sll_real64 :: y

    sll_int32 :: i, j
    sll_int32 :: k1, k2, icell, jcell
    sll_int32 :: ib, jb

    ! Automatic arrays
    sll_real64 :: values(1+max(self%bs1%deg,self%bs2%deg))
    sll_real64 :: work  (1+self%bs2%deg)

    k1 = self%bs1%deg+1
    k2 = self%bs2%deg+1

    icell =  sll_f_find_cell( self%bs1%bsp, x1 )
    jcell =  sll_f_find_cell( self%bs2%bsp, x2 )

    work = 0.0_f64
    do j = 1, k2
      jb = mod( jcell+j-2-self%bs2%offset+self%bs2%n, self%bs2%n ) + 1
      call sll_s_spline_derivatives_at_x( self%bs1%bsp, icell, x1, values(1:k1) )
      do i = 1, k1
        ib = mod( icell+i-2-self%bs1%offset+self%bs1%n, self%bs1%n ) + 1
        work(j) = work(j) + values(i)*self%bcoef(ib,jb)
      end do
    end do

    call sll_s_splines_at_x( self%bs2%bsp, jcell, x2, values(1:k2) )
    y = 0.0_f64
    do j = 1, k2
      y = y + values(j)*work(j)
    end do

  end function sll_f_bspline_2d_eval_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE function sll_f_bspline_2d_eval_deriv_x2( self, x1, x2 ) result( y )

    type(sll_t_bspline_2d), intent(in) :: self
    sll_real64            , intent(in) :: x1
    sll_real64            , intent(in) :: x2
    sll_real64 :: y

    sll_int32 :: i, j
    sll_int32 :: k1, k2, icell, jcell
    sll_int32 :: ib, jb

    ! Automatic arrays
    sll_real64 :: values(1+max(self%bs1%deg,self%bs2%deg))
    sll_real64 :: work  (1+self%bs2%deg)

    k1 = self%bs1%deg+1
    k2 = self%bs2%deg+1

    icell = sll_f_find_cell( self%bs1%bsp, x1 )
    jcell = sll_f_find_cell( self%bs2%bsp, x2 )

    work = 0.0_f64
    do j = 1, k2
      jb = mod( jcell+j-2-self%bs2%offset+self%bs2%n, self%bs2%n ) + 1
      call sll_s_splines_at_x( self%bs1%bsp, icell, x1, values(1:k1) )
      do i = 1, k1
        ib = mod( icell+i-2-self%bs1%offset+self%bs1%n, self%bs1%n ) + 1
        work(j) = work(j) + values(i)*self%bcoef(ib,jb)
      end do
    end do

    call sll_s_spline_derivatives_at_x( self%bs2%bsp, jcell, x2, values(1:k2) )
    y = 0.0_f64
    do j = 1, k2
      y = y + values(j)*work(j)
    end do

  end function sll_f_bspline_2d_eval_deriv_x2

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bspline_2d_eval_array( self, x1, x2, y )

    type(sll_t_bspline_2d), intent(in   ) :: self
    sll_real64            , intent(in   ) :: x1(:,:)
    sll_real64            , intent(in   ) :: x2(:,:)
    sll_real64            , intent(  out) :: y (:,:)

    sll_int32 :: i1, i2
    sll_int32 :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_bspline_2d_eval( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine sll_s_bspline_2d_eval_array

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bspline_2d_eval_array_deriv_x1( self, x1, x2, y )

    type(sll_t_bspline_2d), intent(in   ) :: self
    sll_real64            , intent(in   ) :: x1(:,:)
    sll_real64            , intent(in   ) :: x2(:,:)
    sll_real64            , intent(  out) :: y (:,:)

    sll_int32 :: i1, i2
    sll_int32 :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_bspline_2d_eval_deriv_x1( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine sll_s_bspline_2d_eval_array_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine sll_s_bspline_2d_eval_array_deriv_x2( self, x1, x2, y )

    type(sll_t_bspline_2d), intent(in   ) :: self
    sll_real64            , intent(in   ) :: x1(:,:)
    sll_real64            , intent(in   ) :: x2(:,:)
    sll_real64            , intent(  out) :: y (:,:)

    sll_int32 :: i1, i2
    sll_int32 :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_bspline_2d_eval_deriv_x2( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine sll_s_bspline_2d_eval_array_deriv_x2

  !-----------------------------------------------------------------------------
  subroutine sll_s_bspline_2d_free( self )

    type(sll_t_bspline_2d), intent(inout) :: self

    ! Deallocate local 2D arrays
    deallocate( self%bwork )
    deallocate( self%bcoef )

    ! Free memory of 1D B-splines
    call sll_s_bspline_1d_free( self%bs1 )
    call sll_s_bspline_1d_free( self%bs2 )

  end subroutine sll_s_bspline_2d_free

end module sll_m_bspline_2d
