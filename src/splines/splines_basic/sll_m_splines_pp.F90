!> @ingroup splines
!> @authors
!> Katharina Kormann, IPP
!> Benedikt Perse, IPP
!> @brief Splines in pp form
!> @details Implement splines on a uniform grid, periodic boundaries implemented for degree 0-6 and clamped boundaries for degree 1-3

module sll_m_splines_pp

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_t_spline_pp_1d, &
      sll_t_spline_pp_2d, &
      sll_t_spline_pp_3d, &
      sll_s_spline_pp_b_to_pp_1d, &
      sll_s_spline_pp_b_to_pp_2d, &
      sll_s_spline_pp_b_to_pp_2d_periodic, &
      sll_s_spline_pp_b_to_pp_3d, &
      sll_s_spline_pp_b_to_pp_3d_clamped, &
      sll_s_spline_pp_b_to_pp_3d_clamped_2full, &
      sll_s_spline_pp_horner_m_1d, &
      sll_s_spline_pp_horner_m_2d, &
      sll_s_spline_pp_horner_m_3d, &
      sll_s_spline_pp_horner_primitive_1d, &
      sll_f_spline_pp_horner_1d, &
      sll_f_spline_pp_horner_2d, &
      sll_f_spline_pp_horner_3d, &
      sll_f_spline_pp_horner_3d_d1, &
      sll_f_spline_pp_horner_3d_d2, &
      sll_f_spline_pp_horner_3d_d3, &
      sll_f_spline_pp_horner_derivative_3d, &
      sll_s_spline_pp_init_1d, &
      sll_s_spline_pp_init_2d, &
      sll_s_spline_pp_init_3d, &
      sll_s_spline_pp_free_1d, &
      sll_s_spline_pp_free_2d, &
      sll_s_spline_pp_free_3d, &
      sll_s_spline_pp_pp_to_b_1d, &
      sll_s_spline_evaluate_basis_b_form_1d_clamped, &
      sll_s_spline_evaluate_basis_b_form_1d_periodic, &
      sll_s_spline_evaluate_basis_b_form_3d_clamped, &
      sll_s_spline_evaluate_basis_b_form_3d_periodic, &
      sll_p_boundary_periodic, &
      sll_p_boundary_clamped, &
      sll_p_boundary_clamped_square, &
      sll_p_boundary_clamped_cubic, &
      sll_p_boundary_clamped_clampeddiri, &
      sll_p_boundary_clampeddiri, &
      sll_p_boundary_clampeddiri_clamped

   private
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> Parameter to specify boundary conditions
   sll_int32, parameter :: sll_p_boundary_periodic = 0 !< Parameter specifying periodic boundary conditions
   sll_int32, parameter :: sll_p_boundary_clamped = 1 !< Parameter specifying clamped boundary conditions
   sll_int32, parameter :: sll_p_boundary_clamped_clampeddiri = 2 !< Parameter specifying clamped boundary conditions with the last point left out (for Dirichlet conditions on the right side of the domain)
   sll_int32, parameter :: sll_p_boundary_clampeddiri = 3 !< Parameter specifying clamped boundary conditions with the first and last point left out (for Dirichlet conditions)
   sll_int32, parameter :: sll_p_boundary_clampeddiri_clamped = 4 !< Parameter specifying clamped boundary conditions with the first point left out (for Dirichlet conditions on the left side of the domain)
   sll_int32, parameter :: sll_p_boundary_clamped_square = 5 !< Parameter specifying clamped boundary conditions with square spline for the 0-form, boundary spline of degree 2 is set to zero
   sll_int32, parameter :: sll_p_boundary_clamped_cubic = 6  !< Parameter specifying clamped boundary conditions with a cubic spline for the 0-form, boundary spline of degree 3 is set to zero

   sll_real64, parameter :: inv_2 = 1._f64/2._f64
   sll_real64, parameter :: inv_3 = 1._f64/3._f64
   sll_real64, parameter :: inv_4 = 1._f64/4._f64
   sll_real64, parameter :: inv_6 = 1._f64/6._f64
   sll_real64, parameter :: inv_7 = 1._f64/7._f64
   sll_real64, parameter :: inv_8 = 1._f64/8._f64
   sll_real64, parameter :: inv_12 = 1._f64/12._f64
   sll_real64, parameter :: inv_16 = 1._f64/16._f64
   sll_real64, parameter :: inv_18 = 1._f64/18._f64
   sll_real64, parameter :: inv_20 = 1._f64/20._f64
   sll_real64, parameter :: inv_24 = 1._f64/24._f64
   sll_real64, parameter :: inv_30 = 1._f64/30._f64
   sll_real64, parameter :: inv_36 = 1._f64/36._f64
   sll_real64, parameter :: inv_48 = 1._f64/48._f64
   sll_real64, parameter :: inv_60 = 1._f64/60._f64
   sll_real64, parameter :: inv_72 = 1._f64/72._f64
   sll_real64, parameter :: inv_80 = 1._f64/80._f64
   sll_real64, parameter :: inv_120 = 1._f64/120._f64
   sll_real64, parameter :: inv_240 = 1._f64/240._f64
   sll_real64, parameter :: inv_144 = 1._f64/144._f64
   ! sll_real64, parameter :: inv_600 = 1._f64/600._f64
   sll_real64, parameter :: inv_720 = 1._f64/720._f64
   sll_real64, parameter :: inv_5040 = 1._f64/5040._f64
   !>arbitrary degree 1d spline
   type :: sll_t_spline_pp_1d

      sll_int32 :: degree !< degree of 1d spline
      sll_real64, allocatable :: poly_coeffs(:, :) !< poly_coeffs[j,i] coefficient of x^{degree+1-j} for B-spline in interval degree+2-i function  size= (degree+1, degree+1)
      sll_real64, allocatable :: poly_coeffs_fp(:, :) !< poly_coeffs[j,i] coefficient of x^{deg+1-j} for primitive of B-spline function in interval degree+2-i  size= (degree+1, degree+1)
      sll_real64, allocatable :: poly_coeffs_fpa(:, :) !< poly_coeffs[j,i] coefficient of x^{deg+1-j+1} for primitive of B-spline function in interval degree+2-i  size= (degree+2, degree+1)
      sll_real64, allocatable :: poly_coeffs_fd(:, :) !< poly_coeffs[j,i] coefficient of x^{deg+1-j} for derivative of B-spline function in interval degree+2-i  size= (degree+1, degree+1)
      sll_int32 :: n_cells !< number of gridcells
      sll_int32 :: n_coeffs !< no. of coefficients
      sll_real64, allocatable :: scratch_b(:) !< scratch data for b_to_pp-converting
      sll_real64, allocatable :: scratch_p(:) !< scratch data for b_to_pp-converting

      sll_real64, allocatable :: poly_coeffs_boundary_left(:, :, :) !< poly_coeffs_left[j,i,k] coefficient of x^{degree+1-j} for B-spline in interval degree+2-i function  size= (degree+1, degree+1,degree-1) for kth boundary cell
      sll_real64, allocatable :: poly_coeffs_boundary_right(:, :, :) !< poly_coeffs[j,i,k] coefficient of x^{degree+1-j} for B-spline in interval degree+2-i function  size= (degree+1, degree+1, degree-1) for (N-k)th boundary cell
      sll_real64, allocatable :: poly_coeffs_fd_boundary_left(:, :, :) !< poly_coeffs_left[j,i,k] coefficient of x^{degree+1-j} for B-spline in interval degree+2-i function  size= (degree+1, degree+1,degree-1) for kth boundary cell
      sll_real64, allocatable :: poly_coeffs_fd_boundary_right(:, :, :) !< poly_coeffs[j,i,k] coefficient of x^{degree+1-j} for B-spline in interval degree+2-i function  size= (degree+1, degree+1, degree-1) for (N-k)th boundary cell
      sll_int32 :: boundary_conditions = sll_p_boundary_periodic !< flag for boundary conditions

   end type sll_t_spline_pp_1d
   !>arbitrary degree 2d spline
   type :: sll_t_spline_pp_2d
      type(sll_t_spline_pp_1d) spline1 !< first 1d spline
      type(sll_t_spline_pp_1d) spline2 !< second 1d spline

   end type sll_t_spline_pp_2d

   !>arbitrary degree 3d spline
   type :: sll_t_spline_pp_3d
      type(sll_t_spline_pp_1d) spline1 !< first 1d spline
      type(sll_t_spline_pp_1d) spline2 !< second 1d spline
      type(sll_t_spline_pp_1d) spline3 !< third 1d spline

      sll_real64, allocatable :: scratch_b(:)  !< scratch data for b_to_pp-converting
      sll_real64, allocatable :: scratch_p(:)  !< scratch data for b_to_pp-converting
      sll_real64, allocatable :: scratch_pp(:) !< scratch data for b_to_pp-converting
      sll_real64, allocatable :: scratch_coeffs(:, :) !< scratch data for spline-evaluation

   end type sll_t_spline_pp_3d

contains

   !> Compute b_coeffs from coefficients of the monomials (in \a pp_coeffs)
   subroutine sll_s_spline_pp_pp_to_b_1d(self, n_cells, pp_coeffs, b_coeffs)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(out) :: b_coeffs(0:n_cells - 1)   !< coefficients of spline in B-form
      sll_real64, intent(in):: pp_coeffs(self%degree + 1, 0:n_cells - 1)  !< coefficients of spline in pp-form

      sll_int32 :: i, j, imin

      b_coeffs = 0.0_f64
      do i = 0, n_cells - 1
         imin = i - self%degree
         do j = 0, self%degree
            b_coeffs(modulo(imin + j, n_cells)) = b_coeffs(modulo(imin + j, n_cells)) + &
                                                  sum(pp_coeffs(:, i)*self%poly_coeffs(:, j + 1))
         end do
      end do

   end subroutine sll_s_spline_pp_pp_to_b_1d

   !> Convert 1d spline in B form to spline in pp form for clamped spline
   subroutine sll_s_spline_pp_b_to_pp_1d_clamped(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells + self%degree)   !< coefficients of spline in B-form
      sll_real64, intent(out):: pp_coeffs(self%degree + 1, n_cells)  !< coefficients of spline in pp-form
      sll_int32 :: i, upper

      do i = 1, self%degree - 1
         call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs_boundary_left(:, :, i), &
                                               b_coeffs(i:i + self%degree), pp_coeffs(:, i))
      end do

      upper = n_cells - self%degree + 1
      do i = self%degree, upper
         call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs, &
                                               b_coeffs(i:i + self%degree), pp_coeffs(:, i))
      end do

      do i = upper + 1, n_cells
         call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs_boundary_right(:, :, i - upper), &
                                               b_coeffs(i:i + self%degree), pp_coeffs(:, i))
      end do

   end subroutine sll_s_spline_pp_b_to_pp_1d_clamped

   !> Convert 1d spline in B form to spline in pp form for clamped spline and Dirichlet conditions on the right boundary
   subroutine sll_s_spline_pp_b_to_pp_1d_clamped_clampeddiri(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells + self%degree - 1)   !< coefficients of spline in B-form
      sll_real64, intent(out):: pp_coeffs(self%degree + 1, n_cells)  !< coefficients of spline in pp-form
      sll_int32 :: i, upper

      do i = 1, self%degree - 1
         call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs_boundary_left(:, :, i), &
                                               b_coeffs(i:i + self%degree), pp_coeffs(:, i))
      end do

      upper = n_cells - self%degree + 1
      do i = self%degree, upper
         call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs, &
                                               b_coeffs(i:i + self%degree), pp_coeffs(:, i))
      end do

      do i = upper + 1, n_cells - 1
         call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs_boundary_right(:, :, i - upper), &
                                               b_coeffs(i:i + self%degree), pp_coeffs(:, i))
      end do
      call sll_s_spline_pp_b_to_pp_1d_cella(self%degree, self%poly_coeffs_boundary_right(:, :, i - upper), &
                                            [b_coeffs(i:i + self%degree - 1), 0.0_f64], pp_coeffs(:, i))

   end subroutine sll_s_spline_pp_b_to_pp_1d_clamped_clampeddiri

   !> Convert 1d spline in B form to spline in pp form with periodic boundary conditions
   subroutine sll_s_spline_pp_b_to_pp_1d(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells)   !< coefficients of spline in B-form
      sll_real64, intent(out):: pp_coeffs(self%degree + 1, n_cells)  !< coefficients of spline in pp-form
      sll_int32 :: i

      do i = 1, self%degree
         call sll_s_spline_pp_b_to_pp_1d_cell(self, (/b_coeffs(n_cells - self%degree + i:n_cells), b_coeffs(1:i)/), pp_coeffs(:, i))
      end do

      do i = self%degree + 1, n_cells
         call sll_s_spline_pp_b_to_pp_1d_cell(self, b_coeffs(i - self%degree:i), pp_coeffs(:, i))
      end do

   end subroutine sll_s_spline_pp_b_to_pp_1d

   !> Convert 1d spline in B form in a cell to spline in pp form with periodic boundary conditions
   subroutine sll_s_spline_pp_b_to_pp_1d_cell(self, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_real64, intent(in)  :: b_coeffs(self%degree + 1)   !< coefficients of spline in B-form
      sll_real64, intent(out) :: pp_coeffs(self%degree + 1)  !< coefficients of spline in pp-form
      sll_int32 :: i
      sll_int32 :: j
      sll_int32 :: degp1
      degp1 = self%degree + 1
      pp_coeffs = 0.0_f64
      do i = 1, degp1
         do j = 1, degp1
            pp_coeffs(j) = pp_coeffs(j) + b_coeffs(i)*self%poly_coeffs(j, i)
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_1d_cell

   !> Convert 1d spline in B form in a cell to spline in pp form for specified pp_coefficients of the b_splines
   subroutine sll_s_spline_pp_b_to_pp_1d_cella(degree, pp_b, b_coeffs, pp_coeffs)
      sll_int32, intent(in) :: degree !< spline degree
      sll_real64, intent(in) :: pp_b(degree + 1, degree + 1) !< pp coefficients of the b spline in this interval
      sll_real64, intent(in)  :: b_coeffs(degree + 1)   !< coefficients of spline in B-form
      sll_real64, intent(out) :: pp_coeffs(degree + 1)  !< coefficients of spline in pp-form
      sll_int32 :: i
      sll_int32 :: j
      sll_int32 :: degp1
      degp1 = degree + 1
      pp_coeffs = 0.0_f64
      do i = 1, degp1
         do j = 1, degp1
            pp_coeffs(j) = pp_coeffs(j) + b_coeffs(i)*pp_b(j, i)
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_1d_cella

   !> Convert 2d spline in B form to spline in pp form with periodic boundary conditions
   !> This is a special case of the procedure sll_s_spline_pp_b_to_pp_2d for the double periodic case to avoid the select case statements
   subroutine sll_s_spline_pp_b_to_pp_2d_periodic(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_2d), intent(inout)::  self !< arbitrary degree 2d spline
      sll_int32, intent(in) ::  n_cells(2) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
      sll_real64, intent(out):: pp_coeffs((self%spline1%degree + 1)*(self%spline2%degree + 1), n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
      sll_int32 :: i, j
      sll_int32 :: degree1, degree2
      degree1 = self%spline1%degree
      degree2 = self%spline2%degree

      do j = 1, n_cells(2)
         do i = 1, n_cells(1)
            call sll_s_spline_pp_b_to_pp_2d_cell(self%spline1, self%spline2, n_cells, b_coeffs, pp_coeffs, i, j)
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_2d_periodic

   !> Convert 2d spline in B form to spline in pp form
   subroutine sll_s_spline_pp_b_to_pp_2d(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_2d), intent(inout)::  self !< arbitrary degree 2d spline
      sll_int32, intent(in) ::  n_cells(2) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(:)!(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
      sll_real64, intent(out):: pp_coeffs((self%spline1%degree + 1)*(self%spline2%degree + 1), n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
      sll_int32 :: i, j, l1, l2
      sll_int32 :: degree1, degree2
      sll_real64 :: pp_coeffs_local(self%spline1%degree + 1, self%spline2%degree + 1)
      sll_int32 :: offset1, offset2, cell, n_coeffs(2), index, upper(2)

      degree1 = self%spline1%degree
      degree2 = self%spline2%degree
      n_coeffs(1) = self%spline1%n_coeffs
      n_coeffs(2) = self%spline2%n_coeffs

      select case (self%spline1%boundary_conditions)
      case (sll_p_boundary_periodic)
         offset1 = -degree1
      case (sll_p_boundary_clampeddiri:sll_p_boundary_clamped_cubic)
         offset1 = -1
      case default
         offset1 = 0
      end select

      upper(1) = n_cells(1) - degree1 + 1
      upper(2) = n_cells(2) - degree2 + 1
      do j = 1, n_cells(2)
         do i = 1, n_cells(1)
            cell = (j - 1)*n_cells(1) + i
            do l2 = 0, degree2
               select case (self%spline2%boundary_conditions)
               case (sll_p_boundary_periodic)
                  offset2 = modulo(j - degree2 + l2 - 1, n_coeffs(2))*n_coeffs(1) ! periodic
               case (sll_p_boundary_clampeddiri)
                  offset2 = (j + l2 - 2)*n_coeffs(1)
                  if (j == 1 .and. l2 == 0) then
                     pp_coeffs_local(:, 1) = 0.0_f64
                     !print*, 'a'
                     !continue
                     !goto 100
                     cycle
                  elseif (j == n_cells(2) .and. l2 == degree2) then
                     pp_coeffs_local(:, degree2 + 1) = 0.0_f64
                     !continue
                     !goto 100
                     cycle
                  end if
               case (sll_p_boundary_clamped_clampeddiri)
                  offset2 = (j + l2 - 1)*n_coeffs(1)
                  if (j == n_cells(2) .and. l2 == degree2) then
                     pp_coeffs_local(:, degree2 + 1) = 0.0_f64
                     !continue
                     !goto 100
                     cycle
                  end if
               case (sll_p_boundary_clampeddiri_clamped)
                  offset2 = (j + l2 - 2)*n_coeffs(1)
                  if (j == 1 .and. l2 == 0) then
                     pp_coeffs_local(:, 1) = 0.0_f64
                     !continue
                     !goto 100
                     cycle
                  end if
                  !case ( sll_p_boundary_clampeddiri:sll_p_boundary_clampeddiri_clamped)
                  !   offset2 = (j+l2-2)*n_coeffs(1)
               case default
                  offset2 = (j + l2 - 1)*n_coeffs(1)
               end select

               do l1 = 0, degree1
                  select case (self%spline1%boundary_conditions)
                  case (sll_p_boundary_periodic)
                     index = modulo(i + offset1 + l1 - 1, n_coeffs(1)) + 1
                     !print*, i, j, l1, l2, offset2, index
                     self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                  case (sll_p_boundary_clamped)
                     index = i + offset1 + l1
                     self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                     ! Set to zero for Dirichlet if we are at the boundary
                  case (sll_p_boundary_clamped_square)
                     index = i + offset1 + l1
                     self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                     ! Set to zero for Dirichlet if we are at the boundary
                  case (sll_p_boundary_clamped_cubic)
                     index = i + offset1 + l1
                     self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                     ! Set to zero for Dirichlet if we are at the boundary
                  case (sll_p_boundary_clamped_clampeddiri)
                     if (i == n_cells(1) .and. l1 == degree1) then
                        self%spline1%scratch_b(l1 + 1) = 0.0_f64
                     else
                        index = i + offset1 + l1
                        self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                     end if
                  case (sll_p_boundary_clampeddiri_clamped)
                     if (i == 1 .and. l1 == 0) then
                        self%spline1%scratch_b(l1 + 1) = 0.0_f64
                     else
                        index = i + offset1 + l1
                        self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                     end if
                  case (sll_p_boundary_clampeddiri)
                     if (i == 1 .and. l1 == 0) then
                        self%spline1%scratch_b(l1 + 1) = 0.0_f64
                     elseif (i == n_cells(1) .and. l1 == degree1) then
                        self%spline1%scratch_b(l1 + 1) = 0.0_f64
                     else
                        index = i + offset1 + l1
                        self%spline1%scratch_b(l1 + 1) = b_coeffs(offset2 + index)
                     end if
                  end select
               end do
               ! For clamped splines, we need to use the boundary coefficients
               select case (self%spline1%boundary_conditions)
               case (sll_p_boundary_periodic)
                  call sll_s_spline_pp_b_to_pp_1d_cella(degree1, self%spline1%poly_coeffs, &
                                                        self%spline1%scratch_b, pp_coeffs_local(:, l2 + 1))
               case default
                  if (i > degree1 - 1) then
                     if (i < upper(1) + 1) then
                        call sll_s_spline_pp_b_to_pp_1d_cella(degree1, self%spline1%poly_coeffs, &
                                                              self%spline1%scratch_b, pp_coeffs_local(:, l2 + 1))
                     else
                        call sll_s_spline_pp_b_to_pp_1d_cella(degree1, &
                                                              self%spline1%poly_coeffs_boundary_right(:, :, i - upper(1)), &
                                                              self%spline1%scratch_b, pp_coeffs_local(:, l2 + 1))
                     end if
                  else

                     call sll_s_spline_pp_b_to_pp_1d_cella(degree1, &
                                                           self%spline1%poly_coeffs_boundary_left(:, :, i), &
                                                           self%spline1%scratch_b, pp_coeffs_local(:, l2 + 1))
                  end if

               end select
               !call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, &
               !     self%spline1%scratch_b,pp_coeffs_local(:,l2+1))
               !100          continue
            end do
            do l1 = 0, degree1
               self%spline2%scratch_b = pp_coeffs_local(l1 + 1, :)
               select case (self%spline2%boundary_conditions)
               case (sll_p_boundary_periodic)
                  call sll_s_spline_pp_b_to_pp_1d_cella(degree2, self%spline2%poly_coeffs, &
                                                        self%spline2%scratch_b, self%spline2%scratch_p)
               case default
                  if (j > degree2 - 1) then
                     if (j < upper(2) + 1) then
                        call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                                                              self%spline2%poly_coeffs, &
                                                              self%spline2%scratch_b, self%spline2%scratch_p)
                     else
                        call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                                                              self%spline2%poly_coeffs_boundary_right(:, :, j - upper(2)), &
                                                              self%spline2%scratch_b, self%spline2%scratch_p)
                     end if
                  else

                     call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                                                           self%spline2%poly_coeffs_boundary_left(:, :, j), &
                                                           self%spline2%scratch_b, self%spline2%scratch_p)
                  end if

               end select

               !call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, &
               !     self%spline2%scratch_b,self%spline2%scratch_p)
               do l2 = 0, degree2
                  pp_coeffs(l2*(degree1 + 1) + l1 + 1, cell) = self%spline2%scratch_p(l2 + 1)
               end do
            end do
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_2d

!!$  !> Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions
!!$  subroutine sll_s_spline_pp_b_to_pp_2d_cell_clamped(spline1,spline2,n_cells, b_coeffs, pp_coeffs,i,j)
!!$    type( sll_t_spline_pp_1d), intent(inout)::  spline1 !< arbitrary degree 1d spline
!!$    type( sll_t_spline_pp_1d), intent(inout)::  spline2 !< arbitrary degree 1d spline
!!$    sll_int32, intent(in)    :: n_cells(2) !< number of gridcells
!!$    sll_real64,intent(in)    :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
!!$    sll_real64,intent(inout) :: pp_coeffs((spline1%degree+1)*(spline2%degree+1),n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
!!$
!!$    !> convert b-coefficients to pp-coefficients in first dimension
!!$    do j=1, n_cells(2)
!!$       do i=1, n_cells(1)
!!$          do l=0, degree(2)
!!$             spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+l)*n_cells(1):i+(j-degp2+l)*n_cells(1))
!!$
!!$  end subroutine sll_s_spline_pp_b_to_pp_2d_cell_clamped

   !> Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions
   subroutine sll_s_spline_pp_b_to_pp_2d_cell(spline1, spline2, n_cells, b_coeffs, pp_coeffs, i, j)
      type(sll_t_spline_pp_1d), intent(inout)::  spline1 !< arbitrary degree 1d spline
      type(sll_t_spline_pp_1d), intent(inout)::  spline2 !< arbitrary degree 1d spline
      sll_int32, intent(in)    :: n_cells(2) !< number of gridcells
      sll_real64, intent(in)    :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
      sll_real64, intent(inout) :: pp_coeffs((spline1%degree + 1)*(spline2%degree + 1), n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
      sll_int32, intent(in)    :: i, j !< indices
      sll_int32 :: k, l
      sll_int32 :: degree1, degree2, degp1, degp2
      degree1 = spline1%degree
      degree2 = spline2%degree
      degp1 = degree1 + 1
      degp2 = degree2 + 1
      !> convert b-coefficients in pp-coefficients in first dimension
      if (i > degree1) then
         if (j > degree2) then
            do l = 0, degree2
               spline1%scratch_b = b_coeffs(i - degree1 + (j - degp2 + l)*n_cells(1):i + (j - degp2 + l)*n_cells(1))
      call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1 + l*degp1:degp1*(l + 1), i + n_cells(1)*(j - 1)))
            end do
         else
            !> use of modulo for boundary cells in second dimension
            do l = 0, degree2
             spline1%scratch_b=b_coeffs(i-degree1+modulo(j-degp2+l,n_cells(2))*n_cells(1):i+modulo(j-degp2+l,n_cells(2))*n_cells(1))
      call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1 + l*degp1:degp1*(l + 1), i + n_cells(1)*(j - 1)))
            end do
         end if
      else
         !> use of modulo for boundary cells in both dimensions
         do l = 0, degree2
            do k = 0, degree1
           spline1%scratch_b(k + 1) = b_coeffs(modulo(i - degp1 + k, n_cells(1)) + 1 + modulo(j - degp2 + l, n_cells(2))*n_cells(1))
            end do

      call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1 + l*degp1:degp1*(l + 1), i + n_cells(1)*(j - 1)))
         end do
      end if

      !> convert b-coefficients in pp_coefficients in second dimension
      do k = 1, degp1
         do l = 1, degp2
            spline2%scratch_b(l) = pp_coeffs(k + degp1*(l - 1), i + n_cells(1)*(j - 1))
         end do
         call sll_s_spline_pp_b_to_pp_1d_cell(spline2, spline2%scratch_b, spline2%scratch_p)
         do l = 1, degp2
            pp_coeffs(k + degp1*(l - 1), i + n_cells(1)*(j - 1)) = spline2%scratch_p(l)
         end do
      end do
   end subroutine sll_s_spline_pp_b_to_pp_2d_cell

   !> Convert 3d spline in B form to spline in pp form with periodic boundary conditions
   subroutine sll_s_spline_pp_b_to_pp_3d(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 2d spline
      sll_int32, intent(in) ::  n_cells(3) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells(1)*n_cells(2), n_cells(2))   !< coefficients of spline in B-form
    sll_real64, intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32 :: i, j, k

      do k = 1, n_cells(3)
         do j = 1, n_cells(2)
            do i = 1, n_cells(1)
               call sll_s_spline_pp_b_to_pp_3d_cell(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
            end do
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_3d

   !> Convert 3d spline in B form in a cell to spline in pp form with periodic boundary conditions
   subroutine sll_s_spline_pp_b_to_pp_3d_cell(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 3d spline
      sll_int32, intent(in)    :: n_cells(3) !< number of gridcells
      sll_real64, intent(in)    :: b_coeffs(n_cells(1)*n_cells(2)*n_cells(3))   !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32, intent(in)    :: i, j, k !< indices
      !local variables
      sll_int32 :: l, m, n
      sll_int32 :: degree1, degree2, degree3, degp1, degp2, degp3
      degree1 = self%spline1%degree
      degree2 = self%spline2%degree
      degree3 = self%spline3%degree
      degp1 = degree1 + 1
      degp2 = degree2 + 1
      degp3 = degree3 + 1
      !> convert b-coefficients in pp-coefficients in first dimension
      if (i > degree1) then
         if (j > degree2) then
            if (k > degree3) then
               do m = 0, degree2
                  do n = 0, degree3
                   self%spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+m)*n_cells(1)+(k-degp3+n)*n_cells(1)*n_cells(2):i+(j-degp2+m)*n_cells(1)+(k-degp3+n)*n_cells(1)*n_cells(2))
                   call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2))) 
                  end do
               end do
            else
               !> use of modulo for boundary cells in third dimension
               do m = 0, degree2
                  do n = 0, degree3
                   self%spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+m)*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2):i+(j-degp2+m)*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2))
                   call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
                  end do
               end do
            end if
         else
            !> use of modulo for boundary cells in second and third dimension
            do m = 0, degree2
               do n = 0, degree3
                self%spline1%scratch_b=b_coeffs(i-degree1+modulo(j-degp2+m,n_cells(2))*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2):i+modulo(j-degp2+m,n_cells(2))*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2))
                call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
               end do
            end do
         end if
      else
         !> use of modulo for boundary cells in all three dimensions
         do m = 0, degree2
            do n = 0, degree3
               do l = 0, degree1
                self%spline1%scratch_b(l+1)=b_coeffs(modulo(i-degp1+l,n_cells(1))+1+modulo(j-degp2+m,n_cells(2))*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
            end do
         end do
      end if

      !> convert b-coefficients in pp-coefficients in second dimension
      do l = 1, degp1
         do n = 1, degp3
            do m = 1, degp2
               self%spline2%scratch_b(m) = self%scratch_pp(l + degp1*(m - 1) + degp1*degp2*(n - 1))
            end do
            call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, self%spline2%scratch_b, self%spline2%scratch_p)
            do m = 1, degp2
               self%scratch_pp(l + degp1*(m - 1) + degp1*degp2*(n - 1)) = self%spline2%scratch_p(m)
            end do
         end do
      end do

      !> convert b-coefficients in pp-coefficients in third dimension
      do l = 1, degp1
         do m = 1, degp2
            do n = 1, degp3
               self%spline3%scratch_b(n) = self%scratch_pp(l + degp1*(m - 1) + degp1*degp2*(n - 1))
            end do
            call sll_s_spline_pp_b_to_pp_1d_cell(self%spline3, self%spline3%scratch_b, self%spline3%scratch_p)
            do n = 1, degp3
               pp_coeffs(l+degp1*(m-1)+degp1*degp2*(n-1),i+n_cells(1)*(j-1)+n_cells(1)*n_cells(2)*(k-1))=self%spline3%scratch_p(n)
            end do
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_3d_cell

   !> Convert 3d spline in B form in a cell to spline in pp form with clamped boundary in first direction and periodic boundary conditions in the other directions
   subroutine sll_s_spline_pp_b_to_pp_3d_cella(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 3d spline
      sll_int32, intent(in)    :: n_cells(3) !< number of gridcells
      sll_real64, intent(in)    :: b_coeffs((n_cells(1) + self%spline1%degree)*n_cells(2)*n_cells(3)) !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32, intent(in)    :: i, j, k !< indices
      !local variables
      sll_int32 :: l, m, n
      sll_int32 :: degree1, degree2, degree3, degp1, degp2, degp3, n_dofs
      degree1 = self%spline1%degree
      degree2 = self%spline2%degree
      degree3 = self%spline3%degree
      degp1 = degree1 + 1
      degp2 = degree2 + 1
      degp3 = degree3 + 1
      n_dofs = n_cells(1) + degree1
      !> convert b-coefficients in pp-coefficients in first dimension
      if (i >= degree1 .and. i < n_cells(1) - degree1 + 2) then
         if (j > degree2) then
            if (k > degree3) then
               do m = 0, degree2
                  do n = 0, degree3
                   self%spline1%scratch_b=b_coeffs(i+(j-degp2+m)*n_dofs+(k-degp3+n)*n_dofs*n_cells(2):i+degree1+(j-degp2+m)*n_dofs+(k-degp3+n)*n_dofs*n_cells(2))
                   call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2))) 
                  end do
               end do
            else
               !> use of modulo for boundary cells in third dimension
               do m = 0, degree2
                  do n = 0, degree3
                   self%spline1%scratch_b=b_coeffs(i+(j-degp2+m)*n_dofs+modulo(k-degp3+n,n_cells(3))*n_dofs*n_cells(2):i+degree1+(j-degp2+m)*n_dofs+modulo(k-degp3+n,n_cells(3))*n_dofs*n_cells(2))
                   call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
                  end do
               end do
            end if
         else
            !> use of modulo for boundary cells in second and third dimension
            do m = 0, degree2
               do n = 0, degree3
                self%spline1%scratch_b=b_coeffs(i+modulo(j-degp2+m,n_cells(2))*n_dofs+modulo(k-degp3+n,n_cells(3))*n_dofs*n_cells(2):i+degree1+modulo(j-degp2+m,n_cells(2))*n_dofs+modulo(k-degp3+n,n_cells(3))*n_dofs*n_cells(2))
                call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
               end do
            end do
         end if
      else if (i < degree1) then
         !> use of modulo for boundary cells in all three dimensions
         do m = 0, degree2
            do n = 0, degree3
               do l = 0, degree1
        self%spline1%scratch_b(l+1)=b_coeffs(i+l+modulo(j-degp2+m,n_cells(2))*n_dofs+modulo(k-degp3+n,n_cells(3))*n_dofs*n_cells(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(degree1, self%spline1%poly_coeffs_boundary_left(:,:,i),self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
            end do
         end do
      else if (i >= n_cells(1) - degree1 + 2) then
         !> use of modulo for boundary cells in all three dimensions
         do m = 0, degree2
            do n = 0, degree3
               do l = 0, degree1
        self%spline1%scratch_b(l+1)=b_coeffs(i+l+modulo(j-degp2+m,n_cells(2))*n_dofs+modulo(k-degp3+n,n_cells(3))*n_dofs*n_cells(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(degree1, self%spline1%poly_coeffs_boundary_right(:,:,i-(n_cells(1) - degree1 +1)),self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
            end do
         end do
      end if

      !> convert b-coefficients in pp-coefficients in second dimension
      do l = 1, degp1
         do n = 1, degp3
            do m = 1, degp2
               self%spline2%scratch_b(m) = self%scratch_pp(l + degp1*(m - 1) + degp1*degp2*(n - 1))
            end do
            call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, self%spline2%scratch_b, self%spline2%scratch_p)
            do m = 1, degp2
               self%scratch_pp(l + degp1*(m - 1) + degp1*degp2*(n - 1)) = self%spline2%scratch_p(m)
            end do
         end do
      end do

      !> convert b-coefficients in pp-coefficients in third dimension
      do l = 1, degp1
         do m = 1, degp2
            do n = 1, degp3
               self%spline3%scratch_b(n) = self%scratch_pp(l + degp1*(m - 1) + degp1*degp2*(n - 1))
            end do
            call sll_s_spline_pp_b_to_pp_1d_cell(self%spline3, self%spline3%scratch_b, self%spline3%scratch_p)
            do n = 1, degp3
               pp_coeffs(l+degp1*(m-1)+degp1*degp2*(n-1),i+n_cells(1)*(j-1)+n_cells(1)*n_cells(2)*(k-1))=self%spline3%scratch_p(n)
            end do
         end do
      end do
   end subroutine sll_s_spline_pp_b_to_pp_3d_cella

   !> Convert 3d spline in B form in a cell to spline in pp form with clamped boundary in all three directions
   subroutine sll_s_spline_pp_b_to_pp_3d_cella2f(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 3d spline
      sll_int32, intent(in)    :: n_cells(3) !< number of gridcells
      sll_real64, intent(in)    :: b_coeffs((n_cells(1) + self%spline1%degree)*n_cells(2)*(n_cells(3) + self%spline3%degree)) !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32, intent(in)    :: i, j, k !< indices
      !local variables
      sll_int32 :: l, m, n
      sll_int32 :: deg(3), degp(3), n_dofs(3), upper(3)
      deg(1) = self%spline1%degree
      deg(2) = self%spline2%degree
      deg(3) = self%spline3%degree
      degp = deg + 1
      n_dofs = n_cells + deg
      n_dofs(2) = n_cells(2)
      upper = n_cells - deg + 1
      !> convert b-coefficients in pp-coefficients in first dimension
      if (i >= deg(1) .and. i <= upper(1)) then
         if (j > deg(2)) then
            do n = 0, deg(3)
               do m = 0, deg(2)
                self%spline1%scratch_b=b_coeffs(i+(j-degp(2)+m)*n_dofs(1)+(k-1+n)*n_dofs(1)*n_dofs(2):i+deg(1)+(j-degp(2)+m)*n_dofs(1)+(k-1+n)*n_dofs(1)*n_dofs(2))
                call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2)))) 
               end do
            end do
         else
            !> use of modulo for boundary cells in second dimension
            do n = 0, deg(3)
               do m = 0, deg(2)
                self%spline1%scratch_b=b_coeffs(i+modulo(j-degp(2)+m,n_cells(2))*n_dofs(1)+(k-1+n)*n_dofs(1)*n_dofs(2):i+deg(1)+modulo(j-degp(2)+m,n_cells(2))*n_dofs(1)+(k-1+n)*n_dofs(1)*n_dofs(2))
                call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2)))) 
               end do
            end do
         end if
      else if (i < deg(1)) then
         do n = 0, deg(3)
            do m = 0, deg(2)
               do l = 0, deg(1)
   self%spline1%scratch_b(l + 1) = b_coeffs(i + l + modulo(j - degp(2) + m, n_cells(2))*n_dofs(1) + (k + n - 1)*n_dofs(1)*n_dofs(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(1), self%spline1%poly_coeffs_boundary_left(:,:,i),self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2))))
            end do
         end do
      else if (i > upper(1)) then
         do n = 0, deg(3)
            do m = 0, deg(2)
               do l = 0, deg(1)
   self%spline1%scratch_b(l + 1) = b_coeffs(i + l + modulo(j - degp(2) + m, n_cells(2))*n_dofs(1) + (k + n - 1)*n_dofs(1)*n_dofs(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(1), self%spline1%poly_coeffs_boundary_right(:,:,i-upper(1)),self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2))))
            end do
         end do
      end if

      !> convert b-coefficients in pp-coefficients in second dimension
      do n = 0, deg(3)
         do l = 0, deg(1)
            do m = 0, deg(2)
               self%spline2%scratch_b(m + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
            end do
            call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, self%spline2%scratch_b, self%spline2%scratch_p)
            do m = 0, deg(2)
               self%scratch_pp(1 + l + m*degp(1) + n*degp(1)*degp(2)) = self%spline2%scratch_p(m + 1)
            end do
         end do
      end do

      !> convert b-coefficients in pp-coefficients in third dimension
      if (k >= deg(3) .and. k <= upper(3)) then
         do m = 0, deg(2)
            do l = 0, deg(1)
               do n = 0, deg(3)
                  self%spline3%scratch_b(n + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
               call sll_s_spline_pp_b_to_pp_1d_cell(self%spline3, self%spline3%scratch_b, self%spline3%scratch_p)
               do n = 0, deg(3)
         pp_coeffs(1+l + m*degp(1)+ n*degp(1)*degp(2), i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = self%spline3%scratch_p(n+1)
               end do
            end do
         end do
      else if (k < deg(3)) then
         do m = 0, deg(2)
            do l = 0, deg(1)
               do n = 0, deg(3)
                  self%spline3%scratch_b(n + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(3), self%spline3%poly_coeffs_boundary_left(:,:,k),self%spline3%scratch_b, self%spline3%scratch_p)
               do n = 0, deg(3)
         pp_coeffs(1+l + m*degp(1)+ n*degp(1)*degp(2), i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = self%spline3%scratch_p(n+1)
               end do
            end do
         end do
      else if (k > upper(3)) then
         do m = 0, deg(2)
            do l = 0, deg(1)
               do n = 0, deg(3)
                  self%spline3%scratch_b(n + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(3), self%spline3%poly_coeffs_boundary_right(:,:,k-upper(3)),self%spline3%scratch_b, self%spline3%scratch_p)
               do n = 0, deg(3)
         pp_coeffs(1+l + m*degp(1)+ n*degp(1)*degp(2), i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = self%spline3%scratch_p(n+1)
               end do
            end do
         end do
      end if

   end subroutine sll_s_spline_pp_b_to_pp_3d_cella2f

   !> Convert 3d spline in B form in a cell to spline in pp form with clamped boundary in all three directions
   subroutine sll_s_spline_pp_b_to_pp_3d_cellaf(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 3d spline
      sll_int32, intent(in)    :: n_cells(3) !< number of gridcells
    sll_real64,intent(in)    :: b_coeffs((n_cells(1)+self%spline1%degree)*(n_cells(2)+self%spline2%degree)*(n_cells(3)+self%spline3%degree)) !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32, intent(in)    :: i, j, k !< indices
      !local variables
      sll_int32 :: l, m, n
      sll_int32 :: deg(3), degp(3), n_dofs(3), upper(3)
      deg(1) = self%spline1%degree
      deg(2) = self%spline2%degree
      deg(3) = self%spline3%degree
      degp = deg + 1
      n_dofs = n_cells + deg
      upper = n_cells - deg + 1
      !> convert b-coefficients in pp-coefficients in first dimension
      if (i >= deg(1) .and. i <= upper(1)) then
         do n = 0, deg(3)
            do m = 0, deg(2)
             self%spline1%scratch_b=b_coeffs(i+(j-1+m)*n_dofs(1)+(k-1+n)*n_dofs(1)*n_dofs(2):i+deg(1)+(j-1+m)*n_dofs(1)+(k-1+n)*n_dofs(1)*n_dofs(2))
             call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2)))) 
            end do
         end do
      else if (i < deg(1)) then
         do n = 0, deg(3)
            do m = 0, deg(2)
               do l = 0, deg(1)
                  self%spline1%scratch_b(l + 1) = b_coeffs(i + l + (j + m - 1)*n_dofs(1) + (k + n - 1)*n_dofs(1)*n_dofs(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(1), self%spline1%poly_coeffs_boundary_left(:,:,i),self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2))))
            end do
         end do
      else if (i > upper(1)) then
         do n = 0, deg(3)
            do m = 0, deg(2)
               do l = 0, deg(1)
                  self%spline1%scratch_b(l + 1) = b_coeffs(i + l + (j + m - 1)*n_dofs(1) + (k + n - 1)*n_dofs(1)*n_dofs(2))
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(1), self%spline1%poly_coeffs_boundary_right(:,:,i-upper(1)),self%spline1%scratch_b, self%scratch_pp(1+m*degp(1)+n*degp(1)*degp(2):degp(1)*(1+m+n*degp(2))))
            end do
         end do
      end if

      !> convert b-coefficients in pp-coefficients in second dimension
      if (j >= deg(2) .and. j <= upper(2)) then
         do n = 0, deg(3)
            do l = 0, deg(1)
               do m = 0, deg(2)
                  self%spline2%scratch_b(m + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
               call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, self%spline2%scratch_b, self%spline2%scratch_p)
               do m = 0, deg(2)
                  self%scratch_pp(1 + l + m*degp(1) + n*degp(1)*degp(2)) = self%spline2%scratch_p(m + 1)
               end do
            end do
         end do
      else if (j < deg(2)) then
         do n = 0, deg(3)
            do l = 0, deg(1)
               do m = 0, deg(2)
                  self%spline2%scratch_b(m + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(2), self%spline2%poly_coeffs_boundary_left(:,:,j),self%spline2%scratch_b, self%spline2%scratch_p)
               do m = 0, deg(2)
                  self%scratch_pp(1 + l + m*degp(1) + n*degp(1)*degp(2)) = self%spline2%scratch_p(m + 1)
               end do
            end do
         end do
      else if (j > upper(2)) then
         do n = 0, deg(3)
            do l = 0, deg(1)
               do m = 0, deg(2)
                  self%spline2%scratch_b(m + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(2), self%spline2%poly_coeffs_boundary_right(:,:,j-upper(2)),self%spline2%scratch_b, self%spline2%scratch_p)
               do m = 0, deg(2)
                  self%scratch_pp(1 + l + m*degp(1) + n*degp(1)*degp(2)) = self%spline2%scratch_p(m + 1)
               end do
            end do
         end do
      end if

      !> convert b-coefficients in pp-coefficients in third dimension
      if (k >= deg(3) .and. k <= upper(3)) then
         do m = 0, deg(2)
            do l = 0, deg(1)
               do n = 0, deg(3)
                  self%spline3%scratch_b(n + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
               call sll_s_spline_pp_b_to_pp_1d_cell(self%spline3, self%spline3%scratch_b, self%spline3%scratch_p)
               do n = 0, deg(3)
         pp_coeffs(1+l + m*degp(1)+ n*degp(1)*degp(2), i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = self%spline3%scratch_p(n+1)
               end do
            end do
         end do
      else if (k < deg(3)) then
         do m = 0, deg(2)
            do l = 0, deg(1)
               do n = 0, deg(3)
                  self%spline3%scratch_b(n + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(3), self%spline3%poly_coeffs_boundary_left(:,:,k),self%spline3%scratch_b, self%spline3%scratch_p)
               do n = 0, deg(3)
         pp_coeffs(1+l + m*degp(1)+ n*degp(1)*degp(2), i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = self%spline3%scratch_p(n+1)
               end do
            end do
         end do
      else if (k > upper(3)) then
         do m = 0, deg(2)
            do l = 0, deg(1)
               do n = 0, deg(3)
                  self%spline3%scratch_b(n + 1) = self%scratch_pp(1 + l + degp(1)*m + degp(1)*degp(2)*n)
               end do
             call sll_s_spline_pp_b_to_pp_1d_cella(deg(3), self%spline3%poly_coeffs_boundary_right(:,:,k-upper(3)),self%spline3%scratch_b, self%spline3%scratch_p)
               do n = 0, deg(3)
         pp_coeffs(1+l + m*degp(1)+ n*degp(1)*degp(2), i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = self%spline3%scratch_p(n+1)
               end do
            end do
         end do
      end if

   end subroutine sll_s_spline_pp_b_to_pp_3d_cellaf

   !> Convert 3d spline in B form to spline in pp form for clamped spline
   subroutine sll_s_spline_pp_b_to_pp_3d_clamped(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs((n_cells(1) + self%spline1%degree)*n_cells(2)*n_cells(3))   !< coefficients of spline in B-form
    sll_real64,intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32 :: i, j, k

      do k = 1, n_cells(3)
         do j = 1, n_cells(2)
            do i = 1, n_cells(1)
               call sll_s_spline_pp_b_to_pp_3d_cella(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
            end do
         end do
      end do
   end subroutine sll_s_spline_pp_b_to_pp_3d_clamped

   !> Convert 3d spline in B form to spline in pp form for clamped spline
   subroutine sll_s_spline_pp_b_to_pp_3d_clamped_2full(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs((n_cells(1) + self%spline1%degree)*n_cells(2)*(n_cells(3) + self%spline3%degree))   !< coefficients of spline in B-form
    sll_real64,intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32 :: i, j, k

      do k = 1, n_cells(3)
         do j = 1, n_cells(2)
            do i = 1, n_cells(1)
               call sll_s_spline_pp_b_to_pp_3d_cella2f(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
            end do
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_3d_clamped_2full

   !> Convert 3d spline in B form to spline in pp form for clamped spline
   subroutine sll_s_spline_pp_b_to_pp_3d_clamped_full(self, n_cells, b_coeffs, pp_coeffs)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
    sll_real64,intent(in) :: b_coeffs((n_cells(1)+self%spline1%degree)*(n_cells(2)+self%spline2%degree)*(n_cells(3)+self%spline3%degree))   !< coefficients of spline in B-form
    sll_real64,intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
      sll_int32 :: i, j, k

      do k = 1, n_cells(3)
         do j = 1, n_cells(2)
            do i = 1, n_cells(1)
               call sll_s_spline_pp_b_to_pp_3d_cellaf(self, n_cells, b_coeffs, pp_coeffs, i, j, k)
            end do
         end do
      end do

   end subroutine sll_s_spline_pp_b_to_pp_3d_clamped_full

   subroutine sll_s_spline_evaluate_basis_b_form_1d_clamped(self, n_cells, b_coeffs, val)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells + self%degree) !< coefficients of spline in B-form
      sll_real64, intent(out):: val(:) !< array of values
      !local variables
      sll_real64 :: pp_coeffs(self%degree + 1, n_cells) !< coefficients of sp

      call sll_s_spline_pp_b_to_pp_1d_clamped(self, n_cells, b_coeffs, pp_coeffs)

      call sll_s_spline_evaluate_basis_pp_form_1d(self, n_cells, pp_coeffs, val)

   end subroutine sll_s_spline_evaluate_basis_b_form_1d_clamped

   subroutine sll_s_spline_evaluate_basis_b_form_1d_periodic(self, n_cells, b_coeffs, val)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells) !< coefficients of spline in B-form
      sll_real64, intent(out):: val(:) !< array of values
      !local variables
      sll_real64 :: pp_coeffs(self%degree + 1, n_cells) !< coefficients of sp

      call sll_s_spline_pp_b_to_pp_1d(self, n_cells, b_coeffs, pp_coeffs)

      call sll_s_spline_evaluate_basis_pp_form_1d(self, n_cells, pp_coeffs, val)

   end subroutine sll_s_spline_evaluate_basis_b_form_1d_periodic

   subroutine sll_s_spline_evaluate_basis_pp_form_1d(self, n_cells, pp_coeffs, val)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells !< number of gridcells
      sll_real64, intent(in) :: pp_coeffs(self%degree + 1, n_cells) !< coefficients of sp
      sll_real64, intent(out):: val(:) !< array of values
      !local variables
      sll_int32 :: i

      do i = 1, n_cells
         val(i) = pp_coeffs(self%degree + 1, i)
         !val(i)=sll_f_spline_pp_horner_1d(self%degree, pp_coeffs, 0._f64, i)
      end do
      val(n_cells + 1) = sll_f_spline_pp_horner_1d(self%degree, pp_coeffs, 1._f64, n_cells)

   end subroutine sll_s_spline_evaluate_basis_pp_form_1d

   subroutine sll_s_spline_evaluate_basis_b_form_3d_clamped(self, n_cells, b_coeffs, val)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells(3) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs((n_cells(1) + self%spline1%degree)*n_cells(2)*n_cells(3)) !< coefficients of spline in B-form
      sll_real64, intent(out):: val(:) !< array of values

      call sll_s_spline_pp_b_to_pp_3d_clamped(self, n_cells, b_coeffs, self%scratch_coeffs)

      call sll_s_spline_evaluate_basis_pp_form_3d(self, n_cells, self%scratch_coeffs, val)

   end subroutine sll_s_spline_evaluate_basis_b_form_3d_clamped

   subroutine sll_s_spline_evaluate_basis_b_form_3d_periodic(self, n_cells, b_coeffs, val)
      type(sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells(3) !< number of gridcells
      sll_real64, intent(in) :: b_coeffs(n_cells(1)*n_cells(2)*n_cells(3)) !< coefficients of spline in B-form
      sll_real64, intent(out):: val(:) !< array of values

      call sll_s_spline_pp_b_to_pp_3d(self, n_cells, b_coeffs, self%scratch_coeffs)

      call sll_s_spline_evaluate_basis_pp_form_3d(self, n_cells, self%scratch_coeffs, val)

   end subroutine sll_s_spline_evaluate_basis_b_form_3d_periodic

   subroutine sll_s_spline_evaluate_basis_pp_form_3d(self, n_cells, pp_coeffs, val)
      type(sll_t_spline_pp_3d), intent(in)::  self !< arbitrary degree 1d spline
      sll_int32, intent(in) ::  n_cells(3) !< number of gridcells
    sll_real64,intent(in) :: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3)) !< coefficients of sp
      sll_real64, intent(out):: val(:) !< array of values
      !local variables
      sll_int32 :: i, j, k

      do k = 1, n_cells(3)
         do j = 1, n_cells(2)
            do i = 1, n_cells(1)
               !val(i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2)) = pp_coeffs(self%spline1%degree+1,i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))
             val(i+(j-1)*(n_cells(1)+1)+(k-1)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [0._f64,0._f64,0._f64], [i,j,k], n_cells)
            end do
          val(n_cells(1)+1+(j-1)*(n_cells(1)+1)+(k-1)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [1._f64,0._f64,0._f64], [n_cells(1),j,k], n_cells)
         end do
       val(n_cells(1)+1+n_cells(2)*(n_cells(1)+1)+(k-1)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [1._f64,1._f64,0._f64], [n_cells(1),n_cells(2),k], n_cells)
      end do
    val(n_cells(1)+1+n_cells(2)*(n_cells(1)+1)+n_cells(3)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [1._f64,1._f64,1._f64], [n_cells(1),n_cells(2),n_cells(3)], n_cells)

      do i = 1, n_cells(1)
         do k = 1, n_cells(3)
          val(i+n_cells(2)*(n_cells(1)+1)+(k-1)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [0._f64,1._f64,0._f64], [i,n_cells(2),k], n_cells)
         end do
       val(i+n_cells(2)*(n_cells(1)+1)+n_cells(3)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [0._f64,1._f64,1._f64], [i,n_cells(2),n_cells(3)], n_cells)
      end do
      do j = 1, n_cells(2)
         do i = 1, n_cells(1)
          val(i+(j-1)*(n_cells(1)+1)+n_cells(3)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [0._f64,0._f64,1._f64], [i,j,n_cells(3)], n_cells)
         end do
       val(n_cells(1)+1+(j-1)*(n_cells(1)+1)+n_cells(3)*(n_cells(1)+1)*(n_cells(2)+1)) = sll_f_spline_pp_horner_3d([self%spline1%degree, self%spline2%degree, self%spline3%degree], pp_coeffs, [1._f64,0._f64,1._f64], [n_cells(1),j,n_cells(3)], n_cells)
      end do

   end subroutine sll_s_spline_evaluate_basis_pp_form_3d

   !> Perform a 1d hornerschema on the poly_coeffs
   subroutine sll_s_spline_pp_horner_m_1d(self, val, degree, x)
      type(sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline
      sll_real64, intent(out):: val(:) !< array of values
      sll_int32, intent(in)  :: degree !< degree of the spline
      sll_real64, intent(in) :: x !< point at which we evaluate our spline
      sll_int32 :: i

      do i = 1, degree + 1!size(val)
         val(i) = sll_f_spline_pp_horner_1d(degree, self%poly_coeffs, x, i)
      end do
   end subroutine sll_s_spline_pp_horner_m_1d

   !> Perform two times a 1d hornerschema on the poly_coeffs
   subroutine sll_s_spline_pp_horner_m_2d(self, val, degree, x)
      type(sll_t_spline_pp_2d), intent(in)::  self !< arbitrary degree 2d spline
      sll_real64, intent(out):: val(:, :) !< array of values
      sll_int32, intent(in)  :: degree(2) !< degree of the spline
      sll_real64, intent(in) :: x(2) !< point at which we evaluate our spline
      sll_int32 :: i

      do i = 1, degree(1) + 1!size(val,1)
         val(i, 1) = sll_f_spline_pp_horner_1d(degree(1), self%spline1%poly_coeffs, x(1), i)
      end do
      do i = 1, degree(2) + 1!size(val,2)
         val(i, 2) = sll_f_spline_pp_horner_1d(degree(2), self%spline2%poly_coeffs, x(2), i)
      end do
   end subroutine sll_s_spline_pp_horner_m_2d

   !> Perform three times a 1d hornerschema on the poly_coeffs
   subroutine sll_s_spline_pp_horner_m_3d(self, val, degree, x)
      type(sll_t_spline_pp_3d), intent(in)::  self !< arbitrary degree 3d spline
      sll_real64, intent(out):: val(:, :) !< array of values
      sll_int32, intent(in)  :: degree(3) !< degree of the spline
      sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
      sll_int32 :: i

      do i = 1, degree(1) + 1!size(val,1)
         val(i, 1) = sll_f_spline_pp_horner_1d(degree(1), self%spline1%poly_coeffs, x(1), i)
      end do
      do i = 1, degree(2) + 1!size(val,2)
         val(i, 2) = sll_f_spline_pp_horner_1d(degree(2), self%spline2%poly_coeffs, x(2), i)
      end do
      do i = 1, degree(3) + 1!size(val,3)
         val(i, 3) = sll_f_spline_pp_horner_1d(degree(3), self%spline3%poly_coeffs, x(3), i)
      end do
   end subroutine sll_s_spline_pp_horner_m_3d

   !> Perform a 1d hornerschema on the pp_coeffs evaluate at x
   subroutine sll_s_spline_pp_horner_primitive_1d(val, degree, pp_coeffs, x)
      sll_real64, intent(out):: val(:) !< array of values
      sll_int32, intent(in)  :: degree !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x !< point at which we evaluate our spline
      sll_int32 :: i

      do i = 1, size(val)
         val(i) = sll_f_spline_pp_horner_1d(degree, pp_coeffs, x, i)*x
      end do
   end subroutine sll_s_spline_pp_horner_primitive_1d

   !> Perform a 1d hornerschema on the pp_coeffs at index
   function sll_f_spline_pp_horner_1d(degree, pp_coeffs, x, index) result(res)
      sll_int32, intent(in) :: degree !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x !< point at which we evaluate our spline
      sll_int32, intent(in) :: index !< index of cell in which is x
      sll_real64 :: res !< value of the splinefunction at point x
      sll_int32 :: i

      !print*, pp_coeffs(:,index)
      !print*, x
      res = pp_coeffs(1, index)
      do i = 1, degree
         res = res*x + pp_coeffs(i + 1, index)
      end do
   end function sll_f_spline_pp_horner_1d

   !> Perform a 1d hornerschema on the pp_coeffs at index
   function sll_f_spline_pp_horner_derivative_1d(degree, pp_coeffs, x, index) result(res)
      sll_int32, intent(in) :: degree !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x !< point at which we evaluate our spline
      sll_int32, intent(in) :: index !< index of cell in which is x
      sll_real64 :: res !< value of the splinefunction at point x
      sll_int32 :: i

      !print*, pp_coeffs(:,index)
      !print*, x
      res = pp_coeffs(1, index)*degree
      do i = 1, degree - 1
         res = res*x + pp_coeffs(i + 1, index)*(degree - i)
      end do
   end function sll_f_spline_pp_horner_derivative_1d

   !> Perform a 2d hornerschema on the pp_coeffs at the indices
   function sll_f_spline_pp_horner_2d(degree, pp_coeffs, x, indices, n_cells) result(res)
      sll_int32, intent(in) :: degree(2) !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x(2) !< point at which we evaluate our spline
      sll_int32, intent(in) :: indices(2) !< indices of cell in which is x
      sll_int32, intent(in) :: n_cells(2) !< number of gridcells
      sll_real64 :: res !< value of the splinefunction at point x
      sll_real64 :: pp_coeffs_1d(degree(2) + 1, 1)
      sll_int32  :: i
      !> Perform a 1d hornerschema in the first dimension
      do i = 0, degree(2)
       pp_coeffs_1d(i+1,1)=sll_f_spline_pp_horner_1d(degree(1), pp_coeffs(1+i*(degree(1)+1):(degree(1)+1)*(i+1),1+(indices(2)-1)*n_cells(1):n_cells(1)*indices(2)), x(1),indices(1))
      end do
      !> Perform a 1d hornerschema in the second dimension
      res = sll_f_spline_pp_horner_1d(degree(2), pp_coeffs_1d, x(2), 1)
   end function sll_f_spline_pp_horner_2d

   !> Perform a 3d hornerschema on the pp_coeffs at the indices
   function sll_f_spline_pp_horner_3d(degree, pp_coeffs, x, indices, n_cells) result(res)
      sll_int32, intent(in) :: degree(3) !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
      sll_int32, intent(in) :: indices(3) !< indices of cell in which is x
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_real64 :: res !< value of the splinefunction at point x
      sll_real64 :: pp_coeffs_2d((degree(2) + 1)*(degree(3) + 1), 1)
      sll_int32  :: i, j
      sll_int32  :: degp1, degp2
      degp1 = degree(1) + 1
      degp2 = degree(2) + 1
      !> Perform a 1d hornerschema in the first dimension
      do j = 0, degree(3)
         do i = 0, degree(2)
          pp_coeffs_2d(1+i+j*degp2,1)= sll_f_spline_pp_horner_1d(degree(1),pp_coeffs(1+i*degp1+j*degp1*degp2:degp1+i*degp1+j*degp1*degp2,1+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2)):n_cells(1)*(indices(2)+(indices(3)-1)*n_cells(2))),x(1),indices(1))
         end do
      end do
      !> Perform a 2d hornerschema in the second and third dimension
      res = sll_f_spline_pp_horner_2d(degree(2:3), pp_coeffs_2d, x(2:3), [1, 1], [1, 1])

   end function sll_f_spline_pp_horner_3d

   !> Perform a 3d hornerschema on the pp_coeffs at the indices
   function sll_f_spline_pp_horner_3d_d1(degree, pp_coeffs, x, indices, n_cells) result(res)
      sll_int32, intent(in) :: degree(3) !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
      sll_int32, intent(in) :: indices(3) !< indices of cell in which is x
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_real64 :: res !< value of the splinefunction at point x
      sll_real64 :: pp_coeffs_2d((degree(2) + 1)*(degree(3) + 1), 1)
      sll_real64 :: pp_coeffs_1d(degree(3) + 1, 1)
      sll_int32  :: i, j
      sll_int32  :: degp1, degp2
      degp1 = degree(1) + 1
      degp2 = degree(2) + 1
      !> Perform a 1d hornerschema in the first dimension
      do j = 0, degree(3)
         do i = 0, degree(2)
          pp_coeffs_2d(1+i+j*degp2,1)= sll_f_spline_pp_horner_derivative_1d(degree(1),pp_coeffs(1+i*degp1+j*degp1*degp2:degp1+i*degp1+j*degp1*degp2,indices(1)+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2)):indices(1)+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2))), x(1), 1)
         end do
      end do
      !> Perform a 2d hornerschema in the second and third dimension
      res = sll_f_spline_pp_horner_2d(degree(2:3), pp_coeffs_2d, x(2:3), [1, 1], [1, 1])

      !> Perform a 1d hornerschema in the first dimension
      do i = 0, degree(3)
         pp_coeffs_1d(i+1,1)=sll_f_spline_pp_horner_1d(degree(2), pp_coeffs_2d(1+i*(degree(2)+1):(degree(2)+1)*(i+1),:), x(2), 1)
      end do
      !> Perform a 1d hornerschema in the second dimension
      res = sll_f_spline_pp_horner_1d(degree(3), pp_coeffs_1d, x(3), 1)*n_cells(1)

   end function sll_f_spline_pp_horner_3d_d1

   !> Perform a 3d hornerschema on the pp_coeffs at the indices
   function sll_f_spline_pp_horner_3d_d2(degree, pp_coeffs, x, indices, n_cells) result(res)
      sll_int32, intent(in) :: degree(3) !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
      sll_int32, intent(in) :: indices(3) !< indices of cell in which is x
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_real64 :: res !< value of the splinefunction at point x
      sll_real64 :: pp_coeffs_2d((degree(2) + 1)*(degree(3) + 1), 1)
      sll_real64 :: pp_coeffs_1d(degree(3) + 1, 1)
      sll_int32  :: i, j
      sll_int32  :: degp1, degp2
      degp1 = degree(1) + 1
      degp2 = degree(2) + 1
      !> Perform a 1d hornerschema in the first dimension
      do j = 0, degree(3)
         do i = 0, degree(2)
          pp_coeffs_2d(1+i+j*degp2,1)= sll_f_spline_pp_horner_1d(degree(1),pp_coeffs(1+i*degp1+j*degp1*degp2:degp1+i*degp1+j*degp1*degp2,indices(1)+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2)):indices(1)+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2))), x(1), 1)
         end do
      end do
      !> Perform a 2d hornerschema in the second and third dimension
      res = sll_f_spline_pp_horner_2d(degree(2:3), pp_coeffs_2d, x(2:3), [1, 1], [1, 1])

      !> Perform a 1d hornerschema in the first dimension
      do i = 0, degree(3)
 pp_coeffs_1d(i+1,1)=sll_f_spline_pp_horner_derivative_1d(degree(2), pp_coeffs_2d(1+i*(degree(2)+1):(degree(2)+1)*(i+1),:), x(2), 1)
      end do
      !> Perform a 1d hornerschema in the second dimension
      res = sll_f_spline_pp_horner_1d(degree(3), pp_coeffs_1d, x(3), 1)*n_cells(2)

   end function sll_f_spline_pp_horner_3d_d2

   !> Perform a 3d hornerschema on the pp_coeffs at the indices
   function sll_f_spline_pp_horner_3d_d3(degree, pp_coeffs, x, indices, n_cells) result(res)
      sll_int32, intent(in) :: degree(3) !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
      sll_int32, intent(in) :: indices(3) !< indices of cell in which is x
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_real64 :: res !< value of the splinefunction at point x
      sll_real64 :: pp_coeffs_2d((degree(2) + 1)*(degree(3) + 1), 1)
      sll_real64 :: pp_coeffs_1d(degree(3) + 1, 1)
      sll_int32  :: i, j
      sll_int32  :: degp1, degp2
      degp1 = degree(1) + 1
      degp2 = degree(2) + 1
      !> Perform a 1d hornerschema in the first dimension
      do j = 0, degree(3)
         do i = 0, degree(2)
          pp_coeffs_2d(1+i+j*degp2,1)= sll_f_spline_pp_horner_1d(degree(1),pp_coeffs(1+i*degp1+j*degp1*degp2:degp1+i*degp1+j*degp1*degp2,indices(1)+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2)):indices(1)+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2))), x(1), 1)
         end do
      end do
      !> Perform a 2d hornerschema in the second and third dimension
      res = sll_f_spline_pp_horner_2d(degree(2:3), pp_coeffs_2d, x(2:3), [1, 1], [1, 1])

      !> Perform a 1d hornerschema in the first dimension
      do i = 0, degree(3)
         pp_coeffs_1d(i+1,1)=sll_f_spline_pp_horner_1d(degree(2), pp_coeffs_2d(1+i*(degree(2)+1):(degree(2)+1)*(i+1), :), x(2), 1)
      end do
      !> Perform a 1d hornerschema in the second dimension
      res = sll_f_spline_pp_horner_derivative_1d(degree(3), pp_coeffs_1d, x(3), 1)*n_cells(3)

   end function sll_f_spline_pp_horner_3d_d3

   !> Perform a 3d hornerschema on the pp_coeffs at the indices
   function sll_f_spline_pp_horner_derivative_3d(degree, pp_coeffs, x, indices, n_cells, component) result(res)
      sll_int32, intent(in) :: degree(3) !< degree of the spline
      sll_real64, intent(in) :: pp_coeffs(:, :)  !< coefficients of spline in pp-form
      sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
      sll_int32, intent(in) :: indices(3) !< indices of cell in which is x
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells
      sll_int32, intent(in) :: component
      sll_real64 :: res !< value of the splinefunction at point x
      sll_real64, allocatable :: coeffs(:, :), pp_coeffs_2d(:, :)
      sll_real64 :: exponents(maxval(degree) + 1, 3)
      sll_int32  :: h, i, j
      sll_int32  :: deg(3), degp1, degp2
      degp1 = degree(1) + 1
      degp2 = degree(2) + 1

      deg = degree
      deg(component) = degree(component) - 1
      allocate (coeffs((deg(1) + 1), 1))
      allocate (pp_coeffs_2d((deg(2) + 1)*(deg(3) + 1), 1))

      exponents = 1._f64
      do i = 1, degree(component)
         exponents(i, component) = degree(component) + 1 - i
      end do
      !> Perform a 1d hornerschema in the first dimension
      do j = 0, deg(3)
         do i = 0, deg(2)
            do h = 0, deg(1)
             coeffs(h+1, 1) = pp_coeffs(1+h+i*degp1+j*degp1*degp2, indices(1) + (indices(2)-1)*n_cells(1)+ (indices(3)-1)*n_cells(1)*n_cells(2))*exponents(j+1,3)*exponents(i+1,2)*exponents(h+1,1)
            end do
            pp_coeffs_2d(1 + i + j*(deg(2) + 1), 1) = sll_f_spline_pp_horner_1d(deg(1), coeffs, x(1), 1)
         end do
      end do

      !> Perform a 2d hornerschema in the second and third dimension
      res = sll_f_spline_pp_horner_2d(deg(2:3), pp_coeffs_2d, x(2:3), [1, 1], [1, 1])*n_cells(component)

   end function sll_f_spline_pp_horner_derivative_3d

   !> Initialize \a sll_t_spline_pp_1d object (Set \a poly_coeffs depending on \a  degree)
   subroutine sll_s_spline_pp_init_1d(self, degree, n_cells, boundary)
      type(sll_t_spline_pp_1d), intent(out) ::  self  !< arbitrary degree 1d spline
      sll_int32, intent(in) :: degree !<degree of spline
      sll_int32, intent(in) :: n_cells !< number of gridcells
      sll_int32, intent(in), optional :: boundary
      sll_int32 :: ierr

      if (present(boundary)) then
         self%boundary_conditions = boundary
      end if

      SLL_ASSERT(n_cells >= degree)
      self%degree = degree
      self%n_cells = n_cells
      SLL_ALLOCATE(self%poly_coeffs(degree + 1, degree + 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%poly_coeffs_fp(degree + 1, degree + 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%poly_coeffs_fpa(degree + 2, degree + 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%poly_coeffs_fd(degree, degree + 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%scratch_b(degree + 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%scratch_p(degree + 1), ierr)
      SLL_ASSERT(ierr == 0)

      SLL_ALLOCATE(self%poly_coeffs_boundary_left(degree + 1, degree + 1, degree - 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%poly_coeffs_boundary_right(degree + 1, degree + 1, degree - 1), ierr)
      SLL_ASSERT(ierr == 0)

      SLL_ALLOCATE(self%poly_coeffs_fd_boundary_left(degree, degree + 1, degree - 1), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%poly_coeffs_fd_boundary_right(degree, degree + 1, degree - 1), ierr)
      SLL_ASSERT(ierr == 0)

      select case (self%boundary_conditions)
      case (sll_p_boundary_periodic)
         self%n_coeffs = n_cells
      case (sll_p_boundary_clamped)
         self%n_coeffs = n_cells + degree
      case (sll_p_boundary_clamped_square)
         self%n_coeffs = n_cells + degree
      case (sll_p_boundary_clamped_cubic)
         self%n_coeffs = n_cells + degree
      case (sll_p_boundary_clamped_clampeddiri)
         self%n_coeffs = n_cells + degree - 1
      case (sll_p_boundary_clampeddiri)
         self%n_coeffs = n_cells + degree - 2
      case (sll_p_boundary_clampeddiri_clamped)
         self%n_coeffs = n_cells + degree - 1
      case default
         SLL_ERROR('sll_s_spline_pp_init_1d', 'Wrong boundary condition.')
      end select

      select case (self%degree)
         !poly_coefficients ordered in rows beginning with the splinepart for the last interval(p,p-1,..,1) and in each row with the coefficient for the highest polynomial degree (x^p, x^p-1,...,x^0)
         !for clamped splines in the first and last degree-1 intervals the splineparts are different
      case (0)
         self%poly_coeffs = reshape((/1._f64/), (/1, 1/))
         self%poly_coeffs_fp = reshape((/1._f64/), (/1, 1/))
         self%poly_coeffs_fpa = reshape((/1._f64, 0._f64/), (/2, 1/))
      case (1)
         self%poly_coeffs = reshape((/-1._f64, 1._f64, 1._f64, 0._f64/), (/2, 2/))
         self%poly_coeffs_fp = reshape((/-inv_2, 1._f64, inv_2, 0._f64/), (/2, 2/))
         self%poly_coeffs_fpa = reshape((/-inv_2, 1._f64, inv_2, inv_2, 0._f64, 0._f64/), (/3, 2/))
         self%poly_coeffs_fd = reshape((/-1._f64, 1._f64/), (/1, 2/))

      case (2)
         self%poly_coeffs = reshape((/inv_2, -1._f64, inv_2, &
                                      -1._f64, 1._f64, inv_2, &
                                      inv_2, 0._f64, 0._f64/), (/3, 3/))

         self%poly_coeffs_fp = reshape((/inv_6, -inv_2, inv_2, &
                                         -inv_3, inv_2, inv_2, &
                                         inv_6, 0._f64, 0._f64/), (/3, 3/))

         self%poly_coeffs_fpa = reshape((/inv_6, -inv_2, inv_2, 5.0_f64*inv_6, &
                                          -inv_3, inv_2, inv_2, inv_6, &
                                          inv_6, 0._f64, 0._f64, 0._f64/), (/4, 3/))

         self%poly_coeffs_fd = reshape((/1._f64, -1._f64, &
                                         -2._f64, 1._f64, &
                                         1._f64, 0._f64/), (/2, 3/))

         ! first interval
         self%poly_coeffs_boundary_left(:, 1, 1) = [1._f64, -2._f64, 1._f64]
         self%poly_coeffs_boundary_left(:, 2, 1) = [-1.5_f64, 2._f64, 0._f64]
         self%poly_coeffs_boundary_left(:, 3, 1) = [inv_2, 0._f64, 0._f64]

         self%poly_coeffs_fd_boundary_left(:, 1, 1) = [2._f64, -2._f64]
         self%poly_coeffs_fd_boundary_left(:, 2, 1) = [-3._f64, 2._f64]
         self%poly_coeffs_fd_boundary_left(:, 3, 1) = [1._f64, 0._f64]

         ! last interval
         self%poly_coeffs_boundary_right(:, 1, 1) = [0.5_f64, -1._f64, 0.5_f64]
         self%poly_coeffs_boundary_right(:, 2, 1) = [-1.5_f64, 1._f64, 0.5_f64]
         self%poly_coeffs_boundary_right(:, 3, 1) = [1.0_f64, 0._f64, 0._f64]

         self%poly_coeffs_fd_boundary_right(:, 1, 1) = [1._f64, -1._f64]
         self%poly_coeffs_fd_boundary_right(:, 2, 1) = [-3._f64, 1._f64]
         self%poly_coeffs_fd_boundary_right(:, 3, 1) = [2.0_f64, 0._f64]

         if (self%boundary_conditions == sll_p_boundary_clamped_square) then
            self%poly_coeffs_boundary_left(:, 1, 1) = 0._f64 !perfect conductor boundary condition
            self%poly_coeffs_boundary_right(:, 3, 1) = 0._f64!perfect conductor boundary condition
         end if

      case (3)
         self%poly_coeffs = reshape((/-inv_6, inv_2, -inv_2, inv_6, &
                                      inv_2, -1._f64, 0._f64, 4._f64*inv_6, &
                                      -inv_2, inv_2, inv_2, inv_6, &
                                      inv_6, 0._f64, 0._f64, 0._f64/), (/4, 4/))

         self%poly_coeffs_fp = reshape((/-inv_24, inv_6, -inv_4, inv_6, &
                                         inv_8, -inv_3, 0._f64, 4._f64*inv_6, &
                                         -inv_8, inv_6, inv_4, inv_6, &
                                         inv_24, 0._f64, 0._f64, 0._f64/), (/4, 4/))

         self%poly_coeffs_fpa = reshape((/-inv_24, inv_6, -inv_4, inv_6, 23._f64*inv_24, &
                                          inv_8, -inv_3, 0._f64, 4._f64*inv_6, inv_2, &
                                          -inv_8, inv_6, inv_4, inv_6, inv_24, &
                                          inv_24, 0._f64, 0._f64, 0._f64, 0._f64/), (/5, 4/))

         self%poly_coeffs_fd = reshape((/-inv_2, 1._f64, -inv_2, &
                                         3._f64*inv_2, -2._f64, 0._f64, &
                                         -3._f64*inv_2, 1._f64, inv_2, &
                                         inv_2, 0._f64, 0._f64/), (/3, 4/))

         ! first interval
         self%poly_coeffs_boundary_left(:, 1, 1) = [-1.0_f64, 3.0_f64, -3.0_f64, 1.0_f64]
         self%poly_coeffs_boundary_left(:, 2, 1) = [1.75_f64, -4.5_f64, 3.0_f64, 0._f64]
         self%poly_coeffs_boundary_left(:, 3, 1) = [-11._f64/12._f64, 1.5_f64, 0.0_f64, 0.0_f64]
         self%poly_coeffs_boundary_left(:, 4, 1) = [inv_6, 0._f64, 0._f64, 0._f64]

         self%poly_coeffs_fd_boundary_left(:, 1, 1) = [-3.0_f64, 6.0_f64, -3.0_f64]
         self%poly_coeffs_fd_boundary_left(:, 2, 1) = [5.25_f64, -9._f64, 3.0_f64]
         self%poly_coeffs_fd_boundary_left(:, 3, 1) = [-11._f64/4._f64, 3._f64, 0.0_f64]
         self%poly_coeffs_fd_boundary_left(:, 4, 1) = [inv_2, 0._f64, 0._f64]

         ! second interval
         self%poly_coeffs_boundary_left(:, 1, 2) = [-0.25_f64, 0.75_f64, -0.75_f64, 0.25_f64]
         self%poly_coeffs_boundary_left(:, 2, 2) = [7._f64/12._f64, -1.25_f64, 0.25_f64, 7._f64/12._f64]
         self%poly_coeffs_boundary_left(:, 3, 2) = [-0.5_f64, 0.5_f64, 0.5_f64, inv_6]
         self%poly_coeffs_boundary_left(:, 4, 2) = [inv_6, 0._f64, 0._f64, 0._f64]

         self%poly_coeffs_fd_boundary_left(:, 1, 2) = [-0.75_f64, 1.5_f64, -0.75_f64]
         self%poly_coeffs_fd_boundary_left(:, 2, 2) = [7._f64/4._f64, -2.5_f64, 0.25_f64]
         self%poly_coeffs_fd_boundary_left(:, 3, 2) = [-1.5_f64, 1._f64, 0.5_f64]
         self%poly_coeffs_fd_boundary_left(:, 4, 2) = [inv_2, 0._f64, 0._f64]

         ! second last interval
         self%poly_coeffs_boundary_right(:, 1, 1) = [-inv_6, inv_2, -inv_2, inv_6]
         self%poly_coeffs_boundary_right(:, 2, 1) = [inv_2, -1._f64, 0._f64, 4._f64*inv_6]
         self%poly_coeffs_boundary_right(:, 3, 1) = [-7.0_f64/12.0_f64, 0.5_f64, 0.5_f64, inv_6]
         self%poly_coeffs_boundary_right(:, 4, 1) = [0.25_f64, 0.0_f64, 0.0_f64, 0.0_f64]

         self%poly_coeffs_fd_boundary_right(:, 1, 1) = [-inv_2, 1._f64, -inv_2]
         self%poly_coeffs_fd_boundary_right(:, 2, 1) = [3._f64*inv_2, -2._f64, 0._f64]
         self%poly_coeffs_fd_boundary_right(:, 3, 1) = [-7.0_f64/4.0_f64, 1._f64, 0.5_f64]
         self%poly_coeffs_fd_boundary_right(:, 4, 1) = [0.75_f64, 0.0_f64, 0.0_f64]

         ! last interval
         self%poly_coeffs_boundary_right(:, 1, 2) = [-inv_6, inv_2, -inv_2, inv_6]
         self%poly_coeffs_boundary_right(:, 2, 2) = [11.0_f64/12.0_f64, -1.25_f64, -0.25_f64, 7.0_f64/12.0_f64]
         self%poly_coeffs_boundary_right(:, 3, 2) = [-1.75_f64, 0.75_f64, 0.75_f64, 0.25_f64]
         self%poly_coeffs_boundary_right(:, 4, 2) = [1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64]

         self%poly_coeffs_fd_boundary_right(:, 1, 2) = [-inv_2, 1._f64, -inv_2]
         self%poly_coeffs_fd_boundary_right(:, 2, 2) = [11.0_f64/4.0_f64, -2.5_f64, -0.25_f64]
         self%poly_coeffs_fd_boundary_right(:, 3, 2) = [-5.25_f64, 1.5_f64, 0.75_f64]
         self%poly_coeffs_fd_boundary_right(:, 4, 2) = [3.0_f64, 0.0_f64, 0.0_f64]

         if (self%boundary_conditions == sll_p_boundary_clamped_cubic) then
            self%poly_coeffs_boundary_left(:, 1, 1) = 0._f64 !perfect conductor boundary condition
            self%poly_coeffs_boundary_right(:, 4, 2) = 0._f64 !perfect conductor boundary condition
         end if
      case (4)
         self%poly_coeffs = reshape((/inv_24, -inv_6, inv_4, -inv_6, inv_24 &
                                      , -inv_6, inv_2, -inv_4, -inv_2, 11._f64*inv_24 &
                                      , inv_4, -inv_2, -inv_4, inv_2, 11._f64*inv_24 &
                                      , -inv_6, inv_6, inv_4, inv_6, inv_24 &
                                      , inv_24, 0._f64, 0._f64, 0._f64, 0._f64/), (/5, 5/))

         self%poly_coeffs_fp = reshape((/inv_120, -inv_24, inv_12, -inv_12, inv_24 &
                                         , -inv_30, inv_8, -inv_12, -inv_4, 11._f64*inv_24 &
                                         , inv_20, -inv_8, -inv_12, inv_4, 11._f64*inv_24 &
                                         , -inv_30, inv_24, inv_12, inv_12, inv_24 &
                                         , inv_120, 0._f64, 0._f64, 0._f64, 0._f64/), (/5, 5/))

         self%poly_coeffs_fpa = reshape((/inv_120, -inv_24, inv_12, -inv_12, inv_24, 119._f64*inv_120 &
                                          , -inv_30, inv_8, -inv_12, -inv_4, 11._f64*inv_24, 93._f64*inv_120 &
                                          , inv_20, -inv_8, -inv_12, inv_4, 11._f64*inv_24, 27._f64*inv_120 &
                                          , -inv_30, inv_24, inv_12, inv_12, inv_24, inv_120 &
                                          , inv_120, 0._f64, 0._f64, 0._f64, 0._f64, 0.0_f64/), (/6, 5/))

      case (5)
         self%poly_coeffs = reshape((/-inv_120, inv_24, -inv_12, inv_12, -inv_24, inv_120 &
                                      , inv_24, -inv_6, inv_6, inv_6, -5._f64*inv_12, 26._f64*inv_120 &
                                      , -inv_12, inv_4, 0._f64, -inv_2, 0._f64, 11._f64*inv_20 &
                                      , inv_12, -inv_6, -inv_6, inv_6, 5._f64*inv_12, 26._f64*inv_120 &
                                      , -inv_24, inv_24, inv_12, inv_12, inv_24, inv_120 &
                                      , inv_120, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64/), (/6, 6/))

         self%poly_coeffs_fp = reshape((/-inv_720, inv_120, -inv_48, inv_36, -inv_48, inv_120 &
                                         , inv_144, -inv_30, inv_24, inv_18, -5._f64*inv_24, 26._f64*inv_120 &
                                         , -inv_72, inv_20, 0._f64, -inv_6, 0._f64, 11._f64*inv_20 &
                                         , inv_72, -inv_30, -inv_24, inv_18, 5._f64*inv_24, 26._f64*inv_120 &
                                         , -inv_144, inv_120, inv_48, inv_36, inv_48, inv_120 &
                                         , inv_720, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64/), (/6, 6/))

         self%poly_coeffs_fpa = reshape((/-inv_720, inv_120, -inv_48, inv_36, -inv_48, inv_120, 719._f64*inv_720 &
                                          , inv_144, -inv_30, inv_24, inv_18, -5._f64*inv_24, 26._f64*inv_120, 662._f64*inv_720 &
                                          , -inv_72, inv_20, 0._f64, -inv_6, 0._f64, 11._f64*inv_20, inv_2 &
                                          , inv_72, -inv_30, -inv_24, inv_18, 5._f64*inv_24, 26._f64*inv_120, 58._f64*inv_720 &
                                          , -inv_144, inv_120, inv_48, inv_36, inv_48, inv_120, inv_720 &
                                          , inv_720, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64/), (/7, 6/))
      case (6)
         self%poly_coeffs = reshape((/inv_720, -inv_120, inv_48, -inv_36, inv_48, -inv_120, inv_720 &
                                      , -inv_120, inv_24, -inv_16, -inv_36, 135._f64*inv_720, -150._f64*inv_720, 57._f64*inv_720 &
                                      , inv_48, -inv_12, inv_24, 160._f64*inv_720, -150._f64*inv_720, -inv_3, 302._f64*inv_720 &
                                      , -inv_36, inv_12, inv_24, -160._f64*inv_720, -150._f64*inv_720, inv_3, 302._f64*inv_720 &
                                      , inv_48, -inv_24, -inv_16, inv_36, 135._f64*inv_720, 150._f64*inv_720, 57._f64*inv_720 &
                                      , -inv_120, inv_120, inv_48, inv_36, inv_48, inv_120, inv_720 &
                                      , inv_720, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64/), (/7, 7/))

         self%poly_coeffs_fp = reshape((/inv_5040, -inv_720, inv_240, -inv_144, inv_144, -inv_240, inv_720 &
                                         , -6._f64*inv_5040, inv_144, -inv_80, -inv_144, inv_16, -75._f64*inv_720, 57._f64*inv_720 &
                                         , 15._f64*inv_5040, -inv_72, inv_120, inv_18, -50._f64*inv_720, -inv_6, 302._f64*inv_720 &
                                         , -20._f64*inv_5040, inv_72, inv_120, -inv_18, -50._f64*inv_720, inv_6, 302._f64*inv_720 &
                                         , 15._f64*inv_5040, -inv_144, -inv_80, inv_144, inv_16, 75._f64*inv_720, 57._f64*inv_720 &
                                         , -6.0_f64*inv_5040, inv_720, inv_240, inv_144, inv_144, inv_240, inv_720 &
                                         , inv_5040, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64, 0.0_f64/), (/7, 7/))

         self%poly_coeffs_fpa = reshape((/inv_5040, -inv_720, inv_240, -inv_144, inv_144, -inv_240, inv_720, 5039._f64*inv_5040 &
                     , -6._f64*inv_5040, inv_144, -inv_80, -inv_144, inv_16, -75._f64*inv_720, 57._f64*inv_720, 4919._f64*inv_5040 &
                      , 15._f64*inv_5040, -inv_72, inv_120, inv_18, -50._f64*inv_720, -inv_6, 302._f64*inv_720, 3728._f64*inv_5040 &
                      , -20._f64*inv_5040, inv_72, inv_120, -inv_18, -50._f64*inv_720, inv_6, 302._f64*inv_720, 1312._f64*inv_5040 &
                       , 15._f64*inv_5040, -inv_144, -inv_80, inv_144, inv_16, 75._f64*inv_720, 57._f64*inv_720, 121._f64*inv_5040 &
                                        , -6.0_f64*inv_5040, inv_720, inv_240, inv_144, inv_144, inv_240, inv_720, 1._f64*inv_5040 &
                                          , inv_5040, 0._f64, 0._f64, 0._f64, 0._f64, 0._f64, 0.0_f64, 0._f64/), (/8, 7/))
      case default
         SLL_ERROR('sll_s_spline_pp_init', 'Not implemented')
      end select

      if (self%boundary_conditions > 0 .and. degree > 3) then
         SLL_ERROR('sll_s_spline_pp_init_1d', 'Specified boundary conditions not implemented for given degree.')
      end if

   end subroutine sll_s_spline_pp_init_1d

   !> Initialize \a sll_t_spline_pp_2d object
   subroutine sll_s_spline_pp_init_2d(self, degree, n_cells, boundary)
      type(sll_t_spline_pp_2d) self !< arbitrary degree 2d spline
      sll_int32, intent(in) :: degree(2) !< degrees of the 1d splines
      sll_int32, intent(in) :: n_cells(2) !< number of gridcells in every dimension
      sll_int32, intent(in), optional :: boundary(2) !< boundary conditions

      if (present(boundary)) then
         call sll_s_spline_pp_init_1d(self%spline1, degree(1), n_cells(1), boundary(1))
         call sll_s_spline_pp_init_1d(self%spline2, degree(2), n_cells(2), boundary(2))
      else
         call sll_s_spline_pp_init_1d(self%spline1, degree(1), n_cells(1))
         call sll_s_spline_pp_init_1d(self%spline2, degree(2), n_cells(2))
      end if

   end subroutine sll_s_spline_pp_init_2d

   !> Initialize \a sll_t_spline_pp_3d object
   subroutine sll_s_spline_pp_init_3d(self, degree, n_cells, boundary)
      type(sll_t_spline_pp_3d) self !< arbitrary degree 3d spline
      sll_int32, intent(in) :: degree(3) !< degrees of the 1d splines
      sll_int32, intent(in) :: n_cells(3) !< number of gridcells in every dimension
      sll_int32, intent(in), optional :: boundary(3) !< boundary conditions
      sll_int32 :: ierr

      if (present(boundary)) then
         call sll_s_spline_pp_init_1d(self%spline1, degree(1), n_cells(1), boundary(1))
         call sll_s_spline_pp_init_1d(self%spline2, degree(2), n_cells(2), boundary(2))
         call sll_s_spline_pp_init_1d(self%spline3, degree(3), n_cells(3), boundary(3))
      else
         call sll_s_spline_pp_init_1d(self%spline1, degree(1), n_cells(1))
         call sll_s_spline_pp_init_1d(self%spline2, degree(2), n_cells(2))
         call sll_s_spline_pp_init_1d(self%spline3, degree(3), n_cells(3))
      end if
      SLL_ALLOCATE(self%scratch_b((degree(2) + 1)*(degree(3) + 1)), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%scratch_p((degree(2) + 1)*(degree(3) + 1)), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%scratch_pp((degree(1) + 1)*(degree(2) + 1)*(degree(3) + 1)), ierr)
      SLL_ASSERT(ierr == 0)
      SLL_ALLOCATE(self%scratch_coeffs((degree(1) + 1)*(degree(2) + 1)*(degree(3) + 1), n_cells(1)*n_cells(2)*n_cells(3)), ierr)
      SLL_ASSERT(ierr == 0)

   end subroutine sll_s_spline_pp_init_3d

   !> Destructor 1d
   subroutine sll_s_spline_pp_free_1d(self)
      type(sll_t_spline_pp_1d), intent(inout) ::  self !< arbitrary degree 1d spline

      deallocate (self%poly_coeffs)
      deallocate (self%poly_coeffs_fp)
      deallocate (self%scratch_b)
      deallocate (self%scratch_p)

   end subroutine sll_s_spline_pp_free_1d

   !> Destructor 2d
   subroutine sll_s_spline_pp_free_2d(self)
      type(sll_t_spline_pp_2d), intent(inout) ::  self !< arbitrary degree 2d spline

      call sll_s_spline_pp_free_1d(self%spline1)
      call sll_s_spline_pp_free_1d(self%spline2)

   end subroutine sll_s_spline_pp_free_2d

   !> Destructor 3d
   subroutine sll_s_spline_pp_free_3d(self)
      type(sll_t_spline_pp_3d), intent(inout) ::  self !< arbitrary degree 3d spline

      call sll_s_spline_pp_free_1d(self%spline1)
      call sll_s_spline_pp_free_1d(self%spline2)
      call sll_s_spline_pp_free_1d(self%spline3)
      deallocate (self%scratch_b)
      deallocate (self%scratch_p)
      deallocate (self%scratch_pp)

   end subroutine sll_s_spline_pp_free_3d

end module sll_m_splines_pp
