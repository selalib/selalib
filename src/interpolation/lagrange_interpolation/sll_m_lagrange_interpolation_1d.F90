module sll_m_lagrange_interpolation_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_periodic

   implicit none

   public :: &
      sll_s_compute_lagrange_interpolation_1d, &
      sll_s_cubic_spline_1d_eval_array, &
      sll_f_new_lagrange_interpolation_1d, &
      sll_t_lagrange_interpolation_1d

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type :: sll_t_lagrange_interpolation_1d
      sll_int32                            :: d !half of stencil
      sll_int32                            :: nb_cell
      sll_int32                            :: bc_type
      sll_int32                            :: index_gap
      sll_real64                           :: alpha
      sll_real64                           :: xmin
      sll_real64                           :: xmax
      sll_real64                           :: deta !< \a deta is the grid spacing
      sll_real64, dimension(:), pointer    :: wj
      sll_real64, dimension(:), pointer    :: wj_scale
      sll_real64, dimension(:), pointer    :: data_out !result=p(x) where p is the polynomial of interpolation
      sll_int32                            :: periodic_last !< \a periodic_last indicates if the input data repeats the first point at the end if we have periodic data. It takes the values 0 (not repeated) or 1 (repeated).  Default : 1.
   end type sll_t_lagrange_interpolation_1d

! integer, parameter :: sll_p_periodic = 0, sll_p_hermite = 3

   interface delete
      module procedure delete_lagrange_interpolation_1D
   end interface

contains  !*****************************************************************************

   function sll_f_new_lagrange_interpolation_1d(num_points, xmin, xmax, bc_type, d, periodic_last)
      type(sll_t_lagrange_interpolation_1d), pointer :: sll_f_new_lagrange_interpolation_1d
      sll_int32 ::ierr
      sll_int32, intent(in) :: d, num_points, bc_type
      sll_int32, intent(in), optional :: periodic_last !< \a periodic_last indicates if the input data repeats the first point at the end if we have periodic data. It takes the values 0 (not repeated) or 1 (repeated).  Default : 1.
      sll_real64 :: xmin, xmax

      SLL_ALLOCATE(sll_f_new_lagrange_interpolation_1d, ierr)

      SLL_ALLOCATE(sll_f_new_lagrange_interpolation_1d%wj(2*d), ierr)
      SLL_ALLOCATE(sll_f_new_lagrange_interpolation_1d%wj_scale(2*d), ierr)
      SLL_ALLOCATE(sll_f_new_lagrange_interpolation_1d%data_out(num_points), ierr)
      sll_f_new_lagrange_interpolation_1d%d = d
      sll_f_new_lagrange_interpolation_1d%xmin = xmin
      sll_f_new_lagrange_interpolation_1d%xmax = xmax
      sll_f_new_lagrange_interpolation_1d%nb_cell = num_points - 1
      sll_f_new_lagrange_interpolation_1d%bc_type = bc_type
      sll_f_new_lagrange_interpolation_1d%deta = (xmax - xmin)/real(sll_f_new_lagrange_interpolation_1d%nb_cell, f64)
      if (present(periodic_last)) then
         sll_f_new_lagrange_interpolation_1d%periodic_last = periodic_last
      else
         sll_f_new_lagrange_interpolation_1d%periodic_last = 1
      end if

   end function sll_f_new_lagrange_interpolation_1d

!> This function computes the weights w_j for the barycentric formula
   subroutine sll_s_compute_lagrange_interpolation_1d(lagrange)
      type(sll_t_lagrange_interpolation_1d), pointer :: lagrange !< \a lagrange is the lagrange interpolator object
      sll_int32 :: i, j
      sll_int32, dimension(1:4*lagrange%d - 2) :: table
      sll_real64, dimension(1:2*lagrange%d) :: wj

      do i = 1, 2*lagrange%d - 1
         table(i) = 2*lagrange%d - 1 - (i - 1)
         table(i + 2*lagrange%d - 1) = i
      end do

      wj = 1.0_f64
      do i = 1, lagrange%d
         do j = 1, 2*lagrange%d - 1
            wj(i) = wj(i)*real(table(i + j - 1), f64)
         end do
         wj(i) = ((-1.0_f64)**(lagrange%d + i))*wj(i)
      end do
      do i = 1, lagrange%d
         wj(i + lagrange%d) = -wj(lagrange%d - i + 1)
      end do
      wj = 1.0_f64/wj

      lagrange%wj = wj

   end subroutine sll_s_compute_lagrange_interpolation_1d

!> This function computes the value at the grid points displacement by -alpha
   subroutine sll_s_cubic_spline_1d_eval_array(fi, alpha, lagrange)
      type(sll_t_lagrange_interpolation_1d), pointer :: lagrange !< \a lagrange is the lagrange interpolator object
      sll_real64, intent(in) :: alpha !< \a alpha is the (negative) displacement
      sll_int32 ::i, j, index_gap
      sll_real64 :: sum1, sum2, beta, h
      sll_real64, dimension(1:lagrange%nb_cell + 1), intent(in) :: fi !< \a fi are the values at the grid points

      lagrange%alpha = -alpha

      h = lagrange%deta! grid size
! How many cells do we displace? alpha = index_gap*h+beta*h
      index_gap = floor(lagrange%alpha/h)
      beta = lagrange%alpha/h - real(index_gap, f64)

! Take care of the case where alpha/h is negative and below machine precision
      if (beta == 1.0_f64) then
         beta = 0.0_f64
         index_gap = index_gap + 1
      end if

! If the displacement is precisely a multiple of h, we need to avoid division by zero
      if (beta == 0.0_f64) then
         select case (lagrange%bc_type)
         case (sll_p_periodic)
            do j = 1, lagrange%nb_cell + lagrange%periodic_last
               lagrange%data_out(j) = fi(modulo(index_gap + j - 1, lagrange%nb_cell) + 1); 
            end do
         case (sll_p_hermite)
            do j = 1, lagrange%nb_cell + 1
               lagrange%data_out(j) = fi(min(max(1, index_gap + j), lagrange%nb_cell + 1)); 
            end do
         end select

      else

         sum2 = 0.0_f64
         do j = 1, 2*lagrange%d
            lagrange%wj_scale(j) = lagrange%wj(j)/(beta + real(lagrange%d - j, f64))
            sum2 = sum2 + lagrange%wj_scale(j)
         end do

         select case (lagrange%bc_type)
         case (sll_p_periodic)

            do i = 1, lagrange%nb_cell + lagrange%periodic_last
               sum1 = 0.0_f64
               do j = 1, lagrange%d*2
               sum1 = sum1 + lagrange%wj_scale(j)*fi(modulo(index_gap + (i - 1) + (j - 1) - (lagrange%d - 1), lagrange%nb_cell) + 1)
               end do
               lagrange%data_out(i) = sum1/sum2
            end do

         case (sll_p_hermite)

            do i = 1, lagrange%nb_cell + 1
               sum1 = 0.0_f64
               do j = 1, lagrange%d*2
                  if (index_gap + (i - 1) + (j - 1) - (lagrange%d - 1) < 0) then
                     sum1 = sum1 + lagrange%wj_scale(j)*fi(1)
                  else if (index_gap + (i - 1) + (j - 1) - (lagrange%d - 1) > lagrange%nb_cell) then
                     sum1 = sum1 + lagrange%wj_scale(j)*fi(lagrange%nb_cell + 1)
                  else
                     sum1 = sum1 + lagrange%wj_scale(j)*fi(index_gap + (i - 1) + (j - 1) - (lagrange%d - 1) + 1)
                  end if
               end do
               lagrange%data_out(i) = sum1/sum2
            end do

         case default
            print *, 'ERROR: sll_s_compute_lagrange_interpolation_1d(): not recognized boundary condition'
            STOP
         end select
      end if

   end subroutine sll_s_cubic_spline_1d_eval_array

   subroutine delete_lagrange_interpolation_1D(sll_m_lagrange_interpolation)
      type(sll_t_lagrange_interpolation_1d), pointer :: sll_m_lagrange_interpolation
      sll_int32                    :: ierr
      SLL_ASSERT(associated(sll_m_lagrange_interpolation))
      SLL_DEALLOCATE(sll_m_lagrange_interpolation%wj, ierr)
      SLL_DEALLOCATE(sll_m_lagrange_interpolation, ierr)
   end subroutine delete_lagrange_interpolation_1D

end module sll_m_lagrange_interpolation_1d

