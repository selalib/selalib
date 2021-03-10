!> @ingroup interpolators
!> @brief
!! Interpolator 1d using dg interpolation
!> @author
!! Katharina Kormann, IPP
!! @details
!! Algorithm described in Crouseilles, Mehrenberger, Vecil, ESAIM:Proceedings, 32, pp. 211-230, 2011.
!!
module sll_m_dg_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_gauss_legendre_integration, only: &
      sll_f_gauss_legendre_points, &
      sll_f_gauss_legendre_weights

   use sll_m_gauss_lobatto_integration, only: &
      sll_f_gauss_lobatto_points, &
      sll_f_gauss_lobatto_weights

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   implicit none

   public :: &
      sll_t_dg_interpolator_1d, &
      sll_s_dg_interpolator_1d_init, &
      sll_s_dg_interpolator_1d_free, &
      sll_s_dg_interpolator_1d_interpolate_array_disp_periodic, &
      sll_s_dg_interpolator_1d_interpolate_array_disp_halo, &
      sll_s_dg_interpolator_1d_dg_to_equi_cell, &
      sll_s_dg_interpolator_1d_dg_to_equi_cell_staggered, &
      sll_p_dg_gauss_legendre, &
      sll_p_dg_gauss_lobatto, &
      sll_p_dg_uniform

   private
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32, parameter :: sll_p_dg_gauss_legendre = 101
   sll_int32, parameter :: sll_p_dg_gauss_lobatto = 102
   sll_int32, parameter :: sll_p_dg_uniform = 103

   !> 1d interpolator using dg
   type :: sll_t_dg_interpolator_1d

      sll_int32 :: n_cells !< No. of cells on dg mesh
      sll_int32 :: degree  !< Degree of the dg space
      sll_int32 :: n_total !< Total number of points
      sll_real64 :: cell_size !< Size of each cell

      sll_real64, pointer :: positions(:)
      sll_real64, pointer :: quadrature_points(:)
      sll_real64, pointer :: quadrature_inv_diffs(:, :)
      sll_real64, pointer :: quadrature_weights(:)
      sll_real64, pointer :: quadrature_inv_weights(:)
      sll_real64, pointer :: weights1(:, :)
      sll_real64, pointer :: weights2(:, :)

   end type sll_t_dg_interpolator_1d

contains

   subroutine sll_s_dg_interpolator_1d_init(self, n_cells, degree, x_min, x_max, type)
      class(sll_t_dg_interpolator_1d), intent(inout) :: self
      sll_int32, intent(in) :: n_cells !< No. of cells on dg mesh
      sll_int32, intent(in) :: degree  !< Degree of the dg space
      sll_real64, intent(in) :: x_min  !< Lower bound of the domain
      sll_real64, intent(in) :: x_max  !< Upper bound of the domain
      sll_int32, intent(in) :: type    !< Descriptor telling if Gauss-Lobatto, Gauss-Legendre or uniform points shall be used as points within the cells

      sll_int32 :: ierr, j, i

      self%n_cells = n_cells
      self%degree = degree
      self%n_total = n_cells*degree

      SLL_ALLOCATE(self%positions(self%n_total), ierr)
      SLL_ALLOCATE(self%quadrature_points(self%degree), ierr)
      SLL_ALLOCATE(self%quadrature_weights(self%degree), ierr)
      SLL_ALLOCATE(self%quadrature_inv_weights(self%degree), ierr)
      SLL_ALLOCATE(self%weights1(self%degree, self%degree), ierr)
      SLL_ALLOCATE(self%weights2(self%degree, self%degree), ierr)
      SLL_ALLOCATE(self%quadrature_inv_diffs(self%degree, self%degree), ierr)

      if (type == sll_p_dg_gauss_legendre) then
         ! Gauss-Legendre points and weights for each cell
         self%quadrature_points = sll_f_gauss_legendre_points(self%degree, 0.0_f64, 1.0_f64)
         self%quadrature_weights = sll_f_gauss_legendre_weights(self%degree, 0.0_f64, 1.0_f64)
      elseif (type == sll_p_dg_gauss_lobatto) then
         ! Gauss-Lobatto points and weights for each cell
         self%quadrature_points = sll_f_gauss_lobatto_points(self%degree, 0.0_f64, 1.0_f64)
         self%quadrature_weights = sll_f_gauss_lobatto_weights(self%degree, 0.0_f64, 1.0_f64)
      elseif (type == sll_p_dg_uniform) then
         ! Uniform points and weights for each cell
         do j = 1, self%degree
            self%quadrature_points(j) = real(j - 1, f64)/real(self%degree - 1, f64)
         end do
         self%quadrature_weights = 1.0_f64/(self%degree - 1)
         self%quadrature_weights(1) = 0.5_f64/(self%degree - 1)
         self%quadrature_weights(self%degree) = self%quadrature_weights(1)
      else
         SLL_ERROR('sll_s_dg_interpolator_1d_init', 'Wrong type of points.')
      end if

      ! Compute cell size
      self%cell_size = (x_max - x_min)/n_cells

      do j = 1, n_cells
         self%positions(1 + (j - 1)*degree:j*degree) = x_min + &
                                                       (real(j - 1, f64) + self%quadrature_points)*self%cell_size
      end do

      ! cache the inverse differences of the points which are
      ! used repeatedly during lagrange polynomial evaluation
      do i = 1, self%degree
         do j = 1, self%degree
            if (i /= j) then
               self%quadrature_inv_diffs(j, i) = 1.0_f64/(self%quadrature_points(i) - self%quadrature_points(j))
            else
               ! this value should never be accessed
               self%quadrature_inv_diffs(j, i) = 1.0_f64
            end if
         end do
      end do

      ! cache the inverse weights which are
      ! used repeatedly during interpolation
      do i = 1, self%degree
         self%quadrature_inv_weights(i) = 1.0_f64/self%quadrature_weights(i)
      end do
   end subroutine sll_s_dg_interpolator_1d_init

   subroutine sll_s_dg_interpolator_1d_free(self)
      class(sll_t_dg_interpolator_1d), intent(inout) :: self

      sll_int32 :: ierr

      SLL_DEALLOCATE(self%positions, ierr)
      SLL_DEALLOCATE(self%quadrature_points, ierr)
      SLL_DEALLOCATE(self%quadrature_inv_diffs, ierr)
      SLL_DEALLOCATE(self%quadrature_weights, ierr)
      SLL_DEALLOCATE(self%quadrature_inv_weights, ierr)
      SLL_DEALLOCATE(self%weights1, ierr)
      SLL_DEALLOCATE(self%weights2, ierr)
   end subroutine sll_s_dg_interpolator_1d_free

!DIR$ ATTRIBUTES INLINE :: sll_s_dg_interpolator_1d_interpolate_array_disp_periodic
   subroutine sll_s_dg_interpolator_1d_interpolate_array_disp_periodic(self, num_pts, data, alpha, output_array)
      class(sll_t_dg_interpolator_1d), intent(inout) :: self
      sll_int32, intent(in) :: num_pts
      sll_real64, intent(in) :: alpha
      sll_real64, dimension(:), intent(in) :: data
      sll_real64, dimension(num_pts), intent(out) :: output_array

      sll_int32 :: i, j, k, r, q, ind1, ind2, degree, n_cells
      sll_real64 :: beta, alph
      sll_real64, pointer :: quadrature_points(:)
      sll_real64, pointer :: quadrature_weights(:)
      sll_real64, pointer :: quadrature_inv_weights(:)
      sll_real64, pointer :: weights1(:, :)
      sll_real64, pointer :: weights2(:, :)
      sll_real64, pointer :: quadrature_inv_diffs(:, :)

      SLL_ASSERT(num_pts == self%n_total)

      ! Normalize displacement
      beta = alpha/self%cell_size
      q = floor(beta)
      alph = beta - real(q, f64)

      ! cache nested values to avoid subsequent dereferencing
      degree = self%degree
      n_cells = self%n_cells
      quadrature_points => self%quadrature_points
      quadrature_inv_diffs => self%quadrature_inv_diffs
      quadrature_weights => self%quadrature_weights
      quadrature_inv_weights => self%quadrature_inv_weights
      weights1 => self%weights1
      weights2 => self%weights2

      ! Precompute the weights

!    do j=1, degree
!       do k =1, degree
!          weights1(j,k) = 0.0_f64
!          weights2(j,k) = 0.0_f64
!          do r=1,degree
!!             self%weights1(j,k) = self%weights1(j,k) + self%quadrature_weights(r)*&
!!                  lagrange_poly( self, k,  alph + (1.0_f64-alph)*self%quadrature_points(r)) *&
!!                  lagrange_poly( self, j, (1.0_f64-alph)*self%quadrature_points(r))
!!             self%weights2(j,k) = self%weights2(j,k) + self%quadrature_weights(r)*&
!!                  lagrange_poly( self, k,  alph*self%quadrature_points(r)) *&
!!                  lagrange_poly( self, j, alph*(self%quadrature_points(r)-1.0_f64)+1.0_f64 )
!             weights1(j,k) = weights1(j,k) + quadrature_weights(r)*&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, k,  alph + (1.0_f64-alph)*quadrature_points(r)) *&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, j,         (1.0_f64-alph)*quadrature_points(r))
!             weights2(j,k) = weights2(j,k) + quadrature_weights(r)*&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, k, alph *quadrature_points(r)) *&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, j, alph*(quadrature_points(r)-1.0_f64)+1.0_f64 )
!          end do
!          weights1(j,k) = weights1(j,k) *  (1.0_f64-alph)
!          weights2(j,k) = weights2(j,k) *  alph
!       end do
!    end do

      weights1(:, :) = 0.0_f64
      do j = 1, degree
         do k = 1, degree
            do r = 1, degree
               weights1(k, j) = weights1(k, j) + quadrature_weights(r)* &
              lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, k, alph + (1.0_f64 - alph)*quadrature_points(r))* &
                        lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, j, (1.0_f64 - alph)*quadrature_points(r))
            end do
            weights1(k, j) = weights1(k, j)*(1.0_f64 - alph)
         end do
      end do

      weights2(:, :) = 0.0_f64
      do j = 1, degree
         do k = 1, degree
            do r = 1, degree
               weights2(k, j) = weights2(k, j) + quadrature_weights(r)* &
                                lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, k, alph*quadrature_points(r))* &
              lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, j, alph*(quadrature_points(r) - 1.0_f64) + 1.0_f64)
            end do
            weights2(k, j) = weights2(k, j)*alph
         end do
      end do

      ! Interpolate
      output_array(:) = 0.0_f64
      do i = 1, n_cells
         do j = 1, degree
!          output_array((i-1)*degree+j) = 0.0_f64
            ind1 = modulo(i - 1 + q, self%n_cells)*degree
            ind2 = modulo(i + q, self%n_cells)*degree
            do k = 1, degree
               output_array((i - 1)*degree + j) = output_array((i - 1)*degree + j) + &
                                                  data(ind1 + k)*weights1(k, j) + data(ind2 + k)*weights2(k, j)
            end do
!          output_array((i-1)*degree+j) = output_array((i-1)*degree+j)/quadrature_weights(j)
            output_array((i - 1)*degree + j) = output_array((i - 1)*degree + j)*quadrature_inv_weights(j)
         end do
      end do
   end subroutine sll_s_dg_interpolator_1d_interpolate_array_disp_periodic

   !> Interpolation routine with halo cells
!DIR$ ATTRIBUTES INLINE :: sll_s_dg_interpolator_1d_interpolate_array_disp_halo
   subroutine sll_s_dg_interpolator_1d_interpolate_array_disp_halo(self, data, alpha, weights1, weights2, output_array)
      class(sll_t_dg_interpolator_1d), intent(inout) :: self
      sll_real64, intent(in) :: alpha
      sll_real64, intent(in) :: data(:)
      sll_real64, intent(inout) :: weights1(self%degree, self%degree)
      sll_real64, intent(inout) :: weights2(self%degree, self%degree)
      sll_real64, intent(inout) :: output_array(self%n_total)

      sll_int32 :: i, j, k, r, q, ind1, ind2, num_pts, degree, n_cells
      sll_real64 :: alph
      sll_real64, pointer :: quadrature_points(:)
      sll_real64, pointer :: quadrature_weights(:)
      sll_real64, pointer :: quadrature_inv_weights(:)
      sll_real64, pointer :: quadrature_inv_diffs(:, :)

      ! cache nested values to avoid subsequent dereferencing
      n_cells = self%n_cells
      degree = self%degree
      quadrature_points => self%quadrature_points
      quadrature_weights => self%quadrature_weights
      quadrature_inv_weights => self%quadrature_inv_weights
      quadrature_inv_diffs => self%quadrature_inv_diffs

      ! Normalize displacement
      !beta = alpha/self%cell_size
      !q    = floor( beta )
      !alph = beta - real(q, f64 )
      alph = alpha

      ! --- Precompute the weights

      ! optimizations: split loops, changed storage order of weights (j,k) -> (k,j)

!    weights1(:,:) = 0.0_f64
!    weights2(:,:) = 0.0_f64
!    do j=1, degree
!       do k =1, degree
!!          weights1(j,k) = 0.0_f64
!!          weights2(j,k) = 0.0_f64
!          do r=1,degree
!!             weights1(j,k) = weights1(j,k) + self%quadrature_weights(r)*&
!!                  lagrange_poly( self, k,  alph + (1.0_f64-alph)*self%quadrature_points(r)) *&
!!                  lagrange_poly( self, j, (1.0_f64-alph)*self%quadrature_points(r))
!!             weights2(j,k) = weights2(j,k) + self%quadrature_weights(r)*&
!!                  lagrange_poly( self, k,  alph*self%quadrature_points(r)) *&
!!                  lagrange_poly( self, j, alph*(self%quadrature_points(r)-1.0_f64)+1.0_f64 )
!
!             weights1(k,j) = weights1(k,j) + quadrature_weights(r)*&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, k, alph + (1.0_f64-alph)*quadrature_points(r)) *&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, j,        (1.0_f64-alph)*quadrature_points(r))
!             weights2(k,j) = weights2(k,j) + quadrature_weights(r)*&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, k, alph* quadrature_points(r)) *&
!                  lagrange_poly_opt( quadrature_points, quadrature_inv_diffs, degree, j, alph*(quadrature_points(r)-1.0_f64)+1.0_f64 )
!          end do
!          weights1(k,j) = weights1(k,j) * (1.0_f64-alph)
!          weights2(k,j) = weights2(k,j) * alph
!       end do
!    end do

      weights1(:, :) = 0.0_f64
      do j = 1, degree
         do k = 1, degree
            do r = 1, degree
               weights1(k, j) = weights1(k, j) + quadrature_weights(r)* &
              lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, k, alph + (1.0_f64 - alph)*quadrature_points(r))* &
                        lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, j, (1.0_f64 - alph)*quadrature_points(r))
            end do
            weights1(k, j) = weights1(k, j)*(1.0_f64 - alph)
         end do
      end do

      weights2(:, :) = 0.0_f64
      do j = 1, degree
         do k = 1, degree
            do r = 1, degree
               weights2(k, j) = weights2(k, j) + quadrature_weights(r)* &
                                lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, k, alph*quadrature_points(r))* &
              lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, j, alph*(quadrature_points(r) - 1.0_f64) + 1.0_f64)
            end do
            weights2(k, j) = weights2(k, j)*alph
         end do
      end do

      ! Interpolate
      output_array(:) = 0.0_f64
      do i = 1, n_cells
         do j = 1, degree
            ind1 = (i - 1)*degree
            ind2 = i*degree
            !ind1 = modulo( i-1+q, n_cells )* degree
            !ind2 = modulo( i+q, n_cells )*degree
            !ind1 = i-1
            do k = 1, degree
!             output_array((i-1)*degree+j) = output_array((i-1)*degree+j) + &
!                  data(ind1+k)*weights1(j,k) + data(ind2+k)*weights2(j,k)
               output_array((i - 1)*degree + j) = output_array((i - 1)*degree + j) + &
                                                  data(ind1 + k)*weights1(k, j) + data(ind2 + k)*weights2(k, j)
            end do
            !output_array((i-1)*degree+j) = output_array((i-1)*degree+j)/quadrature_weights(j)
            output_array((i - 1)*degree + j) = output_array((i - 1)*degree + j)*quadrature_inv_weights(j)
         end do
      end do
   end subroutine sll_s_dg_interpolator_1d_interpolate_array_disp_halo

   subroutine sll_s_dg_interpolator_1d_dg_to_equi_cell(self, npoints, output_array)
      type(sll_t_dg_interpolator_1d), intent(in) :: self
      sll_int32, intent(in) :: npoints
      sll_real64, intent(out) :: output_array(:, :)

      sll_int32 :: j, idx, degree
      sll_real64 :: delta
      sll_real64, pointer :: quadrature_points(:)
      sll_real64, pointer :: quadrature_inv_diffs(:, :)

      delta = 1.0_f64/real(npoints - 1, f64); 
      degree = self%degree
      quadrature_points => self%quadrature_points
      quadrature_inv_diffs => self%quadrature_inv_diffs

      do j = 1, npoints
         do idx = 1, degree
!          output_array(idx,j) = lagrange_poly( self, idx, real(j-1,f64)*delta )
            output_array(idx, j) = lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, idx, real(j - 1, f64)*delta)
         end do
      end do
   end subroutine sll_s_dg_interpolator_1d_dg_to_equi_cell

   subroutine sll_s_dg_interpolator_1d_dg_to_equi_cell_staggered(self, npoints, output_array)
      type(sll_t_dg_interpolator_1d), intent(in) :: self
      sll_int32, intent(in) :: npoints
      sll_real64, intent(out) :: output_array(:, :)

      sll_int32 :: j, idx, degree
      sll_real64 :: delta
      sll_real64, pointer :: quadrature_points(:)
      sll_real64, pointer :: quadrature_inv_diffs(:, :)

      delta = 1.0_f64/real(npoints, f64); 
      degree = self%degree
      quadrature_points => self%quadrature_points
      quadrature_inv_diffs => self%quadrature_inv_diffs

      do j = 1, npoints
         do idx = 1, degree
!          output_array(idx,j) = lagrange_poly( self, idx, real(j-1,f64)*delta )
  output_array(idx, j) = lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, idx, (real(j - 1, f64) + 0.5_f64)*delta)
         end do
      end do
   end subroutine sll_s_dg_interpolator_1d_dg_to_equi_cell_staggered

   !> Evaluate Lagrange polynomial, original version.  Suffers from
   !> repeated dereferencing, and repeated difference computation and division.
   !> Kept for documentation purposes.
   function lagrange_poly(self, idx, x) result(r)
      type(sll_t_dg_interpolator_1d), intent(in) :: self
      sll_int32, intent(in) :: idx
      sll_real64, intent(in) :: x
      sll_real64 :: r

      sll_int32 :: j

      r = 1.0_f64
      do j = 1, self%degree
         if (j .ne. idx) then
            r = r*(x - self%quadrature_points(j))/(self%quadrature_points(idx) - self%quadrature_points(j))
         end if
      end do
   end function lagrange_poly

   !> Evaluate Lagrange polynomial, optimized version.
!DIR$ ATTRIBUTES FORCEINLINE :: lagrange_poly_opt
   function lagrange_poly_opt(quadrature_points, quadrature_inv_diffs, degree, idx, x) result(r)
      sll_real64, intent(in) :: quadrature_points(:)
      sll_real64, intent(in) :: quadrature_inv_diffs(:, :)
      sll_int32, intent(in) :: degree
      sll_int32, intent(in) :: idx
      sll_real64, intent(in) :: x

      sll_real64 :: r
      sll_int32 :: j

      r = 1.0_f64
      do j = 1, degree
         if (j .ne. idx) then
            r = r*(x - quadrature_points(j))*quadrature_inv_diffs(j, idx)
         end if
      end do
   end function lagrange_poly_opt

!  function lagrange_poly_opt_bak( quadrature_points, degree, idx, x ) result(r)
!    sll_real64, intent(in) :: quadrature_points(:)
!    sll_int32, intent(in) :: degree
!    sll_int32, intent(in) :: idx
!    sll_real64, intent(in) :: x
!
!    sll_real64 :: r
!    sll_int32 :: j
!
!    r = 1.0_f64
!    do j=1,degree
!       if (j .ne. idx) then
!          r = r * (x-quadrature_points(j))/(quadrature_points(idx)-quadrature_points(j))
!       end if
!    end do
!  end function lagrange_poly_opt_bak

end module sll_m_dg_interpolator_1d
