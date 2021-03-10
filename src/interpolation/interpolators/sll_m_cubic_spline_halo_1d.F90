#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif

! NOTE: Check the initialization of pba_pow(:) below when changing NUM_TERMS.
! Usual value (and maximum useful value in terms of double precision accuracy): 27
!!! #define NUM_TERMS 27
! Temporarily reduced to 15 terms to enable the test case at 16**6 resolution.
#define NUM_TERMS 15

! Note: Uncomment one of the following macros to select how powers of p_b_a are computed.
! fast, pre-computed array
#define PBA_POW(i)   pba_pow(i)
! slow, repeated computation, original implementation
!#define PBA_POW(i)   (-p_b_a)**(i)
!> @ingroup interpolators
!> @author Katharina Kormann, IPP
!> @brief
!> Interpolator 1d using cubic splines on regular mesh with halo cells
!> @details
!> The module provides an optimized implementation of the interpolator function
!> interpolate_array_disp for a domain with halo cells, i.e. we do not have
!> explicit boundary conditions but compute the Lagrange interpolation for a
!> number of cells provided that enough cells around this domain are present
!> such that no boundary conditions need to be imposed.
!> The module also provides a version for periodic boundary conditions.
!> The implementation is based on the algorithms described in Section 5.4.4
!> (Fast local spline interpolation) of Sonnendrücker, Numerical Methods for the
!> Vlasov-Maxwell equations, to appear.
!!
module sll_m_cubic_spline_halo_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_s_cubic_spline_halo_1d_prepare_exchange, &
      sll_s_cubic_spline_halo_1d_finish_boundary_conditions, &
      sll_s_cubic_spline_halo_1d_compute_interpolant, &
      sll_s_cubic_spline_halo_1d_eval_disp, &
      sll_s_cubic_spline_halo_1d_periodic
   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Compile time constants for fast computations
   sll_real64, parameter :: p_a = sqrt((2.0_f64 + sqrt(3.0_f64))/6.0_f64)
   sll_real64, parameter :: p_r_a = 1.0_f64/p_a
   sll_real64, parameter :: p_b = sqrt((2.0_f64 - sqrt(3.0_f64))/6.0_f64)
   sll_real64, parameter :: p_b_a = p_b/p_a
   sll_real64, parameter :: p_sqrt3 = sqrt(3.0_f64)
   sll_real64, parameter :: p_inv_6 = 1._f64/6._f64

!  $ cat pba_pow.py
!  for i in range(0,28): print ("(-p_b_a)**%d," % i),
   sll_real64, dimension(0:27), parameter :: pba_pow = [ &
                                        (-p_b_a)**0, (-p_b_a)**1, (-p_b_a)**2, (-p_b_a)**3, (-p_b_a)**4, (-p_b_a)**5, (-p_b_a)**6, &
                                    (-p_b_a)**7, (-p_b_a)**8, (-p_b_a)**9, (-p_b_a)**10, (-p_b_a)**11, (-p_b_a)**12, (-p_b_a)**13, &
                                 (-p_b_a)**14, (-p_b_a)**15, (-p_b_a)**16, (-p_b_a)**17, (-p_b_a)**18, (-p_b_a)**19, (-p_b_a)**20, &
                                  (-p_b_a)**21, (-p_b_a)**22, (-p_b_a)**23, (-p_b_a)**24, (-p_b_a)**25, (-p_b_a)**26, (-p_b_a)**27 &
                                                       ]

contains

   !> Compute the part of \a d(0) and \a c(num_points+2) that need to be send to neighboring processors
!DIR$ ATTRIBUTES FORCEINLINE :: sll_s_cubic_spline_halo_1d_prepare_exchange
   subroutine sll_s_cubic_spline_halo_1d_prepare_exchange(fdata, si, num_points, d_0, c_np2)
      sll_real64, intent(in) :: fdata(:) !< Local data values
      sll_int32, intent(in) :: si !< Integer part of the shift
      sll_int32, intent(in) :: num_points !< number of local data values
      HALO_DTYPE, intent(out) :: d_0 !< Initialization for recursion algorithm (forward loop)
      HALO_DTYPE, intent(out) :: c_np2 !< Initializatio for recursion algorithm (backward loop)

      sll_int32 :: i, ind_min

      ! Hard assert: In case less points are used buffer overruns do occur that may stay unnoticed!
      SLL_ASSERT_ALWAYS(num_points > NUM_TERMS)

      if (si > 0) then
         d_0 = 0.0_f64
         ind_min = si
      else
         d_0 = fdata(num_points + si)
         ind_min = 1
      end if
      do i = ind_min, NUM_TERMS
         d_0 = d_0 + PBA_POW(i)*fdata(num_points + si - i)
      end do
      if (si < -1) then
         c_np2 = 0.0_f64
         ind_min = -si - 1
      else
         c_np2 = fdata(2 + si)
         ind_min = 1
      end if
      do i = ind_min, NUM_TERMS
         c_np2 = c_np2 + PBA_POW(i)*fdata(2 + si + i)
      end do
      do i = 1, si + 1
         c_np2 = c_np2 + PBA_POW(i)*fdata(2 + si - i)
      end do
   end subroutine sll_s_cubic_spline_halo_1d_prepare_exchange

   !> Complete \a d(0) and \a c(num_points+2) with local data after their values have been received from the neighboring processors
!DIR$ ATTRIBUTES FORCEINLINE :: sll_s_cubic_spline_halo_1d_finish_boundary_conditions
   subroutine sll_s_cubic_spline_halo_1d_finish_boundary_conditions(fdata, si, num_points, d_0, c_np2)
      sll_real64, intent(in) :: fdata(:) !< Local data values
      sll_int32, intent(in) :: si !< Integer part of the shift
      sll_int32, intent(in) :: num_points !< number of local data values
      HALO_DTYPE, intent(out) :: d_0 !< Initialization for recursion algorithm (forward loop)
      HALO_DTYPE, intent(inout) :: c_np2 !< Initializatio for recursion algorithm (backward loop)

      sll_int32 :: i, ind_min

      if (si > 0) then
         d_0 = d_0 + fdata(si)
      end if
      do i = 1, si - 1
         d_0 = d_0 + PBA_POW(i)*fdata(si - i)
      end do
      d_0 = d_0*p_r_a
      if (si < -1) then
         c_np2 = c_np2 + fdata(num_points + 2 + si)
         ind_min = 1
      else
         ind_min = si + 2
      end if
      do i = 1, -si - 2
         c_np2 = c_np2 + PBA_POW(i)*fdata(num_points + 2 + si + i)
      end do
      do i = ind_min, NUM_TERMS
         c_np2 = c_np2 + PBA_POW(i)*fdata(num_points + 2 + si - i)
      end do
      c_np2 = c_np2*p_sqrt3
   end subroutine sll_s_cubic_spline_halo_1d_finish_boundary_conditions

   !> Compute the coefficients of the local interpolating spline (after \a d(0) and \a c(num_points+2) have been computed
!DIR$ ATTRIBUTES FORCEINLINE :: sll_s_cubic_spline_halo_1d_compute_interpolant
   subroutine sll_s_cubic_spline_halo_1d_compute_interpolant(f, num_points, d, coeffs)
      sll_real64, intent(in)   :: f(0:) !< data values including all points to be interpolated at
      sll_int32, intent(in)   :: num_points !< number of local data points
      sll_real64, intent(out)   :: d(0:) !< helper variable for spline coefficients (vales from forward recursion)
      sll_real64, intent(out)   :: coeffs(0:) !< spline coefficients

      sll_int32                         :: i
      sll_int32                         :: np
      logical, save :: first_call = .true.

      if ((first_call) .and. (NUM_TERMS < 27)) then
         write (*, *) "WARNING: sll_s_cubic_spline_halo_1d uses NUM_TERMS=", NUM_TERMS
      end if
      first_call = .false.

      SLL_ASSERT(size(f) .ge. num_points - 1)
      SLL_ASSERT(size(d) .ge. num_points)
      SLL_ASSERT(size(coeffs) .ge. num_points)
      ! Hard assert: In case less points are used buffer overruns do occur that may stay unnoticed!
      SLL_ASSERT_ALWAYS(num_points > NUM_TERMS)

      np = num_points
      d(0) = f(0)

      do i = 1, np + 1
         d(i) = p_r_a*(f(i) - p_b*d(i - 1))
      end do
      coeffs(np + 2) = f(np + 2)!d1*r_a
      ! remaining coefficients:
      do i = np + 1, 0, -1
         coeffs(i) = p_r_a*(d(i) - p_b*coeffs(i + 1))
      end do
   end subroutine sll_s_cubic_spline_halo_1d_compute_interpolant

   !> This function corresponds to the interpolate_array_disp function of the interpolators but the displacement is normalized and between [0,1]
!DIR$ ATTRIBUTES FORCEINLINE :: sll_s_cubic_spline_halo_1d_eval_disp
   subroutine sll_s_cubic_spline_halo_1d_eval_disp(coeffs, alpha, num_points, fout)
      sll_real64, intent(in) :: coeffs(0:) !< Spline coefficients (centered around the cell into which we displace, i.e. integer part of displacement is already build in here)
      sll_real64, intent(in) :: alpha !< Displacement normalized by dx and only remainder of modulo 1   (in [0,1])
      sll_int32, intent(in) :: num_points !< Number of local data points
      sll_real64, intent(out) :: fout(:) !< Interpolated values

      sll_real64 :: calpha, cim1, ci, cip1, cip2, t1, t2, t3, t4
      sll_int32  :: cell

      calpha = 1.0_f64 - alpha

      do cell = 1, num_points
         cim1 = coeffs(cell - 1)
         ci = coeffs(cell)
         cip1 = coeffs(cell + 1)
         cip2 = coeffs(cell + 2)
         t1 = 3.0_f64*ci
         t3 = 3.0_f64*cip1
         t2 = calpha*(calpha*(calpha*(cim1 - t1) + t1) + t1) + ci
         t4 = alpha*(alpha*(alpha*(cip2 - t3) + t3) + t3) + cip1
         fout(cell) = p_inv_6*(t2 + t4)
      end do
   end subroutine sll_s_cubic_spline_halo_1d_eval_disp

!DIR$ ATTRIBUTES FORCEINLINE :: sll_s_cubic_spline_halo_1d_periodic
   subroutine sll_s_cubic_spline_halo_1d_periodic(fin, alpha, num_cells, fout)
      sll_real64, intent(inout) :: fin(0:)
      sll_real64, intent(in) :: alpha
      sll_int32, intent(in) :: num_cells
      sll_real64, intent(out) :: fout(:)

      call sll_s_cubic_spline_halo_1d_compute_periodic(num_cells, fout, fin)
      call sll_s_cubic_spline_halo_1d_eval_disp(fin, alpha, num_cells, fout)
   end subroutine sll_s_cubic_spline_halo_1d_periodic

   !> Compute the coefficients of the local interpolating spline (after \a d(0) and \a c(num_points+2) have been computed
!DIR$ ATTRIBUTES FORCEINLINE :: sll_s_cubic_spline_halo_1d_compute_periodic
   subroutine sll_s_cubic_spline_halo_1d_compute_periodic(num_points, d, f_coeffs)
      sll_int32, intent(in)   :: num_points !< number of local data points
      sll_real64, intent(out)   :: d(0:) !< helper variable for spline coefficients (vales from forward recursion)
      sll_real64, intent(inout)   :: f_coeffs(0:) !< on input: data values including all points to be interpolated at; on output: spline coefficients

      sll_int32                         :: i
      sll_int32                         :: np
      SLL_ASSERT(size(d) .ge. num_points)
      SLL_ASSERT(size(f_coeffs) .ge. num_points)
      ! Hard assert: In case less points are used buffer overruns do occur that may stay unnoticed!
      SLL_ASSERT_ALWAYS(num_points > NUM_TERMS)

      np = num_points
      d(0) = f_coeffs(0)
      do i = 1, NUM_TERMS
         d(0) = d(0) + PBA_POW(i)*f_coeffs(num_points - i)
      end do
      d(0) = d(0)*p_r_a

      do i = 1, np - 1
         d(i) = p_r_a*(f_coeffs(i) - p_b*d(i - 1))
      end do

      f_coeffs(np) = d(np - 1)
      do i = 1, NUM_TERMS
         f_coeffs(np) = f_coeffs(np) + d(i - 1)*PBA_POW(i)
      end do
      f_coeffs(np) = f_coeffs(np)*p_r_a

      ! remaining coefficients:
      do i = np - 1, 1, -1
         f_coeffs(i) = p_r_a*(d(i - 1) - p_b*f_coeffs(i + 1))
      end do

      f_coeffs(0) = f_coeffs(np)
      f_coeffs(np + 1:np + 2) = f_coeffs(1:2)
   end subroutine sll_s_cubic_spline_halo_1d_compute_periodic

end module sll_m_cubic_spline_halo_1d
