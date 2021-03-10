module m_test_spline_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: &
      f64

   use m_analytical_profiles_1d_base, only: &
      t_profile_1d_info, &
      c_analytical_profile_1d

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic, &
      sll_p_hermite, &
      sll_p_greville

   use sll_m_bsplines, only: &
      sll_c_bsplines, &
      sll_s_bsplines_new

   use sll_m_spline_1d, only: &
      sll_t_spline_1d

   use sll_m_spline_interpolator_1d, only: &
      sll_t_spline_interpolator_1d

   use sll_m_timer, only: &
      sll_t_time_mark, &
      sll_s_set_time_mark, &
      sll_f_time_elapsed_between

   implicit none

   public :: &
      t_spline_1d_test_facility

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> Type for running test
   type :: t_spline_1d_test_facility

      class(c_analytical_profile_1d), pointer :: profile_1d

      class(sll_c_bsplines), allocatable :: bsplines
      type(sll_t_spline_1d)              :: spline_1d
      type(sll_t_spline_interpolator_1d) :: spline_interpolator_1d

      real(wp), allocatable ::  tau(:) ! interpolation points
      real(wp), allocatable :: gtau(:) ! profile values at tau

      real(wp) :: time_init
      real(wp) :: time_compute_interpolant
      real(wp) :: time_eval
      real(wp) :: time_eval_array
      real(wp) :: time_eval_deriv
      real(wp) :: time_eval_array_deriv

   contains

      procedure :: init
      procedure :: free
      procedure :: check_equivalence_scalar_array_methods
      procedure :: evaluate_on_1d_grid
      procedure :: evaluate_at_interpolation_points
      procedure :: evaluate_deriv_on_1d_grid

   end type t_spline_1d_test_facility

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   subroutine init(self, profile_1d, degree, ncells, bc_xmin, bc_xmax, breaks)

      class(t_spline_1d_test_facility), intent(out)         :: self
      class(c_analytical_profile_1d), intent(in), target :: profile_1d
      integer, intent(in)         :: degree
      integer, intent(in)         :: ncells
      integer, intent(in)         :: bc_xmin
      integer, intent(in)         :: bc_xmax
      real(wp), optional, intent(in)         :: breaks(:)

      type(t_profile_1d_info) :: info
      integer                 :: nipts
      integer                 :: i, j, s
      logical                 :: periodic

      real(wp), allocatable :: derivs_xmin(:)
      real(wp), allocatable :: derivs_xmax(:)

      type(sll_t_time_mark) :: t0, t1, t2, t3

      ! Store pointer to 1D profile
      self%profile_1d => profile_1d

      ! Extract information about 1D analytical profile
      call self%profile_1d%get_info(info)

      !TODO: add proper checks
      periodic = (bc_xmin == sll_p_periodic)

      call sll_s_set_time_mark(t0)

      ! Create B-splines (uniform or non-uniform depending on input)
      call sll_s_bsplines_new(self%bsplines, degree, periodic, info%xmin, info%xmax, ncells, breaks)

      ! Initialize 1D spline
      call self%spline_1d%init(self%bsplines)

      ! Initialize 1D spline interpolator
      call self%spline_interpolator_1d%init(self%bsplines, bc_xmin, bc_xmax)

      call sll_s_set_time_mark(t1)

      ! Get spline interpolation points => 'tau' array
      call self%spline_interpolator_1d%get_interp_points(self%tau)
      nipts = size(self%tau)

      ! Evaluate analytical profile at interpolation points
      allocate (self%gtau(nipts))
      do i = 1, nipts
         self%gtau(i) = self%profile_1d%eval(self%tau(i))
      end do

      ! If needed, evaluate x derivatives at xmin
      if (bc_xmin == sll_p_hermite) then
         allocate (derivs_xmin(degree/2))
         s = 1 - modulo(degree, 2) ! shift = 1 for even order, 0 for odd order
         do j = 1, degree/2
            derivs_xmin(j) = self%profile_1d%eval(info%xmin, diff=j - s)
         end do
      end if

      ! If needed, evaluate x derivatives at xmax
      if (bc_xmax == sll_p_hermite) then
         allocate (derivs_xmax(degree/2))
         s = 1 - modulo(degree, 2) ! shift = 1 for even order, 0 for odd order
         do j = 1, degree/2
            derivs_xmax(j) = self%profile_1d%eval(info%xmax, diff=j - s)
         end do
      end if

      ! Compute 1D spline that interpolates analytical 1D profile at points above
      call sll_s_set_time_mark(t2)

      ! Hermite - Greville
      if (bc_xmin == sll_p_hermite .and. bc_xmax == sll_p_greville) then

         call self%spline_interpolator_1d%compute_interpolant( &
            self%spline_1d, &
            self%gtau, &
            derivs_xmin=derivs_xmin)

         ! Greville - Hermite
      else if (bc_xmin == sll_p_greville .and. bc_xmax == sll_p_hermite) then

         call self%spline_interpolator_1d%compute_interpolant( &
            self%spline_1d, &
            self%gtau, &
            derivs_xmax=derivs_xmax)

         ! Hermite - Hermite
      else if (bc_xmin == sll_p_hermite .and. bc_xmax == sll_p_hermite) then

         call self%spline_interpolator_1d%compute_interpolant( &
            self%spline_1d, &
            self%gtau, &
            derivs_xmin=derivs_xmin, &
            derivs_xmax=derivs_xmax)

         ! Periodic or Greville - Greville
      else

         call self%spline_interpolator_1d%compute_interpolant( &
            self%spline_1d, &
            self%gtau)

      end if

      call sll_s_set_time_mark(t3)

      ! Deallocate local arrays
      if (allocated(derivs_xmin)) deallocate (derivs_xmin)
      if (allocated(derivs_xmax)) deallocate (derivs_xmax)

      ! Timings (set to -1 values not yet available)
      self%time_init = sll_f_time_elapsed_between(t0, t1)
      self%time_compute_interpolant = sll_f_time_elapsed_between(t2, t3)
      self%time_eval = -1.0_wp
      self%time_eval_array = -1.0_wp
      self%time_eval_deriv = -1.0_wp
      self%time_eval_array_deriv = -1.0_wp

   end subroutine init

   !-----------------------------------------------------------------------------
   subroutine check_equivalence_scalar_array_methods(self, equiv)

      class(t_spline_1d_test_facility), intent(inout) :: self
      logical, intent(out) :: equiv(3)

      type(t_profile_1d_info) :: info
      real(wp), allocatable   :: x(:)
      real(wp), allocatable   :: y(:)
      real(wp), allocatable   :: ya(:)
      real(wp)                :: r
      integer                 :: i
      type(sll_t_time_mark)   :: t0, t1, t2

      ! Choose size of 1D grid
      integer, parameter :: n = 77
      real(wp), parameter :: npts = real(n, wp)

      ! Allocate local arrays
      allocate (x(n))
      allocate (y(n))
      allocate (ya(n))

      ! Get info from profile
      call self%profile_1d%get_info(info)

      ! Generate random grid in x
      do i = 1, n
         call random_number(r)
         x(i) = info%xmin*(1.0_wp - r) + info%xmax*r
      end do

      ! Compare:
      !   . y(i) =  spline_1d % eval( x(i) )   for all i
      !   . call spline_1d % eval_array( x(:), ya(:) )
      call sll_s_set_time_mark(t0)
      do i = 1, n
         y(i) = self%spline_1d%eval(x(i))
      end do
      call sll_s_set_time_mark(t1)
      call self%spline_1d%eval_array(x, ya)
      call sll_s_set_time_mark(t2)

      self%time_eval = sll_f_time_elapsed_between(t0, t1)/npts
      self%time_eval_array = sll_f_time_elapsed_between(t1, t2)/npts
      equiv(1) = all(y == ya)

      ! Compare:
      !   . y(i) = spline_1d % eval_deriv( x(i) )  for all i
      !   . call spline_1d % eval_array_deriv( x(:), ya(:) )
      call sll_s_set_time_mark(t0)
      do i = 1, n
         y(i) = self%spline_1d%eval_deriv(x(i))
      end do
      call sll_s_set_time_mark(t1)
      call self%spline_1d%eval_array_deriv(x, ya)
      call sll_s_set_time_mark(t2)

      self%time_eval_deriv = sll_f_time_elapsed_between(t0, t1)/npts
      self%time_eval_array_deriv = sll_f_time_elapsed_between(t1, t2)/npts
      equiv(2) = all(y == ya)

      ! Deallocate local arrays
      deallocate (x)
      deallocate (y)
      deallocate (ya)

   end subroutine check_equivalence_scalar_array_methods

   !-----------------------------------------------------------------------------
   subroutine evaluate_at_interpolation_points(self, max_norm_error)

      class(t_spline_1d_test_facility), intent(in) :: self
      real(wp), intent(out) :: max_norm_error

      integer  :: i
      real(wp) :: error

      ! Evaluate 1D spline at interpolation points:
      ! interpolation values should be obtained
      max_norm_error = 0.0_wp
      do i = 1, size(self%tau)
         error = self%gtau(i) - self%spline_1d%eval(self%tau(i))
         max_norm_error = max(max_norm_error, abs(error))
      end do

   end subroutine evaluate_at_interpolation_points

   !-----------------------------------------------------------------------------
   subroutine evaluate_on_1d_grid(self, x, max_norm_error)

      class(t_spline_1d_test_facility), intent(inout) :: self
      real(wp), intent(in) :: x(:)
      real(wp), intent(out) :: max_norm_error

      integer               :: i
      integer               :: n
      real(wp), allocatable :: y(:)
      real(wp)              :: error
      type(sll_t_time_mark) :: t0, t1

      n = size(x)

      ! Array of spline values
      allocate (y(n))

      call sll_s_set_time_mark(t0)

      ! Evaluate 1D spline on given grid
      call self%spline_1d%eval_array(x, y)

      call sll_s_set_time_mark(t1)
      self%time_eval_array = sll_f_time_elapsed_between(t0, t1)/real(n, wp)

      ! Compare spline values to analytical profile and compute max norm of error
      max_norm_error = 0.0_wp
      do i = 1, n
         error = self%profile_1d%eval(x(i)) - y(i)
         max_norm_error = max(max_norm_error, abs(error))
      end do

      ! Free local memory
      deallocate (y)

   end subroutine evaluate_on_1d_grid

   !-----------------------------------------------------------------------------
   subroutine evaluate_deriv_on_1d_grid(self, x, max_norm_error_deriv)

      class(t_spline_1d_test_facility), intent(inout) :: self
      real(wp), intent(in) :: x(:)
      real(wp), intent(out) :: max_norm_error_deriv

      integer               :: i
      integer               :: n
      real(wp), allocatable :: y(:)
      real(wp)              :: error
      type(sll_t_time_mark) :: t0, t1

      n = size(x)

      ! Array of spline values
      allocate (y(n))

      ! x-derivative: compare spline to analytical profile
      call sll_s_set_time_mark(t0)

      call self%spline_1d%eval_array_deriv(x, y)

      call sll_s_set_time_mark(t1)
      self%time_eval_array_deriv = sll_f_time_elapsed_between(t0, t1)/real(n, wp)

      max_norm_error_deriv = 0.0_wp
      do i = 1, n
         error = self%profile_1d%eval(x(i), diff=1) - y(i)
         max_norm_error_deriv = max(max_norm_error_deriv, abs(error))
      end do

      ! Free local memory
      deallocate (y)

   end subroutine evaluate_deriv_on_1d_grid

   !-----------------------------------------------------------------------------
   subroutine free(self)

      class(t_spline_1d_test_facility), intent(inout) :: self

      ! Call destructors of member objects
      call self%bsplines%free()
      call self%spline_1d%free()
      call self%spline_interpolator_1d%free()

      ! Free dynamically allocated memory
      deallocate (self%bsplines)
      deallocate (self%tau)
      deallocate (self%gtau)

   end subroutine free

end module m_test_spline_1d
