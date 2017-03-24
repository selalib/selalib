#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif

program test_cubic_spline_halo_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_t_cubic_spline_interpolator_1d

  use sll_m_interpolators_1d_base, only: &
       sll_c_interpolator_1d

  use sll_m_cubic_spline_halo_1d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  sll_int32, parameter  :: n = 64
!  sll_real64, parameter :: tol = 1.d-14 ! tolerance for error
  sll_real64, parameter :: tol = 4.d-9 ! Update: tolerance for error, compatibel with NUM_TERMS=15
  sll_real64            :: error
  logical               :: passed

  type(sll_t_cubic_spline_interpolator_1d), target :: spline

  sll_real64 :: pdata(n+1)
  sll_real64 :: pdata_per(n+3)
  sll_real64 :: pinterp1(n)
  sll_real64 :: pinterp2(n+1)
  sll_real64 :: coord(n+1)

   sll_int32 :: ierr, i

   sll_real64  :: x_min, x_max, delta

   sll_real64 :: alpha, beta
   sll_int32 :: si

   print*, 'Initialize data and point array'
   x_min = 0.0_f64
   x_max = 2.0_f64 * sll_p_pi
   delta = (x_max - x_min ) / real(n,f64)
   do i=1,n+1
      coord(i) = (i-1)*delta
      pdata(i) = f(coord(i))
   end do

   !print*, pdata

   print*, 'Cubic spline interpolation'
   call spline%init &
        (n+1, x_min, x_max, sll_p_periodic, fast_algorithm=.true. )
   !call spline%compute_interpolants(pdata)

   alpha = 0.25_f64
   passed = .true.

   do si = -2, 2
      beta = (real(si,f64)+ alpha)*delta

      call test_halo_version( si, alpha, pdata, pinterp1 )

      call spline%interpolate_array_disp ( n+1, pdata, beta, pinterp2 )

      error = maxval(abs(pinterp1(1:n)-pinterp2(1:n)))
      print*, 'Error for cell displacement', si , ":", error
      if ( error > tol ) passed = .false.
   end do


   pdata_per(1:n) = pdata(1:n)
   call sll_s_cubic_spline_halo_1d_periodic( pdata_per, alpha, n, pinterp1  )
   beta = alpha*delta
   call spline%interpolate_array_disp ( n+1, pdata, beta, pinterp2 )
   error = maxval(abs(pinterp1(1:n)-pinterp2(1:n)))
   print*, 'Error for sll_s_cubic_spline_halo_1d_periodic:' , ":", error
   if ( error > tol ) passed = .false.

   if (passed .eqv. .true. ) then
      print*, 'PASSED.'
   else
      print*, 'FAILED.'
   end if

contains


  subroutine test_halo_version( si, alpha, fdata, fout )
    sll_int32,  intent( in    ) :: si
    sll_real64, intent( in    ) :: alpha
    sll_real64, intent( in    ) :: fdata(n+1)
    sll_real64, intent(   out ) :: fout(n)

    sll_real64 :: coeffs(n+3)
    sll_real64 :: dvec(n+3)
    HALO_DTYPE :: d_0, c_np1
    sll_real64 :: c_np2
    sll_int32 :: i
    sll_real64 :: fin(n+3)

    c_np2 = fdata(modulo(1+si,n)+1)
    do i =1,27
       c_np2 = c_np2 + (sqrt(3.0_f64)-2.0_f64)**i * (fdata(modulo(n+1-i+si, n)+1)+fdata(modulo(1+i+si, n)+1))
    end do

    c_np2 = sqrt(3.0_f64)*c_np2
    c_np1 = 0.0_f64

    call sll_s_cubic_spline_halo_1d_prepare_exchange( fdata(1:n), si, n, d_0, c_np1 )
    call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( fdata(1:n), si, n, d_0, c_np1 )
    !print*,'cn 2', c_np1, c_np2

    fin(1) = d_0
    fin(n+3) = c_np1

    do i =1,n+1
       fin(i+1) = fdata(modulo(si+i-1, n) + 1)
    end do

    call sll_s_cubic_spline_halo_1d_compute( fin, n, dvec, coeffs )
    call sll_s_cubic_spline_halo_1d_advect ( coeffs, alpha, n, fout )
    !print*, 'halo'
    !print*, coeffs
  end subroutine test_halo_version


  function f(x)
    sll_real64 :: x
    sll_real64 :: f

    f = 2.0_f64*(sin(x) + 2.5_f64 + cos(x))
  end function f

end program test_cubic_spline_halo_1d
