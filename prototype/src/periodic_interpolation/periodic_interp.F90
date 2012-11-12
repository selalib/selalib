module periodic_interp_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  implicit none

  sll_real64, parameter    :: pi = 3.1415926535897932385_8
  sll_real64, parameter    :: twopi = 6.2831853071795864769_8
  complex(8), parameter :: ii = complex(0.0_8, 1.0_8)
  type :: periodic_interp_work
     sll_int32          :: N ! number of cells
     sll_int32          :: interpolator ! what interpolator is used
     sll_int32          :: order  ! order of interpolation (not needed for Fourier)
     sll_real64, dimension(:), pointer :: eigenvalues_Minv ! eigenvalues of M matrix
     complex(8), dimension(:), pointer :: eigenvalues_S   ! eigenvalues of shift matrix
     sll_real64, dimension(:), pointer :: knots   ! for B-splines
     sll_real64, dimension(:), pointer :: wsave  ! workspace for fft
     complex(8), dimension(:), pointer :: modes  ! Fourier modes
     complex(8), dimension(:), pointer :: ufft   ! Fourier transform of function
     complex(8), dimension(:), pointer :: uvfft   ! Fourier transform of function
  end type periodic_interp_work

  interface delete
     module procedure delete_periodic_interp_work
  end interface delete

  enum, bind(C)
     enumerator :: TRIGO = 0, SPLINE = 1
  end enum

contains
  subroutine initialize_periodic_interp(this,N,interpolator,order)
    type(periodic_interp_work) :: this 
    sll_int32 :: N ! number of cells
    sll_int32 :: interpolator  ! kind of interpolator
    sll_int32 :: order  ! order of method
    ! local variables
    sll_int32 :: ierr, i, j, k
    sll_int32 :: p ! spline degree
    sll_int32 :: icoarse  ! coarsening factor for stabilization of trigonometric interpolation
    sll_real64, dimension(order) :: biatx
    sll_real64 :: mode, val

    this%N = N
    this%interpolator = interpolator
    !if (present(order)) then 
       this%order = order 
    !end if
    ! Allocate arrays
    SLL_ALLOCATE(this%wsave(15+4*N),ierr)
    SLL_ALLOCATE(this%ufft(N),ierr)
    SLL_ALLOCATE(this%eigenvalues_Minv(N),ierr)
    SLL_ALLOCATE(this%eigenvalues_S(N),ierr)
    SLL_ALLOCATE(this%modes(0:N-1),ierr)

    ! set up spline parameters
    if ((order/2) /= int(order/2.0)) then
       print*, 'initialize_periodic_interp: order of splines needs to be even.', &
            'Order here is: ', Order
       stop
    end if
    p = order - 1  ! spline degree is one less than order
    ! set up knots array for computation of B-splines using De Boor's routines
    allocate(this%knots(2*order),STAT=ierr)
    
    ! Initialize fft
    call zffti(N, this%wsave)

    select case (interpolator)
    case (TRIGO)
       if (order == 0) then ! no coarsening
          this%eigenvalues_Minv = 1.0_8
       else 
          do i=-p,p+1
             this%knots(i+p+1)=real(i*icoarse,8)
          enddo
          ! to be implemented
          this%eigenvalues_Minv = 1.0_8
       end if

    case (SPLINE)
       do i=-p,p+1
          this%knots(i+p+1)=real(i,8)
       enddo
       call bsplvb(this%knots,p+1,1,0.0_8,p+1,biatx)
       do i=1, N
          this%modes(i-1) = exp(ii*twopi*(i-1)/N)
          this%eigenvalues_Minv(i) = biatx((p+1)/2)
          do j = 1,(p+1)/2
             this%eigenvalues_Minv(i) = this%eigenvalues_Minv(i) &
                  + biatx(j+(p+1)/2)*2*cos(j*twopi*(i-1)/N)
          end do
          this%eigenvalues_Minv(i) = 1.0_8 / this%eigenvalues_Minv(i)
       end do
    case default
       print*, 'periodic_interp_module:interpolator ',interpolator, ' not implemented'
       stop
    end select

  end subroutine initialize_periodic_interp

  subroutine delete_periodic_interp_work(this)
    type(periodic_interp_work) :: this 
    sll_int32 :: ierr
    
    SLL_DEALLOCATE(this%wsave, ierr)
    SLL_DEALLOCATE(this%ufft,ierr)
    SLL_DEALLOCATE(this%eigenvalues_Minv,ierr)
    SLL_DEALLOCATE(this%eigenvalues_S,ierr)
    SLL_DEALLOCATE(this%modes,ierr)

  end subroutine delete_periodic_interp_work

  subroutine periodic_interp(this, u_out, u, alpha) 
    ! interpolate function u given at grid points on a periodic grid
    ! at positions j-alpha (alpha is normalized to the cell size)
    type(periodic_interp_work), intent(inout) :: this
    sll_real64, dimension(:), intent(in)   :: u  ! function to be interpolated
    sll_real64, dimension(:), intent(out)   :: u_out  ! result
    sll_real64, intent(in)     :: alpha ! displacement normalized to cell size
    ! local variables
    sll_int32 :: i, j, k, p, ishift, j0, imode
    sll_real64 :: beta, filter, mode
    complex(8) :: int_fact, z
    ! 
    sll_real64, dimension(this%order) :: biatx

    ! Perform FFT of u
    do i=1, this%N
       this%ufft(i) = complex(u(i),8)
    end do
    call zfftf(this%N, this%ufft, this%wsave)

    ! Compute eigenvalues of shift matrix
    select case (this%interpolator)
    case (TRIGO)
       this%eigenvalues_S(1) = 1.0_8
       this%eigenvalues_S(this%N/2+1) = exp(-ii*pi*alpha)
       do k=1, this%N/2-1
          !filter = 0.5_8*(1+tanh(100*(.35_8*this%N/2-k)/(this%N/2))) !F1
          !filter = 0.5_8*(1+tanh(100*(.25_8*this%N/2-k)/(this%N/2))) !F2
          !filter = 0.5_8*(1+tanh(50*(.25_8*this%N/2-k)/(this%N/2))) !F3
          filter = 1.0_8

          this%eigenvalues_S(k+1) = exp(-ii*twopi*k*alpha/this%N) * filter
          this%eigenvalues_S(this%N-k+1) = exp(ii*twopi*k*alpha/this%N) * filter
       end do
       this%ufft = this%ufft*this%eigenvalues_S
    case (SPLINE)
       !compute eigenvalues of splines evaluated at displaced points
       p =  this%order - 1
       ishift = floor (-alpha)
       beta = -ishift - alpha
       call bsplvb(this%knots,p+1,1,beta,p+1,biatx)

       this%eigenvalues_S = (0.0_8, 0.0_8)
       do i=1, this%N
          do  j = -(p-1)/2,(p+1)/2
             imode = modulo((ishift+j)*(i-1),this%N)
             this%eigenvalues_S(i) = this%eigenvalues_S(i) &
                  + biatx(j+(p+1)/2)*this%modes(imode)
             !+ bspline(p,j-(p+1)/2,beta)*this%modes(imode)
          end do
          
       end do
       this%ufft = this%ufft * this%eigenvalues_S * this%eigenvalues_Minv

    case default
       print*, 'periodic_interp_module:interpolator ',this%interpolator, ' not implemented'
       stop
    end select
    
    ! modulate eigenvalues 
    !do k=1, this%N
    !   print*, k, this%eigenvalues_S(k), &
    !        0.5_8*(1+tanh(100*(.25_8*this%N-k)/this%N))!*this%eigenvalues_S(k)
       !print*, k, this%ufft(k), u(k) 
    !end do
    

    !call system_clock(COUNT=clock_end) ! Stop timing
    ! Calculate the elapsed time in seconds:
    !elapsed_time=real(clock_end-clock_start,8)/real(clock_rate,8)
    !print*,'time eigenvalues', elapsed_time

    ! Perform inverse FFT and normalized
    call zfftb(this%N, this%ufft, this%wsave)
    do i=1, this%N
       u_out(i) = real(this%ufft(i),8) / this%N
    end do

  end subroutine periodic_interp

end module periodic_interp_module
