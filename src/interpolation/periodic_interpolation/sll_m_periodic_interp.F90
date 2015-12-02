module sll_m_periodic_interp
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_fft.h"
use sll_m_arbitrary_degree_splines
use sll_m_constants

  implicit none

  sll_int32,  parameter :: TRIGO = 0, SPLINE = 1, LAGRANGE = 2, TRIGO_FFT_SELALIB = 3
  sll_int32,  parameter :: TRIGO_REAL = 4
  sll_comp64, parameter :: ii_64 = (0.0_f64, 1.0_f64)

  type :: periodic_interp_work
     sll_int32           :: N ! number of cells
     sll_int32           :: interpolator ! what interpolator is used
     sll_int32           :: order  ! order of interpolation (not needed for Fourier)
     sll_real64, pointer :: eigenvalues_Minv(:) ! eigenvalues of M matrix
     sll_comp64, pointer :: eigenvalues_S(:)   ! eigenvalues of shift matrix
     sll_real64, pointer :: wsave(:) ! workspace for fft
     sll_comp64, pointer :: modes(:) ! Fourier modes
     sll_comp64, pointer :: ufft (:) ! Fourier transform of function
     sll_real64, pointer :: buf  (:) ! workspace for lagrange interpolation
     sll_int32           :: sizebuf ! size of workspace for lagrange interpolation
     type(sll_fft_plan), pointer :: pinv, pfwd ! type for lagrange_fft_selalib interpolation
   end type periodic_interp_work

  interface delete
     module procedure delete_periodic_interp_work
  end interface delete


contains
  subroutine initialize_periodic_interp( this, N, interpolator, order )
    type(periodic_interp_work), pointer  :: this
    sll_int32,                intent(in) :: N            ! number of cells
    sll_int32,                intent(in) :: interpolator ! interpolation method
    sll_int32,                intent(in) :: order        ! order of method

    ! local variables
    sll_int32  :: ierr, i, j
    !sll_int32  :: k
    sll_int32  :: p ! spline degree
    sll_int32  :: icoarse  ! coarsening factor for stabilization of trigonometric interpolation
    sll_real64 :: biatx(order)
    !sll_real64 :: mode
    !sll_real64 :: val
    !sll_real64, allocatable :: buf(:) !for fft_selalib

    SLL_ALLOCATE( this, ierr )
    this%N = N
    this%interpolator = interpolator
    this%order = order 

    ! Allocate arrays
    SLL_ALLOCATE(this%wsave(15+4*N),ierr)
    SLL_ALLOCATE(this%ufft(N),ierr)
    SLL_ALLOCATE(this%eigenvalues_Minv(N),ierr)
    SLL_ALLOCATE(this%eigenvalues_S(N),ierr)
    SLL_ALLOCATE(this%modes(0:N-1),ierr)

    ! set up spline parameters
    if ((order/2) /= int(order/2.0)) then
       print*, 'initialize_periodic_interp: order of interpolators needs to be even.', &
            'Order here is: ', Order
       stop
    end if
    p = order - 1  ! spline degree is one less than order

    ! Initialize fft
    call zffti(N, this%wsave)

    select case (interpolator)
    case (TRIGO)
       this%buf=>NULL()
       if (order == 0) then ! no coarsening
          this%eigenvalues_Minv = 1.0_f64
       else 
          icoarse = 1
          ! to be implemented
          this%eigenvalues_Minv = 1.0_f64
       end if

    case (SPLINE)
       this%buf=>NULL()
       biatx = uniform_b_splines_at_x(p, 0.0_f64 )
       do i=1, N
          this%modes(i-1) = exp(ii_64*sll_twopi*(i-1)/N)
          this%eigenvalues_Minv(i) = biatx((p+1)/2)
          do j = 1,(p+1)/2
             this%eigenvalues_Minv(i) = this%eigenvalues_Minv(i) &
                  + biatx(j+(p+1)/2)*2*cos(j*sll_twopi*(i-1)/N)
          end do
          this%eigenvalues_Minv(i) = 1.0_f64 / this%eigenvalues_Minv(i)
       end do
    case (LAGRANGE)
      this%sizebuf=3*N+15
      SLL_ALLOCATE(this%buf(0:this%sizebuf-1),ierr)    
       call dffti(N,this%buf(N:3*N+14))
    case (TRIGO_REAL)
      this%sizebuf=2*N+15
      SLL_ALLOCATE(this%buf(1:this%sizebuf),ierr)    
       call dffti(N,this%buf)
    case (TRIGO_FFT_SELALIB)
       this%sizebuf=N
       SLL_ALLOCATE(this%buf(this%sizebuf),ierr)          
       !SLL_ALLOCATE(buf(N),ierr)
       this%pfwd => fft_new_plan(N,this%buf,this%buf,FFT_FORWARD,FFT_NORMALIZE)
       this%pinv => fft_new_plan(N,this%buf,this%buf,FFT_INVERSE)
       SLL_DEALLOCATE_ARRAY(this%buf,ierr)       
    case default
       print*, 'sll_m_periodic_interp:interpolator ',interpolator, ' not implemented'
       stop
    end select

  end subroutine initialize_periodic_interp

  subroutine delete_periodic_interp_work( this )
    type(periodic_interp_work), pointer :: this 

    sll_int32 :: ierr
    
    SLL_DEALLOCATE(this%wsave, ierr)
    SLL_DEALLOCATE(this%ufft,ierr)
    SLL_DEALLOCATE(this%eigenvalues_Minv,ierr)
    SLL_DEALLOCATE(this%eigenvalues_S,ierr)
    SLL_DEALLOCATE(this%modes,ierr)
    if(associated(this%buf))then
      SLL_DEALLOCATE(this%buf,ierr)
    endif
    SLL_DEALLOCATE(this, ierr )
  end subroutine delete_periodic_interp_work

  subroutine periodic_interp( this, u_out, u, alpha )
    ! interpolate function u given at grid points on a periodic grid
    ! at positions j-alpha (alpha is normalized to the cell size)
    type(periodic_interp_work), pointer :: this
    sll_real64,             intent(out) :: u_out(:) ! result
    sll_real64,             intent(in)  :: u(:)     ! function to be interpolated
    sll_real64,             intent(in)  :: alpha    ! displacement normalized to cell size

    ! local variables
    sll_int32  :: i, j, k, p, ishift
    !sll_int32  :: j0
    sll_int32  ::  imode,n
    sll_real64 :: beta
    sll_comp64 :: filter
    !sll_real64 :: mode
    sll_comp64 :: tmp,tmp2
    !sll_comp64 :: int_fact
    !sll_comp64 :: z
    ! 
    sll_real64 :: biatx(this%order)

    SLL_ASSERT(size(u_out)>=this%N)
    SLL_ASSERT(size(u)>=this%N)

    ! Compute eigenvalues of shift matrix
    select case (this%interpolator)
    case (TRIGO)
       ! Perform FFT of u
       do i=1, this%N
         this%ufft(i) = cmplx(u(i),kind=f64)
       end do
       call zfftf(this%N, this%ufft, this%wsave)
       this%eigenvalues_S(1) = (1.0_f64, 0.0_f64)
       this%eigenvalues_S(this%N/2+1) = exp(-ii_64*sll_pi*alpha)
       do k=1, this%N/2-1
          !filter = 0.5_8*(1+tanh(100*(.35_8*this%N/2-k)/(this%N/2))) !F1
          !filter = 0.5_8*(1+tanh(100*(.25_8*this%N/2-k)/(this%N/2))) !F2
          !filter = 0.5_8*(1+tanh(50*(.25_8*this%N/2-k)/(this%N/2))) !F3
          filter = (1.0_f64, 0.0_f64)

          this%eigenvalues_S(k+1) = exp(-ii_64*sll_twopi*k*alpha/this%N) * filter
          this%eigenvalues_S(this%N-k+1) = exp(ii_64*sll_twopi*k*alpha/this%N) * filter
       end do
       this%ufft = this%ufft*this%eigenvalues_S
       ! Perform inverse FFT and normalized
       call zfftb(this%N, this%ufft, this%wsave)
       do i=1, this%N
          u_out(i) = real(this%ufft(i),f64) / real(this%N,f64)
       end do
    case (SPLINE)
       ! Perform FFT of u
       do i=1, this%N
         this%ufft(i) = cmplx(u(i),kind=f64)
       end do
       call zfftf(this%N, this%ufft, this%wsave)

       !compute eigenvalues of splines evaluated at displaced points
       p =  this%order - 1
       ishift = floor (-alpha)
       beta = -ishift - alpha
       biatx = uniform_b_splines_at_x(p, beta )

       this%eigenvalues_S = (0.0_f64, 0.0_f64)
       do i=1, this%N
          do  j = -(p-1)/2,(p+1)/2
             imode = modulo((ishift+j)*(i-1),this%N)
             this%eigenvalues_S(i) = this%eigenvalues_S(i) &
                  + biatx(j+(p+1)/2)*this%modes(imode)
           end do
          
       end do
       this%ufft = this%ufft * this%eigenvalues_S * this%eigenvalues_Minv
       ! Perform inverse FFT and normalized
       call zfftb(this%N, this%ufft, this%wsave)
       do i=1, this%N
          u_out(i) = real(this%ufft(i),f64) / real(this%N,f64)
       end do
    case (LAGRANGE)
       u_out = u
       call fourier1dperlagodd(this%buf,this%sizebuf,u_out,this%N, &
            alpha/this%N,this%order/2 - 1 )
    case (TRIGO_REAL)
       u_out = u
       !call fourier1dperlagodd(this%buf,this%sizebuf,u_out,this%N, &
       !     alpha/this%N,this%order/2 - 1 )
       !print *,this%sizebuf,this%N
       call fourier1dper(this%buf,this%sizebuf,u_out,this%N,alpha/this%N)
            
    case (TRIGO_FFT_SELALIB)
       u_out = u
       call fft_apply_plan(this%pfwd,u_out,u_out)
       n=this%N
       tmp2=-ii_64*2._f64*sll_pi/n*alpha

         GET_MODE0(tmp,u_out)
         !tmp=tmp*exp(tmp2*real(0,f64))
         SET_MODE0(tmp,u_out)
       do i=1,n/2-1
         GET_MODE_LT_N_2(tmp,u_out,i,n)
         tmp=tmp*exp( tmp2*cmplx(i,0.0_f64,f64) )
         SET_MODE_LT_N_2(tmp,u_out,i,n)
       enddo
         GET_MODE_N_2(tmp,u_out,n)
         tmp=tmp*exp( tmp2*cmplx(0.5_f64*n,0.0_f64,f64) )
         SET_MODE_N_2(tmp,u_out,n)

!*** Without macro
      ! do i=0,this%N/2
      !   tmp=fft_get_mode(this%pfwd,u_out,i)
      !   !print *,i,tmp,alpha
      !   !tmp=tmp*exp(-ii*2._f64*sll_pi/this%N*alpha*real(i,f64))
      !   tmp=tmp*exp(tmp2*real(i,f64))
      !   !print *,i,tmp,exp(-ii*2._f64*sll_pi/this%N*alpha*real(i,f64))
      !   call fft_set_mode(this%pfwd,u_out,tmp,i)
      !   !tmp=fft_get_mode(this%pfwd,u_out,i)
      !   !print *,i,tmp
      ! enddo

       call fft_apply_plan(this%pinv,u_out,u_out)        
    case default
       print*, 'sll_m_periodic_interp:interpolator ',this%interpolator, ' not implemented'
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

  end subroutine periodic_interp

 subroutine fourier1dperlagodd( buf, sizebuf, E, N, alpha, d )
    sll_int32,  intent(in)    :: N, sizebuf, d
    sll_real64, intent(inout) :: buf(0:sizebuf-1)
    sll_real64, intent(inout) :: E(1:N)
    sll_real64, intent(in)    :: alpha

    sll_int32  :: i, ix
    sll_real64 :: rea, ima, reb, imb, tmp, x, a

    !localization
    x=alpha
    !if(abs(x)<1.e-15_rk)then
    !  print *,x,N
    !  x=x-real(floor(x),rk)
    !  print *,x
    !  x=x*real(N,rk)
    !  print *,x
    !  ix=floor(x)
    !  x=x-real(ix,rk)    
    !  print *,x,ix
    !  x=0._rk
    !endif
    x=x-real(floor(x),f64)
    !if(x==1)then
    !  x=0._rk
    !endif  
    x=x*real(N,f64)
    ix=floor(x)
    if(ix==N)then
      x=0._f64;ix=0
    endif
    x=x-real(ix,f64)
    do i=0,N-1
      buf(i)=0._f64
    enddo

    a=1._8;
    do i=2,d
      a=a*(x*x-real(i,f64)*real(i,f64))/(real(d,f64)*real(d,f64))
    enddo
    a=a*(x+1._f64)/real(d,f64)
    a=a*(x-real(d,f64)-1._f64)/real(d,f64)
    buf(ix)=a*(x-1._f64)/real(d,f64)
    buf(mod(ix+1,N))=a*x/real(d,f64)
    a=a*x*(x-1._f64)/(real(d,f64)*real(d,f64))
    do i=-d,-1
      buf(mod(i+ix+N,N))=a/((x-real(i,f64))/real(d,f64))
    enddo  
    do i=2,d+1
      buf(mod(i+ix+N,N))=a/((x-real(i,f64))/real(d,f64));
    enddo
    a=1._f64
    do i=-d,d+1 ! If Y is not present then the imaginary component is set to 0.
      buf(mod(i+ix+N,N))=buf(mod(i+ix+N,N))*a
      a=a*real(d,f64)/real(d+i+1,f64)
    enddo
    a=1._f64;
    do i=d+1,-d,-1
      buf(mod(i+ix+N,N))=buf(mod(i+ix+N,N))*a
      a=a*real(d,f64)/real(i-1-d-1,f64)
    enddo

    call dfftf(N,buf(0:N-1),buf(N:3*N+14))

    call dfftf(N,E,buf(N:3*N+14))
    tmp=1._f64/real(N,f64)
    E(1)=E(1)*tmp*buf(0)
    do i=1,(N-2)/2
      rea=E(2*i);ima=E(2*i+1)
      reb=tmp*buf(2*i-1)
      imb=tmp*buf(2*i)
      E(2*i)=rea*reb-ima*imb
      E(2*i+1)=rea*imb+reb*ima
    enddo
    if(mod(N,2)==0) E(N)=E(N)*tmp*buf(N-1)
    call dfftb(N,E,buf(N:3*N+14))
  end subroutine fourier1dperlagodd

  subroutine fourier1dper( coefd, Ncoef, E, N, alpha )
    sll_int32,  intent(in)    :: N, Ncoef
    sll_real64, intent(in)    :: coefd(1:Ncoef)
    sll_real64, intent(inout) :: E(1:N)
    sll_real64, intent(in)    :: alpha

    sll_int32  :: i
    !sll_int32  :: ix
    sll_real64 :: rea, ima, reb, imb, tmp, x
    !localization
    !
    x=-alpha
    !x=alpha
    x=x-real(floor(x),f64)
    !x=x*real(N,f64)
    !ix=floor(x)
    !x=x-real(ix,f64)
    x=x*2._f64*sll_pi

    call dfftf(N,E,coefd)
    !print *,E

    tmp=1._f64/real(N,f64)
    !print *,E

    E(1)=E(1)*tmp
    do i=1,(N-2)/2
       rea=E(2*i);ima=E(2*i+1)
       reb=tmp*cos(real(i,f64)*x)
       imb=tmp*sin(real(i,f64)*x)
       E(2*i)=rea*reb-ima*imb
       E(2*i+1)=rea*imb+reb*ima
    enddo
    if(mod(N,2)==0) E(N)=E(N)*tmp*cos(0.5_f64*real(N,f64)*x)
    call dfftb(N,E,coefd)
  end subroutine fourier1dper

end module sll_m_periodic_interp
