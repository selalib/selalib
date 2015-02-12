module periodic_interp_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_fft.h"
use sll_arbitrary_degree_splines
use sll_constants

  implicit none

  sll_real64, parameter    :: pi = 3.1415926535897932385_8
  sll_real64, parameter    :: twopi = 6.2831853071795864769_8

  integer, parameter  :: TRIGO = 0, SPLINE = 1, LAGRANGE = 2, TRIGO_FFT_SELALIB = 3
  integer, parameter   :: TRIGO_REAL = 4
  complex(8), parameter :: ii_64 = dcmplx(0.0_8, 1.0_8)

  type :: periodic_interp_work
     sll_int32          :: N ! number of cells
     sll_int32          :: interpolator ! what interpolator is used
     sll_int32          :: order  ! order of interpolation (not needed for Fourier)
     sll_real64, dimension(:), pointer :: eigenvalues_Minv ! eigenvalues of M matrix
     complex(8), dimension(:), pointer :: eigenvalues_S   ! eigenvalues of shift matrix
     sll_real64, dimension(:), pointer :: wsave  ! workspace for fft
     complex(8), dimension(:), pointer :: modes  ! Fourier modes
     complex(8), dimension(:), pointer :: ufft   ! Fourier transform of function
     sll_real64, dimension(:), pointer :: buf  ! workspace for lagrange interpolation
     sll_int32          :: sizebuf ! size of workspace for lagrange interpolation
     type(sll_fft_plan), pointer :: pinv,pfwd ! type for lagrange_fft_selalib interpolation
   end type periodic_interp_work

  interface delete
     module procedure delete_periodic_interp_work
  end interface delete


contains
  subroutine initialize_periodic_interp(this,N,interpolator,order)
    type(periodic_interp_work), pointer :: this 
    sll_int32 :: N ! number of cells
    sll_int32 :: interpolator  ! kind of interpolator
    sll_int32 :: order  ! order of method
    ! local variables
    sll_int32 :: ierr, i, j
    !sll_int32 :: k
    sll_int32 :: p ! spline degree
    sll_int32 :: icoarse  ! coarsening factor for stabilization of trigonometric interpolation
    sll_real64, dimension(order) :: biatx
    !sll_real64 :: mode
    !sll_real64 :: val
    !sll_real64, dimension(:), allocatable :: buf !for fft_selalib

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
          this%eigenvalues_Minv = 1.0_8
       else 
          icoarse = 1
          ! to be implemented
          this%eigenvalues_Minv = 1.0_8
       end if

    case (SPLINE)
       this%buf=>NULL()
       biatx = uniform_b_splines_at_x(p, 0.0_8 )
       do i=1, N
          this%modes(i-1) = exp(ii_64*twopi*(i-1)/N)
          this%eigenvalues_Minv(i) = biatx((p+1)/2)
          do j = 1,(p+1)/2
             this%eigenvalues_Minv(i) = this%eigenvalues_Minv(i) &
                  + biatx(j+(p+1)/2)*2*cos(j*twopi*(i-1)/N)
          end do
          this%eigenvalues_Minv(i) = 1.0_8 / this%eigenvalues_Minv(i)
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
       print*, 'periodic_interp_module:interpolator ',interpolator, ' not implemented'
       stop
    end select

  end subroutine initialize_periodic_interp

  subroutine delete_periodic_interp_work(this)
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

  subroutine periodic_interp(this, u_out, u, alpha) 
    ! interpolate function u given at grid points on a periodic grid
    ! at positions j-alpha (alpha is normalized to the cell size)
    type(periodic_interp_work), pointer :: this
    sll_real64, dimension(:), intent(in)   :: u  ! function to be interpolated
    sll_real64, dimension(:), intent(out)   :: u_out  ! result
    sll_real64, intent(in)     :: alpha ! displacement normalized to cell size
    ! local variables
    sll_int32 :: i, j, k, p, ishift
    !sll_int32 :: j0
    sll_int32 ::  imode,n
    sll_real64 :: beta, filter
    !sll_real64 :: mode
    sll_comp64 :: tmp,tmp2
    !complex(8) :: int_fact
    !complex(8) :: z
    ! 
    sll_real64, dimension(this%order) :: biatx
    
    SLL_ASSERT(size(u_out)>=this%N)
    SLL_ASSERT(size(u)>=this%N)
     

    ! Compute eigenvalues of shift matrix
    select case (this%interpolator)
    case (TRIGO)
       ! Perform FFT of u
       do i=1, this%N
         this%ufft(i) = dcmplx(u(i),8)
       end do
       call zfftf(this%N, this%ufft, this%wsave)
       this%eigenvalues_S(1) = 1.0_8
       this%eigenvalues_S(this%N/2+1) = exp(-ii_64*pi*alpha)
       do k=1, this%N/2-1
          !filter = 0.5_8*(1+tanh(100*(.35_8*this%N/2-k)/(this%N/2))) !F1
          !filter = 0.5_8*(1+tanh(100*(.25_8*this%N/2-k)/(this%N/2))) !F2
          !filter = 0.5_8*(1+tanh(50*(.25_8*this%N/2-k)/(this%N/2))) !F3
          filter = 1.0_8

          this%eigenvalues_S(k+1) = exp(-ii_64*twopi*k*alpha/this%N) * filter
          this%eigenvalues_S(this%N-k+1) = exp(ii_64*twopi*k*alpha/this%N) * filter
       end do
       this%ufft = this%ufft*this%eigenvalues_S
       ! Perform inverse FFT and normalized
       call zfftb(this%N, this%ufft, this%wsave)
       do i=1, this%N
          u_out(i) = real(this%ufft(i),8) / this%N
       end do
    case (SPLINE)
       ! Perform FFT of u
       do i=1, this%N
         this%ufft(i) = dcmplx(u(i),8)
       end do
       call zfftf(this%N, this%ufft, this%wsave)

       !compute eigenvalues of splines evaluated at displaced points
       p =  this%order - 1
       ishift = floor (-alpha)
       beta = -ishift - alpha
       biatx = uniform_b_splines_at_x(p, beta )

       this%eigenvalues_S = (0.0_8, 0.0_8)
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
          u_out(i) = real(this%ufft(i),8) / this%N
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
         tmp=tmp*exp(tmp2*real(0,f64))
         SET_MODE0(tmp,u_out)
       do i=1,n/2-1
         GET_MODE_LT_N_2(tmp,u_out,i,n)
         tmp=tmp*exp(tmp2*real(i,f64))
         SET_MODE_LT_N_2(tmp,u_out,i,n)
       enddo
         GET_MODE_N_2(tmp,u_out,n)
         tmp=tmp*exp(tmp2*real(n/2,f64))
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

 

  end subroutine periodic_interp

 subroutine fourier1dperlagodd(buf,sizebuf,E,N,alpha,d)
    integer,intent(in)::N,sizebuf,d
    real(8),dimension(0:sizebuf-1),intent(inout)::buf
    real(8),dimension(1:N),intent(inout)::E
    real(8),intent(in)::alpha
    integer::i,ix
    real(8)::rea,ima,reb,imb,tmp,x,a
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
    x=x-real(floor(x),8)
    !if(x==1)then
    !  x=0._rk
    !endif  
    x=x*real(N,8)
    ix=floor(x)
    if(ix==N)then
      x=0._8;ix=0
    endif
    x=x-real(ix,8)    
    do i=0,N-1
      buf(i)=0._8
    enddo

    a=1._8;
    do i=2,d
      a=a*(x*x-real(i,8)*real(i,8))/(real(d,8)*real(d,8))
    enddo
    a=a*(x+1._8)/real(d,8)
    a=a*(x-real(d,8)-1._8)/real(d,8)
    buf(ix)=a*(x-1._8)/real(d,8)
    buf(mod(ix+1,N))=a*x/real(d,8)
    a=a*x*(x-1._8)/(real(d,8)*real(d,8))  
    do i=-d,-1
      buf(mod(i+ix+N,N))=a/((x-real(i,8))/real(d,8))
    enddo  
    do i=2,d+1
      buf(mod(i+ix+N,N))=a/((x-real(i,8))/real(d,8));
    enddo
    a=1._8;
    do i=-d,d+1
      buf(mod(i+ix+N,N))=buf(mod(i+ix+N,N))*a
      a=a*real(d,8)/real(d+i+1,8)
    enddo
    a=1._8;
    do i=d+1,-d,-1
      buf(mod(i+ix+N,N))=buf(mod(i+ix+N,N))*a
      a=a*real(d,8)/real(i-1-d-1,8)
    enddo

    call dfftf(N,buf(0:N-1),buf(N:3*N+14))

    call dfftf(N,E,buf(N:3*N+14))
    tmp=1._8/real(N,8);            
    E(1)=E(1)*tmp*buf(0)
    do i=1,(N-2)/2
      rea=E(2*i);ima=E(2*i+1)
      reb=tmp*buf(2*i-1);imb=tmp*buf(2*i);
      E(2*i)=rea*reb-ima*imb
      E(2*i+1)=rea*imb+reb*ima
    enddo
    if(mod(N,2)==0)E(N)=E(N)*tmp*buf(N-1)
    call dfftb(N,E,buf(N:3*N+14))
  end subroutine fourier1dperlagodd

  subroutine fourier1dper(coefd,Ncoef,E,N,alpha)
    integer,intent(in)::N,Ncoef
    real(8),dimension(1:Ncoef),intent(in)::coefd
    real(8),dimension(1:N),intent(inout)::E
    real(8),intent(in)::alpha
    integer::i
    !integer::ix
    real(8)::rea,ima,reb,imb,tmp,x
    !localization
    !
    x=-alpha
    !x=alpha
    x=x-real(floor(x),f64)
    !x=x*real(N,8)
    !ix=floor(x)
    !x=x-real(ix,8)    
    x=x*2._f64*sll_pi

    
    
    
    call dfftf(N,E,coefd)

    !print *,E

    

    tmp=1._f64/real(N,f64);
    !print *,E

    E(1)=E(1)*tmp
    do i=1,(N-2)/2
       rea=E(2*i);ima=E(2*i+1)
       reb=tmp*cos(real(i,f64)*x);imb=tmp*sin(real(i,f64)*x);
       E(2*i)=rea*reb-ima*imb
       E(2*i+1)=rea*imb+reb*ima
    enddo
    if(mod(N,2)==0)E(N)=E(N)*tmp*cos(0.5_f64*real(N,f64)*x)
    call dfftb(N,E,coefd)
  end subroutine fourier1dper
end module periodic_interp_module
