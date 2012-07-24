module poisson_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use polar_kind
  use sll_fft
!!$  use fftpack_module
  use sll_tridiagonal
  use numeric_constants
  implicit none

contains

  !>subroutine poisson_solve_polar(ftab,rmin,dr,Nr,Ntheta,pfwd,pinv,phitab)
  !>poisson solver for polar system
  !>ftab : distribution function, size (nr+1)X(ntheta+1)
  !>phitab : solution of laplacien phi = -f
  !>phitab(1,:) and phitab(nr+1,:) must be known as boudary condition
  !>rmin : radius of the hole
  !>dr : size of r step
  !>Nr and Ntheta : number Step in directions r and theta
  !>pfwd and pinv : initialization object for FFt forward and inverse
  !>initialization must be done outside the solver
  subroutine poisson_solve_polar(adv)

    type(polar_vp_data), intent(inout), pointer :: adv

    sll_real64 :: rmin,dr
    sll_int32 :: nr, ntheta

    sll_real64 :: r, ind_k
    sll_int32::i,k,err
    sll_real64, dimension(:), pointer :: buf

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    rmin=adv%data%rmin
    dr=adv%data%dr

    SLL_ALLOCATE(buf(2*ntheta+15),err)

    adv%f_fft=adv%f

    call dffti(ntheta,buf)
    do i=1,nr+1
!!$       call fft_apply_plan(pfwd,ffttab(i,:),ffttab(i,:))
       call dfftf(ntheta,adv%f_fft(i,1:ntheta),buf)
    end do
    adv%f_fft=adv%f_fft/real(ntheta,f64)

    ! poisson solver
    do k=0,ntheta-1
       ind_k=real(floor(real(k+1,f64)/2.0_f64),f64)
       !PRINT*,"k=",ind_k
!!$    do k=0,ntheta-1
!!$       ind_k=real(k,f64)
       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          adv%a(3*i)=1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
          adv%a(3*i-1)=-2.0_f64/dr**2-(ind_k/r)**2
          adv%a(3*i-2)=1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
!!$          fk(i)=fft_get_mode(pfwd,ffttab(i,:),k)
       enddo
       adv%a(1)=0.0_f64
       adv%a(3*nr+3)=0.0_f64
       adv%a(2)=1.0_f64
       adv%a(3*nr+2)=1.0_f64

       call setup_cyclic_tridiag(adv%a,nr+1,adv%cts,adv%ipiv)
       call solve_cyclic_tridiag(adv%cts,adv%ipiv,adv%f_fft(:,k+1),nr+1,adv%phi(:,k+1))

!!$       do i=1,nr+1
!!$          call fft_set_mode(pinv,phitab(i,:),phik(i),k)
!!$       end do
    end do

    ! FFT INVERSE
    do i=1,Nr+1
!!$       call fft_apply_plan(pinv,phitab(i,1:ntheta),phitab(i,1:ntheta))
       call dfftb(ntheta,adv%phi(i,1:ntheta),buf)
    end do

    adv%phi(:,ntheta+1)=adv%phi(:,1)

  end subroutine poisson_solve_polar


  !>subroutine derivate_fft(nr,ntheta,phi,derivated)
  !>this routine is used to make derivation using the fft
  !>in the VP probleme it is used to compute derivation in direction theta
  !>nr and ntheta : number of points in direction r and theta
  !>phi : field we want to derivate in direction theta
  !>derivate : grad(phi)
  !>
  !>this routine is writen for compute_grad_field, so derivated(1,:,:) = d_r(phi)
  subroutine derivate_fft(adv)

    implicit none

    type(polar_vp_data), intent(in), pointer :: adv

    sll_int32 :: nr, ntheta
    sll_real64, dimension(:), pointer :: buf
    sll_real64 :: temp
    sll_int32 :: i,j,k,err

    nr=adv%data%nr
    ntheta=adv%data%ntheta

    !SLL_ALLOCATE(phi_copie(nr+1,ntheta),err)
    SLL_ALLOCATE(buf(2*ntheta+15),err)

    !!!!can't be done with fftpack r2r

  end subroutine derivate_fft

end module poisson_polar
