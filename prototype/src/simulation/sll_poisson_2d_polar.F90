module poisson_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use polar_kind
  use sll_fft
  use sll_tridiagonal
  use numeric_constants
  implicit none

contains

  !>subroutine poisson_solve_polar(adv)
  !>poisson solver for polar system : Laplacian(phi)=-f
  !>adv : polar_vp_data object, all datas and needed objet are inside
  !>initialization must be done outside the solver
  subroutine poisson_solve_polar(adv)

    type(polar_vp_data), intent(inout), pointer :: adv

    sll_real64 :: rmin,dr
    sll_int32 :: nr, ntheta

    sll_real64 :: r
    sll_int32::i,k,err, ind_k
    !sll_real64, dimension(:), pointer :: buf

    nr=adv%data%nr
    ntheta=adv%data%ntheta
    rmin=adv%data%rmin
    dr=adv%data%dr

!!$    SLL_ALLOCATE(buf(2*ntheta+15),err)

    adv%f_fft=adv%f

!!$    call dffti(ntheta,buf)
    do i=1,nr+1
       call fft_apply_plan(adv%pfwd,adv%f_fft(i,:),adv%f_fft(i,:))
!!$       call dfftf(ntheta,adv%f_fft(i,1:ntheta),buf)
    end do
    !adv%f_fft=adv%f_fft/real(ntheta,f64)

   ! poisson solver
    do ind_k=0,ntheta/2
       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          adv%a(3*i)=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
          adv%a(3*i-1)=2.0_f64/dr**2+(real(ind_k,f64)/r)**2
          adv%a(3*i-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
          adv%fk(i)=fft_get_mode(adv%pfwd,adv%f_fft(i,1:ntheta),ind_k)
       enddo
       adv%a(1)=0.0_f64
       adv%a(3*nr+3)=0.0_f64
       !adv%a(2)=1.0_f64
       !adv%a(3*nr+2)=1.0_f64

       call setup_cyclic_tridiag(adv%a,nr+1,adv%cts,adv%ipiv)
!!$       call solve_cyclic_tridiag(adv%cts,adv%ipiv,adv%f_fft(:,k+1),nr+1,adv%phi(:,k+1))
       call solve_cyclic_tridiag(adv%cts,adv%ipiv,adv%fk,nr+1,adv%phik)

       do i=1,nr+1
          call fft_set_mode(adv%pinv,adv%phi(i,1:ntheta),adv%phik(i),ind_k)
       end do
    end do

    ! FFT INVERSE
    do i=1,Nr+1
       call fft_apply_plan(adv%pinv,adv%phi(i,1:ntheta),adv%phi(i,1:ntheta))
!!$       call dfftb(ntheta,adv%phi(i,1:ntheta),buf)
    end do

    adv%phi(1,:)=0.0_f64
    adv%phi(:,ntheta+1)=adv%phi(:,1)

  end subroutine poisson_solve_polar


!!$  !>subroutine derivate_fft(adv)
!!$  !>this routine is used to make derivation using the fft
!!$  !>in the CG probleme it is used to compute derivation in direction theta
!!$  !>adv : polar_vp_data object
!!$  !>WARNING :  if your program works using fftpack real to real fft, then this routine can not be used
!!$  subroutine derivate_fft(adv)
!!$
!!$    implicit none
!!$
!!$    type(polar_vp_data), intent(in), pointer :: adv
!!$
!!$    sll_int32 :: nr, ntheta
!!$    sll_real64, dimension(:), pointer :: buf
!!$    sll_real64 :: temp
!!$    sll_int32 :: i,j,k,err
!!$
!!$    nr=adv%data%nr
!!$    ntheta=adv%data%ntheta
!!$
!!$    !SLL_ALLOCATE(phi_copie(nr+1,ntheta),err)
!!$    SLL_ALLOCATE(buf(2*ntheta+15),err)
!!$
!!$    !!!!can't be done with fftpack r2r
!!$
!!$  end subroutine derivate_fft

end module poisson_polar
