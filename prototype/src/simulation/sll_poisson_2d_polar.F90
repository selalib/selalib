module poisson_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

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
  subroutine poisson_solve_polar(ftab,rmin,dr,Nr,Ntheta,pfwd,pinv,phitab)

    sll_real64, intent(in) :: rmin,dr
    sll_int32, intent(in) :: nr, ntheta
    sll_real64, dimension(:,:), intent(in), pointer :: ftab
    type(sll_fft_plan),intent(in), pointer ::pfwd, pinv
    sll_real64, dimension(:,:), intent(inout), pointer :: phitab

    sll_real64 :: r, ind_k
    sll_int32::i,k,err
    sll_real64, dimension(:,:), pointer :: ffttab
    ! for the tridiag solver
    sll_real64, dimension(:), pointer :: cts
    sll_real64, dimension(:), pointer :: a
    sll_int32, dimension(:), pointer :: ipiv

    !pour faire la fft avec fftpack
    sll_real64, dimension(:), pointer :: buf

    SLL_ALLOCATE(ffttab(nr+1,ntheta),err)
    SLL_ALLOCATE(a(3*(nr+1)),err)
    SLL_ALLOCATE(cts(7*(nr+1)),err)
    SLL_ALLOCATE(ipiv(nr+1),err)

    !travail avec fftpack
    SLL_ALLOCATE(buf(2*ntheta+15),err)

    ! copy of ftab
    ! we work with ffttab not to modify ftab
    ffttab=-ftab(:,1:ntheta)
!!$    ffttab(1,:)=phitab(1,1:ntheta)
!!$    ffttab(nr+1,:)=phitab(nr+1,1:ntheta)

    do i=1,nr+1
       call fft_apply_plan(pfwd,ffttab(i,:),ffttab(i,:))
    end do
!!$    !travail avec fftpack
!!$    call dffti(ntheta,buf)
!!$    do i=1,nr+1
!!$       call dfftf(ntheta,ffttab(i,:),buf)
!!$    end do
!!$    ffttab=ffttab/real(ntheta,f64)

    ! poisson solver
    do k=1,ntheta
       if (k<=ntheta) then
          ind_k=real(k-1,f64)
       else
          ind_k=real(ntheta-(k-1),f64)
       end if
!!$       ind_k=real(k-1,f64)
!!$
       if (ind_k==0.0_f64) then
          phitab(:,k)=0.0_f64
       else
!!$          ! build the matrix
          do i=1,Nr+1
             r=rmin+real(i-1,f64)*dr
             a(3*i)=1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
             a(3*i-1)=-2.0_f64/dr**2-(ind_k/r)**2
             a(3*i-2)=1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
          enddo
          a(1)=0.0_f64
          a(3*nr+3)=0
!!$          !Dirichlet conditions
!!$          a(2)=0.0_f64
!!$          a(4)=0.0_f64
!!$          a(3*nr+2)=0.0_f64
!!$          a(3*nr)=0.0_f64
!!$          !Neuman condition
!!$          a(1)=0.0_f64
!!$          a(2)=real(k-1,f64)/dr**2+1.0_f64/(2.0_f64*rmin*dr)
!!$          a(3)=0.0_f64
!!$          !Dirichlet condition
!!$          a(3*(nr+1))=0.0_f64
!!$          a(3*(nr+1)-1)=1.0_f64
!!$          a(3*(nr+1)-2)=0.0_f64

          call setup_cyclic_tridiag(a,nr+1,cts,ipiv)
          call solve_cyclic_tridiag(cts,ipiv,ffttab(:,k),nr+1,phitab(:,k))

       end if
    end do

    ! FFT INVERSE
!!$    do i=1,Nr+1
!!$       call fft_apply_plan(pinv,phitab(i,1:ntheta),phitab(i,1:ntheta))
!!$    end do
    do i=1,nr+1
!       call dfftb(ntheta,phitab(i,1:ntheta))
    end do

    phitab(:,ntheta+1)=phitab(:,1)

    !phitab=phitab*real(ntheta,f64)

    SLL_DEALLOCATE(ffttab,err)
    SLL_DEALLOCATE(cts,err)
    SLL_DEALLOCATE(ipiv,err)
    SLL_DEALLOCATE(a,err)

  end subroutine poisson_solve_polar

end module poisson_polar
