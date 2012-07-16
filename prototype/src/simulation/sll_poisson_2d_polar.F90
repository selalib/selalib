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
    sll_real64, dimension(:), pointer :: buf
    sll_comp64, dimension(:), pointer :: fk,phik
    ! for the tridiag solver
    sll_real64, dimension(:), pointer :: cts
    sll_real64, dimension(:), pointer :: a
    sll_int32, dimension(:), pointer :: ipiv

    SLL_ALLOCATE(ffttab(nr+1,ntheta),err)
    SLL_ALLOCATE(buf(2*ntheta+15),err)
    SLL_ALLOCATE(a(3*(nr+1)),err)
    SLL_ALLOCATE(cts(7*(nr+1)),err)
    SLL_ALLOCATE(ipiv(nr+1),err)
!!$    SLL_ALLOCATE(fk(nr+1),err)
!!$    SLL_ALLOCATE(phik(nr+1),err)

    ! copy of ftab
    ! we work with ffttab not to modify ftab
    ffttab=-ftab(:,1:ntheta)
    ffttab(1,:)=0.0_f64
    ffttab(nr+1,:)=0.0_f64

    call dffti(ntheta,buf)
    do i=1,nr+1
!!$       call fft_apply_plan(pfwd,ffttab(i,:),ffttab(i,:))
       call dfftf(ntheta,ffttab(i,:),buf)
    end do
    ffttab=ffttab/real(ntheta,f64)

    ! poisson solver
    do k=0,ntheta-1
       ind_k=real(floor(real(k+1,f64)/2.0_f64),f64)
       !PRINT*,"k=",ind_k
!!$    do k=0,ntheta-1
!!$       ind_k=real(k,f64)
       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          a(3*i)=1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
          a(3*i-1)=-2.0_f64/dr**2-(ind_k/r)**2
          a(3*i-2)=1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
!!$          fk(i)=fft_get_mode(pfwd,ffttab(i,:),k)
       enddo
       a(1)=0.0_f64
       a(3*nr+3)=0.0_f64
       a(2)=1.0_f64
       a(3*nr+2)=1.0_f64

!!$       call setup_cyclic_tridiag(a,nr+1,cts,ipiv)
!!$       call solve_cyclic_tridiag(cts,ipiv,fk,nr+1,phik)
       call setup_cyclic_tridiag(a,nr+1,cts,ipiv)
       call solve_cyclic_tridiag(cts,ipiv,ffttab(:,k+1),nr+1,phitab(:,k+1))

!!$       do i=1,nr+1
!!$          call fft_set_mode(pinv,phitab(i,:),phik(i),k)
!!$       end do
    end do

    ! FFT INVERSE
    do i=1,Nr+1
!!$       call fft_apply_plan(pinv,phitab(i,1:ntheta),phitab(i,1:ntheta))
       call dfftb(ntheta,phitab(i,1:ntheta),buf)
    end do

    phitab(:,ntheta+1)=phitab(:,1)

    SLL_DEALLOCATE(ffttab,err)
    SLL_DEALLOCATE(cts,err)
    SLL_DEALLOCATE(ipiv,err)
    SLL_DEALLOCATE(a,err)
    SLL_DEALLOCATE(buf,err)

  end subroutine poisson_solve_polar


  !>subroutine derivate_fft(nr,ntheta,phi,derivated)
  !>this routine is used to make derivation using the fft
  !>in the VP probleme it is used to compute derivation in direction theta
  !>nr and ntheta : number of points in direction r and theta
  !>phi : field we want to derivate in direction theta
  !>derivate : grad(phi)
  !>
  !>this routine is writen for compute_grad_field, so derivated(1,:,:) = d_r(phi)
  subroutine derivate_fft(nr,ntheta,phi,derivated)

    implicit none

    sll_int32, intent(in) :: nr, ntheta
    sll_real64, dimension(:,:), intent(in), pointer :: phi
    sll_real64, dimension(:,:,:), intent(inout), pointer :: derivated

    sll_real64,dimension(:,:), pointer :: phi_copie
    sll_real64, dimension(:), pointer :: buf
    sll_real64 :: temp
    sll_int32 :: i,j,k,err

    !SLL_ALLOCATE(phi_copie(nr+1,ntheta),err)
    SLL_ALLOCATE(buf(2*ntheta+15),err)

    !phi_copie=phi(:,1:ntheta)
    derivated(2,:,:)=phi
    call dffti(ntheta,buf)
    do i=1,nr+1
       !call dfftf(ntheta,phi_copie(i,:),buf)
       call dfftf(ntheta,derivated(2,i,1:ntheta),buf)
    end do

    do i=1,nr+1
       do j=1,ntheta/2
          temp=derivated(2,i,2*j-1)*real(j,f64)
          derivated(2,i,2*j-1)=-derivated(2,i,2*(j-1)+2)*real(j,f64)
          derivated(2,i,2*j)=temp
       end do
    end do

    do i=1,nr+1
       call dfftb(ntheta,derivated(2,i,:),buf)
    end do
    !derivated(2,:,1:ntheta)=phi_copie
    derivated(2,:,ntheta+1)=derivated(2,:,1)

  end subroutine derivate_fft

end module poisson_polar
