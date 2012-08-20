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

  !>subroutine poisson_solve_polar(data,f,phi,f_fft,fk,phik,a,cts,ipiv,pfwd,pinv)
  !>poisson solver for polar system : -Laplacian(phi)=f
  !>data : polar_data object, contains data about the domain
  !>f : distribution function, size (nr+1)*(ntheta+1)
  !>phi : unknown field, size (nr+1)*(ntheta+1)
  !>f_fft : copy of f not to overwrite f
  !>fk and phik : complex vectors, size ntheta/2+1
  !>a : matrix for Laplacian in polar coordinates
  !>cts and ipiv : vectors for the tridiagonal solver, see the section about the tridiagonal solver
  !>pfwd and pinv : sll_fft_plan object for fft forward and inverse
  !>initialization must be done outside the solver
  subroutine poisson_solve_polar(data,f,phi,f_fft,fk,phik,a,cts,ipiv,pfwd,pinv)

    type(polar_data), intent(in), pointer :: data
    type(sll_fft_plan), intent(in), pointer :: pfwd, pinv
    sll_real64, dimension(data%nr+1,data%ntheta+1), intent(inout) :: f,phi,f_fft
    sll_real64, dimension(:), intent(inout) :: cts,a
    sll_int32, dimension(:), intent(inout) :: ipiv
    sll_comp64, dimension(:), intent(inout) :: fk,phik

    sll_real64 :: rmin,dr
    sll_int32 :: nr, ntheta

    sll_real64 :: r
    sll_int32::i,k,err, ind_k

    nr=data%nr
    ntheta=data%ntheta
    rmin=data%rmin
    dr=data%dr

    f_fft=f

    do i=1,nr+1
       call fft_apply_plan(pfwd,f_fft(i,1:ntheta),f_fft(i,1:ntheta))
    end do

   ! poisson solver
    do ind_k=0,ntheta/2
       do i=1,nr+1
          r=rmin+real(i-1,f64)*dr
          a(3*i)=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
          a(3*i-1)=2.0_f64/dr**2+(real(ind_k,f64)/r)**2
          a(3*i-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
          fk(i)=fft_get_mode(pfwd,f_fft(i,1:ntheta),ind_k)
       enddo
       a(1)=0.0_f64
       a(3*nr+3)=0.0_f64
       !a(2)=1.0_f64
       !a(3*nr+2)=1.0_f64

       call setup_cyclic_tridiag(a,nr+1,cts,ipiv)
       call solve_cyclic_tridiag(cts,ipiv,fk,nr+1,phik)

       do i=1,nr+1
          call fft_set_mode(pinv,phi(i,1:ntheta),phik(i),ind_k)
       end do
    end do

    ! FFT INVERSE
    do i=1,Nr+1
       call fft_apply_plan(pinv,phi(i,1:ntheta),phi(i,1:ntheta))
    end do

    phi(1,:)=0.0_f64
    phi(nr+1,:)=0.0_f64
    phi(:,ntheta+1)=phi(:,1)

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
