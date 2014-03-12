module poisson2dpp_seq

use used_precision
use geometry_module
use fft_module
#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_constants


implicit none
public :: new, dealloc, solve, transposexy, transposeyx

type, public:: poisson2dpp
   complex(wp), dimension(:,:), pointer :: rhot, ext, eyt
   type(fftclass) :: fftx, ffty
   type(geometry) :: geomx
   logical :: transposed
   sll_int32 :: jstartx,jendx,istartk,iendk
end type poisson2dpp

interface new
   module procedure new_poisson2dpp
end interface
interface dealloc
   module procedure dealloc_poisson2dpp
end interface
interface solve
   module procedure solve_poisson2dpp
end interface

contains

subroutine new_poisson2dpp(this,rho,geomx,iflag)

   type(poisson2dpp),intent(out) :: this
   sll_real64, dimension(:,:), intent(in) :: rho
   type(geometry),intent(in)  :: geomx
   sll_int32, intent(out) :: iflag

   sll_int32 :: ierr ! indicateur d'erreur
   sll_int32 :: nxh1

   ! indicateur d'erreur
   iflag = 0
   ! on commence par utiliser la fonction rho(x,y)
   this%transposed=.false.

   this%jstartx = 1
   this%jendx = geomx%ny

   ! la taille totale de la zone en kx est nxh1
   ! ATTENTION : les tableaux concernees par ce decoupage sont rhot,
   ! ext et eyt. Ce sont des complexes. Pour cette raison leur taille est
   ! la moitie de celle des tableaux reels correspondants.
   nxh1 = geomx%nx/2
   this%istartk = 1  ! cas sequentiel
   this%iendk = nxh1+1 ! cas sequentiel

   ! initialisation de la geometrie
   this%geomx=geomx
   ! allocation memoire
   SLL_ALLOCATE(this%rhot(this%geomx%ny,this%istartk:this%iendk),ierr)
   SLL_ALLOCATE(this%ext(this%geomx%ny,this%istartk:this%iendk),ierr)
   SLL_ALLOCATE(this%eyt(this%geomx%ny,this%istartk:this%iendk),ierr)
   ! initialisation des fft (le nombre de points pris en compte est n)
   call initfft(this%fftx, rho, this%geomx%nx)
   call initfft(this%ffty, this%rhot, this%geomx%ny)

end subroutine new_poisson2dpp

subroutine dealloc_poisson2dpp(this)

   type(poisson2dpp),intent(out) :: this
   sll_int32 :: ierr

   SLL_DEALLOCATE(this%rhot,ierr)
   SLL_DEALLOCATE(this%ext,ierr)
   SLL_DEALLOCATE(this%eyt,ierr)

end subroutine dealloc_poisson2dpp

subroutine solve_poisson2dpp(this,ex,ey,rho,nrj)

   type(poisson2dpp) :: this
   sll_real64, dimension(:,:) :: ex,ey,rho
   sll_real64, optional :: nrj
   sll_int32 :: ik,jk,i,j ! indices de boucle
   sll_int32 :: istart
   sll_real64 :: kx0,ky0, kx,ky, kx2, k2
   sll_real64, dimension(size(rho,1),size(rho,2)) :: phi

    ! faire une FFT dans la direction x de rho
   call fft(this%fftx,rho(:,this%jstartx:this%jendx))

   ! transposition de rho stockee dans rhot
   call transposexy(this,rho)

   ! faire une FFT dans la direction y de rhot
   call fft(this%ffty,this%rhot(:,this%istartk:this%iendk))

   ! calcul de la transformee de Fourier de E a partir de celle de rho
   kx0=2*sll_pi/(this%geomx%nx*this%geomx%dx)
   ky0=2*sll_pi/(this%geomx%ny*this%geomx%dy)

   ! pour jk=1
   jk = 1    
   this%ext(1,1)=0._8
   this%eyt(1,1)=0._8

   istart = 2

   do ik=istart,this%iendk
      kx= (ik-1)*kx0
      k2 = kx*kx
      this%ext(jk,ik)=-dcmplx(zero,kx/k2)*this%rhot(jk,ik)
      this%eyt(jk,ik)= 0.
      this%rhot(jk,ik)= this%rhot(jk,ik)/k2
   end do

    ! pour jk > 1
   do ik=this%istartk,this%iendk
      kx= (ik-1)*kx0
      kx2 = kx*kx
      do jk =  2, this%geomx%ny/2+1
         ky= (jk-1)*ky0
         k2= kx2 +ky*ky

         this%ext(jk,ik)=-dcmplx(zero,kx/k2)*this%rhot(jk,ik)
         this%eyt(jk,ik)=-dcmplx(zero,ky/k2)*this%rhot(jk,ik)
         this%rhot(jk,ik)= this%rhot(jk,ik)/k2
      end do
      do jk = this%geomx%ny/2+2 , this%geomx%ny     
         ky= (jk-1-this%geomx%ny)*ky0
         k2= kx2 +ky*ky

         this%ext(jk,ik)=-dcmplx(zero,kx/k2)*this%rhot(jk,ik)
         this%eyt(jk,ik)=-dcmplx(zero,ky/k2)*this%rhot(jk,ik)
         this%rhot(jk,ik)= this%rhot(jk,ik)/k2
      end do
   end do

    call fftinv(this%ffty,this%ext(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%eyt(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%rhot(:,this%istartk:this%iendk))

    ! transposition de E
    !call transposeyx(this,ex,ey)
    call transposeyx(this,ex,ey,phi) 

    call fftinv(this%fftx,ex(:,this%jstartx:this%jendx))
    call fftinv(this%fftx,ey(:,this%jstartx:this%jendx))
    call fftinv(this%fftx,phi(:,this%jstartx:this%jendx))

    rho=phi

    if (present(nrj)) then 
       nrj=0._wp
       do i=1,this%geomx%nx
          do j=1,this%geomx%ny
             nrj=nrj+ex(i,j)*ex(i,j)+ey(i,j)*ey(i,j)          
          enddo
       enddo
   
       nrj=nrj*this%geomx%dx*this%geomx%dy
       if (nrj>1.e-30) then 
          nrj=0.5_wp*log(nrj)
       else
          nrj=-10**9
       endif
   end if

end subroutine solve_poisson2dpp

subroutine transposexy(this,rho)

   type(poisson2dpp) :: this
   sll_real64, dimension(:,:) :: rho
   sll_int32 :: i,j ! indices de boucle
   sll_int32 :: nxh1

   nxh1 = this%geomx%nx/2
   do i=max(this%istartk,2),this%iendk-1
      do j=1,this%geomx%ny
         this%rhot(j,i)=dcmplx(rho(2*i-2,j),rho(2*i-1,j))
      end do
   end do
   if (this%istartk.eq.1) then
      do j=1,this%geomx%ny
         this%rhot(j,1) = dcmplx(rho(1,j),0.0)
      end do
   end if
   if (this%iendk.eq.nxh1+1) then
      do j=1,this%geomx%ny
         this%rhot(j,nxh1+1) = dcmplx(rho(2*nxh1,j),0.0)
      end do
   end if
    
end subroutine transposexy

subroutine transposeyx(this,ex,ey,rho)

   type(poisson2dpp) :: this
   sll_real64, dimension(:,:) :: ex,ey,rho
   sll_int32 :: i,j ! indices de boucle
   sll_int32 :: nxh1

   ! transposition
   nxh1 = this%geomx%nx/2
   do j=1,this%geomx%ny
      if (this%istartk.eq.1) then
         ex(1,j)= real (this%ext(j,1))
         ey(1,j)= real (this%eyt(j,1))
         rho(1,j)= real (this%rhot(j,1))
      end if
      if (this%iendk.eq.nxh1+1) then
         ex(this%geomx%nx,j)= real (this%ext(j,nxh1+1))
         ey(this%geomx%nx,j)= real (this%eyt(j,nxh1+1))
         rho(this%geomx%nx,j)= real (this%rhot(j,nxh1+1))
      end if
      do i=max(this%istartk,2),this%iendk-1
         ex(2*i-2,j)=real(this%ext(j,i))
         ex(2*i-1,j)=imag(this%ext(j,i))
         ey(2*i-2,j)=real(this%eyt(j,i))
         ey(2*i-1,j)=imag(this%eyt(j,i))
         rho(2*i-2,j)=real(this%rhot(j,i))
         rho(2*i-1,j)=imag(this%rhot(j,i))
      end do
   end do

end subroutine transposeyx

end module poisson2dpp_seq
