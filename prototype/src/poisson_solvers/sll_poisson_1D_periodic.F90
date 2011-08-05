module sll_poisson_1D_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
  use numeric_constants
  use geometry1d_module
  use fft1d_module

  implicit none
  private
  public :: new, solve

  type, public                               :: poisson1dp
     sll_real64, dimension(:), pointer       :: rho, ex
     type(fft1dclass)                        :: fftx
     sll_int32                               :: ncx
  end type poisson1dp

  interface new
     module procedure new_poisson1dp
  end interface

  interface solve
     module procedure solve_poisson1dp !, solve_poisson1dp_axisymetrique
  end interface

contains
  subroutine new_poisson1dp(this,ncx,iflag)
    type(poisson1dp),intent(out)             :: this
    sll_int32,intent(in)                     :: ncx
    sll_int32, intent(out)                   :: iflag 

    ! Indicateur d'erreur
    iflag = 0
    ! Initialisation de la geometrie
    this%ncx=ncx
    ! Initialisation des fft (le nombre de points pris en compte est n
    call initdoubfft(this%fftx, this%ncx)
    ! Allocate auxiliary arrays for fft
    SLL_ALLOCATE(this%rho(ncx+2),iflag)
    SLL_ALLOCATE(this%ex(ncx+2),iflag)

  end subroutine new_poisson1dp

  subroutine solve_poisson1dp(this,ex,rho)
    type(poisson1dp)                          :: this
    type(field_1D_vec1), pointer              :: ex,rho
    sll_int32                                 :: i, ik
    sll_int32                                 :: nxh1
    sll_real64                                :: mesh_length, kx0, kx, k2
    !sll_real64, dimension(this%geomx%nx+1)   :: rhoaux

    ! copy rho into auxiliary array for fftpack
    this%rho(1:this%ncx+1) = FIELD_DATA(rho) 

    ! Faire une FFT dans la direction x de rho
    call fft(this%fftx, this%rho)

    ! Calcul de la transformee de Fourier de E a partir de celle de rho
    mesh_length = GET_FIELD_NCELLS_X1( rho ) * GET_FIELD_DELTA_X1( rho )
    kx0=2*sll_pi/(mesh_length)
    nxh1 = (this%ncx-2)/2 

    ! La moyenne de Ex est nulle donc les composantes de Fourier 
    ! correspondant a k=0 sont nulles
    this%ex(1) = 0.

    ! Calcul des autres composantes de Fourier
    do ik=1,nxh1
       kx= ik*kx0
       k2 = kx*kx
       this%ex(2*ik) = kx/k2*this%rho(2*ik+1)
       this%ex(2*ik+1) = -kx/k2*this%rho(2*ik)
       this%rho(2*ik) = 1/k2*this%rho(2*ik)
       this%rho(2*ik+1) = 1/k2*this%rho(2*ik+1)
    end do

    this%ex(this%ncx)= 0.          ! because Im(rho_n/2)=0
    this%rho(this%ncx) = ((GET_FIELD_DELTA_X1( rho )/sll_pi)**2)*GET_FIELD_NCELLS_X1( rho )

 
    ! Faire une FFT inverse  dans la direction x de E
    call fftinv(this%fftx,this%ex)
    call fftinv(this%fftx,this%rho)

    ! Copy local Ex into field data structure
    do i=1, this%ncx
       FIELD_1D_AT_I(ex,i) = this%ex(i)
    end do
    ! complete last term by periodicity
    FIELD_1D_AT_I(ex,this%ncx+1) = this%ex(1)

  end subroutine solve_poisson1dp

!!$  subroutine solve_poisson1dp_axisymetrique(ex,rho,geomx)
!!$    sll_real64, dimension(:)                   :: ex, rho
!!$    type(geometry1d),intent(in)              :: geomx
!!$    sll_int32                                  :: i, j         ! indices de boucle
!!$    ! Variables de travail
!!$    sll_real64                                 :: pas, integrale, xi, xim1
!!$
!!$    if(mod(geomx%nx,2)==0) then
!!$       write(0,*) 'nx must be odd in axisymmetric mode.'
!!$       stop
!!$    endif
!!$
!!$    integrale=0
!!$    ex=0        
!!$
!!$    pas=geomx%dx
!!$    do i=2+(geomx%nx-1)/2,geomx%nx
!!$       xi=geomx%x0+(i-1)*pas
!!$       xim1=geomx%x0+(i-2)*pas
!!$       integrale=integrale+geomx%dx*(xim1*rho(i-1)+xi*rho(i))/2.
!!$       ex(i) = integrale/xi
!!$       ex(geomx%nx-i+1) = -ex(i)
!!$    end do
!!$
!!$  end subroutine solve_poisson1dp_axisymetrique



end module sll_poisson_1D_periodic
