module sll_poisson_1d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_1d.h"
  use numeric_constants
  use fft1d_module

  implicit none
  private
  public :: new, solve

  type, public :: poisson_1d_periodic
     sll_real64, dimension(:), pointer       :: rhs
     sll_real64, dimension(:), pointer       :: phi
     type(fft1dclass)                        :: fft
     sll_int32                               :: n_cells
  end type poisson_1d_periodic

  interface new
     module procedure new_poisson_1d_periodic
  end interface

  interface solve
     module procedure solve_poisson_1d_periodic 
  end interface

contains

  subroutine new_poisson_1d_periodic(this,n_cells,error)
    type(poisson_1d_periodic),intent(out) :: this
    sll_int32,intent(in)                  :: n_cells
    sll_int32, intent(out)                :: error 

    ! Indicateur d'erreur
    error = 0
    ! Initialisation de la geometrie
    this%n_cells=n_cells
    ! Initialisation des fft (le nombre de points pris en compte est n
    call initdoubfft(this%fft, this%n_cells)
    ! Allocate auxiliary arrays for fft
    SLL_ALLOCATE(this%rhs(n_cells+2),error)
    SLL_ALLOCATE(this%phi(n_cells+2),error)

  end subroutine new_poisson_1d_periodic

  ! Even though the variable is called phi. This routine computes the electric
  ! field from the charge density rho
  subroutine solve_poisson_1d_periodic(this,phi,rhs)
    type(poisson_1d_periodic)                 :: this
    type(scalar_field_1d)                     :: phi
    type(scalar_field_1d)                     :: rhs
    sll_int32                                 :: i, ik
    sll_int32                                 :: nxh1
    sll_real64                                :: kx0, kx, k2

    ! Check that phi and rhs are both associated to the 
    ! same mesh with the right number of cells 
    ! that has been initialized in new_poisson_1d_periodic
    SLL_ASSERT(associated(phi%mesh,target=rhs%mesh))
    SLL_ASSERT(rhs%mesh%nc_eta1 == this%n_cells )

    ! copy rhs into auxiliary array for fftpack
    this%rhs(1:this%n_cells+1) = FIELD_DATA(rhs) 

    ! Faire une FFT dans la direction x de rhs
    call fft(this%fft, this%rhs)

    ! Calcul de la transformee de Fourier de E a partir de celle de rhs
    kx0  = 2*sll_pi/(rhs%mesh%x1_at_node(this%n_cells+1)  &
         - rhs%mesh%x1_at_node(1))
    print*, 'kx0 ', this%n_cells, rhs%mesh%x1_at_node(1), rhs%mesh%x1_at_node(this%n_cells+1)
    nxh1 = (this%n_cells-2)/2 

    ! La moyenne de Ex est nulle donc les composantes de Fourier 
    ! correspondant a k=0 sont nulles
    this%phi(1) = 0.

    ! Calcul des autres composantes de Fourier
    do ik=1,nxh1
       kx= ik*kx0
       k2 = kx*kx
       this%phi(2*ik) = kx/k2*this%rhs(2*ik+1)
       this%phi(2*ik+1) = -kx/k2*this%rhs(2*ik)
       this%rhs(2*ik) = 1/k2*this%rhs(2*ik)
       this%rhs(2*ik+1) = 1/k2*this%rhs(2*ik+1)
    end do

    this%phi(this%n_cells)= 0.          ! because Im(rhs/2)=0
    !PN this%rhs(this%n_cells) = ((GET_FIELD_DELTA_X1( rhs )/sll_pi)**2)*GET_FIELD_NCELLS_X1( rhs )

 
    ! Faire une FFT inverse  dans la direction x de E
    call fftinv(this%fft,this%phi)
    call fftinv(this%fft,this%rhs)

    ! Copy local Ex into field data structure
    do i=1, this%n_cells
       FIELD_1D_AT_I(phi,i) = this%phi(i)
    end do
    ! complete last term by periodicity
    FIELD_1D_AT_I(phi,this%n_cells+1) = this%phi(1)

  end subroutine solve_poisson_1d_periodic

!!$  subroutine solve_poisson1dp_axisymetrique(phi,rhs,geomx)
!!$    sll_real64, dimension(:)                   :: phi, rhs
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
!!$    phi=0        
!!$
!!$    pas=geomx%dx
!!$    do i=2+(geomx%nx-1)/2,geomx%nx
!!$       xi=geomx%x0+(i-1)*pas
!!$       xim1=geomx%x0+(i-2)*pas
!!$       integrale=integrale+geomx%dx*(xim1*rhs(i-1)+xi*rhs(i))/2.
!!$       phi(i) = integrale/xi
!!$       phi(geomx%nx-i+1) = -phi(i)
!!$    end do
!!$
!!$  end subroutine solve_poisson1dp_axisymetrique



end module sll_poisson_1d_periodic
