module sll_poisson_1d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
  use numeric_constants
 !PN use geometry1d_module
  use fft1d_module

  implicit none
  private
  public :: new, solve

  type, public :: poisson_1d_periodic
     sll_real64, dimension(:), pointer       :: rhs
     sll_real64, dimension(:), pointer       :: e_field
     type(fft1dclass)                        :: fft
     sll_int32                               :: n_cells
  end type poisson_1d_periodic

  interface new
     module procedure new_poisson_1d_periodic
  end interface

  interface solve
     module procedure solve_poisson_1d_periodic !, solve_poisson1dp_axisymetrique
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
    SLL_ALLOCATE(this%e_field(n_cells+2),error)

  end subroutine new_poisson_1d_periodic

  subroutine solve_poisson_1d_periodic(this,e_field,rhs)
    type(poisson_1d_periodic)                 :: this
    type(field_1D_vec1), pointer              :: e_field
    type(field_1D_vec1), pointer              :: rhs
    sll_int32                                 :: i, ik
    sll_int32                                 :: nxh1
    sll_real64                                :: kx0, kx, k2

    ! Check that e_field and rhs are both associated to the 
    ! same mesh with the right number of cells 
    ! that has been initialized in new_poisson_1d_periodic
    SLL_ASSERT(associated(e_field%descriptor,target=rhs%descriptor))
    SLL_ASSERT(rhs%descriptor%nc_eta1 == this%n_cells )

    ! copy rhs into auxiliary array for fftpack
    this%rhs(1:this%n_cells+1) = FIELD_DATA(rhs) 

    ! Faire une FFT dans la direction x de rhs
    call fft(this%fft, this%rhs)

    ! Calcul de la transformee de Fourier de E a partir de celle de rhs
    kx0  = 2*sll_pi/(this%n_cells * rhs%descriptor%delta_eta1 )
    nxh1 = (this%n_cells-2)/2 

    ! La moyenne de Ex est nulle donc les composantes de Fourier 
    ! correspondant a k=0 sont nulles
    this%e_field(1) = 0.

    ! Calcul des autres composantes de Fourier
    do ik=1,nxh1
       kx= ik*kx0
       k2 = kx*kx
       this%e_field(2*ik) = kx/k2*this%rhs(2*ik+1)
       this%e_field(2*ik+1) = -kx/k2*this%rhs(2*ik)
       this%rhs(2*ik) = 1/k2*this%rhs(2*ik)
       this%rhs(2*ik+1) = 1/k2*this%rhs(2*ik+1)
    end do

    this%e_field(this%n_cells)= 0.          ! because Im(rhs/2)=0
    !PN this%rhs(this%n_cells) = ((GET_FIELD_DELTA_X1( rhs )/sll_pi)**2)*GET_FIELD_NCELLS_X1( rhs )

 
    ! Faire une FFT inverse  dans la direction x de E
    call fftinv(this%fft,this%e_field)
    call fftinv(this%fft,this%rhs)

    ! Copy local Ex into field data structure
    do i=1, this%n_cells
       FIELD_1D_AT_I(e_field,i) = this%e_field(i)
    end do
    ! complete last term by periodicity
    FIELD_1D_AT_I(e_field,this%n_cells+1) = this%e_field(1)

  end subroutine solve_poisson_1d_periodic

!!$  subroutine solve_poisson1dp_axisymetrique(e_field,rhs,geomx)
!!$    sll_real64, dimension(:)                   :: e_field, rhs
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
!!$    e_field=0        
!!$
!!$    pas=geomx%dx
!!$    do i=2+(geomx%nx-1)/2,geomx%nx
!!$       xi=geomx%x0+(i-1)*pas
!!$       xim1=geomx%x0+(i-2)*pas
!!$       integrale=integrale+geomx%dx*(xim1*rhs(i-1)+xi*rhs(i))/2.
!!$       e_field(i) = integrale/xi
!!$       e_field(geomx%nx-i+1) = -e_field(i)
!!$    end do
!!$
!!$  end subroutine solve_poisson1dp_axisymetrique



end module sll_poisson_1d_periodic
