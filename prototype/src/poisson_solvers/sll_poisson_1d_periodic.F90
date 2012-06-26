module sll_poisson_1d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use numeric_constants
  use fft1d_module

  implicit none
  private
  public :: new, solve

  type, public :: poisson_1d_periodic
     sll_real64, dimension(:), pointer       :: rhs
     sll_real64, dimension(:), pointer       :: e_field
     type(fft1dclass)                        :: fft
     sll_int32                               :: nc_eta1
     sll_real64                              :: eta1_min
     sll_real64                              :: eta1_max
  contains
     procedure :: initialize => new_poisson_1d_periodic
     procedure :: solve => solve_poisson_1d_periodic
  end type poisson_1d_periodic

  interface new
     module procedure new_poisson_1d_periodic
  end interface

  interface solve
     module procedure solve_poisson_1d_periodic 
  end interface

contains

  subroutine new_poisson_1d_periodic(this,eta1_min,eta1_max,nc_eta1,error)
    class(poisson_1d_periodic),intent(out) :: this
    sll_int32,intent(in)                   :: nc_eta1
    sll_int32, intent(out)                 :: error 
    sll_real64, intent(in)                 :: eta1_min
    sll_real64, intent(in)                 :: eta1_max

    ! Indicateur d'erreur
    error = 0
    ! Initialisation de la geometrie
    this%nc_eta1=nc_eta1
    this%eta1_min=eta1_min
    this%eta1_max=eta1_max
    ! Initialisation des fft (le nombre de points pris en compte est n
    call initdoubfft(this%fft, this%nc_eta1)
    ! Allocate auxiliary arrays for fft
    SLL_ALLOCATE(this%rhs(nc_eta1+2),error)
    SLL_ALLOCATE(this%e_field(nc_eta1+2),error)

  end subroutine new_poisson_1d_periodic


  subroutine solve_poisson_1d_periodic(this, efield, rhs)
    class(poisson_1d_periodic)                :: this
    sll_real64, dimension(:), intent(out)     :: efield
    sll_real64, dimension(:), intent(in)      :: rhs
    sll_int32                                 :: i, ik
    sll_int32                                 :: nxh1
    sll_real64                                :: kx0, kx, k2

    ! Check that e_field and rhs are both associated to the 
    ! same mesh with the right number of cells 
    ! that has been initialized in new_poisson_1d_periodic
    !SLL_ASSERT(associated(efield%mesh,target=rhs%mesh))
    !SLL_ASSERT(rhs%mesh%nc_eta1 == this%nc_eta1 )
    SLL_ASSERT(size(efield)==this%nc_eta1+1)
    SLL_ASSERT(size(rhs)==this%nc_eta1+1)

    ! copy rhs into auxiliary array for fftpack
    this%rhs = rhs 

    ! Faire une FFT dans la direction x de rhs
    call fft(this%fft, this%rhs)

    ! Calcul de la transformee de Fourier de E a partir de celle de rhs
    !kx0  = 2*sll_pi/(rhs%mesh%x1_at_node(this%nc_eta1+1)  &
         !- rhs%mesh%x1_at_node(1))
    kx0  = 2*sll_pi/(this%eta1_max-this%eta1_min)
    !print*,' delta_eta1 = ', this%delta_eta1 
    !print*,' nc_eta1*delta_eta1 = ', this%delta_eta1 * this%nc_eta1
    !print*,rhs%mesh%x1_at_node(this%nc_eta1+1) - rhs%mesh%x1_at_node(1)
    !print*, 'kx0 ', this%nc_eta1, rhs%mesh%x1_at_node(1), rhs%mesh%x1_at_node(this%nc_eta1+1)
    nxh1 = (this%nc_eta1-2)/2 

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

    this%e_field(this%nc_eta1)= 0.          ! because Im(rhs/2)=0
    !PN this%rhs(this%nc_eta1) = ((GET_FIELD_DELTA_X1( rhs )/sll_pi)**2)*GET_FIELD_NCELLS_X1( rhs )

 
    ! Faire une FFT inverse  dans la direction x de E
    call fftinv(this%fft,this%e_field)
    call fftinv(this%fft,this%rhs)

    ! Copy local Ex into field data structure
    do i=1, this%nc_eta1
       efield(i) = this%e_field(i)
    end do
    ! complete last term by periodicity
    efield(this%nc_eta1+1) = this%e_field(1)

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
