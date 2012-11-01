module sll_poisson_1d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use numeric_constants

  implicit none
  private
  public :: new, solve

  type, public :: poisson_1d_periodic
     sll_int32                         :: nc_eta1
     sll_real64                        :: eta1_min
     sll_real64                        :: eta1_max
     sll_real64, dimension(:), pointer :: wsave
     sll_real64, dimension(:), pointer :: work
#ifndef STDF95
  contains
     procedure :: initialize => new_poisson_1d_periodic
     procedure :: solve => solve_poisson_1d_periodic
#endif
  end type poisson_1d_periodic

  interface new
     module procedure new_poisson_1d_periodic
  end interface

  interface solve
     module procedure solve_poisson_1d_periodic 
  end interface

contains

  subroutine new_poisson_1d_periodic(this,eta1_min,eta1_max,nc_eta1,error)
#ifdef STDF95
    type(poisson_1d_periodic),intent(out) :: this
#else
    class(poisson_1d_periodic),intent(out) :: this
#endif
    sll_int32,intent(in)                   :: nc_eta1
    sll_int32, intent(out)                 :: error 
    sll_real64, intent(in)                 :: eta1_min
    sll_real64, intent(in)                 :: eta1_max

    error = 0
    ! geometry
    this%nc_eta1  = nc_eta1
    this%eta1_min = eta1_min
    this%eta1_max = eta1_max

    SLL_ALLOCATE(this%wsave(2*this%nc_eta1+15),error)

    call dffti(this%nc_eta1,this%wsave)

    ! Allocate auxiliary arrays for fft in order to keep rhs unchanged
    SLL_ALLOCATE(this%work(nc_eta1+1),error)

  end subroutine new_poisson_1d_periodic


  subroutine solve_poisson_1d_periodic(this, field, rhs)
#ifdef STDF95
    type(poisson_1d_periodic),intent(inout) :: this
#else
    class(poisson_1d_periodic),intent(inout) :: this
#endif
    sll_real64, dimension(:), intent(out)     :: field
    sll_real64, dimension(:), intent(in)      :: rhs
    sll_int32                                 :: ik
    sll_real64                                :: kx0, kx, k2

    ! Check that field and rhs are both associated to the 
    ! same mesh with the right number of cells 
    ! that has been initialized in new_poisson_1d_periodic
    SLL_ASSERT(size(field)==this%nc_eta1+1)
    SLL_ASSERT(size(rhs)==this%nc_eta1+1)

    ! copy rhs into auxiliary array for fftpack
    ! in order to keep rhs unchanged
    this%work = rhs 

    ! Forward FFT 
    call dfftf( this%nc_eta1, this%work, this%wsave)

    this%work = this%work /this%nc_eta1      ! normalize FFT

    kx0  = 2_f64*sll_pi/(this%eta1_max-this%eta1_min)

    ! La moyenne de Ex est nulle donc les composantes de Fourier 
    ! correspondant a k=0 sont nulles
    field(1) = 0.

    ! Calcul des autres composantes de Fourier
    do ik=1,(this%nc_eta1-2)/2 
       kx= ik*kx0
       k2 = kx*kx
       field(2*ik)       = kx/k2*this%work(2*ik+1)
       field(2*ik+1)     = -kx/k2*this%work(2*ik)
       this%work(2*ik)   = 1/k2*this%work(2*ik)
       this%work(2*ik+1) = 1/k2*this%work(2*ik+1)
    end do

    field(this%nc_eta1)= 0.          ! because Im(rhs/2)=0
 
    ! Backward FFT 

    call dfftb( this%nc_eta1, field,  this%wsave )

    ! complete last term by periodicity
    field(this%nc_eta1+1) = field(1)

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
