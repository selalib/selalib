!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> Module to solve Poisson equation on one dimensional mesh using FFT
!> transform.
module sll_poisson_1d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_constants

  implicit none
  private
  public :: initialize, new, solve

  !> Solver data structure
  type, public :: poisson_1d_periodic
     sll_int32                         :: nc_eta1 !< number of cells
     sll_real64                        :: eta1_min !< left corner
     sll_real64                        :: eta1_max !< right corner
     sll_real64, dimension(:), pointer :: wsave !< array used by fftpack
     sll_real64, dimension(:), pointer :: work  !< array used by fftpack
  end type poisson_1d_periodic

  !> Create a new poisson solver on 1d mesh
  interface new
     module procedure new_poisson_1d_periodic
  end interface

  !> Initialize a new poisson solver on 1d mesh
  interface initialize
     module procedure initialize_poisson_1d_periodic
  end interface

  !> Solve the Poisson equation on 1d mesh and compute the potential
  interface solve
     module procedure solve_poisson_1d_periodic 
  end interface

contains

  !> Create a new solver
  !> @return
  function new_poisson_1d_periodic(eta1_min,eta1_max,nc_eta1,error) &
     result(this)
     type(poisson_1d_periodic),pointer :: this     !< Solver data structure
     sll_int32,intent(in)              :: nc_eta1  !< number of cells
     sll_int32, intent(out)            :: error    !< error code
     sll_real64, intent(in)            :: eta1_min !< left corner
     sll_real64, intent(in)            :: eta1_max !< right corner

     SLL_ALLOCATE(this, error)
     call initialize_poisson_1d_periodic(this,eta1_min,eta1_max,nc_eta1,error)

  end function new_poisson_1d_periodic 
  
  !> Initialize the solver
  subroutine initialize_poisson_1d_periodic(this,eta1_min,eta1_max,nc_eta1,error)

    type(poisson_1d_periodic),intent(out) :: this     !< Solver data structure
    sll_int32,intent(in)                  :: nc_eta1  !< number of cells
    sll_int32, intent(out)                :: error    !< error code
    sll_real64, intent(in)                :: eta1_min !< left corner
    sll_real64, intent(in)                :: eta1_max !< right corner

    error = 0
    ! geometry
    this%nc_eta1  = nc_eta1
    this%eta1_min = eta1_min
    this%eta1_max = eta1_max

    SLL_ALLOCATE(this%wsave(2*this%nc_eta1+15),error)

    call dffti(this%nc_eta1,this%wsave)

    ! Allocate auxiliary arrays for fft in order to keep rhs unchanged
    SLL_ALLOCATE(this%work(nc_eta1+1),error)

  end subroutine initialize_poisson_1d_periodic

  subroutine solve_poisson_1d_periodic(this, field, rhs)

    type(poisson_1d_periodic),intent(inout) :: this
    sll_real64, dimension(:), intent(out)   :: field
    sll_real64, dimension(:), intent(in)    :: rhs
    sll_int32                               :: ik
    sll_real64                              :: kx0, kx, k2

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
