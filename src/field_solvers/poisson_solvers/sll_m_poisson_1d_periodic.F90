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

!> @ingroup poisson_solvers
!> Module to sll_o_solve Poisson equation on one dimensional mesh using FFT
!> transform.
module sll_m_poisson_1d_periodic

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi
  use sll_m_poisson_1d_base, only: &
    sll_c_poisson_1d_base

  implicit none

  public :: &
    sll_o_initialize, &
    sll_o_new, &
    sll_t_poisson_1d_periodic, &
    sll_o_solve,               &
    sll_f_new_poisson_1d_periodic

  private

  !> Solver data structure
  type :: sll_t_poisson_1d_periodic
     sll_int32                         :: nc_eta1 !< number of cells
     sll_real64                        :: eta1_min !< left corner
     sll_real64                        :: eta1_max !< right corner
     sll_real64, dimension(:), pointer :: wsave !< array used by fftpack
     sll_real64, dimension(:), pointer :: work  !< array used by fftpack
  end type sll_t_poisson_1d_periodic

  type,extends(sll_c_poisson_1d_base) :: sll_c_poisson_1d_periodic  
  
    type(sll_t_poisson_1d_periodic), pointer :: poiss
  
  contains
    procedure, pass(poisson) :: sll_o_initialize => &
      initialize_poisson_1d_periodic_wrapper
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_1d_periodic
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_1d_periodic
      
  end type sll_c_poisson_1d_periodic 

  !> Create a sll_o_new poisson solver on 1d mesh
  interface sll_o_new
     module procedure new_poisson_1d_periodic
  end interface

  !> sll_o_initialize a sll_o_new poisson solver on 1d mesh
  interface sll_o_initialize
     module procedure initialize_poisson_1d_periodic
  end interface

  !> sll_o_solve the Poisson equation on 1d mesh and compute the potential
  interface sll_o_solve
     module procedure solve_poisson_1d_periodic 
  end interface

contains

  !> Create a sll_o_new solver
  !> @return
  function new_poisson_1d_periodic(eta1_min,eta1_max,nc_eta1,error) &
     result(this)
     type(sll_t_poisson_1d_periodic),pointer :: this     !< Solver data structure
     sll_int32,intent(in)              :: nc_eta1  !< number of cells
     sll_int32, intent(out)            :: error    !< error code
     sll_real64, intent(in)            :: eta1_min !< left corner
     sll_real64, intent(in)            :: eta1_max !< right corner

     SLL_ALLOCATE(this, error)
     call initialize_poisson_1d_periodic(this,eta1_min,eta1_max,nc_eta1,error)

  end function new_poisson_1d_periodic 
  
  !> sll_o_initialize the solver
  subroutine initialize_poisson_1d_periodic(this,eta1_min,eta1_max,nc_eta1,error)

    type(sll_t_poisson_1d_periodic),intent(out) :: this     !< Solver data structure
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

    type(sll_t_poisson_1d_periodic),intent(inout) :: this
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

    kx0  = 2_f64*sll_p_pi/(this%eta1_max-this%eta1_min)

    ! La moyenne de Ex est nulle donc les composantes de Fourier 
    ! correspondant a k=0 sont nulles
    field(1) = 0.0_f64

    ! Calcul des autres composantes de Fourier
    do ik=1,(this%nc_eta1-2)/2 
       kx= ik*kx0
       k2 = kx*kx
       field(2*ik)       = kx/k2*this%work(2*ik+1)
       field(2*ik+1)     = -kx/k2*this%work(2*ik)
       this%work(2*ik)   = 1/k2*this%work(2*ik)
       this%work(2*ik+1) = 1/k2*this%work(2*ik+1)
    end do

    field(this%nc_eta1)= 0.0_f64          ! because Im(rhs/2)=0
 
    ! Backward FFT 

    call dfftb( this%nc_eta1, field,  this%wsave )

    ! complete last term by periodicity
    field(this%nc_eta1+1) = field(1)

  end subroutine solve_poisson_1d_periodic

  function sll_f_new_poisson_1d_periodic( &
    eta1_min, &
    eta1_max, &
    nc_eta1) &
    result(poisson)
      
    type(sll_c_poisson_1d_periodic),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_1d_periodic_wrapper( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1)    
  end function sll_f_new_poisson_1d_periodic
  
  subroutine initialize_poisson_1d_periodic_wrapper( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1)
    class(sll_c_poisson_1d_periodic) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32 :: ierr

    
    poisson%poiss => sll_o_new( &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      ierr)
    
  end subroutine initialize_poisson_1d_periodic_wrapper
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_1d_periodic( poisson, phi, rho )
    class(sll_c_poisson_1d_periodic), target :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: phi
    
    print *,'#compute_phi_from_rho_1d_periodic'
    print *,'#not implemented yet'
    phi = 0._f64
    if(.not.(associated(poisson%poiss)))then
      print *,'#poisson%poiss not associated'
    endif
    print *,maxval(rho)  
    stop
    
  end subroutine compute_phi_from_rho_1d_periodic

  ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
  subroutine compute_E_from_rho_1d_periodic( poisson, E, rho )
    class(sll_c_poisson_1d_periodic) :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: E
      
    call sll_o_solve(poisson%poiss, E, rho)
           
  end subroutine compute_E_from_rho_1d_periodic
  
end module sll_m_poisson_1d_periodic
