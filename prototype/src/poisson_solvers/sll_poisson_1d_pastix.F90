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

!> @brief
!> Solve Poisson equation 1D using Pastix solver
!> @details
!> More about Pastix http://pastix.gforge.inria.fr
module sll_poisson_1d_pastix

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_pastix.h"

  implicit none
  private
  public :: initialize, new, solve

  !> Solver data structure
  type, public :: poisson_1d_pastix
     sll_int32                         :: nc_eta1 !< number of cells
     sll_real64                        :: eta1_min !< left corner
     sll_real64                        :: eta1_max !< right corner
     type(pastix_solver)               :: pastix
  end type poisson_1d_pastix

  !> Create a new poisson solver on 1d mesh
  interface new
     module procedure new_poisson_1d_pastix
  end interface

  !> Initialize a new poisson solver on 1d mesh
  interface initialize
     module procedure initialize_poisson_1d_pastix
  end interface

  !> Solve the Poisson equation on 1d mesh and compute the potential
  interface solve
     module procedure solve_poisson_1d_pastix 
  end interface

contains

  !> Create a new solver
  function new_poisson_1d_pastix(eta1_min,eta1_max,nc_eta1,error) &
     result(this)
     type(poisson_1d_pastix),pointer :: this     !< Solver data structure
     sll_int32,intent(in)              :: nc_eta1  !< number of cells
     sll_int32, intent(out)            :: error    !< error code
     sll_real64, intent(in)            :: eta1_min !< left corner
     sll_real64, intent(in)            :: eta1_max !< right corner

     SLL_ALLOCATE(this, error)
     call initialize_poisson_1d_pastix(this,eta1_min,eta1_max,nc_eta1,error)

  end function new_poisson_1d_pastix 
  
  !> Initialize the solver
  subroutine initialize_poisson_1d_pastix(this,eta1_min,eta1_max,nc_eta1,error)

    type(poisson_1d_pastix),intent(out) :: this     !< Solver data structure
    sll_int32,intent(in)                :: nc_eta1  !< number of cells
    sll_int32, intent(out)              :: error    !< error code
    sll_real64, intent(in)              :: eta1_min !< left corner
    sll_real64, intent(in)              :: eta1_max !< right corner
    sll_int32                           :: i
    sll_int32                           :: j
    sll_int32                           :: nnzeros

    error = 0
    ! geometry
    this%nc_eta1  = nc_eta1
    this%eta1_min = eta1_min
    this%eta1_max = eta1_max

    nnzeros = 3*(nc_eta1-1) + 2
    call initialize(this%pastix, nc_eta1-1, nnzeros)

    j=1
    do i = 1, nc_eta1-1
       this%pastix%colptr(i) = j
       this%pastix%row(j)    = i
       this%pastix%avals(j) = 2
       j=j+1
       if (i /= nc_eta1-1) then
          this%pastix%row(j)   = i+1
          this%pastix%avals(j) = -1.
          j=j+1
       end if
    end do
    this%pastix%colptr(nc_eta1) = j

    call factorize(this%pastix)

  end subroutine initialize_poisson_1d_pastix

  subroutine solve_poisson_1d_pastix(this, field, rhs)

    type(poisson_1d_pastix),intent(inout)   :: this
    sll_real64, dimension(:), intent(inout) :: field
    sll_real64, dimension(:), intent(inout) :: rhs

    ! Check that field and rhs are both associated to the 
    ! same mesh with the right number of cells 
    ! that has been initialized in new_poisson_1d_pastix
    SLL_ASSERT(size(field)==this%nc_eta1+1)
    SLL_ASSERT(size(rhs)==this%nc_eta1+1)
    field = rhs
    call solve(this%pastix, field(2:this%nc_eta1))

  end subroutine solve_poisson_1d_pastix

end module sll_poisson_1d_pastix
