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
module sll_m_poisson_1d_hmf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi, sll_p_i1
  use sll_m_poisson_1d_base, only: &
    sll_c_poisson_1d_base
  use sll_m_fft, only: sll_t_fft, &
    sll_s_fft_init_c2r_1d,        &
    sll_s_fft_init_r2c_1d,        &
    sll_s_fft_exec_c2r_1d,        &
    sll_s_fft_exec_r2c_1d,        &
    sll_s_fft_free

  implicit none

  public :: &
    sll_o_initialize, &
    sll_o_new, &
    sll_t_poisson_1d_hmf, &
    sll_o_solve,               &
    sll_f_new_poisson_1d_hmf

  private

  !> Solver data structure
  type :: sll_t_poisson_1d_hmf
     sll_int32                         :: nc_eta1   !< number of cells
     sll_real64                        :: eta1_min  !< left corner
     sll_real64                        :: eta1_max  !< right corner
     sll_real64             , pointer  :: tmp(:)    !<  array to store rho
     sll_comp64             , pointer  :: rhok(:)   !< fft(rho)
     type(sll_t_fft)                   :: fw        !< forward fft plan
     type(sll_t_fft)                   :: bw        !< backward fft plan
  end type sll_t_poisson_1d_hmf

  type,extends(sll_c_poisson_1d_base) :: sll_c_poisson_1d_hmf  
  
    type(sll_t_poisson_1d_hmf), pointer :: poiss
  
  contains
    procedure, pass(poisson) :: sll_o_initialize => &
      initialize_poisson_1d_hmf_wrapper
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_1d_hmf
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_1d_hmf
      
  end type sll_c_poisson_1d_hmf 

  !> Create a sll_o_new poisson solver on 1d mesh
  interface sll_o_new
     module procedure new_poisson_1d_hmf
  end interface

  !> sll_o_initialize a sll_o_new poisson solver on 1d mesh
  interface sll_o_initialize
     module procedure initialize_poisson_1d_hmf
  end interface

  !> sll_o_solve the Poisson equation on 1d mesh and compute the potential
  interface sll_o_solve
     module procedure solve_poisson_1d_hmf 
  end interface

contains

  !> Create a sll_o_new solver
  !> @return
  function new_poisson_1d_hmf(eta1_min,eta1_max,nc_eta1,error) &
     result(self)
     type(sll_t_poisson_1d_hmf),pointer :: self     !< Solver data structure
     sll_int32,intent(in)              :: nc_eta1  !< number of cells
     sll_int32, intent(out)            :: error    !< error code
     sll_real64, intent(in)            :: eta1_min !< left corner
     sll_real64, intent(in)            :: eta1_max !< right corner

     SLL_ALLOCATE(self, error)
     call initialize_poisson_1d_hmf(self,eta1_min,eta1_max,nc_eta1,error)

  end function new_poisson_1d_hmf 
  
  !> sll_o_initialize the solver
  subroutine initialize_poisson_1d_hmf(self,eta1_min,eta1_max,nc_eta1,error)

    type(sll_t_poisson_1d_hmf),intent(out) :: self     !< Solver data structure
    sll_int32,intent(in)                  :: nc_eta1  !< number of cells
    sll_int32, intent(out)                :: error    !< error code
    sll_real64, intent(in)                :: eta1_min !< left corner
    sll_real64, intent(in)                :: eta1_max !< right corner

    sll_int32 :: i

    error = 0
    ! geometry
    self%nc_eta1  = nc_eta1
    self%eta1_min = eta1_min
    self%eta1_max = eta1_max

    SLL_ALLOCATE(self%rhok(1:nc_eta1),error)
    SLL_ALLOCATE(self%tmp(1:nc_eta1),error)

    call sll_s_fft_init_r2c_1d(self%fw,nc_eta1,self%tmp,self%rhok)
    call sll_s_fft_init_c2r_1d(self%bw,nc_eta1,self%rhok,self%tmp)

  end subroutine initialize_poisson_1d_hmf

  subroutine solve_poisson_1d_hmf(self, field, rhs)

    type(sll_t_poisson_1d_hmf),intent(inout) :: self
    sll_real64, dimension(:), intent(out)   :: field
    sll_real64, dimension(:), intent(in)    :: rhs

    SLL_ASSERT(size(rhs)==self%nc_eta1+1)

    self%tmp = rhs(1:self%nc_eta1)
    call sll_s_fft_exec_r2c_1d(self%fw, self%tmp, self%rhok)

    self%rhok(2)  = self%rhok(2) * sll_p_pi * sll_p_i1 ! We keep only one mode
    self%rhok(1)  = 0.0_f64
    self%rhok(3:) = 0.0_f64

    call sll_s_fft_exec_c2r_1d(self%bw, self%rhok, self%tmp)

    field(1:self%nc_eta1) = self%tmp / self%nc_eta1
    ! complete last term by periodicity
    field(self%nc_eta1+1) = field(1)

  end subroutine solve_poisson_1d_hmf

  function sll_f_new_poisson_1d_hmf( &
    eta1_min, &
    eta1_max, &
    nc_eta1) &
    result(poisson)
      
    type(sll_c_poisson_1d_hmf),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_1d_hmf_wrapper( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1)    
  end function sll_f_new_poisson_1d_hmf
  
  subroutine initialize_poisson_1d_hmf_wrapper( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1)
    class(sll_c_poisson_1d_hmf) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32 :: ierr

    
    poisson%poiss => sll_o_new( &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      ierr)
    
  end subroutine initialize_poisson_1d_hmf_wrapper
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_1d_hmf( poisson, phi, rho )
    class(sll_c_poisson_1d_hmf), target :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: phi
    
    print *,'#compute_phi_from_rho_1d_hmf'
    print *,'#not implemented yet'
    phi = 0._f64
    if(.not.(associated(poisson%poiss)))then
      print *,'#poisson%poiss not associated'
    endif
    print *,maxval(rho)  
    stop
    
  end subroutine compute_phi_from_rho_1d_hmf

  ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
  subroutine compute_E_from_rho_1d_hmf( poisson, E, rho )
    class(sll_c_poisson_1d_hmf) :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: E
      
    call sll_o_solve(poisson%poiss, E, rho)
           
  end subroutine compute_E_from_rho_1d_hmf
  
end module sll_m_poisson_1d_hmf
