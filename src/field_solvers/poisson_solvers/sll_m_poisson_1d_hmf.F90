!**************************************************************
!  Copyright INRIA
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
!> @author Pierre Navaro
!> @brief
!> Module to solve Poisson equation for the HMF model 
!> @details
!> You could find more information in
!>
!> "On numerical Landau damping for splitting methods applied to the Vlasov-HMF model"
!>
!> by Erwan Faou, Romain Horsin and Frédéric Rousset.
!> https://arxiv.org/abs/1510.06555
module sll_m_poisson_1d_hmf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: sll_p_pi, sll_p_i1, sll_p_i0
  use sll_m_poisson_1d_base, only: sll_c_poisson_1d_base
  use sll_m_fft, only: sll_t_fft, &
    sll_s_fft_init_c2r_1d,        &
    sll_s_fft_init_r2c_1d,        &
    sll_s_fft_exec_c2r_1d,        &
    sll_s_fft_exec_r2c_1d,        &
    sll_s_fft_free

  implicit none

  public ::                     &
    sll_s_init_poisson_1d_hmf,  &
    sll_t_poisson_1d_hmf,       &
    sll_s_solve_poisson_1d_hmf, &
    sll_s_free_poisson_1d_hmf

  private

  !> Implementation of the poisson 1d solver for the Vlasov-HMF model
  type, extends(sll_c_poisson_1d_base) :: sll_t_poisson_1d_hmf  
  
     sll_int32                         :: nc_eta1   !< number of cells
     sll_real64                        :: eta1_min  !< left corner
     sll_real64                        :: eta1_max  !< right corner
     sll_real64             , pointer  :: tmp(:)    !<  array to store rho
     sll_comp64             , pointer  :: rhok(:)   !< fft(rho)
     type(sll_t_fft)                   :: fw        !< forward fft plan
     type(sll_t_fft)                   :: bw        !< backward fft plan
  
  contains

    procedure, pass(self)    :: init => sll_s_init_poisson_1d_hmf
    procedure, pass(self)    :: free => sll_s_free_poisson_1d_hmf
    procedure, pass(poisson) :: compute_phi_from_rho => compute_phi_from_rho_1d_hmf
    procedure, pass(poisson) :: compute_e_from_rho => compute_E_from_rho_1d_hmf
      
  end type sll_t_poisson_1d_hmf 

contains

  !> Initialize the poisson 1d hmf solver
  subroutine sll_s_init_poisson_1d_hmf(self,eta1_min,eta1_max,nc_eta1,error)

    class(sll_t_poisson_1d_hmf), intent(out)  :: self     !< Solver structure
    sll_int32,intent(in)                      :: nc_eta1  !< number of cells
    sll_int32, intent(out)                    :: error    !< error code
    sll_real64, intent(in)                    :: eta1_min !< left corner
    sll_real64, intent(in)                    :: eta1_max !< right corner

    error = 0
    ! geometry
    self%nc_eta1  = nc_eta1
    self%eta1_min = eta1_min
    self%eta1_max = eta1_max

    SLL_ALLOCATE(self%rhok(1:nc_eta1),error)
    SLL_ALLOCATE(self%tmp(1:nc_eta1),error)

    call sll_s_fft_init_r2c_1d(self%fw,nc_eta1,self%tmp,self%rhok)
    call sll_s_fft_init_c2r_1d(self%bw,nc_eta1,self%rhok,self%tmp)

  end subroutine sll_s_init_poisson_1d_hmf

  !> Solve the 1d equation for the Vlasov-HMF model
  subroutine sll_s_solve_poisson_1d_hmf(self, field, rhs)

    class(sll_t_poisson_1d_hmf),intent(inout) :: self
    sll_real64, dimension(:), intent(out)     :: field
    sll_real64, dimension(:), intent(in)      :: rhs

    SLL_ASSERT(size(rhs)==self%nc_eta1+1)

    self%tmp = rhs(1:self%nc_eta1)
    call sll_s_fft_exec_r2c_1d(self%fw, self%tmp, self%rhok)

    self%rhok(2)  = self%rhok(2) * sll_p_pi * sll_p_i1 ! We keep only one mode
    self%rhok(1)  = sll_p_i0
    self%rhok(3:) = sll_p_i0

    call sll_s_fft_exec_c2r_1d(self%bw, self%rhok, self%tmp)

    field(1:self%nc_eta1) = self%tmp / real(self%nc_eta1,f64)
    ! complete last term by periodicity
    field(self%nc_eta1+1) = field(1)

  end subroutine sll_s_solve_poisson_1d_hmf

  !> solves \f[ -\Delta \phi = \rho \f] in 2d
  subroutine compute_phi_from_rho_1d_hmf( poisson, phi, rho )
    class(sll_t_poisson_1d_hmf), target :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: phi
    
    stop ' compute phi from rho is not implemented '
    
  end subroutine compute_phi_from_rho_1d_hmf

  !> solves \f[ E = -\nabla \phi \f] with \f[ -\Delta \phi = \rho \f] in 2d 
  subroutine compute_E_from_rho_1d_hmf( poisson, E, rho )
    class(sll_t_poisson_1d_hmf) :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: E
      
    call sll_s_solve_poisson_1d_hmf(poisson, E, rho)
           
  end subroutine compute_E_from_rho_1d_hmf

  !> Free the poisson 1d hmf solver
  subroutine sll_s_free_poisson_1d_hmf(self, error)

    class(sll_t_poisson_1d_hmf) :: self
    sll_int32, intent(inout)    :: error    !< error code

    deallocate(self%rhok,stat=error)
    deallocate(self%tmp,stat=error)

  end subroutine sll_s_free_poisson_1d_hmf
  
end module sll_m_poisson_1d_hmf
