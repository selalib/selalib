!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations with coordinate transformation in 3D
!> The linear systems are solved using PLAF
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_preconditioner_fft
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"


  use sll_m_linear_solver_spline_mass_fft, only : &
       sll_t_linear_solver_spline_mass_fft

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_compute_mass_eig

  use sll_m_linear_solver_block, only : &
       sll_t_linear_solver_block

  implicit none

  public :: &
       sll_t_preconditioner_fft

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: sll_t_preconditioner_fft
     sll_real64 :: delta_x(3)       !< cell size
     sll_int32  :: n_dofs0(3)        !< number of Degrees of Freedom
     sll_int32  :: n_dofs1(3)        !< number of Degrees of Freedom 
     sll_int32  :: s_deg_0(3)       !< spline degree 0-forms
     sll_int32  :: s_deg_1(3)       !< spline degree 1-forms
     type(sll_t_linear_solver_spline_mass_fft) :: inverse_mass1_1d(3)  !< Fourier solver for 1-form mass matrix
     type(sll_t_linear_solver_spline_mass_fft) :: inverse_mass2_1d(3)  !< Fourier solver for 2-form mass matrix
     type(sll_t_linear_solver_block)           :: inverse_mass1_3d !< block matrix solver
     type(sll_t_linear_solver_block)           :: inverse_mass2_3d !< block matrix solver
     sll_real64, allocatable :: mass_line0_1(:) !< massline
     sll_real64, allocatable :: mass_line0_2(:) !< massline
     sll_real64, allocatable :: mass_line0_3(:) !< massline
     sll_real64, allocatable :: mass_line1_1(:) !< massline
     sll_real64, allocatable :: mass_line1_2(:) !< massline
     sll_real64, allocatable :: mass_line1_3(:) !< massline

   contains

     procedure :: &
          init => init_3d_trafo
     procedure :: &
          free => free_3d_trafo
  

  end type sll_t_preconditioner_fft

contains

 
  subroutine init_3d_trafo( self, Lx, n_cells, s_deg_0, clamped )
    class(sll_t_preconditioner_fft), intent( inout ) :: self      !< fft preconditioner
    sll_real64,                      intent( in    ) :: Lx(3)     !< domain length
    sll_int32,                       intent( in    ) :: n_cells(3) !< number of degrees of freedom (here number of cells and grid points)
    sll_int32,                       intent( in    ) :: s_deg_0(3)   !< highest spline degree
    logical, optional :: clamped !< logical for clamped boundary conditions
    ! local variables
    sll_int32  :: j
    
    sll_real64, allocatable :: eig_values_mass_0_1(:)
    sll_real64, allocatable :: eig_values_mass_0_2(:)
    sll_real64, allocatable :: eig_values_mass_0_3(:)
    sll_real64, allocatable :: eig_values_mass_1_1(:)
    sll_real64, allocatable :: eig_values_mass_1_2(:)
    sll_real64, allocatable :: eig_values_mass_1_3(:)


    self%delta_x = Lx / real(n_cells,f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1
    self%n_dofs0 = n_cells
    self%n_dofs1 = n_cells
       
    if( present(clamped) )then
       if(clamped) then
          self%n_dofs0(1) = n_cells(1)+self%s_deg_0(1)
          self%n_dofs1(1) = n_cells(1)+self%s_deg_1(1)
       end if
    end if
    ! Allocate scratch data
    allocate( self%mass_line0_1(s_deg_0(1)+1) )
    allocate( self%mass_line1_1(s_deg_0(1)) )
    allocate( self%mass_line0_2(s_deg_0(2)+1) )
    allocate( self%mass_line1_2(s_deg_0(2)) )
    allocate( self%mass_line0_3(s_deg_0(3)+1) )
    allocate( self%mass_line1_3(s_deg_0(3)) )
    call sll_s_spline_fem_mass_line ( self%s_deg_0(1), self%mass_line0_1 )
    call sll_s_spline_fem_mass_line ( self%s_deg_1(1), self%mass_line1_1 )

    call sll_s_spline_fem_mass_line ( self%s_deg_0(2), self%mass_line0_2 )
    call sll_s_spline_fem_mass_line ( self%s_deg_1(2), self%mass_line1_2 )

    call sll_s_spline_fem_mass_line ( self%s_deg_0(3), self%mass_line0_3 )
    call sll_s_spline_fem_mass_line ( self%s_deg_1(3), self%mass_line1_3 )

    self%mass_line0_1 = self%delta_x(1) * self%mass_line0_1 
    self%mass_line1_1 = self%delta_x(1) * self%mass_line1_1  

    self%mass_line0_2 = self%delta_x(2) * self%mass_line0_2
    self%mass_line1_2 = self%delta_x(2) * self%mass_line1_2 

    self%mass_line0_3 = self%delta_x(3) * self%mass_line0_3
    self%mass_line1_3 = self%delta_x(3) * self%mass_line1_3 

    allocate( eig_values_mass_0_1( self%n_dofs0(1) ) )
    allocate( eig_values_mass_0_2( self%n_dofs0(2) ) )
    allocate( eig_values_mass_0_3( self%n_dofs0(3) ) )
    allocate( eig_values_mass_1_1( self%n_dofs1(1) ) )
    allocate( eig_values_mass_1_2( self%n_dofs1(2) ) )
    allocate( eig_values_mass_1_3( self%n_dofs1(3) ) )
    
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs0(1), self%s_deg_0(1), self%mass_line0_1, &
         eig_values_mass_0_1 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs0(2), self%s_deg_0(2), self%mass_line0_2, &
         eig_values_mass_0_2 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs0(3), self%s_deg_0(3), self%mass_line0_3, &
         eig_values_mass_0_3 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs1(1), self%s_deg_1(1), self%mass_line1_1, &
         eig_values_mass_1_1 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs1(2), self%s_deg_1(2), self%mass_line1_2, &
         eig_values_mass_1_2 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs1(3), self%s_deg_1(3), self%mass_line1_3, &
         eig_values_mass_1_3 )

    eig_values_mass_0_1 = 1._f64/eig_values_mass_0_1
    eig_values_mass_0_2 = 1._f64/eig_values_mass_0_2
    eig_values_mass_0_3 = 1._f64/eig_values_mass_0_3
    eig_values_mass_1_1 = 1._f64/eig_values_mass_1_1
    eig_values_mass_1_2 = 1._f64/eig_values_mass_1_2
    eig_values_mass_1_3 = 1._f64/eig_values_mass_1_3
        
    call self%inverse_mass1_1d(1)%create( self%n_dofs1, eig_values_mass_1_1, eig_values_mass_0_2, eig_values_mass_0_3 )
    call self%inverse_mass1_1d(2)%create( self%n_dofs0, eig_values_mass_0_1, eig_values_mass_1_2, eig_values_mass_0_3 )
    call self%inverse_mass1_1d(3)%create( self%n_dofs0, eig_values_mass_0_1, eig_values_mass_0_2, eig_values_mass_1_3 )

    call self%inverse_mass2_1d(1)%create( self%n_dofs0, eig_values_mass_0_1, eig_values_mass_1_2, eig_values_mass_1_3 )
    call self%inverse_mass2_1d(2)%create( self%n_dofs1, eig_values_mass_1_1, eig_values_mass_0_2, eig_values_mass_1_3 )
    call self%inverse_mass2_1d(3)%create( self%n_dofs1, eig_values_mass_1_1, eig_values_mass_1_2, eig_values_mass_0_3 )

    call self%inverse_mass1_3d%create( 3, 3 )
    call self%inverse_mass2_3d%create( 3, 3 )
    do j=1,3
       call self%inverse_mass1_3d%set( j, j, self%inverse_mass1_1d(j) )
       call self%inverse_mass2_3d%set( j, j, self%inverse_mass2_1d(j) )
    end do


  end subroutine init_3d_trafo



  subroutine free_3d_trafo( self )
    class(sll_t_preconditioner_fft) :: self !< fft preconditioner
    !local variable
    sll_int32 :: j

    do j=1, 3
       call self%inverse_mass1_1d(j)%free()
       call self%inverse_mass2_1d(j)%free()
    end do

    call self%inverse_mass1_3d%free()
    call self%inverse_mass2_3d%free()

  end subroutine free_3d_trafo


end module sll_m_preconditioner_fft
