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
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> Yaman Güçlü (yaman.guclu@gmail.com)
!> @brief 
!> the plan is to adapt former 
!> simulation_4d_drift_kinetic_field_aligned_polar
!> to curvilinear case
!> for further comments, see sll_docs/fcisl/fcisl_report.pdf

module sll_m_sim_bsl_dk_3d1v_curv_field_aligned

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"

  use sll_m_cartesian_meshes, only: &
    sll_cartesian_mesh_1d, &
    sll_cartesian_mesh_2d

  use sll_m_collective, only: &
    sll_get_collective_rank, &
    sll_world_collective

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_coordinate_transformation_2d_base

  use sll_m_dk_curv_mesh, only: &
    init_dk_curv_mesh

  use sll_m_sim_base, only: &
    sll_simulation_base_class

  use sll_m_utilities, only: &
    sll_new_file_id

  implicit none

  public :: &
    sll_t_sim_sl_dk_3d1v_curv_field_aligned

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_simulation_base_class) :: &
    sll_t_sim_sl_dk_3d1v_curv_field_aligned

    !mesh
    type(sll_cartesian_mesh_2d), pointer :: m_x1x2
    type(sll_cartesian_mesh_1d), pointer :: m_x3
    type(sll_cartesian_mesh_1d), pointer :: m_x4
    class(sll_coordinate_transformation_2d_base), pointer :: transformation



    contains

      procedure :: init_from_file => init_sim_sl_dk_3d1v_curv_field_aligned
      procedure :: run            => run_sim_sl_dk_3d1v_curv_field_aligned


  end type sll_t_sim_sl_dk_3d1v_curv_field_aligned


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !----------------------------------------------------------------------------
  !> Initialize simulation from input file
  subroutine init_sim_sl_dk_3d1v_curv_field_aligned( sim, filename )
    class(sll_t_sim_sl_dk_3d1v_curv_field_aligned), &
                      intent(inout) :: sim
    character(len=*), intent(in)    :: filename

    call init_dk_curv_mesh( &
      filename, &
      sim%m_x1x2, &
      sim%m_x3, &
      sim%m_x4, &
      sim%transformation, &
      sll_get_collective_rank(sll_world_collective))
    
  end subroutine init_sim_sl_dk_3d1v_curv_field_aligned

  !----------------------------------------------------------------------------
  !> Run simulation
  subroutine run_sim_sl_dk_3d1v_curv_field_aligned(sim)
    class(sll_t_sim_sl_dk_3d1v_curv_field_aligned), intent(inout) :: sim
  end subroutine run_sim_sl_dk_3d1v_curv_field_aligned
  
  

end module sll_m_sim_bsl_dk_3d1v_curv_field_aligned
