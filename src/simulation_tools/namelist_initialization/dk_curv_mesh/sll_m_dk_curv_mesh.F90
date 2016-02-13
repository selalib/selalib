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
!> @brief 
!> initialization from dk_curv_mesh namelist 
!> this module is then used in sl_dk_3d1v_curv_field_aligned
!> we hope to be able to factorize namelists initializations
!> to other simulations
!> and prevent from too long initializations of a simulation

module sll_m_dk_curv_mesh
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_1d, &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_2d

  use sll_m_common_coordinate_transformations, only: &
    sll_f_sharped_geo_jac11, &
    sll_f_sharped_geo_jac12, &
    sll_f_sharped_geo_jac21, &
    sll_f_sharped_geo_jac22, &
    sll_f_sharped_geo_x1, &
    sll_f_sharped_geo_x2, &
    sll_f_polar_jac11, &
    sll_f_polar_jac12, &
    sll_f_polar_jac21, &
    sll_f_polar_jac22, &
    sll_f_polar_shear_jac11, &
    sll_f_polar_shear_jac12, &
    sll_f_polar_shear_jac21, &
    sll_f_polar_shear_jac22, &
    sll_f_polar_shear_x1, &
    sll_f_polar_shear_x2, &
    sll_f_polar_x1, &
    sll_f_polar_x2

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_utilities, only: &
    sll_s_new_file_id

  implicit none

  public :: &
    sll_s_init_dk_curv_mesh

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



contains
  !----------------------------------------------------------------------------
  !> compute num_params_mesh
  function compute_num_params_mesh(mesh_case) result(res)
    character(len=*), intent(in) :: mesh_case
    sll_int32 :: res
    character(len=256) :: err_msg

    select case (mesh_case)
      case ("SLL_POLAR_MESH")
        res = 1
      case ("SLL_POLAR_SHEAR_MESH")
        res = 4
      case ("SLL_D_SHAPED_MESH")
        res = 9
      case default
       err_msg = 'bad mesh_case: '//trim(mesh_case) 
       SLL_ERROR( 'compute_num_params_mesh', trim( err_msg ))
    end select
    
  end function compute_num_params_mesh

  !----------------------------------------------------------------------------
  !> compute params_mesh

  subroutine compute_params_mesh( &
    mesh_case, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    mesh_alpha, &
    params_mesh, &
    num_params_mesh )
    character(len=*), intent(in) :: mesh_case
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_real64, intent(in) :: mesh_alpha(:)
    sll_real64, intent(out) :: params_mesh(:)
    sll_int32, intent(in) :: num_params_mesh
    
    character(len=256) :: err_msg
    
    if(size(params_mesh)<num_params_mesh)then
      err_msg = 'bad size for params_mesh: '
      SLL_ERROR( 'compute_params_mesh', trim( err_msg ))      
    endif

    params_mesh(1:num_params_mesh) = 0._f64

    select case (mesh_case)
      case ("SLL_POLAR_MESH")
      case ("SLL_POLAR_SHEAR_MESH")
        params_mesh(1) = mesh_alpha(1)
        params_mesh(2) = mesh_alpha(2)        
      case ("SLL_D_SHAPED_MESH")
        params_mesh(1:5) = mesh_alpha(1:5)
        params_mesh(6) = eta1_min
        params_mesh(7) = eta1_max
        params_mesh(8) = eta2_min
        params_mesh(9) = eta2_max
      case default
       err_msg = 'bad mesh_case: '//trim(mesh_case) 
       SLL_ERROR( 'compute_params_mesh', trim( err_msg ))
    end select
    
  end subroutine compute_params_mesh   

  !----------------------------------------------------------------------------
  !> select transformation

  subroutine select_transformation( &
    mesh_case, &
    mesh_2d, &
    params_mesh, &
    num_params_mesh, &
    transformation)
    character(len=*), intent(in) :: mesh_case
    type(sll_t_cartesian_mesh_2d), pointer :: mesh_2d
    sll_real64, intent(out) :: params_mesh(:)
    sll_int32, intent(in) :: num_params_mesh
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation

    character(len=256) :: err_msg

    if(size(params_mesh)<num_params_mesh)then
      err_msg = 'bad size for params_mesh: '
      SLL_ERROR( 'select_transformation', trim( err_msg ))      
    endif


    select case (mesh_case)
      case ("SLL_POLAR_MESH") 
        transformation => sll_f_new_coordinate_transformation_2d_analytic( &
          "analytic_polar_transformation", &
          mesh_2d, &
          sll_f_polar_x1, &
          sll_f_polar_x2, &
          sll_f_polar_jac11, &
          sll_f_polar_jac12, &
          sll_f_polar_jac21, &
          sll_f_polar_jac22, &
          params_mesh  )     
      case ("SLL_POLAR_SHEAR_MESH") 
        transformation => sll_f_new_coordinate_transformation_2d_analytic( &
          "analytic_polar_shear_transformation", &
          mesh_2d, &
          sll_f_polar_shear_x1, &
          sll_f_polar_shear_x2, &
          sll_f_polar_shear_jac11, &
          sll_f_polar_shear_jac12, &
          sll_f_polar_shear_jac21, &
          sll_f_polar_shear_jac22, &
          params_mesh )     
      case ("SLL_D_SHAPED_MESH")
        transformation => sll_f_new_coordinate_transformation_2d_analytic( &
          "analytic_D_SHAPED_transformation", &
          mesh_2d, &
          sll_f_sharped_geo_x1, &
          sll_f_sharped_geo_x2, &
          sll_f_sharped_geo_jac11, &
          sll_f_sharped_geo_jac12, &
          sll_f_sharped_geo_jac21, &
          sll_f_sharped_geo_jac22, &
          params_mesh )    
      case default
        err_msg = 'bad mesh_case: '//trim(mesh_case) 
        SLL_ERROR( 'select_transformation', trim( err_msg ))
    end select  
  
  end subroutine select_transformation


  !----------------------------------------------------------------------------
  !> Initialize simulation from input file
  subroutine sll_s_init_dk_curv_mesh( &
    filename, &
    m_x1x2, &
    m_x3, &
    m_x4, &
    transformation, &
    proc_id)
    character(len=*), intent(in)    :: filename
    type(sll_t_cartesian_mesh_2d), pointer :: m_x1x2
    type(sll_t_cartesian_mesh_1d), pointer :: m_x3
    type(sll_t_cartesian_mesh_1d), pointer :: m_x4
    class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
    sll_int32, intent(in) :: proc_id
    
    !--> mesh
    character(len=256) :: mesh_case !< poloidal plane
    sll_int32 :: num_cells_eta1
    sll_int32 :: num_cells_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: mesh_alpha1
    sll_real64 :: mesh_alpha2
    sll_real64 :: mesh_alpha3
    sll_real64 :: mesh_alpha4
    sll_real64 :: mesh_alpha5
    sll_real64 :: mesh_alpha6
    sll_int32  :: num_cells_x3  !< z
    sll_int32  :: num_cells_x4  !< v
    sll_real64 :: z_min
    sll_real64 :: z_max
    sll_real64 :: v_min
    sll_real64 :: v_max
    
    !--> local variables
    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    sll_real64, allocatable :: params_mesh(:)
    sll_int32 :: num_params_mesh
    sll_real64 :: mesh_alpha(6)
    

    namelist /mesh/ &
      mesh_case, &
      num_cells_eta1, &
      num_cells_eta2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      mesh_alpha1, &
      mesh_alpha2, &
      mesh_alpha3, &
      mesh_alpha4, &
      mesh_alpha5, &
      mesh_alpha6, &
      num_cells_x3, &
      num_cells_x4, &
      z_min, &
      z_max, &
      v_min, &
      v_max

    !> first default values
    
    mesh_case = "SLL_POLAR_MESH"
    num_cells_eta1 = 32
    num_cells_eta2 = 64
    eta1_min = 0.1_f64
    eta1_max = 14.5_f64
    eta2_min = 0._f64
    eta2_max = 6.283185307179586477_f64
    mesh_alpha1 = 0._f64
    mesh_alpha2 = 0._f64
    mesh_alpha3 = 0._f64
    mesh_alpha4 = 0._f64
    mesh_alpha5 = 0._f64
    mesh_alpha6 = 0._f64
    num_cells_x3 = 32
    num_cells_x4 = 64
    z_min = 0.0_f64
    z_max = 1506.759067_f64
    v_min = -7.32_f64
    v_max = 7.32_f64
    
    call sll_s_new_file_id(namelist_id, ierr)
    open(unit = namelist_id, file=trim(filename)//'.nml',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       err_msg = 'failed to open file '//trim(filename)//'.nml' 
       SLL_ERROR( 'sll_s_init_dk_curv_mesh', trim( err_msg ))
    end if
    read(namelist_id, mesh) 


    if(proc_id==0)then
      print *,'#Mesh:'
      print *,'#mesh_case=',trim(mesh_case)
      print *,'#num_cells_eta1=',num_cells_eta1
      print *,'#num_cells_eta2=',num_cells_eta2
      print *,'#eta1_min=',eta1_min
      print *,'#eta1_max=',eta1_max
      print *,'#eta2_min=',eta2_min
      print *,'#eta2_max=',eta2_max
      print *,'#mesh_alpha1=',mesh_alpha1
      print *,'#mesh_alpha2=',mesh_alpha2
      print *,'#mesh_alpha3=',mesh_alpha3
      print *,'#mesh_alpha4=',mesh_alpha4
      print *,'#mesh_alpha5=',mesh_alpha5
      print *,'#mesh_alpha6=',mesh_alpha6
      print *,'#num_cells_x3=',num_cells_x3
      print *,'#num_cells_x4=',num_cells_x4
      print *,'#z_min=',z_min
      print *,'#z_max=',z_max
      print *,'#v_min=',v_min
      print *,'#v_max=',v_max
    endif    

    m_x1x2 => sll_f_new_cartesian_mesh_2d( &
      num_cells_eta1, &
      num_cells_eta2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max ) 
    m_x3 => sll_f_new_cartesian_mesh_1d(num_cells_x3,eta_min=z_min,eta_max=z_max)
    m_x4 => sll_f_new_cartesian_mesh_1d(num_cells_x4,eta_min=v_min,eta_max=v_max)
    
    num_params_mesh =  compute_num_params_mesh(mesh_case)
    SLL_ALLOCATE(params_mesh(num_params_mesh),ierr)    
    mesh_alpha(1:3) = (/mesh_alpha1,mesh_alpha2,mesh_alpha3/)
    mesh_alpha(4:6) = (/mesh_alpha4,mesh_alpha5,mesh_alpha6/)
    
    call compute_params_mesh( &
      mesh_case, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      mesh_alpha, &
      params_mesh, &
      num_params_mesh )
    
    call select_transformation( &
      mesh_case, &
      m_x1x2, &
      params_mesh, &
      num_params_mesh, &
      transformation)


    
    
  end subroutine sll_s_init_dk_curv_mesh




end module sll_m_dk_curv_mesh 
