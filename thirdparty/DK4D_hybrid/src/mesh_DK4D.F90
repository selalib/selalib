!===========================================================================
!> Initialization of the mesh for 4D drift-kinetic hybrid
!>  simulation
!>
!> \date 2014-08-19
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module mesh_DK4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use input_DK4D_module
  use sll_m_cartesian_meshes
  use sll_m_coordinate_transformation_2d_base
  use sll_m_common_coordinate_transformations
  use sll_m_coordinate_transformations_2d

  implicit none

  !===========================================================================
  !> Physical mesh in (x,y,eta3,vpar)
  !---------------------------------------------------------------------------
  type, public :: mesh_DK4D_t

    !> Geometry type
    sll_int32 :: geometry_type

    !> Major radius R0
    sll_real64 :: major_radius

    !> 2D mesh in (eta1,eta2) direction
    type(sll_t_cartesian_mesh_2d), pointer :: eta1_eta2_mesh2d
    sll_real64, dimension(:) , pointer :: eta1_grid
    sll_real64, dimension(:) , pointer :: eta2_grid

    !> 1D mesh in eta3 direction
    type(sll_t_cartesian_mesh_1d), pointer :: eta3_mesh1d
    sll_real64, dimension(:) , pointer :: eta3_grid

    !> 1D mesh in vpar direction
    type(sll_t_cartesian_mesh_1d), pointer :: vpar_mesh1d
    sll_real64, dimension(:) , pointer :: vpar_grid

    !> Coordinate transformations F
    class(sll_c_coordinate_transformation_2d_base), pointer :: &
        transf_eta1eta2_xy

    !> 2D (x,y) mesh
    sll_real64, dimension(:,:), pointer :: xgrid_2d
    sll_real64, dimension(:,:), pointer :: ygrid_2d
    !> rgrid_2d = sqrt(xgrid_2d**2+ygrid_2d**2) 
    sll_real64, dimension(:,:), pointer :: rgrid_2d

  end type mesh_DK4D_t
  !---------------------------------------------------------------------------

  sll_int32, parameter, public :: CIRCULAR_GEOMETRY = 1

contains

  !===========================================================================
  !> Mesh: Allocation
  !---------------------------------------------------------------------------
  subroutine new_mesh_DK4D( mesh4d, &
      num_cells_x1, num_cells_x2, num_cells_x3, num_cells_x4, &
      eta1_min, eta1_max, &
      eta2_min, eta2_max, &
      eta3_min, eta3_max, &
      vpar_min, vpar_max, &
      scheme_type )

    type(mesh_DK4D_t) , intent(inout) :: mesh4d
    sll_int32 , intent(in) :: num_cells_x1
    sll_int32 , intent(in) :: num_cells_x2
    sll_int32 , intent(in) :: num_cells_x3
    sll_int32 , intent(in) :: num_cells_x4
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_real64, intent(in) :: eta3_min
    sll_real64, intent(in) :: eta3_max
    sll_real64, intent(in) :: vpar_min
    sll_real64, intent(in) :: vpar_max
    sll_int32 , intent(in) :: scheme_type

    sll_int32 :: Nx, Ny
    sll_int32 :: ierr

    !*** Initialization of the geometry type ***
    mesh4d%geometry_type = CIRCULAR_GEOMETRY

    !*** Initialization 2D mesh in (eta1,eta2) direction and ***
    !***  initialization of the coordinate transformation    ***
    !***  from (eta1,eta2) to (x,y)                          ***
    select case (scheme_type)

    case ( BSL_HYBRID )
      !--> 2D mesh in (eta1,eta2)
      mesh4d%eta1_eta2_mesh2d => sll_f_new_cartesian_mesh_2d( &
          num_cells_x1, &
          num_cells_x2, &
          eta1_min=0._f64, &
          eta1_max=1._f64, &
          eta2_min=0._f64, &
          eta2_max=1._f64 )
      !--> coordinate transformation from (eta1,eta2) to (x,y)
      mesh4d%transf_eta1eta2_xy => sll_f_new_coordinate_transformation_2d_analytic( &
          "polar_transformation", &
          mesh4d%eta1_eta2_mesh2d, &
          sll_f_x1_polar_f, &
          sll_f_x2_polar_f, &
          sll_f_deriv_x1_polar_f_eta1, &
          sll_f_deriv_x1_polar_f_eta2, &
          sll_f_deriv_x2_polar_f_eta1, &
          sll_f_deriv_x2_polar_f_eta2, &
          (/eta1_min,eta1_max/) )
      
    case ( BSL_POLAR )
      !--> 2D mesh in (eta1,eta2)
      mesh4d%eta1_eta2_mesh2d => sll_f_new_cartesian_mesh_2d( &
          num_cells_x1, &
          num_cells_x2, &
          eta1_min=eta1_min, &
          eta1_max=eta1_max, &
          eta2_min=eta2_min, &
          eta2_max=eta2_max )
      !--> coordinate transformation from (eta1,eta2) to (x,y)
      mesh4d%transf_eta1eta2_xy => sll_f_new_coordinate_transformation_2d_analytic( &
          "polar_transformation", &
          mesh4d%eta1_eta2_mesh2d, &
          sll_f_polar_x1, &
          sll_f_polar_x2, &
          sll_f_deriv_x1_polar_f_eta1, &
          sll_f_deriv_x1_polar_f_eta2, &
          sll_f_deriv_x2_polar_f_eta1, &
          sll_f_deriv_x2_polar_f_eta2, &
          (/0._f64,1._f64/) )

    case default
      print*, ' Scheme type = ', scheme_type, ' not treated '
      stop
      call exit(-1)

    end select

    SLL_ALLOCATE( mesh4d%eta1_grid(num_cells_x1+1), ierr )
    SLL_ALLOCATE( mesh4d%eta2_grid(num_cells_x2+1), ierr )

    !--> Writing of the transformation in the output file 
    !-->  'polar_transformation.xmf'
    call  mesh4d%transf_eta1eta2_xy%write_to_file(0)

    !*** 2D mesh in (x,y) direction ***
    Nx = num_cells_x1 + 1
    Ny = num_cells_x2 + 1
    SLL_ALLOCATE( mesh4d%xgrid_2d(Nx,Ny), ierr )
    SLL_ALLOCATE( mesh4d%ygrid_2d(Nx,Ny), ierr )
    SLL_ALLOCATE( mesh4d%rgrid_2d(Nx,Ny), ierr )

    !*** 1D mesh in eta3 direction ***
    mesh4d%eta3_mesh1d => sll_f_new_cartesian_mesh_1d( &
        num_cells_x3, &
        eta_min=eta3_min, &
        eta_max=eta3_max )
    SLL_ALLOCATE( mesh4d%eta3_grid(num_cells_x3+1), ierr )

    !*** 1D mesh in vpar direction ***
    mesh4d%vpar_mesh1d => sll_f_new_cartesian_mesh_1d( &
        num_cells_x4, &
        eta_min=vpar_min, &
        eta_max=vpar_max )
    SLL_ALLOCATE( mesh4d%vpar_grid(num_cells_x4+1), ierr )
    
  end subroutine new_mesh_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Mesh: Deallocation
  !---------------------------------------------------------------------------
  subroutine delete_mesh_DK4D( mesh4d )

    type(mesh_DK4D_t), intent(inout) :: mesh4d
    
    sll_int32 :: ierr

    call sll_o_delete( mesh4d%eta1_eta2_mesh2d )
    SLL_DEALLOCATE( mesh4d%eta1_grid, ierr )
    SLL_DEALLOCATE( mesh4d%eta2_grid, ierr )

    SLL_DEALLOCATE( mesh4d%xgrid_2d, ierr )
    SLL_DEALLOCATE( mesh4d%ygrid_2d, ierr )
    SLL_DEALLOCATE( mesh4d%rgrid_2d, ierr )

    call sll_o_delete( mesh4d%eta3_mesh1d )
    SLL_DEALLOCATE( mesh4d%eta3_grid, ierr )

    call sll_o_delete( mesh4d%vpar_mesh1d )
    SLL_DEALLOCATE( mesh4d%vpar_grid, ierr )
    
  end subroutine delete_mesh_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Mesh: Initialization
  !---------------------------------------------------------------------------
  subroutine init_mesh_DK4D( mesh4d )

    use sll_m_constants

    type(mesh_DK4D_t), intent(inout) :: mesh4d

    sll_int32 :: ieta1, ieta2, ieta3, ivpar
    sll_int32 :: Neta1, Neta2, Neta3, Nvpar
    sll_real64 :: x_point, y_point
    sll_real64 :: Lz

    Neta1 = size(mesh4d%eta1_grid)
    Neta2 = size(mesh4d%eta2_grid)
    Neta3 = size(mesh4d%eta3_grid)
    Nvpar = size(mesh4d%vpar_grid)

    !-> Initialization of the grid in eta1 direction
    do ieta1 = 1,Neta1
      mesh4d%eta1_grid(ieta1) = mesh4d%eta1_eta2_mesh2d%eta1_min + &
          (ieta1-1)*mesh4d%eta1_eta2_mesh2d%delta_eta1
    end do

    !-> Initialization of the grid in eta2 direction
    do ieta2 = 1,Neta2
      mesh4d%eta2_grid(ieta2) = mesh4d%eta1_eta2_mesh2d%eta2_min + &
          (ieta2-1)*mesh4d%eta1_eta2_mesh2d%delta_eta2
    end do

    !-> Initialization of the grid in eta3 direction
    do ieta3 = 1,Neta3
      mesh4d%eta3_grid(ieta3) = mesh4d%eta3_mesh1d%eta_min + &
          (ieta3-1)*mesh4d%eta3_mesh1d%delta_eta
    end do

    !-> Initialization of the grid in vpar direction
    do ivpar = 1,Nvpar
      mesh4d%vpar_grid(ivpar) = mesh4d%vpar_mesh1d%eta_min + &
          (ivpar-1)*mesh4d%vpar_mesh1d%delta_eta
    end do

    !-> Initialization of the grid in (x,y) direction depending 
    !->  of the transformation
    do ieta2 = 1,Neta2
      do ieta1 = 1,Neta1
        x_point = mesh4d%transf_eta1eta2_xy%x1_at_node(ieta1,ieta2)
        y_point = mesh4d%transf_eta1eta2_xy%x2_at_node(ieta1,ieta2)
        mesh4d%xgrid_2d(ieta1,ieta2) = x_point            
        mesh4d%ygrid_2d(ieta1,ieta2) = y_point
        mesh4d%rgrid_2d(ieta1,ieta2) = sqrt(x_point**2 + y_point**2)
      end do
    end do

    !--> Initialisation of R0 = major_radius
    Lz = abs(mesh4d%eta3_mesh1d%eta_max - &
        mesh4d%eta3_mesh1d%eta_min)
    mesh4d%major_radius = Lz/(2._f64*sll_p_pi)

  end subroutine init_mesh_DK4D


  !===========================================================================
  !> Mesh: Print in HDF5
  !---------------------------------------------------------------------------
  subroutine print_mesh_DK4D( mesh4d, filename )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
         sll_o_hdf5_ser_write_array, &
         sll_s_hdf5_ser_file_close, &
         sll_t_hdf5_ser_handle

    type(mesh_DK4D_t), intent(in) :: mesh4d
    character(len=*) , intent(in) :: filename

    sll_int32 :: my_rank

    !--> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle

    my_rank = sll_f_get_collective_rank(sll_v_world_collective)
    if (my_rank.eq.0) then
      call sll_s_hdf5_ser_file_create( trim(filename), handle, file_err )
      call sll_o_hdf5_ser_write_array( handle, mesh4d%eta1_grid, &
          'eta1_grid', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%eta2_grid, &
          'eta2_grid', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%eta3_grid, &
          'eta3_grid', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%vpar_grid, &
          'vpar_grid', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%xgrid_2d, &
          'xgrid_2d', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%ygrid_2d, &
          'ygrid_2d', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%rgrid_2d, &
          'rgrid_2d', file_err)
      call sll_s_hdf5_ser_file_close( handle, file_err )
    end if

  end subroutine print_mesh_DK4D
  !---------------------------------------------------------------------------

end module mesh_DK4D_module
!---------------------------------------------------------------------------
