!===========================================================================
!> Initialization of the mesh for 4D Vlasov-Poisson hybrid
!>  simulation
!>
!> \date 2015-02-26
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module mesh_VP4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use input_VP4D_module
  use sll_m_cartesian_meshes
  use sll_m_coordinate_transformation_2d_base
  use sll_m_common_coordinate_transformations
  use sll_m_coordinate_transformations_2d
  use sll_m_coordinate_transformations_2d_nurbs
  use sll_m_utilities, only : sll_s_new_file_id

  implicit none

  !===========================================================================
  !> Physical mesh in (x,y,eta3,vpar)
  !---------------------------------------------------------------------------
  type, public :: mesh_VP4D_t

    !> Geometry type
    sll_int32 :: mapping_type
    sll_int32 :: analytical_formula

    !> 2D mesh in (eta1,eta2) direction
    type(sll_t_cartesian_mesh_2d), pointer :: eta1_eta2_mesh2d
    sll_real64, dimension(:) , pointer :: eta1_grid
    sll_real64, dimension(:) , pointer :: eta2_grid
    sll_real64 :: eta1_min, eta1_max, eta1_length
    sll_real64 :: eta2_min, eta2_max, eta2_length

    !> 2D (x,y) physical mesh
    sll_real64, dimension(:,:), pointer :: xgrid_2d
    sll_real64, dimension(:,:), pointer :: ygrid_2d

    !> 2D mesh in (vx,vy) direction
    type(sll_t_cartesian_mesh_2d), pointer :: vx_vy_mesh2d
    sll_real64, dimension(:)   , pointer :: vx_grid
    sll_real64, dimension(:)   , pointer :: vy_grid
    sll_real64 :: vx_min, vx_max, vx_length
    sll_real64 :: vy_min, vy_max, vy_length

    !> Coordinate transformations F
    character(len=80) :: mesh_filename
    class(sll_c_coordinate_transformation_2d_base), pointer :: &
        transf_eta1eta2_xy

    !> Jacobian matrix and its inverse
    sll_real64, dimension(:,:,:,:), pointer :: Jacobian_matrix
    sll_real64, dimension(:,:,:,:), pointer :: inv_Jacobian_matrix

    !> Jacobian = determinant of the Jacobian matrix
    sll_real64, dimension(:,:), pointer :: Jacobian

    !> Used for integration in all directions
    sll_real64, dimension(:), pointer :: coef_int_deta1
    sll_real64, dimension(:), pointer :: coef_int_deta2
    sll_real64, dimension(:), pointer :: coef_int_dvx
    sll_real64, dimension(:), pointer :: coef_int_dvy

    !> Mapping file name (only used when the mesh is defined with CAID)
    character(len=80) :: mapping_filename

  end type mesh_VP4D_t
  !---------------------------------------------------------------------------

  sll_real64, private, parameter :: unused_ = 0._f64

contains

  !===========================================================================
  !> Mesh: Allocation
  !---------------------------------------------------------------------------
  subroutine new_mesh_VP4D( mesh4d, &
      num_cells_x1, num_cells_x2, num_cells_x3, num_cells_x4, &
      eta1_min, eta1_max, &
      eta2_min, eta2_max, &
      vx_min, vx_max, &
      vy_min, vy_max, &
      mapping_type, &
      analytical_formula, &
      mapping_filename, &
      colella_coeff )

    type(mesh_VP4D_t), intent(inout) :: mesh4d
    sll_int32        , intent(in) :: num_cells_x1
    sll_int32        , intent(in) :: num_cells_x2
    sll_int32        , intent(in) :: num_cells_x3
    sll_int32        , intent(in) :: num_cells_x4
    sll_real64       , intent(in) :: eta1_min
    sll_real64       , intent(in) :: eta1_max
    sll_real64       , intent(in) :: eta2_min
    sll_real64       , intent(in) :: eta2_max
    sll_real64       , intent(in) :: vx_min
    sll_real64       , intent(in) :: vx_max
    sll_real64       , intent(in) :: vy_min
    sll_real64       , intent(in) :: vy_max
    sll_int32        , intent(in) :: mapping_type
    sll_int32        , intent(in) :: analytical_formula
    character(len=80), intent(in) :: mapping_filename
    sll_real64       , intent(in) :: colella_coeff

    !--> Local variables
    sll_int32  :: Nx, Ny, Neta1, Neta2
    sll_int32  :: Nvx, Nvy
    sll_int32  :: ierr

    !*** Initialization of the mesh and the coordinate transformation***
    mesh4d%mapping_type       = mapping_type
    mesh4d%analytical_formula = analytical_formula

    select case ( mesh4d%mapping_type )
    case ( ANALYTICAL )
      mesh4d%eta1_eta2_mesh2d => sll_f_new_cartesian_mesh_2d( &
          num_cells_x1, &
          num_cells_x2, &
          eta1_min=eta1_min, &
          eta1_max=eta1_max, &
          eta2_min=eta2_min, &
          eta2_max=eta2_max )
      
      select case ( mesh4d%analytical_formula )
        !--> Initialization 2D mesh in (eta1,eta2) direction 

      case ( IDENTITY )
        print*,'====> Analytical formula: IDENTITY mapping'
        !--> Initialization of the coordinate transformation
        !-->   from (eta1,eta2) to (x,y)
        mesh4d%transf_eta1eta2_xy => &
            sll_f_new_coordinate_transformation_2d_analytic( &
            "identity_transformation", &
            mesh4d%eta1_eta2_mesh2d, &
            sll_f_identity_x1, &
            sll_f_identity_x2, &
            sll_f_identity_jac11, &
            sll_f_identity_jac12, &
            sll_f_identity_jac21, &
            sll_f_identity_jac22, &
            (/0.05_f64/) )
        !\todo: Understand why 0.05 ?
      
      case ( COLELLA )
        print*,'====> Analytical formula: COLELLA mapping'
        mesh4d%transf_eta1eta2_xy => &
            sll_f_new_coordinate_transformation_2d_analytic( &
            "alpha_0.1", &
            mesh4d%eta1_eta2_mesh2d, &
            sll_f_sinprod_x1, &
            sll_f_sinprod_x2, &
            sll_f_sinprod_jac11, &
            sll_f_sinprod_jac12, &
            sll_f_sinprod_jac21, &
            sll_f_sinprod_jac22, &
            (/ colella_coeff,colella_coeff,eta1_max,eta2_max /) )
      
      case default
        print*, '====> Analytical formula ', &
            mesh4d%analytical_formula, ' not treated'
      end select

    case ( CAID_FILE )
      print*,'====>  CAID mapping using file: ',mapping_filename
      mesh4d%transf_eta1eta2_xy => &
          sll_f_new_nurbs_2d_transformation_from_file(trim(mapping_filename))
      mesh4d%eta1_eta2_mesh2d => mesh4d%transf_eta1eta2_xy%get_cartesian_mesh()
      write(mesh4d%mapping_filename,'(A)') mapping_filename
      call CAID_collocpoints_modif_( mesh4d )
     
    case default
      print*, ' Mapping type = ', mapping_type, ' not treated'
      stop
      call exit(-1)
    end select

    !--> Writing of the transformation in the output file 
    !-->  '<mapping>_transformation.xmf'
    call  mesh4d%transf_eta1eta2_xy%write_to_file(sll_p_io_mtv)

    Neta1  = mesh4d%eta1_eta2_mesh2d%num_cells1+1
    Neta2  = mesh4d%eta1_eta2_mesh2d%num_cells2+1
    Nx     = Neta1
    Ny     = Neta2

    !*** Initialization 2D mesh in (x,y) direction ***
    SLL_ALLOCATE( mesh4d%eta1_grid(Nx), ierr )
    SLL_ALLOCATE( mesh4d%eta2_grid(Ny), ierr )
    SLL_ALLOCATE( mesh4d%xgrid_2d(Nx,Ny), ierr )
    SLL_ALLOCATE( mesh4d%ygrid_2d(Nx,Ny), ierr )

    !*** 2D mesh in (vx,vy) direction ***
    mesh4d%vx_vy_mesh2d => sll_f_new_cartesian_mesh_2d( &
        num_cells_x3, &
        num_cells_x4, &
        eta1_min=vx_min, &
        eta1_max=vx_max, &
        eta2_min=vy_min, &
        eta2_max=vy_max )
    Nvx = mesh4d%vx_vy_mesh2d%num_cells1+1
    Nvy = mesh4d%vx_vy_mesh2d%num_cells2+1
    
    SLL_ALLOCATE( mesh4d%vx_grid(Nvx), ierr )
    SLL_ALLOCATE( mesh4d%vy_grid(Nvy), ierr )

    !*** Allocation for the Jacobian of the transformation ***
    SLL_CLEAR_ALLOCATE( mesh4d%Jacobian_matrix(Neta1,Neta2,2,2), ierr )
    SLL_CLEAR_ALLOCATE( mesh4d%inv_Jacobian_matrix(Neta1,Neta2,2,2), ierr )
    SLL_CLEAR_ALLOCATE( mesh4d%Jacobian(Neta1,Neta2), ierr )

    !*** Allocation for the coefficients for integral computation ***
    SLL_ALLOCATE( mesh4d%coef_int_deta1(Nx), ierr )
    SLL_ALLOCATE( mesh4d%coef_int_deta2(Ny), ierr )
    SLL_ALLOCATE( mesh4d%coef_int_dvx(Nvx), ierr )
    SLL_ALLOCATE( mesh4d%coef_int_dvy(Nvy), ierr )

  end subroutine new_mesh_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Mesh: Deallocation
  !---------------------------------------------------------------------------
  subroutine delete_mesh_VP4D( mesh4d )

    type(mesh_VP4D_t), intent(inout) :: mesh4d
    
    sll_int32 :: ierr

    call sll_o_delete( mesh4d%eta1_eta2_mesh2d )
    call sll_o_delete( mesh4d%vx_vy_mesh2d )

    SLL_DEALLOCATE( mesh4d%xgrid_2d, ierr )
    SLL_DEALLOCATE( mesh4d%ygrid_2d, ierr )
    
    SLL_DEALLOCATE( mesh4d%eta1_grid, ierr )
    SLL_DEALLOCATE( mesh4d%eta2_grid, ierr )
    SLL_DEALLOCATE( mesh4d%vx_grid, ierr )
    SLL_DEALLOCATE( mesh4d%vy_grid, ierr )

    SLL_DEALLOCATE( mesh4d%Jacobian_matrix, ierr )
    SLL_DEALLOCATE( mesh4d%inv_Jacobian_matrix, ierr )
    SLL_DEALLOCATE( mesh4d%Jacobian, ierr )

    SLL_DEALLOCATE( mesh4d%coef_int_deta1, ierr )
    SLL_DEALLOCATE( mesh4d%coef_int_deta2, ierr )
    SLL_DEALLOCATE( mesh4d%coef_int_dvx, ierr )
    SLL_DEALLOCATE( mesh4d%coef_int_dvy, ierr )

  end subroutine delete_mesh_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Mesh: Initialization
  !---------------------------------------------------------------------------
  subroutine init_mesh_VP4D( mesh4d )

    use sll_m_constants

    type(mesh_VP4D_t), intent(inout) :: mesh4d

    sll_int32 :: ieta1, ieta2, ivx, ivy
    sll_int32 :: Neta1, Neta2, Nvx, Nvy
    sll_real64 :: x_point, y_point
    sll_real64 :: delta_eta1_tmp, delta_eta2_tmp
    sll_real64 :: delta_vx_tmp, delta_vy_tmp

    Neta1 = size(mesh4d%eta1_grid)
    Neta2 = size(mesh4d%eta2_grid)
    Nvx   = size(mesh4d%vx_grid)
    Nvy   = size(mesh4d%vy_grid)

    !-> Initialization of the grid in eta1 direction
    do ieta1 = 1,Neta1
      mesh4d%eta1_grid(ieta1) = mesh4d%eta1_eta2_mesh2d%eta1_min + &
          (ieta1-1)*mesh4d%eta1_eta2_mesh2d%delta_eta1
    end do
    mesh4d%eta1_min    = mesh4d%eta1_grid(1)
    mesh4d%eta1_max    = mesh4d%eta1_grid(Neta1)
    mesh4d%eta1_length = abs(mesh4d%eta1_max-mesh4d%eta1_min)

    !-> Initialization of the grid in eta2 direction
    do ieta2 = 1,Neta2
      mesh4d%eta2_grid(ieta2) = mesh4d%eta1_eta2_mesh2d%eta2_min + &
          (ieta2-1)*mesh4d%eta1_eta2_mesh2d%delta_eta2
    end do
    mesh4d%eta2_min    = mesh4d%eta2_grid(1)
    mesh4d%eta2_max    = mesh4d%eta2_grid(Neta2)
    mesh4d%eta2_length = abs(mesh4d%eta2_max-mesh4d%eta2_min)

    !-> Initialization of the grid in vx direction
    do ivx = 1,Nvx
      mesh4d%vx_grid(ivx) = mesh4d%vx_vy_mesh2d%eta1_min + &
          (ivx-1)*mesh4d%vx_vy_mesh2d%delta_eta1
    end do
    mesh4d%vx_min    = mesh4d%vx_grid(1)
    mesh4d%vx_max    = mesh4d%vx_grid(Nvx)
    mesh4d%vx_length = abs(mesh4d%vx_max-mesh4d%vx_min)

    !-> Initialization of the grid in vy direction
    do ivy = 1,Nvy
      mesh4d%vy_grid(ivy) = mesh4d%vx_vy_mesh2d%eta1_min + &
          (ivy-1)*mesh4d%vx_vy_mesh2d%delta_eta2
    end do
    mesh4d%vy_min    = mesh4d%vy_grid(1)
    mesh4d%vy_max    = mesh4d%vy_grid(Nvy)
    mesh4d%vy_length = abs(mesh4d%vy_max-mesh4d%vy_min)

    !-> Initialization of the grid in (x,y) direction depending 
    !->  of the transformation
    do ieta2 = 1,Neta2
      do ieta1 = 1,Neta1
        x_point = mesh4d%transf_eta1eta2_xy%x1_at_node(ieta1,ieta2)
        y_point = mesh4d%transf_eta1eta2_xy%x2_at_node(ieta1,ieta2)
        mesh4d%xgrid_2d(ieta1,ieta2) = x_point            
        mesh4d%ygrid_2d(ieta1,ieta2) = y_point
      end do
    end do

    !-> Computation of the Jacobian matrix, its inverse and its determinant 
    call computeJacobian_mesh_VP4D_( mesh4d )

    !-> Initialization of the coefficients for integral computation
    delta_eta1_tmp = mesh4d%eta1_eta2_mesh2d%delta_eta1
    delta_eta2_tmp = mesh4d%eta1_eta2_mesh2d%delta_eta2
    delta_vx_tmp   = mesh4d%vx_vy_mesh2d%delta_eta1
    delta_vy_tmp   = mesh4d%vx_vy_mesh2d%delta_eta2
    do ieta1 = 2,Neta1-1
      mesh4d%coef_int_deta1(ieta1) = delta_eta1_tmp
    end do
    mesh4d%coef_int_deta1(1)     = 0.5_f64*delta_eta1_tmp
    mesh4d%coef_int_deta1(Neta1) = 0.5_f64*delta_eta1_tmp

    do ieta2 = 2,Neta2-1
      mesh4d%coef_int_deta2(ieta2) = delta_eta2_tmp
    end do
    mesh4d%coef_int_deta2(1)     = 0.5_f64*delta_eta2_tmp
    mesh4d%coef_int_deta2(Neta2) = 0.5_f64*delta_eta2_tmp

    do ivx = 2,Nvx-1
      mesh4d%coef_int_dvx(ivx) = delta_vx_tmp
    end do
    mesh4d%coef_int_dvx(1)   = 0.5_f64*delta_vx_tmp
    mesh4d%coef_int_dvx(Nvx) = 0.5_f64*delta_vx_tmp

    do ivy = 2,Nvy-1
      mesh4d%coef_int_dvy(ivy) = delta_vy_tmp
    end do
    mesh4d%coef_int_dvy(1)   = 0.5_f64*delta_vy_tmp
    mesh4d%coef_int_dvy(Nvy) = 0.5_f64*delta_vy_tmp

  end subroutine init_mesh_VP4D


  !===========================================================================
  !> Caid Mesh modification (temporary routine which must be improved)
  !---------------------------------------------------------------------------
  subroutine CAID_collocpoints_modif_( mesh4d )

    use sll_m_interpolators_2d_base
    use sll_m_coordinate_transformations_2d_nurbs
    
    type(mesh_VP4D_t), intent(inout) :: mesh4d

    !---> For control point number
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_int32 :: IO_stat
    sll_int32 :: input_file_id
    !---> For modification of the control point positions
    sll_int32 :: icoef1, icoef2
    sll_int32 :: Ncoef1, Ncoef2
    sll_int32 :: ierr
    sll_real64, dimension(:,:), pointer :: control_pts1_2d
    sll_real64, dimension(:,:), pointer :: control_pts2_2d

    class(sll_c_interpolator_2d), pointer :: x1_interp_tmp
    class(sll_c_interpolator_2d), pointer :: x2_interp_tmp
    class(sll_c_coordinate_transformation_2d_base), pointer :: transf_tmp

    namelist /shape/ num_pts1, num_pts2 

    !*** Reading of the number of control points in each directions ***
    call sll_s_new_file_id( input_file_id, ierr )
    open(unit=input_file_id, &
        file=trim(mesh4d%mapping_filename), STATUS="OLD", IOStat=IO_stat)
    read( input_file_id, shape )
    Ncoef1 = num_pts1
    Ncoef2 = num_pts2
    close( input_file_id )
    print*, "Ncoef1 = ", Ncoef1
    print*, "Ncoef2 = ", Ncoef2

    !*** Modification of the control point values ***
    SLL_ALLOCATE( control_pts1_2d(Ncoef1,Ncoef2), ierr )
    SLL_ALLOCATE( control_pts2_2d(Ncoef1,Ncoef2), ierr )

    transf_tmp => mesh4d%transf_eta1eta2_xy
    select type (transf_tmp)
    type is (sll_t_coordinate_transformation_2d_nurbs)
      x1_interp_tmp => transf_tmp%x1_interp
      x2_interp_tmp => transf_tmp%x2_interp
      control_pts1_2d = x1_interp_tmp%get_coefficients()
      control_pts2_2d = x2_interp_tmp%get_coefficients()
      do icoef2 = 1,Ncoef2
        do icoef1 = 1,Ncoef1
          control_pts1_2d(icoef1,icoef2) = &
              control_pts1_2d(icoef1,icoef2) + 1._f64
          control_pts2_2d(icoef1,icoef2) = &
              control_pts2_2d(icoef1,icoef2) + 1._f64
        end do
      end do
      call x1_interp_tmp%set_coefficients( &
        coeffs_2d=control_pts1_2d, &
        coeff2d_size1=Ncoef1, &
        coeff2d_size2=Ncoef2 )
      call x2_interp_tmp%set_coefficients( &
        coeffs_2d=control_pts2_2d, &
        coeff2d_size1=Ncoef1, &
        coeff2d_size2=Ncoef2 )
    end select

    SLL_DEALLOCATE( control_pts1_2d, ierr )
    SLL_DEALLOCATE( control_pts2_2d, ierr )
    
  end subroutine CAID_collocpoints_modif_
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the Jacobian matrix, its inverse and its determinant
  !---------------------------------------------------------------------------
  subroutine computeJacobian_mesh_VP4D_( mesh4d )

    type(mesh_VP4D_t), intent(inout) :: mesh4d

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: eta1_tmp, eta2_tmp
    
    Neta1 = size(mesh4d%eta1_grid)
    Neta2 = size(mesh4d%eta2_grid)

    do ieta2 = 1,Neta2      
      eta2_tmp = mesh4d%eta2_grid(ieta2)
      do ieta1 = 1,Neta1
        eta1_tmp = mesh4d%eta1_grid(ieta1)
        mesh4d%Jacobian_matrix(ieta1,ieta2,:,:)     =  &
            mesh4d%transf_eta1eta2_xy%jacobian_matrix(eta1_tmp,eta2_tmp)
        mesh4d%Jacobian(ieta1,ieta2)                =  &
            mesh4d%Jacobian_matrix(ieta1,ieta2,1,1) * &
            mesh4d%Jacobian_matrix(ieta1,ieta2,2,2) - &
            mesh4d%Jacobian_matrix(ieta1,ieta2,1,2) * &
            mesh4d%Jacobian_matrix(ieta1,ieta2,2,1)
        mesh4d%inv_Jacobian_matrix(ieta1,ieta2,:,:) =  &
            mesh4d%transf_eta1eta2_xy%inverse_jacobian_matrix( &
            eta1_tmp,eta2_tmp)
       end do
    end do

  end subroutine computeJacobian_mesh_VP4D_
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Mesh: Print in HDF5
  !---------------------------------------------------------------------------
  subroutine print_mesh_VP4D( mesh4d, filename )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
         sll_o_hdf5_ser_write_array, &
         sll_s_hdf5_ser_file_close, &
         sll_t_hdf5_ser_handle

    type(mesh_VP4D_t), intent(in) :: mesh4d
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
      call sll_o_hdf5_ser_write_array( handle, mesh4d%vx_grid, &
          'vx_grid', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%vy_grid, &
          'vy_grid', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%xgrid_2d, &
          'xgrid_2d', file_err)
      call sll_o_hdf5_ser_write_array( handle, mesh4d%ygrid_2d, &
          'ygrid_2d', file_err)

      call sll_s_hdf5_ser_file_close( handle, file_err )
    end if

  end subroutine print_mesh_VP4D
  !---------------------------------------------------------------------------

end module mesh_VP4D_module
!---------------------------------------------------------------------------
