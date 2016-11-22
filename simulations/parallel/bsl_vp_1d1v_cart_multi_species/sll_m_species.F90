module sll_m_species

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_base, only: &
    sll_t_advection_1d_base_ptr

  use sll_m_binary_io, only: &
    sll_s_binary_file_close, &
    sll_s_binary_file_create, &
    sll_s_binary_read_array_0d, &
    sll_s_binary_read_array_2d, &
    sll_s_binary_write_array_0d, &
    sll_s_binary_write_array_2d

  use sll_m_buffer_loader_utilities, only: &
    sll_s_compute_displacements_array_2d, &
    sll_s_load_buffer_2d, &
    sll_f_receive_counts_array_2d

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_s_collective_gatherv_real64, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_common_array_initializers, only: &
    sll_f_beam_initializer_2d, &
    sll_f_bump_on_tail_initializer_2d, &
    sll_f_landau_initializer_2d, &
    sll_i_scalar_initializer_2d, &
    sll_f_two_stream_instability_initializer_2d

  use sll_m_fft, only: &
    sll_s_fft_exec_r2r_1d, &
    sll_f_fft_get_mode_r2c_1d, &
    sll_t_fft

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_1d

  use hdf5, only: hid_t
  use sll_m_hdf5_io_serial, only: &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close, &
    sll_o_hdf5_ser_write_array

  use sll_m_parallel_array_initializer, only: &
    sll_o_2d_parallel_array_initializer_cartesian

  use sll_m_primitives, only: &
    sll_s_function_to_primitive, &
    sll_s_primitive_to_function

  use sll_m_remapper, only: &
    sll_o_apply_remap_2d, &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_2d_real64

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_m_xdmf, only: &
    sll_s_xdmf_close, &
    sll_o_xdmf_open, &
    sll_o_xdmf_write_array

  use sll_mpi, only: &
    mpi_sum

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_thread_num

#endif
  implicit none

  public :: &
    sll_s_advection_x1, &
    sll_s_advection_x2, &
    sll_s_change_initial_function_species, &
    sll_s_compute_rho, &
    sll_s_diagnostics, &
    sll_s_initialize_species, &
    sll_s_read_restart_file, &
    sll_s_set_initial_function, &
    sll_p_advective, &
    sll_p_conservative, &
    sll_t_species, &
    sll_s_write_f, &
    sll_s_write_restart_file

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type sll_t_species
  character                                                :: label
  sll_int32                                                :: num_bloc_x2
  sll_int32                                                :: num_dof_x2
  type(sll_t_cartesian_mesh_2d),                   pointer :: mesh2d
  procedure(sll_i_scalar_initializer_2d), nopass,  pointer :: init_func
  sll_real64, dimension(:),                        pointer :: params
  sll_real64                                               :: nrj0
  sll_real64, dimension(:),                        pointer :: x1_array
  sll_real64, dimension(:),                        pointer :: x2_array
  sll_int32,  dimension(:),                        pointer :: bloc_index_x2
  sll_real64                                               :: kx
  sll_real64                                               :: eps
  sll_real64, dimension(:),                        pointer :: integration_weight
  sll_real64                                               :: factor_x1
  type(sll_t_advection_1d_base_ptr), dimension(:), pointer :: advect_x1
  type(sll_t_advection_1d_base_ptr), dimension(:), pointer :: advect_x2
  sll_int32                                                :: advection_form_x2
  sll_real64                                               :: alpha
  type(sll_t_remap_plan_2d_real64),                pointer :: remap_plan_x1_x2
  type(sll_t_remap_plan_2d_real64),                pointer :: remap_plan_x2_x1
  type(sll_t_layout_2d),                           pointer :: layout_x1
  type(sll_t_layout_2d),                           pointer :: layout_x2
  sll_real64, dimension(:,:),                      pointer :: f_x1
  sll_real64, dimension(:,:),                      pointer :: f_x2
  sll_real64, dimension(:,:),                      pointer :: f_x1_init
  sll_real64, dimension(:),                        pointer :: rho
  sll_real64, dimension(:),                        pointer :: rho_loc
  sll_real64, dimension(:),                        pointer :: x2_array_unit
  sll_real64, dimension(:),                        pointer :: x2_array_middle
  sll_real64, dimension(:),                        pointer :: node_positions_x2
  sll_real64, dimension(:,:),                      pointer :: f_visu 
  sll_real64, dimension(:),                        pointer :: f_visu_buf1d
  sll_real64, dimension(:),                        pointer :: f_x1_buf1d
  sll_real64, dimension(:),                        pointer :: f_hat_x2_loc
  sll_real64, dimension(:),                        pointer :: f_hat_x2
  sll_real64, dimension(:,:),                      pointer :: f1d_omp_in
  sll_real64, dimension(:,:),                      pointer :: f1d_omp_out

  sll_real64                                               :: mass           
  sll_real64                                               :: momentum       
  sll_real64                                               :: l1norm         
  sll_real64                                               :: l2norm         
  sll_real64                                               :: kinetic_energy 

end type sll_t_species

integer, parameter :: sll_p_advective = 0
integer, parameter :: sll_p_conservative = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine sll_s_set_initial_function(sp &
, initial_function_case            & 
, kmode                            &
, eps                              &
, sigma                            &
, v0                               &
, factor1                          &
, alpha_gaussian)       

type(sll_t_species),   intent(inout) :: sp

character(len=256), intent(in) :: initial_function_case
sll_real64        , intent(in) :: kmode
sll_real64        , intent(in) :: eps
sll_real64        , intent(in) :: sigma
sll_real64        , intent(in) :: v0
sll_real64        , intent(in) :: factor1
sll_real64        , intent(in) :: alpha_gaussian

sll_int32 :: ierr

sp%nrj0 = 0._f64
sp%kx   = kmode
sp%eps  = eps
select case (initial_function_case)
case ("SLL_LANDAU")
  sp%init_func => sll_f_landau_initializer_2d
  SLL_ALLOCATE(sp%params(4),ierr)
  sp%params(1) = kmode
  sp%params(2) = eps
  sp%params(3) = v0
  sp%params(4) = sigma
  sp%nrj0 = 0._f64  !compute the right value
  !(0.5_f64*eps*sll_p_pi)**2/(kmode_x1*kmode_x2) &
  !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
  !for the moment
  sp%kx  = kmode
  sp%eps = eps
case ("SLL_BUMP_ON_TAIL")
  sp%init_func => sll_f_bump_on_tail_initializer_2d
  SLL_ALLOCATE(sp%params(2),ierr)
  sp%params(1) = kmode
  sp%params(2) = eps
  sp%nrj0 = 0._f64  !compute the right value
  !(0.5_f64*eps*sll_p_pi)**2/(kmode_x1*kmode_x2) &
  !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
  !for the moment
  sp%kx = kmode
  sp%eps = eps
case ("SLL_TWO_STREAM_INSTABILITY")
  sp%init_func => sll_f_two_stream_instability_initializer_2d
  SLL_ALLOCATE(sp%params(4),ierr)
  sp%params(1) = kmode
  sp%params(2) = eps
  sp%params(3) = sigma
  sp%params(4) = factor1
  sp%nrj0 = 0._f64  !compute the right value
  !(0.5_f64*eps*sll_p_pi)**2/(kmode_x1*kmode_x2) &
  !*(1._f64/kmode_x1**2+1._f64/kmode_x2**2)
  !for the moment
  sp%kx  = kmode
  sp%eps = eps
case ("SLL_BEAM")  
  sp%init_func => sll_f_beam_initializer_2d
  SLL_ALLOCATE(sp%params(1),ierr)
  sp%params(1) = alpha_gaussian
case default
  print *,'#init_func_case not implemented'
  stop
end select

end subroutine sll_s_set_initial_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sll_s_initialize_species(sp, label, nb_mode)

type(sll_t_species), intent(inout) :: sp
character,     intent(in)    :: label
sll_int32,     intent(in)    :: nb_mode

sll_int32  :: nc_x1, nc_x2
sll_int32  :: np_x1, np_x2
sll_int32  :: ierr
sll_int32  :: psize
sll_int32  :: prank
sll_int32  :: nproc_x1
sll_int32  :: nproc_x2
logical    :: mpi_master
sll_int32  :: global_indices(2)
sll_int32  :: local_size_x1
sll_int32  :: local_size_x2
sll_int32  :: i

sp%label = label
prank = sll_f_get_collective_rank( sll_v_world_collective )
psize = sll_f_get_collective_size( sll_v_world_collective )
mpi_master = merge(.true., .false., prank == 0)

nc_x1 = sp%mesh2d%num_cells1
np_x1 = sp%mesh2d%num_cells1+1
nc_x2 = sp%mesh2d%num_cells2
np_x2 = sp%mesh2d%num_cells2+1

if(mpi_master) then
  SLL_ALLOCATE(sp%f_visu(np_x1,sp%num_dof_x2),ierr)
  SLL_ALLOCATE(sp%f_visu_buf1d(np_x1*sp%num_dof_x2),ierr)
else
  SLL_ALLOCATE(sp%f_visu(1:1,1:1),ierr)
  SLL_ALLOCATE(sp%f_visu_buf1d(1:1),ierr)
endif

sp%layout_x1 => sll_f_new_layout_2d( sll_v_world_collective )
sp%layout_x2 => sll_f_new_layout_2d( sll_v_world_collective )

nproc_x1 = psize
nproc_x2 = 1
call sll_o_initialize_layout_with_distributed_array( &
     np_x1, sp%num_dof_x2, nproc_x1, nproc_x2, sp%layout_x2 )
call sll_o_initialize_layout_with_distributed_array( &
     np_x1, sp%num_dof_x2, nproc_x2, nproc_x1, sp%layout_x1 )

call sll_o_compute_local_sizes( sp%layout_x2, local_size_x1, local_size_x2 )

SLL_ALLOCATE(sp%f_x2(local_size_x1,local_size_x2),ierr)

call sll_o_compute_local_sizes( sp%layout_x1, local_size_x1, local_size_x2 )

global_indices(1:2) = sll_o_local_to_global( sp%layout_x1, (/1, 1/) )
SLL_ALLOCATE(sp%f_x1(local_size_x1,local_size_x2),ierr)    
SLL_ALLOCATE(sp%f_x1_init(local_size_x1,local_size_x2),ierr)    
SLL_ALLOCATE(sp%f_x1_buf1d(local_size_x1*local_size_x2),ierr)    

sp%remap_plan_x1_x2 => sll_o_new_remap_plan(sp%layout_x1, sp%layout_x2, sp%f_x1)
sp%remap_plan_x2_x1 => sll_o_new_remap_plan(sp%layout_x2, sp%layout_x1, sp%f_x2)

SLL_ALLOCATE(sp%rho(np_x1),ierr)
SLL_ALLOCATE(sp%rho_loc(np_x1),ierr)


SLL_ALLOCATE(sp%x2_array_unit(np_x2),ierr)
SLL_ALLOCATE(sp%x2_array_middle(np_x2),ierr)
SLL_ALLOCATE(sp%node_positions_x2(sp%num_dof_x2),ierr)
SLL_ALLOCATE(sp%f_hat_x2_loc(nb_mode+1),ierr)
SLL_ALLOCATE(sp%f_hat_x2(nb_mode+1),ierr)

sp%x2_array_unit(1:np_x2) = &
  (sp%x2_array(1:np_x2)-sp%x2_array(1))/(sp%x2_array(np_x2)-sp%x2_array(1))
do i = 1, np_x2-1
   sp%x2_array_middle(i) = 0.5_f64*(sp%x2_array(i)+sp%x2_array(i+1))
end do
sp%x2_array_middle(np_x2) = &
  sp%x2_array_middle(1)+sp%x2_array(np_x2)-sp%x2_array(1)

select case (sp%advection_form_x2)
case (sll_p_advective)
  sp%node_positions_x2(1:sp%num_dof_x2) = sp%x2_array(1:sp%num_dof_x2)
case (sll_p_conservative)
  sp%node_positions_x2(1:sp%num_dof_x2) = sp%x2_array_middle(1:sp%num_dof_x2)
case default
  print *,'#advection_form_x2=',sp%advection_form_x2
  print *,'#not implemented'
  stop
end select  

call sll_o_2d_parallel_array_initializer_cartesian( &
  sp%layout_x1,                                   &
  sp%x1_array,                                    &
  sp%node_positions_x2,                           &
  sp%f_x1,                                        &
  sp%init_func,                                   &
  sp%params)

call sll_o_2d_parallel_array_initializer_cartesian( &
  sp%layout_x1,                                   &
  sp%x1_array,                                    &
  sp%node_positions_x2,                           &
  sp%f_x1_init,                                   &
  sll_f_landau_initializer_2d,                      &
  [sp%params(1),0._f64,sp%params(3),sp%params(4)])

end subroutine sll_s_initialize_species

subroutine sll_s_change_initial_function_species( &
  sp,                                       &
  init_func,                                &
  params,                                   &
  num_params)

type(sll_t_species), intent(inout) :: sp

procedure(sll_i_scalar_initializer_2d), pointer    :: init_func
sll_int32,                            intent(in) :: num_params
sll_real64, dimension(num_params),    intent(in) :: params

sll_int32 :: ierr
sll_int32 :: i

sp%init_func => init_func

if(associated(sp%params))then
  SLL_DEALLOCATE(sp%params,ierr)
endif

if(num_params<1)then
  print *,'#num_params should be >=1 in sll_s_change_initial_function_species'
  stop
endif
SLL_ALLOCATE(sp%params(num_params),ierr)
if(size(params)<num_params)then
  print *,'#size of params is not good in sll_s_change_initial_function_species'
  stop
endif

do i=1,num_params
  sp%params(i) = params(i)
enddo

end subroutine sll_s_change_initial_function_species



subroutine sll_s_compute_rho(sp)
type(sll_t_species), intent(inout) :: sp

sll_int32 :: np_x1
sll_int32 :: np_x2
sll_int32 :: global_indices(2)
sll_int32 :: local_size_x1
sll_int32 :: local_size_x2
sll_int32 :: i, ig

np_x1 = sp%mesh2d%num_cells1+1
np_x2 = sp%mesh2d%num_cells2+1

call sll_o_compute_local_sizes(sp%layout_x1,local_size_x1,local_size_x2)
global_indices = sll_o_local_to_global( sp%layout_x1, (/1, 1/) )

sp%rho_loc = 0.0_f64
ig = global_indices(2)-1
do i=1,np_x1
  sp%rho_loc(i)=sp%rho_loc(i) &
    +sum(sp%f_x1(i,1:local_size_x2) &
    *sp%integration_weight(1+ig:local_size_x2+ig))
end do

call sll_o_collective_allreduce( sll_v_world_collective, &
                               sp%rho_loc,           &
                               np_x1,                &
                               MPI_SUM,              &
                               sp%rho )

end subroutine sll_s_compute_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_diagnostics(sp, pfwd, buf_fft, nb_mode)
type(sll_t_species), intent(inout) :: sp
sll_int32,     intent(in)    :: nb_mode
type(sll_t_fft)              :: pfwd
sll_real64                   :: buf_fft(:)

sll_real64  :: tmp(5), tmp_loc(5)
sll_int32   :: global_indices(2)
sll_int32   :: local_size_x1, local_size_x2   
sll_int32   :: i, ig, k

sp%mass            = 0._f64
sp%momentum        = 0._f64
sp%l1norm          = 0._f64
sp%l2norm          = 0._f64
sp%kinetic_energy  = 0._f64

tmp_loc            = 0._f64

call sll_o_compute_local_sizes(sp%layout_x1, local_size_x1, local_size_x2)
global_indices = sll_o_local_to_global( sp%layout_x1, (/1, 1/) )
ig = global_indices(2)-1               

do i = 1, sp%mesh2d%num_cells1
  tmp_loc(1)= tmp_loc(1)+sum(sp%f_x1(i,1:local_size_x2) &
    *sp%integration_weight(1+ig:local_size_x2+ig))
  tmp_loc(2)= tmp_loc(2)+sum(abs(sp%f_x1(i,1:local_size_x2)) &
    *sp%integration_weight(1+ig:local_size_x2+ig))
  tmp_loc(3)= tmp_loc(3)+sum((sp%f_x1(i,1:local_size_x2))**2 &
    *sp%integration_weight(1+ig:local_size_x2+ig))
  tmp_loc(4)= tmp_loc(4) +sum(sp%f_x1(i,1:local_size_x2) &
    *sp%x2_array(global_indices(2)-1+1:global_indices(2)-1+local_size_x2) &
    *sp%integration_weight(1+ig:local_size_x2+ig))          
  tmp_loc(5)= tmp_loc(5)+sum(sp%f_x1(i,1:local_size_x2) &
    *sp%x2_array(global_indices(2)-1+1:global_indices(2)-1+local_size_x2)**2 &
    *sp%integration_weight(1+ig:local_size_x2+ig) )          
end do

call sll_o_collective_allreduce(sll_v_world_collective, tmp_loc, 5, MPI_SUM, tmp)

sp%mass           = tmp(1)  * sp%mesh2d%delta_eta1 
sp%l1norm         = tmp(2)  * sp%mesh2d%delta_eta1 
sp%l2norm         = tmp(3)  * sp%mesh2d%delta_eta1 
sp%momentum       = tmp(4)  * sp%mesh2d%delta_eta1 
sp%kinetic_energy = 0.5_f64 * tmp(5) * sp%mesh2d%delta_eta1 

sp%f_hat_x2_loc(1:nb_mode+1) = 0._f64
do i=1,local_size_x2
  buf_fft = sp%f_x1(1:sp%mesh2d%num_cells1,i)
  call sll_s_fft_exec_r2r_1d(pfwd,buf_fft,buf_fft)
  do k=0,nb_mode
     sp%f_hat_x2_loc(k+1) = sp%f_hat_x2_loc(k+1)   &
       +abs(sll_f_fft_get_mode_r2c_1d(pfwd,buf_fft,k))**2 &
       *sp%integration_weight(ig+i)
  enddo
enddo

call sll_o_collective_allreduce( sll_v_world_collective, sp%f_hat_x2_loc, &
  nb_mode+1, MPI_SUM, sp%f_hat_x2 )

end subroutine sll_s_diagnostics

subroutine plot_f_cartesian( iplot,             &
                             f,                 &
                             node_positions_x1, &
                             nnodes_x1,         &
                             node_positions_x2, &
                             nnodes_x2,         &
                             array_name,        &
                             spec_name,         &
                             time)


sll_real64, dimension(:),   intent(in) :: node_positions_x1
sll_real64, dimension(:),   intent(in) :: node_positions_x2    
character(len=*),           intent(in) :: array_name !< field name
character(len=*),           intent(in) :: spec_name !< name of the sll_t_species
sll_int32,                  intent(in) :: nnodes_x1
sll_int32,                  intent(in) :: nnodes_x2
sll_int32,                  intent(in) :: iplot
sll_real64, dimension(:,:), intent(in) :: f

sll_int32                               :: file_id
integer(hid_t)                          :: hfile_id
sll_int32                               :: error
sll_real64, dimension(:,:), allocatable :: x1
sll_real64, dimension(:,:), allocatable :: x2
sll_int32                               :: i
sll_int32                               :: j
character(len=4)                        :: cplot
sll_real64                              :: time

if (iplot == 1) then

  SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
  SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
  do j = 1,nnodes_x2
    do i = 1,nnodes_x1
      x1(i,j) = node_positions_x1(i) !x1_min+real(i-1,f32)*dx1
      x2(i,j) = node_positions_x2(j) !x2_min+real(j-1,f32)*dx2
    end do
  end do
#ifndef NOHDF5
  call sll_s_hdf5_ser_file_create("cartesian_mesh_"//trim(spec_name)//"-x1.h5", &
    hfile_id,error)
  call sll_o_hdf5_ser_write_array(hfile_id,x1,"/x1",error)
  call sll_s_hdf5_ser_file_close(hfile_id, error)
  call sll_s_hdf5_ser_file_create("cartesian_mesh_"//trim(spec_name)//"-x2.h5", &
    hfile_id,error)
  call sll_o_hdf5_ser_write_array(hfile_id,x2,"/x2",error)
  call sll_s_hdf5_ser_file_close(hfile_id, error)
#endif

  deallocate(x1)
  deallocate(x2)

end if

call sll_s_int2string(iplot,cplot)
call sll_o_xdmf_open(trim(array_name)//cplot//".xmf", &
  "cartesian_mesh_"//trim(spec_name), &
  nnodes_x1,nnodes_x2,file_id,error)
write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"
call sll_o_xdmf_write_array(trim(array_name)//cplot,f,"values", &
  error,file_id,"Node")
call sll_s_xdmf_close(file_id,error)

end subroutine plot_f_cartesian

subroutine sll_s_write_f( sp, iplot, array_name, time)

  type(sll_t_species),    intent(in) :: sp
  sll_int32,        intent(in) :: iplot
  character(len=*), intent(in) :: array_name
  sll_real64,       intent(in) :: time
  sll_int32                    :: i
  sll_int32                    :: ierr
  sll_int32                    :: nc_x1
  sll_int32                    :: nc_x2
  sll_int32                    :: psize
  sll_int32                    :: prank
  logical                      :: mpi_master
  sll_int32                    :: local_size_x1
  sll_int32                    :: local_size_x2
  sll_int32,       allocatable :: collective_displs(:)
  sll_int32,       allocatable :: collective_recvcnts(:)

  nc_x1 = sp%mesh2d%num_cells1
  nc_x2 = sp%mesh2d%num_cells2

  prank = sll_f_get_collective_rank( sll_v_world_collective )
  psize = sll_f_get_collective_size( sll_v_world_collective )
  mpi_master = merge(.true., .false., prank == 0)

  SLL_ALLOCATE(collective_displs(psize),ierr)
  SLL_ALLOCATE(collective_recvcnts(psize),ierr)

  call sll_o_compute_local_sizes(sp%layout_x1, local_size_x1, local_size_x2)
  call sll_s_compute_displacements_array_2d(sp%layout_x1, psize, collective_displs )
  collective_recvcnts = sll_f_receive_counts_array_2d(sp%layout_x1, psize )

  call sll_s_load_buffer_2d( sp%layout_x1, sp%f_x1-sp%f_x1_init, sp%f_x1_buf1d )
  call sll_s_collective_gatherv_real64(    &
    sll_v_world_collective,                &
    sp%f_x1_buf1d,                       &
    local_size_x1*local_size_x2,         &
    collective_recvcnts,                 &
    collective_displs,                   &
    0,                                   &
    sp%f_visu_buf1d )

  sp%f_visu = reshape(sp%f_visu_buf1d, shape(sp%f_visu))

  if (mpi_master) then
    do i=1,sp%num_dof_x2
      sp%f_visu_buf1d(i) = sum(sp%f_visu(1:nc_x1,i))*sp%mesh2d%delta_eta1
    enddo

    call sll_o_gnuplot_1d(                     &
      sp%f_visu_buf1d(1:sp%num_dof_x2),      &
      sp%node_positions_x2(1:sp%num_dof_x2), &
      'int_'//array_name,                    &
      iplot )

    call plot_f_cartesian(       &
        iplot,                   &
        sp%f_visu,               &
        sp%x1_array,             &
        nc_x1+1,                 &
        sp%node_positions_x2,    &
        sp%num_dof_x2,           &
        "delta"//array_name,     &
        sp%label,                &
        time )        
  end if

  call sll_s_load_buffer_2d( sp%layout_x1, sp%f_x1, sp%f_x1_buf1d )
  call sll_s_collective_gatherv_real64(   &
    sll_v_world_collective,               &
    sp%f_x1_buf1d,                      &
    local_size_x1*local_size_x2,        &
    collective_recvcnts,                &
    collective_displs,                  &
    0,                                  &
    sp%f_visu_buf1d )

  sp%f_visu = reshape(sp%f_visu_buf1d, shape(sp%f_visu))


  if (mpi_master) then

    do i=1,sp%num_dof_x2
      sp%f_visu_buf1d(i) = sum(sp%f_visu(1:nc_x1,i))*sp%mesh2d%delta_eta1
    enddo

    call sll_o_gnuplot_1d( &
      sp%f_visu_buf1d(1:sp%num_dof_x2),      &
      sp%node_positions_x2(1:sp%num_dof_x2), &
      'int'//array_name//'dx',               &
      iplot )

    call plot_f_cartesian(       &
        iplot,                   &
        sp%f_visu,               &
        sp%x1_array,             &
        nc_x1+1,                 &
        sp%node_positions_x2,    &
        sp%num_dof_x2,           &
        array_name,              &
        sp%label,                &
        time )        

  end if

  print *,prank,'#max'//array_name//" = ",maxval(sp%f_visu), minval(sp%f_visu)

end subroutine sll_s_write_f

subroutine sll_s_write_restart_file(sp, iplot, time)

  type(sll_t_species),    intent(in) :: sp
  sll_int32,        intent(in) :: iplot
  sll_real64,       intent(in) :: time
  character(len=4)             :: cplot
  character(len=4)             :: cproc
  sll_int32                    :: prank
  sll_int32                    :: ierr
  sll_int32                    :: restart_id
  sll_int32                    :: local_size_x1
  sll_int32                    :: local_size_x2

  prank = sll_f_get_collective_rank( sll_v_world_collective )
  call sll_s_int2string(prank,cproc)
  call sll_s_int2string(iplot,cplot) 
  call sll_s_binary_file_create( &
    'f'//'_plot_'//cplot//'_proc_'//cproc//sp%label//'.rst', restart_id, ierr )
  call sll_s_binary_write_array_0d(restart_id,time,ierr)
  call sll_o_compute_local_sizes(sp%layout_x1, local_size_x1, local_size_x2)
  call sll_s_binary_write_array_2d(restart_id, &
    sp%f_x1(1:local_size_x1,1:local_size_x2),ierr)
  call sll_s_binary_file_close(restart_id,ierr)    

end subroutine sll_s_write_restart_file

subroutine sll_s_read_restart_file(restart_file, sp, time)

  character(len=*), intent(in)    :: restart_file
  type(sll_t_species),    intent(in)    :: sp
  sll_real64,       intent(inout) :: time
  character(len=4)                :: cproc
  sll_int32                       :: prank
  sll_int32                       :: ierr
  sll_int32                       :: restart_id
  sll_int32                       :: local_size_x1
  sll_int32                       :: local_size_x2
  logical                         :: file_exists

  prank = sll_f_get_collective_rank( sll_v_world_collective )
  call sll_s_int2string(prank,cproc)
  call sll_o_compute_local_sizes(sp%layout_x1, local_size_x1, local_size_x2)

  print*, restart_file

  if (restart_file == "no_restart_file") return

  inquire(file=restart_file//'_proc_'//cproc//sp%label//'.rst', &
          exist=file_exists)

  if(.not.(file_exists))then
    print *,'#file ',restart_file//sp%label//'_proc_'//cproc//sp%label//'.rst'
    print *,'does not exist'
    stop
  end if

  open(unit=restart_id, &
       file=restart_file//'_proc_'//cproc//sp%label//'.rst', &
       access="stream", &
       form='unformatted', &
       iostat=ierr)      

  if( ierr /= 0 ) then
    print *, 'ERROR while opening file ', &
    restart_file//'_proc_'//cproc//sp%label//'.rst', &
   '. Called from run_vp2d_cartesian().'
    stop
  end if
  print *, &
    '#read restart file '//restart_file//sp%label//'_proc_'//cproc//'.rst'      
  call sll_s_binary_read_array_0d(restart_id,time,ierr)
  call sll_s_binary_read_array_2d(restart_id, &
    sp%f_x1(1:local_size_x1,1:local_size_x2),ierr)
  call sll_s_binary_file_close(restart_id,ierr)

end subroutine sll_s_read_restart_file

subroutine sll_s_advection_x1(sp, split_step, dt)
type(sll_t_species) , intent(inout) :: sp
sll_real64    , intent(in)    :: split_step
sll_real64    , intent(in)    :: dt

sll_real64    :: alpha_omp
sll_int32     :: local_size_x1
sll_int32     :: local_size_x2
sll_int32     :: np_x1
sll_int32     :: i_omp
sll_int32     :: ig_omp
sll_int32     :: tid
sll_int32     :: global_indices(2)

np_x1 = sp%mesh2d%num_cells1+1
call sll_o_compute_local_sizes( sp%layout_x1, local_size_x1, local_size_x2 )
global_indices(1:2) = sll_o_local_to_global(sp%layout_x1, (/1, 1/) )
tid = 1          
!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(tid,i_omp,ig_omp,alpha_omp) 
!$ tid = omp_get_thread_num()+1
!$OMP DO
do i_omp = 1, local_size_x2
  ig_omp = i_omp+global_indices(2)-1
  alpha_omp = sp%factor_x1*sp%node_positions_x2(ig_omp) * split_step
  sp%f1d_omp_in(1:np_x1,tid) = sp%f_x1(1:np_x1,i_omp)
  call sp%advect_x1(tid)%ptr%advect_1d_constant(&
       alpha_omp,                               &
       dt,                                      &
       sp%f1d_omp_in(1:np_x1,tid),              &
       sp%f1d_omp_out(1:np_x1,tid))
  
  sp%f_x1(1:np_x1,i_omp)=sp%f1d_omp_out(1:np_x1,tid)

end do
!$OMP END DO
!$OMP END PARALLEL
    
end subroutine sll_s_advection_x1

subroutine sll_s_advection_x2(sp, factor, efield, split_step, dt)
type(sll_t_species) , intent(inout) :: sp
sll_real64    , intent(in)    :: factor
sll_real64    , intent(in)    :: efield(:)
sll_real64    , intent(in)    :: split_step
sll_real64    , intent(in)    :: dt

sll_real64    :: alpha_omp
sll_int32     :: local_size_x1
sll_int32     :: local_size_x2
sll_int32     :: nc_x2
sll_int32     :: i_omp
sll_int32     :: ig_omp
sll_int32     :: tid
sll_int32     :: global_indices(2)
sll_real64    :: mean_omp

call sll_o_apply_remap_2d( sp%remap_plan_x1_x2, sp%f_x1, sp%f_x2 )

nc_x2 = sp%mesh2d%num_cells2
call sll_o_compute_local_sizes(sp%layout_x2, local_size_x1, local_size_x2)
global_indices(1:2) = sll_o_local_to_global( sp%layout_x2, (/1, 1/) )

tid         = 1

!$OMP PARALLEL &
!$OMP PRIVATE(tid,ig_omp,i_omp,alpha_omp,mean_omp) &
!$OMP FIRSTPRIVATE(dt, split_step,nc_x2)
!$ tid = omp_get_thread_num()+1
!$OMP DO 
do i_omp = 1,local_size_x1

  ig_omp=i_omp+global_indices(1)-1
  alpha_omp = factor * efield(ig_omp) * split_step
  sp%f1d_omp_in(1:sp%num_dof_x2,tid) = sp%f_x2(i_omp,1:sp%num_dof_x2)

  if (sp%advection_form_x2==sll_p_conservative) then
    call sll_s_function_to_primitive(sp%f1d_omp_in(:,tid),sp%x2_array_unit, &
      nc_x2,mean_omp)
  endif

  call sp%advect_x2(tid)%ptr%advect_1d_constant(  &
           alpha_omp,                             &
           dt,                                    &
           sp%f1d_omp_in(1:sp%num_dof_x2,tid),    &
           sp%f1d_omp_out(1:sp%num_dof_x2,tid))

  if (sp%advection_form_x2==sll_p_conservative) then
    call sll_s_primitive_to_function(sp%f1d_omp_out(:,tid),sp%x2_array_unit, &
      nc_x2,mean_omp)
  endif

  sp%f_x2(i_omp,1:sp%num_dof_x2) = sp%f1d_omp_out(1:sp%num_dof_x2,tid)

end do
!$OMP END DO
!$OMP END PARALLEL

call sll_o_apply_remap_2d( sp%remap_plan_x2_x1, sp%f_x2, sp%f_x1 )
    
end subroutine sll_s_advection_x2

end module sll_m_species
