program unit_test_initializers_4d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_4d, &
    sll_o_display, &
    operator(*)

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_common_coordinate_transformations, only: &
    sll_f_identity_jac11, &
    sll_f_identity_jac12, &
    sll_f_identity_jac21, &
    sll_f_identity_jac22, &
    sll_f_identity_x1, &
    sll_f_identity_x2

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_4d, &
    sll_o_local_to_global, &
    sll_f_new_layout_4d, &
    sll_o_view_lims

  use m_parallel_array_output, only: &
    write_mesh_4d, &
    write_xmf_file

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define MPI_MASTER 0

type(sll_t_layout_4d),                              pointer :: layout
type(sll_t_cartesian_mesh_2d),                  pointer :: mx
type(sll_t_cartesian_mesh_2d),                  pointer :: mv
type(sll_t_cartesian_mesh_4d),                  pointer :: mesh_4d
class(sll_c_coordinate_transformation_2d_base), pointer :: tx
sll_real64, dimension(2)                              :: params_identity
sll_real64, dimension(:,:,:,:), allocatable           :: f


sll_int32  :: prank
sll_int32  :: comm
sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
sll_int64  :: psize
sll_int32  :: i,j,k,l,error
sll_real64 :: vx,vy,v2,x,y
sll_real64 :: kx, ky, eps = 0.1_f64
sll_int32  :: gi, gj, gk, gl
sll_int32  :: global_indices(4)

call sll_s_boot_collective()

prank = sll_f_get_collective_rank(sll_v_world_collective)
psize = int(sll_f_get_collective_size(sll_v_world_collective),i64)
comm  = sll_v_world_collective%comm

params_identity(:) = (/0.0_f64, 0.0_f64/) ! for identity this can be whatever
! initialize the logical meshes
mx => sll_f_new_cartesian_mesh_2d(63,63)
mv => sll_f_new_cartesian_mesh_2d(63,63)

! initialize the transformation
tx => sll_f_new_coordinate_transformation_2d_analytic( &
     "identity_transformation",                  &
     mx,                                         &
     sll_f_identity_x1,                                &
     sll_f_identity_x2,                                &
     sll_f_identity_jac11,                             &
     sll_f_identity_jac12,                             &
     sll_f_identity_jac21,                             &
     sll_f_identity_jac22,                             &
     params_identity )

! initialize the array

! THIS TEST NEEDS TO BE FINISHED IN A COMPLETE WAY. FOR NOW, THE MODULE
! WILL BE PARTIALLY TESTED ON A SIMULATION...

mesh_4d => mx * mv

if (prank == MPI_MASTER) call sll_o_display(mesh_4d)

call write_mesh_4d(mesh_4d)

!Layout for plotting
layout => sll_f_new_layout_4d( sll_v_world_collective )

call sll_o_initialize_layout_with_distributed_array( mx%num_cells1+1,    &
                                               mx%num_cells2+1,    &
                                               mv%num_cells1+1,    &
                                               mv%num_cells2+1,    &
                                               1,1,1,int(psize,4), &
                                               layout)

if ( prank == MPI_MASTER ) call sll_o_view_lims( layout )
flush( output_unit )

call sll_o_compute_local_sizes(layout,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
SLL_CLEAR_ALLOCATE(f(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l), error)

kx  = 2_f64*sll_p_pi/(mx%eta1_max-mx%eta1_min)
ky  = 2_f64*sll_p_pi/(mx%eta2_max-mx%eta2_min)
    
do l=1,loc_sz_l 
do k=1,loc_sz_k
do j=1,loc_sz_j
do i=1,loc_sz_i

   global_indices = sll_o_local_to_global(layout,(/i,j,k,l/)) 
   gi = global_indices(1)
   gj = global_indices(2)
   gk = global_indices(3)
   gl = global_indices(4)

   x  = mesh_4d%eta1_min+(gi-1)*mesh_4d%delta_eta1
   y  = mesh_4d%eta2_min+(gj-1)*mesh_4d%delta_eta2
   vx = mesh_4d%eta3_min+(gk-1)*mesh_4d%delta_eta3
   vy = mesh_4d%eta4_min+(gl-1)*mesh_4D%delta_eta4

   v2 = vx*vx+vy*vy

   f(i,j,k,l) = landau_cos_prod(eps,kx, ky, x, y, v2)

end do
end do
end do
end do

call write_xmf_file(mesh_4d, f, layout, iplot=1)

print *, 'PASSED'

contains

  function landau_cos_prod(eps, kx, ky, x, y, v2)
    sll_real64 :: landau_cos_prod
    sll_real64, intent(in) :: x, kx
    sll_real64, intent(in) :: y, ky
    sll_real64, intent(in) :: eps, v2

    landau_cos_prod = (1._f64+eps*cos(kx*x)*cos(ky*y))/(2*sll_p_pi)*exp(-0.5_f64*v2)

  end function landau_cos_prod

end program unit_test_initializers_4d
