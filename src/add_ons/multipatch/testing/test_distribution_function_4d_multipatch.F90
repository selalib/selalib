program unit_test_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_common_array_initializers, only: &
    sll_i_scalar_initializer_4d, &
    sll_f_landau_initializer_4d

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_coordinate_transformation_multipatch, only: &
    sll_t_coordinate_transformation_multipatch_2d

  use sll_m_distribution_function_4d_multipatch, only: &
    sll_s_compute_charge_density_multipatch, &
    sll_t_distribution_function_4d_multipatch, &
    sll_f_new_distribution_function_4d_multipatch, &
    sll_o_delete

  use sll_m_scalar_field_2d_multipatch, only: &
    sll_f_new_scalar_field_multipatch_2d, &
    sll_t_scalar_field_multipatch_2d, &
    sll_o_delete

  use sll_m_utilities, only: &
    sll_f_is_even



  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_coordinate_transformation_multipatch_2d)      :: t_mp
  class(sll_t_scalar_field_multipatch_2d), pointer         :: rho_mp
  type(sll_t_distribution_function_4d_multipatch), pointer :: f_mp
  sll_int32 :: npts_x3
  sll_int32 :: npts_x4
  sll_int32 :: power2
  sll_int32 :: world_size
  sll_int32 :: rank
  sll_int32 :: nproc_factor1
  sll_int32 :: nproc_factor2
  sll_real64, dimension(5) :: landau_params
  type(sll_t_cartesian_mesh_2d), pointer :: meshv
  procedure(sll_i_scalar_initializer_4d), pointer        :: init_function


  call sll_s_boot_collective()
  rank = sll_f_get_collective_rank(sll_v_world_collective)
  world_size = sll_f_get_collective_size(sll_v_world_collective)
  power2 = int(log(real(world_size))/log(2.0))

  ! special case N = 1, so power2 = 0
  if(power2 == 0) then
     nproc_factor1 = 1
     nproc_factor2 = 1
  end if
    
  if(sll_f_is_even(power2)) then
     nproc_factor1 = 2**(power2/2)
     nproc_factor2 = 2**(power2/2)
  else 
     nproc_factor1 = 2**((power2-1)/2)
     nproc_factor2 = 2**((power2+1)/2)
  end if
  
  npts_x3 = 32
  npts_x4 = 32

  meshv => sll_f_new_cartesian_mesh_2d(npts_x3, npts_x4, &
       eta1_min=-9.0_f64, eta1_max=9.0_f64, &
       eta2_min=-9.0_f64, eta2_max=9.0_f64)

  ! THE FOLLOWING VALUES MUST BE  CONSISTENT WITH THE LIMITS OF THE PHYSICAL 
  ! DOMAIN AS DEFINED BY THE COORDINATE TRANSFORMATION
  landau_params(1) = 0.0_f64
  landau_params(2) = 4.0*sll_p_pi
  landau_params(3) = 0.0_f64
  landau_params(4) = 1.0_f64
  landau_params(5) = 0.05_f64

!  call mp%read_from_file("identity_mp_info.nml")
  call t_mp%read_from_file("square_4p_n10")

  rho_mp => sll_f_new_scalar_field_multipatch_2d("rho_unit_test_df_multipatch", t_mp)
  call rho_mp%allocate_memory()

  ! The spatial dimensions of each patch are defined by the coordinate
  ! transformation, the velocity dimensions need to  be passed explicitly.
  f_mp => sll_f_new_distribution_function_4d_multipatch( sll_v_world_collective, &
       t_mp, meshv, nproc_factor1, nproc_factor2 )

  init_function => sll_f_landau_initializer_4d
  call f_mp%initialize( init_function, landau_params )

  call sll_s_compute_charge_density_multipatch( f_mp, rho_mp )

  call rho_mp%update_interpolation_coefficients()

  call f_mp%set_to_sequential_x1x2()
  call f_mp%set_to_sequential_x3x4()

  if(rank == 0) then
     call rho_mp%write_to_file(1) ! iplot must be > 0 [YG - 06.10.2015]
  end if

  call sll_o_delete(f_mp)
  call sll_o_delete(rho_mp)
  call t_mp%delete()


  print *, "PASSED"
  call sll_s_halt_collective()

end program unit_test_2d
