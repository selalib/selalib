program unit_test_2d
#include "sll_working_precision.h"
  use sll_coordinate_transformation_multipatch_module
  use sll_distribution_function_4d_multipatch_module
  use sll_collective
  use sll_common_array_initializers_module
  implicit none

  type(sll_coordinate_transformation_multipatch_2d)      :: t_mp
  class(sll_scalar_field_multipatch_2d), pointer         :: rho_mp
  type(sll_distribution_function_4d_multipatch), pointer :: f_mp
  sll_int32 :: npts_x3
  sll_int32 :: npts_x4
  sll_int32 :: power2
  sll_int32 :: world_size
  sll_int32 :: rank
  sll_int32 :: nproc_factor1
  sll_int32 :: nproc_factor2
  sll_real64, dimension(5) :: landau_params
  type(sll_logical_mesh_2d), pointer :: meshv

  call sll_boot_collective()
  rank = sll_get_collective_rank(sll_world_collective)
  world_size = sll_get_collective_size(sll_world_collective)
  power2 = int(log(real(world_size))/log(2.0))

  ! special case N = 1, so power2 = 0
  if(power2 == 0) then
     nproc_factor1 = 1
     nproc_factor2 = 1
  end if
    
  if(is_even(power2)) then
     nproc_factor1 = 2**(power2/2)
     nproc_factor2 = 2**(power2/2)
  else 
     nproc_factor1 = 2**((power2-1)/2)
     nproc_factor2 = 2**((power2+1)/2)
  end if
  
  npts_x3 = 32
  npts_x4 = 32

  meshv => new_logical_mesh_2d(npts_x3, npts_x4, &
       eta1_min=-9.0_f64, eta1_max=9.0_f64, &
       eta2_min=-9.0_f64, eta2_max=9.0_f64)

  ! THE FOLLOWING VALUES MUST BE  CONSISTENT WITH THE LIMITS OF THE PHYSICAL 
  ! DOMAIN AS DEFINED BY THE COORDINATE TRANSFORMATION
  landau_params(1) = 0.0_f64
  landau_params(2) = 4.0*sll_pi
  landau_params(3) = 0.0_f64
  landau_params(4) = 1.0_f64
  landau_params(5) = 0.05_f64

!  call mp%read_from_file("identity_mp_info.nml")
  call t_mp%read_from_file("square_4p_n10")

  rho_mp => new_scalar_field_multipatch_2d("rho_unit_test_df_multipatch", t_mp)
  call rho_mp%allocate_memory()

  ! The spatial dimensions of each patch are defined by the coordinate
  ! transformation, the velocity dimensions need to  be passed explicitly.
  f_mp => sll_new_distribution_function_4d_multipatch( sll_world_collective, &
       t_mp, meshv, nproc_factor1, nproc_factor2 )

  call f_mp%initialize( sll_landau_initializer_4d, landau_params )

  call compute_charge_density_multipatch( f_mp, rho_mp )

  call rho_mp%update_interpolation_coefficients()

  if(rank == 0) then
     call rho_mp%write_to_file(0)
  end if

  call f_mp%delete()
!  call mp%delete()

  print *, "PASSED"
  call sll_halt_collective()

end program unit_test_2d
