program test_charge_to_density

#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_accumulators
use sll_m_charge_to_density
use sll_m_cartesian_meshes

implicit none

type(sll_t_cartesian_mesh_2d)                      :: mesh_2d
type(sll_t_charge_accumulator_2d),         pointer :: charge
type(sll_t_charge_accumulator_2d_cs),      pointer :: charge_cs
sll_real64, dimension(:,:),                pointer :: rho
sll_real64, dimension(:,:),                pointer :: e_1, e_2
type(sll_t_electric_field_accumulator),    pointer :: e_accumulator
type(sll_t_electric_field_accumulator_cs), pointer :: e_accumulator_cs

allocate(rho(10,10))
allocate(e_1(10,10))
allocate(e_2(10,10))

call sll_s_cartesian_mesh_2d_init( mesh_2d, 10, 10, &
                                   0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64 )

call sll_s_charge_accumulator_2d_init(charge, mesh_2d )
call sll_s_convert_charge_to_rho_2d_per_per( charge, rho )
call sll_s_accumulate_field( e_1, e_2, e_accumulator)

call sll_s_charge_accumulator_2d_cs_init(charge_cs, mesh_2d )
call sll_s_convert_charge_to_rho_2d_per_per_cs( charge_cs, rho )
call sll_s_accumulate_field_cs( e_1, e_2, e_accumulator_cs)

end program test_charge_to_density
