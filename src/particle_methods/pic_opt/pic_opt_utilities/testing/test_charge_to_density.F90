program test_charge_to_density

#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_accumulators
   use sll_m_charge_to_density
   use sll_m_cartesian_meshes

   implicit none

   type(sll_t_cartesian_mesh_2d)             :: mesh_2d
   type(sll_t_charge_accumulator_2d)         :: charge
   type(sll_t_charge_accumulator_2d_cs)      :: charge_cs
   sll_real64, dimension(:, :), allocatable   :: rho
   sll_real64, dimension(:, :), allocatable   :: e_1, e_2
   type(sll_t_electric_field_accumulator)    :: e_accumulator
   type(sll_t_electric_field_accumulator_cs) :: e_accumulator_cs

   sll_int32, parameter :: num_cells1 = 10
   sll_int32, parameter :: num_cells2 = 10

   sll_real64 :: ref = 12100_f64, s

   allocate (rho(num_cells1 + 1, num_cells2 + 1))
   allocate (e_1(num_cells1 + 1, num_cells2 + 1))
   allocate (e_2(num_cells1 + 1, num_cells2 + 1))

   call sll_s_cartesian_mesh_2d_init(mesh_2d, &
                                     num_cells1, num_cells2, 0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64)

   call sll_s_charge_accumulator_2d_init(charge, mesh_2d)
   charge%q_acc%q_sw = 1.0_f64
   call sll_s_convert_charge_to_rho_2d_per_per(charge, rho)
   call sll_s_field_accumulator_2d_init(e_accumulator, mesh_2d)
   call sll_s_accumulate_field(e_1, e_2, e_accumulator)

   s = sum(rho)
   if (abs(s - ref) > 1d-7) stop 'FAILED'
   s = sum(e_1)
   if (abs(s) > 1d-7) stop 'FAILED'
   s = sum(e_2)
   if (abs(s) > 1d-7) stop 'FAILED'

   call sll_s_charge_accumulator_2d_cs_init(charge_cs, mesh_2d)
   charge_cs%q_acc%q_ij = 1.0_f64
   call sll_s_convert_charge_to_rho_2d_per_per_cs(charge_cs, rho)
   call sll_s_field_accumulator_cs_2d_init(e_accumulator_cs, mesh_2d)
   call sll_s_accumulate_field_cs(e_1, e_2, e_accumulator_cs)

   s = sum(rho)
   if (abs(s - ref) > 1d-7) stop 'FAILED'
   s = sum(e_1)
   if (abs(s) > 1d-7) stop 'FAILED'
   s = sum(e_2)
   if (abs(s) > 1d-7) stop 'FAILED'

   print *, 'PASSED'
end program test_charge_to_density
