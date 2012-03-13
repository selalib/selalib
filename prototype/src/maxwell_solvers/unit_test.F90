program testMaxwell
  !-------------------------------------------------------------------
  !  test 1D Maxwell solver based on FFT
  !-------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
!#include "sll_maxwell_solvers.h"
use numeric_constants

use sll_maxwell_2d

implicit none

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

sll_int32  :: nc_eta1, nc_eta2
sll_int32  :: error

type(maxwell_2d),          pointer :: maxwell_TE
type(geometry_2D),         pointer :: geom
type(mesh_descriptor_2D),  pointer :: mesh
type(field_2D_vec3),       pointer :: ExEyHz

sll_int32                          :: i, j
sll_real64                         :: omega
sll_real64                         :: time
sll_int32                          :: istep, nstep
sll_real64                         :: err_l2
sll_real64                         :: dt

sll_real64, dimension(:,:), allocatable :: bz

sll_real64, parameter              :: c = 1.0_f64
sll_real64, parameter              :: cfl = 0.5_f64
sll_int32,  parameter              :: mode = 2

character(len=4)                   :: counter

eta1_min = .0_f64; eta1_max = 1.0_f64
eta2_min = .0_f64; eta2_max = 1.0_f64

geom => new_geometry_2D ('cartesian')

nc_eta1 = 127; nc_eta2 = 127

mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
        PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

delta_eta1 = mesh%delta_eta1
delta_eta2 = mesh%delta_eta2

call write_mesh_2D(mesh)

ExEyHz      => new_field_2D_vec3(mesh)

SLL_ASSERT(associated(exeyhz))

maxwell_TE  => new(ExEyHz)

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2)) / c
nstep = 100

time  = 0.

omega = c * sqrt( (mode*sll_pi/(nc_eta1*delta_eta1))**2   &
      &          +(mode*sll_pi/(nc_eta2*delta_eta2))**2)

SLL_ALLOCATE(bz(nc_eta1+1,nc_eta2+1), error)

do istep = 1, nstep !*** Loop over time

   if (istep == 1) then
      ExEyHz%data%v1 = 0.0_f64
      ExEyHz%data%v2 = 0.0_f64
   end if

   time = time + 0.5_f64*dt

   do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      bz(i,j) =   - cos(mode*sll_pi*(i-0.5_f64)/nc_eta1)    &
                  * cos(mode*sll_pi*(j-0.5_f64)/nc_eta2)    &
                  * cos(omega*time)
   end do  
   end do  

   if (istep == 1) ExEyHz%data(:,:)%v3 = bz

   call solve(maxwell_TE, dt)

   time = time + 0.5_f64*dt

   err_l2 = maxval(abs(maxwell_TE%fields%data%v3 - bz))
   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g12.3,' sec')",advance="no") time
   write(*,"(' erreur L2 = ',g10.5)") sqrt(err_l2)

   call int2string(istep, counter)
   call write_vec1d(ExEyHz%data%v3,mesh%nc_eta1+1,mesh%nc_eta2+1,"hz"//counter,"mesh",0)

end do ! next time step

call delete(maxwell_TE)
call delete_field_2D_vec3( ExEyHz )

end program testMaxwell

