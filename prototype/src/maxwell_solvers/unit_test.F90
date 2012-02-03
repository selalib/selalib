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


sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
sll_int32   :: nc_eta1, nc_eta2
sll_int32   :: error

type(maxwell_2d), pointer :: maxwell_TE
type(geometry_2D),         pointer :: geom
type(mesh_descriptor_2D),  pointer :: mesh
type(field_2D_vec3),       pointer :: ExEyHz
sll_real64                         :: x1, x2
sll_int32                          :: mode
sll_int32                          :: i, j

eta1_min = .0_f64; eta1_max = 1.0_f64
eta2_min = .0_f64; eta2_max = 1.0_f64

geom => new_geometry_2D ('cartesian')

nc_eta1 = 127; nc_eta2 = 127

mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
        PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

call write_mesh_2D(mesh)

ExEyHz       => new_field_2D_vec3(mesh)

maxwell_TE  => new_maxwell_2d(ExEyHz)

mode = 2
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      x1 = eta1_min+(i-1)*mesh%delta_eta1
      x2 = eta2_min+(j-1)*mesh%delta_eta2
      ExEyHz%data(i,j)%v1=1
      ExEyHz%data(i,j)%v2=0
      ExEyHz%data(i,j)%v3=1

   end do
end do

call solve_maxwell_2d(maxwell_TE,error)

write(*,*) " Ex Error = " , 0

call delete_maxwell_2d(maxwell_TE)
call delete_field_2D_vec3( ExEyHz )


end program testMaxwell

