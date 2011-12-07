program testPoisson
  !-------------------------------------------------------------------
  !  test 1D Poisson solver based on FFT
  !-------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
use numeric_constants

!use geometry1d_module
use sll_poisson_1D_periodic
use sll_poisson_2D_periodic

implicit none

type (mesh_descriptor_1D), pointer    :: geomx ! 1D mesh
type (poisson1dp)         :: poiss1dp !champ electrique
type(field_1D_vec1), pointer :: ex, rho, ex_exact

sll_int32      :: iflag, mode, i, ncx 
sll_real64     :: xmin, xmax, dx
sll_real64     :: eta1_max, eta1_min, eta2_max, eta2_min
sll_int32      :: nc_eta1, nc_eta2
sll_real64     :: delta_eta1, delta_eta2

type(geometry_2D), pointer :: geom
type(mesh_descriptor_2D), pointer :: mesh
type(field_2D_vec2), pointer :: u, u_exact
type(field_2D_vec1), pointer :: rhs

!PN!
!PN!!Solveur de Poisson1D, commente car ne parche pas
!PN!
!PN!logical, parameter :: per = .true.
!PN! initialisation of 1D periodic mesh 
!PN!xmin = 0.0; xmax = 2*sll_pi;
!PN!ncx = 128
!PN!
!PN!geomx    => new_mesh_descriptor_1D( xmin, xmax, ncx )
!PN!rho      => new_field_1D_vec1( geomx )
!PN!ex       => new_field_1D_vec1( geomx )
!PN!ex_exact => new_field_1D_vec1( geomx )
!PN!
!PN!! intialize Poisson solver
!PN!call new(poiss1dp,ncx,iflag) 
!PN!
!PN!! set rho
!PN!mode = 7
!PN!dx = GET_FIELD_DELTA_X1( rho )
!PN!do i=1,ncx+1
!PN!   FIELD_1D_AT_I(rho,i) =  mode**2*sin(mode*(i-1)*dx)
!PN!   FIELD_1D_AT_I(ex_exact,i) = -mode*cos(mode*(i-1)*dx)
!PN!end do
!PN!! compute electric field
!PN!call solve(poiss1dp,ex,rho)
!PN!     
!PN!! check solution
!PN!print*,'mode=',mode,'   error=',maxval(abs(FIELD_DATA(ex)-FIELD_DATA(ex_exact)))


eta1_min =  -8.0_f64; eta1_max =  8.0_f64
eta2_min =  -8.0_f64; eta2_max =  8.0_f64 

geom => new_geometry_2D ('cartesian')

nc_eta1 = 100; nc_eta2 = 100

mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
        PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

rhs     => new_field_2D_vec1(mesh)
u       => new_field_2D_vec2(mesh)
u_exact => new_field_2D_vec2(mesh)

end program testPoisson

