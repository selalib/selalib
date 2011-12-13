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
!use sll_poisson_1D_periodic
use sll_poisson_2D_periodic

implicit none

!PN!type (mesh_descriptor_1D), pointer    :: geomx ! 1D mesh
!PN!type (poisson1dp)         :: poiss1dp !champ electrique
!PN!type(field_1D_vec1), pointer :: ex, rho, ex_exact
!PN!sll_int32   :: iflag, mode, i, ncx 
!PN!sll_real64  :: xmin, xmax, dx

sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
sll_int32   :: nc_eta1, nc_eta2
sll_int32   :: error

type(poisson2dp)                  :: poisson_obj
type(geometry_2D),        pointer :: geom
type(mesh_descriptor_2D), pointer :: mesh
type(field_2D_vec2),      pointer :: exy, exy_exact
type(field_2D_vec1),      pointer :: rho, phi
sll_real64                        :: x1, x2
sll_int32                         :: mode

sll_int32                         :: i, j

eta1_min =  .0_f64; eta1_max =  2.0_f64*sll_pi
eta2_min =  .0_f64; eta2_max =  2.0_f64*sll_pi

geom => new_geometry_2D ('cartesian')

nc_eta1 = 127; nc_eta2 = 127

mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
        PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

call write_mesh_2D(mesh)

rho       => new_field_2D_vec1(mesh)
phi       => new_field_2D_vec1(mesh)
exy       => new_field_2D_vec2(mesh)
exy_exact => new_field_2D_vec2(mesh)

call new(poisson_obj,rho%data,mesh,error)

mode = 2
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      x1 = eta1_min+(i-1)*mesh%delta_eta1
      x2 = eta2_min+(j-1)*mesh%delta_eta2
      rho%data(i,j) = -2_f64 * mode**3 * sin(mode*x1)*cos(mode*x2)
      exy%data(i,j)%v1 =  mode**2*cos(mode*x1)*cos(mode*x2)
      exy%data(i,j)%v2 = -mode**2*sin(mode*x1)*sin(mode*x2)
   end do
end do

call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho0","mesh",0)
call write_vec2d(exy%data%v1, exy%data%v2,mesh%nc_eta1+1,mesh%nc_eta2+1,"exy0","mesh",0)

exy_exact%data%v1 = exy%data%v1
exy_exact%data%v2 = exy%data%v2

exy%data%v1 = 0.0
exy%data%v2 = 0.0

call solve(poisson_obj,exy%data%v1,exy%data%v2,rho%data,phi%data,error)

call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho1","mesh",0)
call write_vec2d(exy%data%v1, exy%data%v2,mesh%nc_eta1+1,mesh%nc_eta2+1,"exy1","mesh",0)

write(*,*) " Ex Error = " , maxval(abs(exy_exact%data%v1-exy%data%v1))
write(*,*) " Ey Error = " , maxval(abs(exy_exact%data%v2-exy%data%v2))

call free(poisson_obj)

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



end program testPoisson

