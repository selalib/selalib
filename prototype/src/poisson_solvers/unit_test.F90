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

  implicit none

  type (mesh_descriptor_1D), pointer    :: geomx ! 1D mesh
  type (poisson1dp)         :: poiss1dp !champ electrique
  type(field_1D_vec1), pointer :: ex, rho, ex_exact

  sll_int32      :: iflag, mode, i, ncx 
  sll_real64     :: xmin, xmax, dx
  !logical, parameter :: per = .true.
  ! initialisation of 1D periodic mesh 
  xmin = 0.0; xmax = 2*sll_pi;
  ncx = 128

  geomx    => new_mesh_descriptor_1D( xmin, xmax, ncx )
  rho      => new_field_1D_vec1( geomx )
  ex       => new_field_1D_vec1( geomx )
  ex_exact => new_field_1D_vec1( geomx )

!  sll_real64, dimension(:), pointer    :: rho  ! charge density
!  sll_real64, dimension(:), pointer    :: ex ! electric field
!  sll_real64, dimension(:), pointer    :: ex_exact ! exact electric field
!    call new(geomx,x0,x1,nx,per,iflag)
!  if (iflag.ne.0) stop 'erreur dans l initialisation de geomx'

  ! allocation 
!  SLL_ALLOCATE(rho(geomx%nx+2),iflag)
!  SLL_ALLOCATE(ex(geomx%nx+2),iflag)
!  SLL_ALLOCATE(ex_exact(geomx%nx),iflag)

  ! intialize Poisson solver
  call new(poiss1dp,ncx,iflag) 

  ! set rho
  mode = 7
  dx = GET_FIELD_DELTA_X1( rho )
  do i=1,ncx+1
     FIELD_1D_AT_I(rho,i) =  mode**2*sin(mode*(i-1)*dx)
     FIELD_1D_AT_I(ex_exact,i) = -mode*cos(mode*(i-1)*dx)
  end do
  ! compute electric field
  call solve(poiss1dp,ex,rho)
     
  ! check solution
  print*,'mode=',mode,'   error=',maxval(abs(FIELD_DATA(ex)-FIELD_DATA(ex_exact)))

end program testPoisson

