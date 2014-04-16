!> @brief 
!> Selalib 2D poisson solver for cartesian coordinates, Dirichlet BC's
!   
!> @authors                    
!> Nhung PHAM 
!> Edwin CHACON-GOLCHER
!                                  
!> @details
!> Add in here all sequential poisson solvers that use Dirichlet boundary
!> conditions.
!**************************************************************************

module sll_poisson_2d_dirichlet_cartesian

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
!  use sll_fft
  use sll_constants
   implicit none

  !> Structure to store data from Poisson solver. This
  !> solver is parallel on structured cartesian mesh. 
  type poisson_2d_dirichlet_cartesian
     sll_int32                           :: ncx    !< number of cells  
     sll_int32                           :: ncy    !< number of cells
     sll_real64                          :: Lx     !< domain length 
     sll_real64                          :: Ly     !< domain length
  end type poisson_2d_dirichlet_cartesian

  interface sll_delete
     module procedure delete_poisson_2d_dirichlet_cartesian
  end interface sll_delete

contains

  function new_poisson_2d_dirichlet_cartesian_plan( &
    ncx, &            
    ncy, &            
    Lx, &    
    Ly ) result(plan)

    type (poisson_2d_dirichlet_cartesian), pointer :: plan
    sll_int32                        :: ncx          !< number of cells in x
    sll_int32                        :: ncy          !< number of cells in y
    sll_real64                       :: Lx           !< length x
    sll_real64                       :: Ly           !< length y
    sll_int32     :: ierr
    SLL_ALLOCATE(plan, ierr)

    plan%ncx = ncx
    plan%ncy = ncy
    plan%Lx  = Lx
    plan%Ly  = Ly

  end function new_poisson_2d_dirichlet_cartesian_plan




  !> Note that the equation that is solved is: \f$ \Delta \phi = \rho \f$
  !> Thus the user is responsible for giving the proper sign to the source term.
  subroutine solve_poisson_2d_dirichlet_cartesian(plan, rho, phi)
    type (poisson_2d_dirichlet_cartesian), pointer :: plan !< self object
    sll_real64, dimension(:,:)        :: rho      !< charge density
    sll_real64, dimension(:,:)        :: phi      !< electric potential
    sll_int32                         :: ncx      !< global size
    sll_int32                         :: ncy      !< global size

    sll_int32                         :: i, j
    sll_int32                         :: ii, jj,p,k,nsky,l
    sll_real64                        :: kx, ky
    !sll_real64, dimension(:,:)          :: vkgs
    !sll_int32, dimen                  :: indi,indj


  end subroutine solve_poisson_2d_dirichlet_cartesian


!> Delete the Poisson solver object
  subroutine delete_poisson_2d_dirichlet_cartesian(plan)
    type (poisson_2d_dirichlet_cartesian), pointer :: plan
    sll_int32                                           :: ierr

    if( .not. associated(plan) ) then
       print *, 'ERROR, delete_poisson_2d_dirichlet_cartesian_plan(): ', &
            'passed plan is not associated.'
       STOP
    end if

    SLL_DEALLOCATE(plan, ierr)
  end subroutine delete_poisson_2d_dirichlet_cartesian



end module sll_poisson_2d_dirichlet_cartesian
