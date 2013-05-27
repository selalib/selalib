module sll_mudpack

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none
sll_int32, private :: i, j, k

!> Fishpack solver cartesian 2d
type, public :: mudpack_2d

   sll_int32 :: iprm(16),mgopt(4)
   sll_real64 :: fprm(6)

   sll_int32 :: intl
   sll_int32 :: nxa
   sll_int32 :: nxb
   sll_int32 :: nyc
   sll_int32 :: nyd
   sll_int32 :: ixp
   sll_int32 :: jyq
   sll_int32 :: iex
   sll_int32 :: jey
   sll_int32 :: nx
   sll_int32 :: ny
   sll_int32 :: iguess
   sll_int32 :: maxcy
   sll_int32 :: method
   sll_int32 :: nwork
   sll_int32 :: lwrkqd
   sll_int32 :: itero

   sll_real64 :: xa
   sll_real64 :: xb
   sll_real64 :: yc
   sll_real64 :: yd
   sll_real64 :: tolmax
   sll_real64 :: relmax

end type mudpack_2d

enum, bind(C)
   enumerator :: CARTESIAN_2D = 0, POLAR_2D = 1, CARTESIAN_3D = 2
   enumerator :: PERIODIC=0, DIRICHLET=1
   enumerator :: NEUMANN_RIGHT=2, NEUMANN_LEFT=3, NEUMANN=4
end enum

contains

  !> Initialize the mudpack solver 2d.
  subroutine initialize(this, geometry,                    &
                        eta1_min,eta1_max,nc_eta1,bc_eta1, &
                        eta2_min,eta2_max,nc_eta2,bc_eta2)

    class(mudpack_2d),intent(out) :: this      !< Fishpack solver
    sll_int32, intent(in)         :: nc_eta1   !< x number of cells
    sll_int32, intent(in)         :: nc_eta2   !< y number of cells
    sll_real64, intent(in)        :: eta1_min  !< left side of the domain
    sll_real64, intent(in)        :: eta1_max  !< right side of the domain
    sll_real64, intent(in)        :: eta2_min  !< bottom side of the domain
    sll_real64, intent(in)        :: eta2_max  !< top side of the domain
    sll_int32                     :: geometry  !< polar or cartesian
    sll_int32                     :: bc_eta1   !< boundary condition parameter
    sll_int32                     :: bc_eta2   !< boundary condition parameter
    !sll_real64                   :: elmbda

!    this%geometry = geometry
!    ! Indicateur d'erreur
!    this%error = 0
!    ! Initialisation de la geometrie
!    this%nc_eta1=nc_eta1
!    this%nc_eta2=nc_eta2
!
!    this%eta1_min = eta1_min
!    this%eta2_min = eta2_min
!    this%eta1_max = eta1_max
!    this%eta2_max = eta2_max
!
!    this%bc_eta1 = bc_eta1
!    this%bc_eta2 = bc_eta2
!
!    if (bc_eta1 /= PERIODIC) then
!       SLL_ALLOCATE(this%bda(nc_eta2+1),this%error)
!       SLL_ALLOCATE(this%bdb(nc_eta2+1),this%error)
!    end if
!    if (bc_eta2 /= PERIODIC) then
!       SLL_ALLOCATE(this%bdc(nc_eta1+1),this%error)
!       SLL_ALLOCATE(this%bdd(nc_eta1+1),this%error)
!    end if

  end subroutine initialize

end module sll_mudpack
