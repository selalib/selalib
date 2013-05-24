module sll_mudpack

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none
integer, private :: i, j, k

!> Fishpack solver cartesian 2d
type, public :: mudpack_2d

   integer :: iprm(16),mgopt(4)
   real(8) :: fprm(6)
   integer :: intl
   integer :: nxa
   integer :: nxb
   integer :: nyc
   integer :: nyd
   integer :: ixp
   integer :: jyq
   integer :: iex
   integer :: jey
   integer :: nx
   integer :: ny
   integer :: iguess
   integer :: maxcy
   integer :: method
   integer :: nwork
   integer :: lwrkqd
   integer :: itero

   real(8) :: xa
   real(8) :: xb
   real(8) :: yc
   real(8) :: yd
   real(8) :: tolmax
   real(8) :: relmax

contains

   procedure :: initialize => new_2d
   procedure :: solve => solve_2d

end type mudpack_2d

enum, bind(C)
   enumerator :: CARTESIAN_2D = 0, POLAR_2D = 1, CARTESIAN_3D = 2
   enumerator :: PERIODIC=0, DIRICHLET=1
   enumerator :: NEUMANN_RIGHT=2, NEUMANN_LEFT=3, NEUMANN=4
end enum

contains

  !> Initialize the mudpack solver 2d.
  subroutine new_2d(this, geometry,            &
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

  end subroutine new_2d

  !> Solve routine for mudpack 2d solver
  subroutine solve_2d(this, field)

     implicit none
     class(mudpack_2d),intent(in)     :: this
     sll_real64, dimension(:,:)        :: field
     sll_real64                        :: w
     sll_int32                         :: idimf
     sll_int32                         :: ldimf
     sll_int32                         :: mdimf
  
!     !     auxiliary quantities.
!  
!     idimf = this%nc_eta1+1
!     ldimf = this%nc_eta1+1
!     mdimf = this%nc_eta2+1
!
!     if ( this%geometry == CARTESIAN_2D ) then
!  
!
!     else if ( this%geometry == POLAR_2D ) then
!
!
!     else
!
!        write(*,*) " This kind of geometry is not implemented "
!
!     end if
  
  end subroutine solve_2d

  
end module sll_mudpack
