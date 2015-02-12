!> @brief 
!> Interface module to use fishpack libarary
!> @details
!> You have to download and install fishpack 
!> http://www2.cisl.ucar.edu/resources/legacy/fishpack
module sll_fishpack

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none
integer, private :: i, j, k

!> Fishpack solver cartesian 2d
type, public :: fishpack_2d
   sll_int32                               :: nc_eta1  !< number of cells
   sll_int32                               :: nc_eta2  !< number of cells
   sll_real64                              :: eta1_min !< geometry parameter
   sll_real64                              :: eta1_max !< geometry parameter
   sll_real64                              :: eta2_min !< geometry parameter
   sll_real64                              :: eta2_max !< geometry parameter
   sll_real64, allocatable, dimension(:)   :: bda      !< boundary condition parameter
   sll_real64, allocatable, dimension(:)   :: bdb      !< boundary condition parameter
   sll_real64, allocatable, dimension(:)   :: bdc      !< boundary condition parameter
   sll_real64, allocatable, dimension(:)   :: bdd      !< boundary condition parameter
   sll_int32                               :: bc_eta1  !< boundary condition parameter
   sll_int32                               :: bc_eta2  !< boundary condition parameter
   sll_real64                              :: elmbda   !< fishpack parameter
   sll_real64                              :: pertrb   !< fishpack parameter
   sll_int32                               :: error    !< error code
   sll_int32                               :: geometry !< cartesian or polar
contains
   !> create the solver
   procedure :: initialize => new_2d   
   !> compute the potential
   procedure :: solve => solve_2d      
end type fishpack_2d

!> Fishpack solver cartesian 3d
type, public :: fishpack_3d
   sll_int32                               :: nc_eta1  !< number of cells
   sll_int32                               :: nc_eta2  !< number of cells
   sll_int32                               :: nc_eta3  !< number of cells
   sll_real64                              :: eta1_min !< geometry parameter
   sll_real64                              :: eta1_max !< geometry parameter
   sll_real64                              :: eta2_min !< geometry parameter
   sll_real64                              :: eta2_max !< geometry parameter
   sll_real64                              :: eta3_min !< geometry parameter
   sll_real64                              :: eta3_max !< geometry parameter
   sll_real64, allocatable, dimension(:,:) :: bda      !< boundary condition parameter
   sll_real64, allocatable, dimension(:,:) :: bdb      !< boundary condition parameter
   sll_real64, allocatable, dimension(:,:) :: bdc      !< boundary condition parameter
   sll_real64, allocatable, dimension(:,:) :: bdd      !< boundary condition parameter
   sll_real64, allocatable, dimension(:,:) :: bde      !< boundary condition parameter
   sll_real64, allocatable, dimension(:,:) :: bdf      !< boundary condition parameter
   sll_int32                               :: bc_eta1  !< boundary condition parameter
   sll_int32                               :: bc_eta2  !< boundary condition parameter
   sll_int32                               :: bc_eta3  !< boundary condition parameter
   sll_real64                              :: elmbda   !< fishpack parameter
   sll_real64                              :: pertrb   !< fishpack parameter
   sll_int32                               :: error    !< fishpack parameter
   sll_int32                               :: geometry !< cartesian or cylindrical
contains
   !> create the solver
   procedure :: initialize => new_3d    
   !> compute the potential
   procedure :: solve => solve_3d 
end type fishpack_3d

enum, bind(C)
   enumerator :: CARTESIAN_2D = 0, POLAR_2D = 1, CARTESIAN_3D = 2
   enumerator :: PERIODIC=0, DIRICHLET=1
   enumerator :: NEUMANN_RIGHT=2, NEUMANN_LEFT=3, NEUMANN=4
end enum

contains

  !> Initialize the fishpack solver 2d.
  subroutine new_2d(this, geometry,            &
                    eta1_min,eta1_max,nc_eta1,bc_eta1, &
                    eta2_min,eta2_max,nc_eta2,bc_eta2)

    class(fishpack_2d),intent(out) :: this      !< Fishpack solver
    sll_int32, intent(in)          :: nc_eta1   !< x number of cells
    sll_int32, intent(in)          :: nc_eta2   !< y number of cells
    sll_real64, intent(in)         :: eta1_min  !< left side of the domain
    sll_real64, intent(in)         :: eta1_max  !< right side of the domain
    sll_real64, intent(in)         :: eta2_min  !< bottom side of the domain
    sll_real64, intent(in)         :: eta2_max  !< top side of the domain
    sll_int32                      :: geometry  !< polar or cartesian
    sll_int32                      :: bc_eta1   !< boundary condition parameter
    sll_int32                      :: bc_eta2   !< boundary condition parameter
    !sll_real64                    :: elmbda

    this%geometry = geometry
    ! Indicateur d'erreur
    this%error = 0
    ! Initialisation de la geometrie
    this%nc_eta1=nc_eta1
    this%nc_eta2=nc_eta2

    this%eta1_min = eta1_min
    this%eta2_min = eta2_min
    this%eta1_max = eta1_max
    this%eta2_max = eta2_max

    this%bc_eta1 = bc_eta1
    this%bc_eta2 = bc_eta2

    if (bc_eta1 /= PERIODIC) then
       SLL_ALLOCATE(this%bda(nc_eta2+1),this%error)
       SLL_ALLOCATE(this%bdb(nc_eta2+1),this%error)
    end if
    if (bc_eta2 /= PERIODIC) then
       SLL_ALLOCATE(this%bdc(nc_eta1+1),this%error)
       SLL_ALLOCATE(this%bdd(nc_eta1+1),this%error)
    end if

  end subroutine new_2d

  !> Initialize the fishpack solver 3d.
  subroutine new_3d(this, geometry,            &
                    eta1_min,eta1_max,nc_eta1,bc_eta1, &
                    eta2_min,eta2_max,nc_eta2,bc_eta2, &
                    eta3_min,eta3_max,nc_eta3,bc_eta3  )

    class(fishpack_3d),intent(out) :: this
    sll_int32, intent(in)          :: nc_eta1
    sll_int32, intent(in)          :: nc_eta2
    sll_int32, intent(in)          :: nc_eta3
    sll_int32, intent(in)          :: bc_eta1
    sll_int32, intent(in)          :: bc_eta2
    sll_int32, intent(in)          :: bc_eta3
    sll_real64, intent(in)         :: eta1_min
    sll_real64, intent(in)         :: eta2_min
    sll_real64, intent(in)         :: eta3_min
    sll_real64, intent(in)         :: eta1_max
    sll_real64, intent(in)         :: eta2_max
    sll_real64, intent(in)         :: eta3_max
    sll_int32,  intent(in)         :: geometry

    this%geometry = geometry
    ! Indicateur d'erreur
    this%error = 0
    ! Initialisation de la geometrie
    this%nc_eta1=nc_eta1
    this%nc_eta2=nc_eta2
    this%nc_eta3=nc_eta3

    this%eta1_min = eta1_min
    this%eta2_min = eta2_min
    this%eta3_min = eta3_min
    this%eta1_max = eta1_max
    this%eta2_max = eta2_max
    this%eta3_max = eta3_max

    this%bc_eta1 = bc_eta1
    this%bc_eta2 = bc_eta2
    this%bc_eta3 = bc_eta3

    if (bc_eta1 /= PERIODIC) then
       SLL_ALLOCATE(this%bda(nc_eta2+1,nc_eta3+1),this%error)
       SLL_ALLOCATE(this%bdb(nc_eta2+1,nc_eta3+1),this%error)
    end if
    if (bc_eta2 /= PERIODIC) then
       SLL_ALLOCATE(this%bdc(nc_eta1+1,nc_eta3+1),this%error)
       SLL_ALLOCATE(this%bdd(nc_eta1+1,nc_eta3+1),this%error)
    end if
    if (bc_eta3 /= PERIODIC) then
       SLL_ALLOCATE(this%bde(nc_eta1+1,nc_eta2+1),this%error)
       SLL_ALLOCATE(this%bdf(nc_eta1+1,nc_eta2+1),this%error)
    end if

  end subroutine new_3d


  !> Solve routine for fishpack 2d solver
  subroutine solve_2d(this, field)

     implicit none
     class(fishpack_2d),intent(in)     :: this
     sll_real64, dimension(:,:)        :: field
     sll_real64                        :: w
     sll_int32                         :: idimf
     sll_int32                         :: ldimf
     sll_int32                         :: mdimf
  
  
     !     auxiliary quantities.
  
     idimf = this%nc_eta1+1
     ldimf = this%nc_eta1+1
     mdimf = this%nc_eta2+1

     if ( this%geometry == CARTESIAN_2D ) then
  
        call hwscrt (this%eta1_min, this%eta1_max, this%nc_eta1, &
                     this%bc_eta1, this%bda, this%bdb, &
                     this%eta2_min, this%eta2_max, this%nc_eta2, &
                     this%bc_eta2, this%bdc, this%bdd, &
                     this%elmbda, field, idimf, this%pertrb, this%error)

     else if ( this%geometry == POLAR_2D ) then


        call hwsplr (this%eta1_min, this%eta1_max, this%nc_eta1, &
                     this%bc_eta1, this%bda, this%bdb, &
                     this%eta2_min, this%eta2_max, this%nc_eta2, &
                     this%bc_eta2, this%bdc, this%bdd, &
                     this%elmbda, field, idimf, this%pertrb, this%error, w)

     else

        write(*,*) " This kind of geometry is not implemented "

     end if
  
  end subroutine solve_2d

  !> Solve routine for fishpack 3d solver
  subroutine solve_3d(this, field)

     implicit none
     class(fishpack_3d),intent(in)     :: this
     sll_real64, dimension(:,:,:)      :: field
     sll_real64                        :: w
     sll_int32                         :: ldimf
     sll_int32                         :: mdimf
  
     !     auxiliary quantities.
     ldimf = this%nc_eta1+1
     mdimf = this%nc_eta2+1

     if ( this%geometry == CARTESIAN_3D) then

        call hw3crt (this%eta1_min,this%eta1_max,this%nc_eta1, &
                     this%bc_eta1,this%bda,this%bdb, &
                     this%eta2_min,this%eta2_max,this%nc_eta2, &
                     this%bc_eta2,this%bdc,this%bdd, &
                     this%eta3_min,this%eta3_max,this%nc_eta3, &
                     this%bc_eta3,this%bde,this%bdf,this%elmbda, &
                     ldimf,mdimf,field, this%pertrb,this%error,w)

     else

        write(*,*) " This kind of geometry is not implemented "

     end if
  
  end subroutine solve_3d
  
end module sll_fishpack
