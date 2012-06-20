module fishpack

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none
integer, private :: i, j, k

type, public :: poisson_fishpack
   sll_int32                             :: nc_eta1
   sll_int32                             :: nc_eta2
   sll_int32                             :: nc_eta3
   sll_real64                            :: eta1_min, eta1_max
   sll_real64                            :: eta2_min, eta2_max
   sll_real64                            :: eta3_min, eta3_max
   sll_real64, allocatable, dimension(:) :: bda
   sll_real64, allocatable, dimension(:) :: bdb
   sll_real64, allocatable, dimension(:) :: bdc
   sll_real64, allocatable, dimension(:) :: bdd
   sll_int32                             :: mbdcnd, nbdcnd
   sll_real64                            :: elmbda, pertrb
   sll_int32                             :: error
   sll_int32                             :: geometry
end type poisson_fishpack

enum, bind(C)
   enumerator :: CARTESIAN_2D = 0, POLAR_2D = 1
end enum


contains


  subroutine new_poisson_2d_fishpack(this, geometry,            &
                                     eta1_min,eta1_max,nc_eta1, &
                                     eta2_min,eta2_max,nc_eta2)

    type(poisson_fishpack),intent(out) :: this
    sll_int32, intent(in)              :: nc_eta1
    sll_int32, intent(in)              :: nc_eta2
    sll_real64, intent(in)             :: eta1_min, eta1_max
    sll_real64, intent(in)             :: eta2_min, eta2_max
    sll_int32                          :: geometry
    !sll_int32                          :: mbdcnd, nbdcnd
    !sll_real64                         :: elmbda

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

    SLL_ALLOCATE(this%bda(nc_eta2+1),this%error)
    SLL_ALLOCATE(this%bdb(nc_eta2+1),this%error)
    SLL_ALLOCATE(this%bdc(nc_eta1+1),this%error)
    SLL_ALLOCATE(this%bdd(nc_eta1+1),this%error)

  end subroutine new_poisson_2d_fishpack

  subroutine new_poisson_3d_fishpack(this, geometry,            &
                                     eta1_min,eta1_max,nc_eta1, &
                                     eta2_min,eta2_max,nc_eta2, &
                                     eta3_min,eta3_max,nc_eta3  )

    type(poisson_fishpack),intent(out) :: this
    sll_int32, intent(in)              :: nc_eta1
    sll_int32, intent(in)              :: nc_eta2
    sll_int32, intent(in)              :: nc_eta3
    sll_real64, intent(in)             :: eta1_min, eta1_max
    sll_real64, intent(in)             :: eta2_min, eta2_max
    sll_real64, intent(in)             :: eta3_min, eta3_max
    sll_int32,  intent(in)             :: geometry

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

  end subroutine new_poisson_3d_fishpack


  subroutine solve_poisson_fishpack(this, field)

     implicit none
     type(poisson_fishpack),intent(in) :: this
     sll_real64, dimension(:,:)        :: field
     sll_real64                        :: w
     sll_int32                         :: idimf
  
  
     !     auxiliary quantities.
  
     idimf = this%nc_eta1+1

     if ( this%geometry == CARTESIAN_2D ) then
  
        call hwscrt (this%eta1_min, this%eta1_max, this%nc_eta1, &
                     this%mbdcnd, this%bda, this%bdb, &
                     this%eta2_min, this%eta2_max, this%nc_eta2, &
                     this%nbdcnd, this%bdc, this%bdd, &
                     this%elmbda, field, idimf, this%pertrb, this%error)

     else if ( this%geometry == POLAR_2D ) then


        call hwsplr (this%eta1_min, this%eta1_max, this%nc_eta1, &
                     this%mbdcnd, this%bda, this%bdb, &
                     this%eta2_min, this%eta2_max, this%nc_eta2, &
                     this%nbdcnd, this%bdc, this%bdd, &
                     this%elmbda, field, idimf, this%pertrb, this%error, w)
     else
 
        write(*,*) " This kind of geometry is not implemented "

     end if
  
  end subroutine solve_poisson_fishpack
  
end module fishpack
