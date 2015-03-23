!> @ingroup descriptors
!> @brief Describes different global flags throughout the library
!> @details
!> The intent of this module is to provide a single, library-wide definition
!> of the names used to describe different properties which are applicabale
!> in a global sense
!> <b> How to use-it </b>

module sll_descriptors
#include "sll_working_precision.h"

  implicit none
  
   !
 !SLL_LANDAU_SUM
 !SLL_LANDAU_PROD
 !SLL_TWOSTREAM
 !SLL_BUMPONTAIL
 !SLL_KEEN

  
 
  !>Vlasov - Poisson simulation descriptors
  sll_int32, parameter :: SLL_LANDAU_DIAG    =1
  !>Landau damping
  sll_int32, parameter :: SLL_LANDAU_SUM     =2
  !>Landau damping
  sll_int32, parameter :: SLL_LANDAU_PROD    =3
  !>Two stream
  sll_int32, parameter :: SLL_TWOSTREAM      =4
  !>Bump on tail
  sll_int32, parameter :: SLL_BUMPONTAIL     =5
  
  !>Key - Poisson simulation descriptors
  character(len=*), parameter :: &
	  sll_vp_descriptor_key(1:5) = &
         (/"SLL_LANDAU_DIAG    ",&
           "SLL_LANDAU_SUM     ",&
           "SLL_LANDAU_PROD    ",&
           "SLL_TWOSTREAM      ",&
           "SLL_BUMPONTAIL     "/)
  
   
!   type, abstract :: sll_descriptor
!     sll_int32 :: id
! !     character(len=*) :: name
!   end type sll_descriptor

! 
! !>Vlasov - Poisson simulation descriptors
! type :: sll_vp_descriptor
!      sll_int32, private :: id_priv
!      character(len=32), private :: name_priv
!      contains
! !      procedure, pass(this) :: id=>
! !      procedure, pass(this) :: names=>
! end type sll_vp_descriptor
!
!   type(sll_vp_descriptor), parameter :: SLL_LANDAU_DIAG=sll_vp_descriptor(1,"SLL_LANDAU_DIAG")

 
 
!  
!  contains
!  
! 
!  function sll_vp_descriptor( id, name ) result(desc)
!   sll_int32, intent(in) :: id
!   character(len=*), intent(in) :: name
!   type(sll_vp_descriptor) :: desc
!   
!   desc%id_priv=id
!   desc%name_priv=trim(name)
!  end function sll_vp_descriptor
!   
  
 

end module sll_descriptors
