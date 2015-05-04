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
  
 
  !>Vlasov - Poisson simulation descriptors
  type sll_vlasovpoisson_sim
    sll_int32                  :: id
    character(len=32), private :: pname
  contains
    procedure, pass(self)      :: name=>name_sll_vlasovpoisson_sim
    procedure, pass(self)      :: parse=>parse_sll_vlasovpoisson_sim
  end type sll_vlasovpoisson_sim

  type(sll_vlasovpoisson_sim), parameter :: SLL_LANDAU_DIAG = sll_vlasovpoisson_sim(1,"SLL_LANDAU_DIAG")
  type(sll_vlasovpoisson_sim), parameter :: SLL_LANDAU_SUM  = sll_vlasovpoisson_sim(2,'SLL_LANDAU_SUM' )
  type(sll_vlasovpoisson_sim), parameter :: SLL_LANDAU_PROD = sll_vlasovpoisson_sim(3,'SLL_LANDAU_PROD' )
  type(sll_vlasovpoisson_sim), parameter :: SLL_TWOSTREAM   = sll_vlasovpoisson_sim(4,'SLL_TWOSTREAM' )
  type(sll_vlasovpoisson_sim), parameter :: SLL_BUMPONTAIL  = sll_vlasovpoisson_sim(5,'SLL_BUMPONTAIL' )

  interface operator(.eq.)
    module procedure sll_vlasovpoisson_sim_compare
  end interface
  
!==============================================================================
contains
!==============================================================================
!   function id_sll_vlasovpoisson_sim( self ) result( r )
!     class( sll_vlasovpoisson_sim ), intent( in ) :: self
!     sll_int32 :: r
!     r = self%pid
!   end function id_sll_vlasovpoisson_sim

  !----------------------------------------------------------------------------
  pure function name_sll_vlasovpoisson_sim( self ) result( r )
    class( sll_vlasovpoisson_sim ), intent( in ) :: self
    character(len=len(self%pname))     :: r
    r = self%pname
  end function name_sll_vlasovpoisson_sim
  !----------------------------------------------------------------------------
  
   pure function sll_vlasovpoisson_sim_compare(bc1,bc2) result(compare)
  type(sll_vlasovpoisson_sim), intent(in) :: bc1, bc2
  logical :: compare

   compare = bc1%id .eq. bc2%id

 end function sll_vlasovpoisson_sim_compare
  
  
  
  subroutine parse_sll_vlasovpoisson_sim(self, str)
      class( sll_vlasovpoisson_sim ), intent( inout ) :: self
      character(len=*), intent(in) :: str
      character(len=len(str)) :: strc
      
      !Remove blanks left and right
      strc=adjustl(trim(str))
      
      if (strc .eq. SLL_LANDAU_DIAG%name()) then
        self%id=SLL_LANDAU_DIAG%id
        self%pname=SLL_LANDAU_DIAG%pname
      elseif (strc .eq. SLL_LANDAU_SUM%name()) then
        self%id=SLL_LANDAU_SUM%id
        self%pname=SLL_LANDAU_SUM%pname
      elseif (strc .eq. SLL_LANDAU_PROD%name()) then
        self%id=SLL_LANDAU_PROD%id
        self%pname=SLL_LANDAU_PROD%pname
      elseif (strc .eq. SLL_TWOSTREAM%name()) then
        self%id=SLL_TWOSTREAM%id
        self%pname=SLL_TWOSTREAM%pname
      elseif (strc .eq. SLL_BUMPONTAIL%name()) then
        self%id=SLL_BUMPONTAIL%id
        self%pname=SLL_BUMPONTAIL%pname
!       elseif (str .eq. %name()) then
!         self%pid= %pid
!         self%pname=%pname
      endif
!        SELECT CASE (trim(str))
!              CASE(SLL_LANDAU_DIAG%name())
!                self=SLL_LANDAU_DIAG
!             CASE(SLL_LANDAU_SUM%name())
!                self=SLL_LANDAU_SUM
!             CASE(SLL_LANDAU_PROD%name())
!                self=SLL_LANDAU_PROD
!             CASE DEFAULT
!               print *, "Type not found"
!        END SELECT
!   
  end subroutine parse_sll_vlasovpoisson_sim
  
end module sll_descriptors




!==============================================================================
  
  
!   sll_int32, parameter :: SLL_LANDAU_DIAG    =1
!   !>Landau damping
!   sll_int32, parameter :: SLL_LANDAU_SUM     =2
!   !>Landau damping
!   sll_int32, parameter :: SLL_LANDAU_PROD    =3
!   !>Two stream
!   sll_int32, parameter :: SLL_TWOSTREAM      =4
!   !>Bump on tail
!   sll_int32, parameter :: SLL_BUMPONTAIL     =5
!   
!   !>Key - Poisson simulation descriptors
!   character(len=*), parameter :: &
! 	  sll_vp_descriptor_key(1:5) = &
!          (/"SLL_LANDAU_DIAG    ",&
!            "SLL_LANDAU_SUM     ",&
!            "SLL_LANDAU_PROD    ",&
!            "SLL_TWOSTREAM      ",&
!            "SLL_BUMPONTAIL     "/)
!   
  
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
! end module sll_descriptors
