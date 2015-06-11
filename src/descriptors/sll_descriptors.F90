!This file was generated using the macro descriptors.m4
!Usage: m4 descriptors.m4 > sll_descriptors.F90
!> @ingroup descriptors
!> @brief Describes different global flags throughout the library
!> @details
!> The intent of this module is to provide a single, library-wide definition
!> of the names used to describe different properties which are applicabale
!> in a global sense
!> <b> How to use-it </b>









module sll_descriptors
#include "sll_working_precision.h"
type sll_vlasovpoisson_sim 
    sll_int32                  :: id
    character(len=32), private :: pname
  contains
    procedure,  pass(self)      :: name=>name_vlasovpoisson_sim 
    procedure,  pass(self)      :: parse=>parse_vlasovpoisson_sim
  end type sll_vlasovpoisson_sim
  
  interface operator(.eq.)
    module procedure sll_vlasovpoisson_sim_compare
  end interface 
 
type(sll_vlasovpoisson_sim), parameter :: SLL_LANDAU_DIAG = sll_vlasovpoisson_sim(1,"SLL_LANDAU_DIAG") 
type(sll_vlasovpoisson_sim), parameter :: SLL_LANDAU_SUM = sll_vlasovpoisson_sim(2,"SLL_LANDAU_SUM") 
type(sll_vlasovpoisson_sim), parameter :: SLL_LANDAU_PROD = sll_vlasovpoisson_sim(3,"SLL_LANDAU_PROD") 
type(sll_vlasovpoisson_sim), parameter :: SLL_TWOSTREAM = sll_vlasovpoisson_sim(4,"SLL_TWOSTREAM") 
type(sll_vlasovpoisson_sim), parameter :: SLL_BUMPONTAIL = sll_vlasovpoisson_sim(5,"SLL_BUMPONTAIL") 

type sll_field_eqn 
    sll_int32                  :: id
    character(len=32), private :: pname
  contains
    procedure,  pass(self)      :: name=>name_field_eqn 
    procedure,  pass(self)      :: parse=>parse_field_eqn
  end type sll_field_eqn
  
  interface operator(.eq.)
    module procedure sll_field_eqn_compare
  end interface 
 
type(sll_field_eqn), parameter :: POISSON = sll_field_eqn(1,"POISSON") 
type(sll_field_eqn), parameter :: AMPERE = sll_field_eqn(2,"AMPERE") 
type(sll_field_eqn), parameter :: MAXWELL = sll_field_eqn(3,"MAXWELL") 
type(sll_field_eqn), parameter :: ADIABATIC_WO_ZONAL = sll_field_eqn(4,"ADIABATIC_WO_ZONAL") 
type(sll_field_eqn), parameter :: ADIABATIC = sll_field_eqn(5,"ADIABATIC") 
type(sll_field_eqn), parameter :: QN_ADIABATIC = sll_field_eqn(6,"QN_ADIABATIC") 


contains


 !----------------------------------------------------------------------
 !------------------- vlasovpoisson_sim -----------------------------------------
 pure function name_vlasovpoisson_sim( self ) result( r )
   class( sll_vlasovpoisson_sim ), intent( in ) :: self
   character(len=len ( self%pname ) )     :: r
   r = self%pname
 end function name_vlasovpoisson_sim  
  
 pure function sll_vlasovpoisson_sim_compare(bc1,bc2) result(compare)
  type(sll_vlasovpoisson_sim), intent(in) :: bc1, bc2
  logical :: compare
   compare = bc1%id .eq. bc2%id
 end function sll_vlasovpoisson_sim_compare

 subroutine parse_vlasovpoisson_sim(self, str)
    class( sll_vlasovpoisson_sim ), intent( inout ) :: self
    character(len=*), intent(in) :: str
    character(len= len ( str ) ) :: strc
      
     !Remove blanks left and right
     strc=adjustl(trim(str))
     
    if (strc .eq. SLL_LANDAU_DIAG%name()) then  
        self%id=SLL_LANDAU_DIAG%id
        self%pname=SLL_LANDAU_DIAG%pname
      elseif(strc .eq. SLL_LANDAU_SUM%name()) then  
        self%id=SLL_LANDAU_SUM%id
        self%pname=SLL_LANDAU_SUM%pname
      elseif(strc .eq. SLL_LANDAU_PROD%name()) then  
        self%id=SLL_LANDAU_PROD%id
        self%pname=SLL_LANDAU_PROD%pname
      elseif(strc .eq. SLL_TWOSTREAM%name()) then  
        self%id=SLL_TWOSTREAM%id
        self%pname=SLL_TWOSTREAM%pname
      elseif(strc .eq. SLL_BUMPONTAIL%name()) then  
        self%id=SLL_BUMPONTAIL%id
        self%pname=SLL_BUMPONTAIL%pname
      elseif  (.true.) then
    endif
    
 end subroutine parse_vlasovpoisson_sim
 !----------------------------------------------------------------------
      

 !----------------------------------------------------------------------
 !------------------- field_eqn -----------------------------------------
 pure function name_field_eqn( self ) result( r )
   class( sll_field_eqn ), intent( in ) :: self
   character(len=len ( self%pname ) )     :: r
   r = self%pname
 end function name_field_eqn  
  
 pure function sll_field_eqn_compare(bc1,bc2) result(compare)
  type(sll_field_eqn), intent(in) :: bc1, bc2
  logical :: compare
   compare = bc1%id .eq. bc2%id
 end function sll_field_eqn_compare

 subroutine parse_field_eqn(self, str)
    class( sll_field_eqn ), intent( inout ) :: self
    character(len=*), intent(in) :: str
    character(len= len ( str ) ) :: strc
      
     !Remove blanks left and right
     strc=adjustl(trim(str))
     
    if (strc .eq. POISSON%name()) then  
        self%id=POISSON%id
        self%pname=POISSON%pname
      elseif(strc .eq. AMPERE%name()) then  
        self%id=AMPERE%id
        self%pname=AMPERE%pname
      elseif(strc .eq. MAXWELL%name()) then  
        self%id=MAXWELL%id
        self%pname=MAXWELL%pname
      elseif(strc .eq. ADIABATIC_WO_ZONAL%name()) then  
        self%id=ADIABATIC_WO_ZONAL%id
        self%pname=ADIABATIC_WO_ZONAL%pname
      elseif(strc .eq. ADIABATIC%name()) then  
        self%id=ADIABATIC%id
        self%pname=ADIABATIC%pname
      elseif(strc .eq. QN_ADIABATIC%name()) then  
        self%id=QN_ADIABATIC%id
        self%pname=QN_ADIABATIC%pname
      elseif  (.true.) then
    endif
    
 end subroutine parse_field_eqn
 !----------------------------------------------------------------------
      
end module sll_descriptors
