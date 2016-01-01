!This file was generated using the macro descriptors.m4
!Usage: m4 descriptors.m4 > sll_m_descriptors.F90
!> @ingroup descriptors
!> @brief Describes different global flags throughout the library
!> @details
!> The intent of this module is to provide a single, library-wide definition
!> of the names used to describe different properties which are applicabale
!> in a global sense
!> <b> How to use-it </b>








module sll_m_descriptors
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_p_landau_diag, &
    sll_p_landau_prod, &
    sll_p_landau_sum, &
    sll_t_vlasovpoisson_sim

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
type sll_t_vlasovpoisson_sim 
    sll_int32                  :: id
    character(len=32), private :: pname
  contains
    procedure,  pass(self)      :: name=>name_vlasovpoisson_sim 
    procedure,  pass(self)      :: parse=>parse_vlasovpoisson_sim
  end type sll_t_vlasovpoisson_sim
  
  interface operator(.eq.)
    module procedure sll_vlasovpoisson_sim_compare
  end interface 
 
type(sll_t_vlasovpoisson_sim), parameter :: sll_p_landau_diag = sll_t_vlasovpoisson_sim(1,"sll_p_landau_diag") 
type(sll_t_vlasovpoisson_sim), parameter :: sll_p_landau_sum = sll_t_vlasovpoisson_sim(2,"sll_p_landau_sum") 
type(sll_t_vlasovpoisson_sim), parameter :: sll_p_landau_prod = sll_t_vlasovpoisson_sim(3,"sll_p_landau_prod") 
type(sll_t_vlasovpoisson_sim), parameter :: SLL_TWOSTREAM = sll_t_vlasovpoisson_sim(4,"SLL_TWOSTREAM") 
type(sll_t_vlasovpoisson_sim), parameter :: SLL_BUMPONTAIL = sll_t_vlasovpoisson_sim(5,"SLL_BUMPONTAIL") 

type sll_boundary 
    sll_int32                  :: id
    character(len=32), private :: pname
  contains
    procedure,  pass(self)      :: name=>name_boundary 
    procedure,  pass(self)      :: parse=>parse_boundary
  end type sll_boundary
  
  interface operator(.eq.)
    module procedure sll_boundary_compare
  end interface 
 
type(sll_boundary), parameter :: OPEN = sll_boundary(1,"OPEN") 
type(sll_boundary), parameter :: CLOSED = sll_boundary(2,"CLOSED") 


contains


 !----------------------------------------------------------------------
 !------------------- vlasovpoisson_sim -----------------------------------------
 pure function name_vlasovpoisson_sim( self ) result( r )
   class( sll_t_vlasovpoisson_sim ), intent( in ) :: self
   character(len=10)     :: r
   r = self%pname(1:10)
 end function name_vlasovpoisson_sim  
  
 pure function sll_vlasovpoisson_sim_compare(bc1,bc2) result(compare)
  type(sll_t_vlasovpoisson_sim), intent(in) :: bc1, bc2
  logical :: compare
   compare = bc1%id .eq. bc2%id
 end function sll_vlasovpoisson_sim_compare

 subroutine parse_vlasovpoisson_sim(self, str)
    class( sll_t_vlasovpoisson_sim ), intent( inout ) :: self
    character(len=*), intent(in) :: str
    character(len=3) :: strc
      
     !Remove blanks left and right
     strc=adjustl(trim(str))
     
    if (strc .eq. sll_p_landau_diag%name()) then  
        self%id=sll_p_landau_diag%id
        self%pname=sll_p_landau_diag%pname
      elseif(strc .eq. sll_p_landau_sum%name()) then  
        self%id=sll_p_landau_sum%id
        self%pname=sll_p_landau_sum%pname
      elseif(strc .eq. sll_p_landau_prod%name()) then  
        self%id=sll_p_landau_prod%id
        self%pname=sll_p_landau_prod%pname
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
 !------------------- boundary -----------------------------------------
 pure function name_boundary( self ) result( r )
   class( sll_boundary ), intent( in ) :: self
   character(len=10)     :: r
   r = self%pname(1:10)
 end function name_boundary  
  
 pure function sll_boundary_compare(bc1,bc2) result(compare)
  type(sll_boundary), intent(in) :: bc1, bc2
  logical :: compare
   compare = bc1%id .eq. bc2%id
 end function sll_boundary_compare

 subroutine parse_boundary(self, str)
    class( sll_boundary ), intent( inout ) :: self
    character(len=*), intent(in) :: str
    character(len=3) :: strc
      
     !Remove blanks left and right
     strc=adjustl(trim(str))
     
    if (strc .eq. OPEN%name()) then  
        self%id=OPEN%id
        self%pname=OPEN%pname
      elseif(strc .eq. CLOSED%name()) then  
        self%id=CLOSED%id
        self%pname=CLOSED%pname
      elseif  (.true.) then
    endif
    
 end subroutine parse_boundary
 !----------------------------------------------------------------------
      

end module sll_m_descriptors
