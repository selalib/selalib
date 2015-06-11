dnl Comment for the Generated file
!This file was generated using the macro descriptors.m4
!Usage: m4 descriptors.m4 > sll_descriptors.F90
!> @ingroup descriptors
!> @brief Describes different global flags throughout the library
!> @details
!> The intent of this module is to provide a single, library-wide definition
!> of the names used to describe different properties which are applicabale
!> in a global sense
!> <b> How to use-it </b>


define(`foreach', `pushdef(`$1')_foreach($@)popdef(`$1')')
define(`_arg1', `$1')
define(`_foreach', `ifelse(`$2', `()', `',
  `define(`$1', _arg1$2)$3`'$0(`$1', (shift$2), `$3')')')
divert`'dnl
dnl
dnl
dnl
dnl
dnl
dnl####### Type definition ###########
define(TYPE_DEF, `type sll_$1 
    sll_int32                  :: id
    character(len=32) private :: pname
  contains
    procedure,  pass(self)      :: name=>name_$1 
    procedure,  pass(self)      :: parse=>parse_$1
  end type sll_$1
  
  interface operator(.eq.)
    module procedure sll_$1_compare
  end interface' )
dnl  
dnl  
dnl##### member functions ##########
dnl
dnl ##### define counter ##############
define(`count', `_$0`'define(`_$0',incr(_$0))')dnl
dnl #################
define(DESCRIPTOR_PARAM_DEFINITIONS,` 
define(`_count', `1')dnl
foreach(`descriptorparamarg', $2, 
`type(sll_$1), parameter :: descriptorparamarg = sll_$1(count,"descriptorparamarg") 
')')dnl
dnl
dnl
define(MEMBER_FUNCTIONS, `
 !----------------------------------------------------------------------
 !------------------- $1 -----------------------------------------
 pure function name_$1( self ) result( r )
   class( sll_$1 ), intent( in ) :: self
   character(len=len(self%pname))     :: r
   r = self%pname
 end function name_$1  
  
 pure function sll_$1_compare(bc1,bc2) result(compare)
  type(sll_$1), intent(in) :: bc1, bc2
  logical :: compare
   compare = bc1%id .eq. bc2%id
 end function sll_$1_compare

 subroutine parse_$1(self, str)
    class( sll_$1 ), intent( inout ) :: self
    character(len=*), intent(in) :: str
    character(len=len(str)) :: strc
      
     !Remove blanks left and right
     strc=adjustl(trim(str))
     
    if foreach(`typedescriptor', $2, 
    `(strc .eq. typedescriptor%name()) then  
        self%id=typedescriptor%id
        self%pname=typedescriptor%pname
      elseif')dnl
  (.true.) then
    endif
    
 end subroutine parse_$1
 !----------------------------------------------------------------------
      ')dnl      
dnl ##### DESCRIPTOR HEADER ########
define(DESCRIPTOR_HEADER,`TYPE_DEF($1)
DESCRIPTOR_PARAM_DEFINITIONS($1,$2)')
dnl
dnl ##### DESCRIPTOR BODY ########
define(DESCRIPTOR_BODY,`MEMBER_FUNCTIONS($1,$2)')
dnl
dnl
dnl
dnl      
dnl 
dnl
dnl###########################################################
dnl###### ASSEMBLE FORTRAN SOURCE ############################
dnl###########################################################
dnl
dnl  
dnl  
dnl ############## HEADER #########################  
module sll_descriptors
#include "sll_working_precision.h"
dnl
dnl
dnl
DESCRIPTOR_HEADER(vlasovpoisson_sim,(SLL_LANDAU_DIAG,SLL_LANDAU_SUM,SLL_LANDAU_PROD,SLL_TWOSTREAM,SLL_BUMPONTAIL))
DESCRIPTOR_HEADER(boundary,(OPEN,CLOSED))

contains

DESCRIPTOR_BODY(vlasovpoisson_sim,(SLL_LANDAU_DIAG,SLL_LANDAU_SUM,SLL_LANDAU_PROD,SLL_TWOSTREAM,SLL_BUMPONTAIL))
DESCRIPTOR_BODY(boundary,(OPEN,CLOSED))

dnl
end module sll_descriptors
dnl
dnl ############## SOURCE #########################
