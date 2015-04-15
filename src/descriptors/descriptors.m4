
define(DESCRIPTOR_NAME, vlasovpoisson_sim )dnl 
define(DESCRIPTOR_TAGS, SLL_LANDAU_DIAG SLL_LANDAU_DIAG SLL_LANDAU_DIAG SLL_LANDAU_DIAG  )dnl 



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
  

  
  
  
dnl##### member functions ##########

define(DESCRIPTOR_PARAM_DEFINITIONS, 
foreach(`descriptorparamarg', $2 , 
`type(sll_vlasovpoisson_sim), parameter :: descriptorparamarg = sll_$1(1,"descriptorparamarg") ')dnl
)


define(MEMBER_FUNCTIONS, `$1

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
       elseif (str .eq. %name()) then
       self%pid= %pid
       self%pname=%pname
     endif
      ')

  

dnl###########################################################
dnl###### ASSEMBLE FORTRAN SOURCE ############################
dnl###########################################################

  
  
dnl ############## HEADER #########################  
TYPE_DEF(vlasovpoisson_sim)
DESCRIPTOR_PARAM_DEFINITIONS(vlasovpoisson_sim, (SLL_LANDAU_DIAG, SLL_LANDAU_DIAG, SLL_LANDAU_DIAG) )
 

dnl MEMBER_FUNCTIONS(vlasovpoisson_sim   )





contains



dnl ############## SOURCE #########################



