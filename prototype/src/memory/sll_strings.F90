module sll_strings
implicit none

interface operator (.conc.)
   module procedure sll_concatenate_strings
end interface

contains

function sll_concatenate_strings(s1, s2)
  character(len=*), intent(in) :: s1
  character(len=*), intent(in) :: s2
  character(len=(len(s1)+len(s2))) :: sll_concatenate_strings
  sll_concatenate_strings = s1//s2
end function sll_concatenate_strings



end module sll_strings

