module sll_m_qsort_partition
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  implicit none

  public :: &
    sll_s_qsortc

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

recursive subroutine sll_s_qsortc(A)

  integer, intent(inout), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
    call Partition(A, iq)
    call sll_s_qsortc(A(:iq-1))
    call sll_s_qsortc(A(iq:))
  endif

end subroutine sll_s_qsortc

subroutine partition(A, marker)

  integer, intent(inout), dimension(:) :: A
  integer, intent(out)                 :: marker
  integer                              :: i
  integer                              :: j
  integer                              :: temp
  integer                              :: x

  x = A(1)
  i = 0
  j = size(A) + 1

  do
    j = j-1
    do
      if (A(j) <= x) exit
      j = j-1
    end do
    i = i+1
    do
      if (A(i) >= x) exit
      i = i+1
    end do
    if (i < j) then
      ! exchange A(i) and A(j)
      temp = A(i)
      A(i) = A(j)
      A(j) = temp
    elseif (i == j) then
      marker = i+1
      return
    else
      marker = i
      return
    endif
  end do

end subroutine partition

end module sll_m_qsort_partition
