program test_sprng
implicit none

#include "sprng_f.h"

integer       :: streamnum, nstreams, seed,junk
SPRNG_POINTER :: stream
integer       :: gtype

!print *,'Type in a generator type (integers: 0,1,2,3,4,5):  '
!read *, gtype

do gtype = 0, 5

  streamnum = 0
  nstreams = 1
  seed = 985456376
  
  stream = init_sprng(gtype,streamnum,nstreams,seed,SPRNG_DEFAULT)
  junk = print_sprng(stream)
  !print *, 'Printing information about new stream'
  
  call sub1(stream)
  
  junk = free_sprng(stream)

end do

end

subroutine sub1(stream)
      
#include "sprng_f.h"
SPRNG_POINTER :: stream
real(8) :: rn

print *, 'Printing 3 double precision numbers in [0,1): '
do i = 1, 3
  rn = sprng(stream)
  write(*, "(i6,2x,f19.16)") i, rn
end do

end subroutine sub1
