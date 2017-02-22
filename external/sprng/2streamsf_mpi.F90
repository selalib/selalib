!
!              Demonstrates use of shared and non-shared streams    
! Each process has two streams.  One stream is common to all the    
! processes. The other stream is different on each processor. 
!
! Uncomment the following line to get the interface with pointer checking
! #define CHECK_POINTERS

program twostreamsf_mpi
use, intrinsic :: ISO_C_BINDING
use mpi
implicit none

#include "sprng_f.h"

integer       :: streamnum,commNum, nstreams, seed
SPRNG_POINTER :: stream, commonStream
real(8)       :: rn
integer       :: i
integer       :: myid, nprocs, ierror
integer       :: junk
integer       :: gtype       

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)

streamnum = myid       !This stream is different on each proces
commNum   = nprocs     !This stream is common to all processes
nstreams  = nprocs+1   !extra stream is common to all processes
seed      = 985456376

!--- node 0 is reading in a generator type
if (myid .eq. 0) then
  print *, 'Available generators; use corresponding numeral:'
  print *, '   lfg     --- 0 '
  print *, '   lcg     --- 1 '
  print *, '   lcg64   --- 2 '
  print *, '   cmrg    --- 3 '
  print *, '   mlfg    --- 4 '
  print *, '   pmlcg   --- 5 '
  print *,'Type in a generator type (integers: 0,1,2,3,4,5):  '
!PN set gtype hardly
!  read *, gtype
  print*,'gtype=', gtype
   gtype = 0
endif

call MPI_BCAST(gtype,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

!  This stream is different on each process
stream = init_sprng(gtype,streamnum,nstreams,seed,SPRNG_DEFAULT)
write(6, 44) myid
44 format("Process", i2, ": Print information about new stream")
junk = print_sprng(stream)

!  This stream is identical on each process
commonStream = init_sprng(gtype,commNum,nstreams,seed,SPRNG_DEFAULT)
write (6, 55) myid
55 format ("Process", i2,": This stream is identical on all processes")
junk = print_sprng(commonStream)

do i = 1, 2
  rn = sprng(stream)
  write(6, 66) myid, i, rn
end do

do i = 1, 2
  rn = sprng(commonStream)
  write(6, 77) myid, i, rn
end do

66 format("Process", i2,", random number (distinct stream)",i2,": ", f8.6)
77 format("Process", i2, ", random number (shared stream)",i2,": ", f8.6)
junk = free_sprng(stream)
junk = free_sprng(commonStream)

call MPI_FINALIZE(ierror)

end
