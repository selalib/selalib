subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)

implicit none 

integer :: n, numprocs, myid, s, e
integer :: nlocal
integer :: deficit
!------------------------------------------------------------------------
!  From the MPE library
!  This file contains a routine for producing a decomposition of a 1-d 
!  array when given a number of processors.  It may be used in "direct" 
!  product decomposition.  The values returned assume a "global" domain 
!  in [1:n]
!
! Code      : tmgd3
! Called in : main
! Calls     : --
!------------------------------------------------------------------------
nlocal  = n / numprocs
s	      = myid * nlocal + 1
deficit = mod(n,numprocs)
s	      = s + min(myid,deficit)
if (myid .lt. deficit) then
    nlocal = nlocal + 1
endif
e = s + nlocal - 1
if (e .gt. n .or. myid .eq. numprocs-1) e = n
return
end
