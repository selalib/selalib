      subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
# include "compdir.inc"
      integer n, numprocs, myid, s, e
      integer nlocal
      integer deficit
c------------------------------------------------------------------------
c  From the MPE library
c  This file contains a routine for producing a decomposition of a 1-d 
c  array when given a number of processors.  It may be used in "direct" 
c  product decomposition.  The values returned assume a "global" domain 
c  in [1:n]
c
c Code      : tmgd2
c Called in : main
c Calls     : --
c------------------------------------------------------------------------
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
