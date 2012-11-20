! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123;  Default seeds x,y,z,w.
!  Set your own seeds with statement i=kisset(ix,iy,iz,iw).
module marsaglia
implicit none
private
public :: kiss, kisset
   integer :: x=123456789, y=362436069, z=521288629, w=916191069

contains

   FUNCTION kiss ()
      integer :: kiss
      x = 69069 * x + 1327217885
      y = m (m (m (y, 13), - 17), 5)
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss = x + y + ishft (z, 16) + w
   contains
      function m(k, n)
         integer :: m, k, n
         m = ieor (k, ishft (k, n) )
      end function m
   END FUNCTION kiss
   function kisset (ix, iy, iz, iw)
      integer :: kisset, ix, iy, iz, iw
      x = ix
      y = iy
      z = iz
      w = iw
      kisset = 1
   end function kisset

end module marsaglia

!     PROGRAM test
!     use marsaglia
!     PRINT *, kiss ()
!     PRINT *, kisset (1, 2, 3, 4)
!     PRINT *, kiss ()
!     END PROGRAM test
