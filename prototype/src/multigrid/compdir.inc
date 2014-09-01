      implicit none

#if cdebug

#if NBLOCKGR

      integer nisend,nirecv,nreduce,nallreduce,nalltoall,nwait
      integer nwaitall
      common/comsgr/nisend(2,3),nirecv(2,3),nreduce,nallreduce,
     &              nalltoall,nwait,nwaitall

#else

      integer nsendrecv,nreduce,nallreduce,nalltoall
      common/comsgr/nsendrecv(2,3),nreduce,nallreduce,nalltoall

#endif

      integer nisendfr,nirecvfr,nwaitallfr
      common/comsfr/nisendfr,nirecvfr,nwaitallfr
      double precision timing
      integer nsteptiming
      logical nocterr
      common/mpitiming/timing(100),nsteptiming,nocterr

#endif
