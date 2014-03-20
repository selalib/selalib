# define WMGD 0
# define REALN real(8)
      implicit none

# define NBLOCKGR 1
# define cdebug  0
# if cdebug
# if NBLOCKGR
      integer nisend,nirecv,nreduce,nallreduce,nalltoall,nwait,
     1        nwaitall
      common/comsgr/nisend(2,3),nirecv(2,3),nreduce,nallreduce,
     1              nalltoall,nwait,nwaitall
# else
      integer nsendrecv,nreduce,nallreduce,nalltoall
      common/comsgr/nsendrecv(2,3),nreduce,nallreduce,nalltoall
# endif
      integer nisendfr,nirecvfr,nwaitallfr
      common/comsfr/nisendfr,nirecvfr,nwaitallfr
      double precision timing
      integer nsteptiming
      logical nocterr
      common/mpitiming/timing(100),nsteptiming,nocterr
# endif

