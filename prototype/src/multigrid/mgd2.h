
# define WMGD 0
# define REALN sll_real64
# define NBLOCKGR 1
# define cdebug  0
# define IOUT 6
# if cdebug
# if NBLOCKGR
sll_int32  :: nisend,nirecv,nreduce,nallreduce,nalltoall,nwait,nwaitall
common/comsgr/nisend(2,3),nirecv(2,3),nreduce,nallreduce, &
              nalltoall,nwait,nwaitall
# else
sll_int32  :: nsendrecv,nreduce,nallreduce,nalltoall
common/comsgr/nsendrecv(2,3),nreduce,nallreduce,nalltoall
# endif
sll_int32  :: nisendfr,nirecvfr,nwaitallfr
sll_real64 :: timing
sll_int32  :: nsteptiming
logical    :: nocterr
common/comsfr/nisendfr,nirecvfr,nwaitallfr
common/mpitiming/timing(100),nsteptiming,nocterr
# endif
