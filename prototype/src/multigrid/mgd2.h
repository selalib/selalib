implicit none

# define WMGD 0
# define double_precision 1

# define NBLOCKGR 1

# if DEBUG
# if NBLOCKGR
sll_int32  :: nisend,nirecv,nreduce,nallreduce,nalltoall,nwait,nwaitall
common/comsgr/nisend(2,3),nirecv(2,3),nreduce,nallreduce, &
              nalltoall,nwait,nwaitall
# else
sll_int32  :: nsendrecv,nreduce,nallreduce,nalltoall
common/comsgr/nsendrecv(2,3),nreduce,nallreduce,nalltoall
# endif
sll_int32  :: nisendfr,nirecvfr,nwaitallfr
common/comsfr/nisendfr,nirecvfr,nwaitallfr
sll_real64 :: timing
sll_int32  :: nsteptiming
logical    :: nocterr
common/mpitiming/timing(100),nsteptiming,nocterr
# endif
