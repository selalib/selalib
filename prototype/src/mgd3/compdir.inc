# define WMGD 0
# define double_precision 1
# if double_precision
# define REALN real*8
      implicit double precision(a-h,o-z)
# else
# define REALN real*4
      implicit real*4(a-h,o-z)
# endif
# define NBLOCKGR 1
# define cdebug  0
# if cdebug
# if NBLOCKGR
      integer nisend,nirecv,reduce,nallreduce,nalltoall,nwaitall
      common/comsgr/nisend(3,3),nirecv(3,3),nreduce,nallreduce,
     1              nalltoall,nwaitall
      integer nisendb,nirecvb,reduceb,nallreduceb,nalltoallb,nwaitallb
      common/comsgrb/nisendb(3,3),nirecvb(3,3),nreduceb,nallreduceb,
     1               nalltoallb,nwaitallb
      double precision tisend,tirecv,treduce,tallreduce,talltoall,
     1                 twaitall
      common/mpicomsgr/tisend(3,3),tirecv(3,3),treduce,tallreduce,
     1                 talltoall,twaitall
      double precision tisendb,tirecvb,treduceb,tallreduceb,talltoallb,
     1                 twaitallb
      common/mpicomsgrb/tisendb(3,3),tirecvb(3,3),treduceb,tallreduceb,
     1                  talltoallb,twaitallb
# else
      integer nsendrecv,nreduce,nallreduce,nalltoall
      common/comsgr/nsendrecv(3,3),nreduce,nallreduce,nalltoall
      integer nsendrecvb,nreduceb,nallreduceb,nalltoallb
      common/comsgrb/nsendrecvb(3,3),nreduceb,nallreduceb,nalltoallb
      double precision tsendrecv,treduce,tallreduce,talltoall
      common/mpicomsgr/tsendrecv(3,3),treduce,tallreduce,talltoall
      double precision tsendrecvb,treduceb,tallreduceb,talltoallb
      common/mpicomsgrb/tsendrecvb(3,3),treduceb,tallreduceb,talltoallb
# endif
# if NBLOCKFR
      integer nisendfr,nirecvfr,nwaitallfr
      common/comsfr/nisendfr,nirecvfr,nwaitallfr
      integer nisendfrb,nirecvfrb,nwaitallfrb
      common/comsfrb/nisendfrb,nirecvfrb,nwaitallfrb
      double precision tisendfr,tirecvfr,twaitallfr
      common/mpicomsfr/tisendfr,tirecvfr,twaitallfr
      double precision tisendfrb,tirecvfrb,twaitallfrb
      common/mpicomsfrb/tisendfrb,tirecvfrb,twaitallfrb
# else
      integer nsendfr,nrecvfr
      common/comsfr/nsendfr,nrecvfr
      integer nsendfrb,nrecvfrb
      common/comsfrb/nsendfrb,nrecvfrb
      double precision tsendfr,trecvfr
      common/mpicomsfr/tsendfr,trecvfr
      double precision tsendfrb,trecvfrb
      common/mpicomsfrb/tsendfrb,trecvfrb
# endif
      double precision timing
      integer nsteptiming
      logical nocterr
      common/mpitiming/timing(150),nsteptiming,nocterr
      double precision timingb
      integer nsteptimingb
      logical nocterrb
      common/mpitiming/timingb(150),nsteptimingb,nocterrb
# endif
