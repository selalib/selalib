
#if ! defined (USE_R3_INFO) && (defined (USE_R3_INFO_MIX) || \
    defined (USE_R3_INFO_MPI) || defined (USE_R3_INFO_SMA) || \
    defined (USE_R3_INFO_OMP) )
#define USE_R3_INFO
#endif

! This headder file defines all r3_info parameters and common blocks
! warning need to be -i8/-i4 compatible.

! Written : r. richter, sgi, oct 2006

      INTEGER (4), PARAMETER :: r3_info_mx = 256
      INTEGER (4), PARAMETER :: r3_part_mx = 10
      INTEGER (4), PARAMETER :: r3_stak_mx = 64

! ...

      INTEGER (4)        :: r3_mype, r3_npes, r3_root, r3_comm
      COMMON /r3_info_a/    r3_mype, r3_npes, r3_root, r3_comm

#if defined (USE_R3_INFO_OMP)
!$OMP THREADPRIVATE (/r3_info_a/)
#endif

#if defined (USE_R3_INFO_MIX)
      INTEGER (4)        :: r3_mytask, r3_ntasks
      COMMON /r3_info_b/    r3_mytask, r3_ntasks

!$OMP THREADPRIVATE (/r3_info_b/)
#endif

#if defined (USE_R3_INFO)
      INTEGER (4)        :: r3_info_n, r3_part_n, r3_stak_n
      COMMON /r3_info_1/    r3_info_n, r3_part_n, r3_stak_n

      INTEGER (4)        :: r3_info_stak(r3_stak_mx)
      COMMON /r3_info_2/    r3_info_stak

      INTEGER (8)        :: r3_info_count(-6:r3_info_mx,r3_part_mx)
      COMMON /r3_info_3/    r3_info_count

      REAL (8)           :: r3_info_start(2,-6:r3_info_mx,r3_part_mx)
      COMMON /r3_info_4/    r3_info_start

      REAL (8)           :: r3_info_sncds(2,-6:r3_info_mx,r3_part_mx)
      COMMON /r3_info_5/    r3_info_sncds

      CHARACTER (len=80) :: r3_info_name(-6:r3_info_mx)
      CHARACTER (len=80) :: r3_part_name(2,-4:r3_part_mx)
      COMMON /r3_info_6/    r3_info_name, r3_part_name

      INTEGER (4)        :: histx_start, histx_stop
      COMMON /r3_histx/     histx_start, histx_stop

#if defined (USE_R3_INFO_OMP)
!$OMP THREADPRIVATE (/r3_info_1/, /r3_info_2/, /r3_info_3/, /r3_info_4/)
!$OMP THREADPRIVATE (/r3_info_5/, /r3_info_6/, /r3_histx/ )
#endif

#endif

      INTEGER (4), SAVE  :: r3_info_index_0    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_1    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_2    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_3    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_4    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_5    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_6    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_7    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_8    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_9    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_a    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_b    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_c    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_d    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_e    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_f    = -1_4
      INTEGER (4), SAVE  :: r3_info_index_x(2) = (/ -1_4, -1_4 /)

#if defined (USE_R3_INFO_OMP)
!$OMP THREADPRIVATE (r3_info_index_0, r3_info_index_1, r3_info_index_2)
!$OMP THREADPRIVATE (r3_info_index_3, r3_info_index_4, r3_info_index_5)
!$OMP THREADPRIVATE (r3_info_index_6, r3_info_index_7, r3_info_index_8)
!$OMP THREADPRIVATE (r3_info_index_9, r3_info_index_a, r3_info_index_b)
!$OMP THREADPRIVATE (r3_info_index_c, r3_info_index_d, r3_info_index_e)
!$OMP THREADPRIVATE (r3_info_index_f, r3_info_index_x)
#endif

