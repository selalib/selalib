#if ! defined (USE_R3_INFO) && (defined (USE_R3_INFO_MIX) || \
defined (USE_R3_INFO_MPI) || defined (USE_R3_INFO_SMA) || \
defined (USE_R3_INFO_OMP) )
#define USE_R3_INFO
#endif
      
SUBROUTINE r3_info_init ()
  
  ! Written : r. richter, sgi, oct 2006
  use read_write_module  
  IMPLICIT NONE
  
#include "r3_info.h"
  
  ! Functions
  INTEGER (4) :: omp_get_thread_num, omp_get_num_threads
  LOGICAL (4) :: omp_in_parallel
  
  
#if defined (USE_R3_INFO_MIX)
  IF (omp_in_parallel()) &
    STOP 'CALL r3_info_init () in parrallel region is not allowed'
  
!$OMP PARALLEL  
  r3_mytask = omp_get_thread_num()
  r3_ntasks = omp_get_num_threads()
!$OMP END PARALLEL
  
  CALL r3_info_init_x ()
  
#elif defined (USE_R3_INFO_OMP)
  IF (omp_in_parallel()) &
    STOP 'CALL r3_info_init () in parrallel region is not allowed'
  
!$OMP PARALLEL
  CALL r3_info_init_x ()
!$OMP END PARALLEL
  
#else
  CALL r3_info_init_x ()
#endif
  
  RETURN
END SUBROUTINE r3_info_init
      
SUBROUTINE r3_info_init_x ()
  
  ! Written : r. richter, sgi, oct 2006
  
  IMPLICIT NONE
  
#include <mpif.h>
#include "r3_info.h"
  
  ! Local parameters
  
  INTEGER (4)        :: info
  CHARACTER (len=80) :: chvalue
  
  ! Functions
  
  INTEGER (4)        :: my_pe, num_pes, omp_get_thread_num, &
    omp_get_num_threads
  
  
  ! Initialyze communication parameters
  
  r3_comm = -1
  
#if defined (USE_R3_INFO_MIX) || defined (USE_R3_INFO_MPI)
  r3_comm = mpi_comm_world
  CALL mpi_comm_rank (r3_comm, r3_mype, info)
  CALL mpi_comm_size (r3_comm, r3_npes, info)
#elif defined (USE_R3_INFO_SMA)
  r3_mype = my_pe()
  r3_npes = num_pes()
#elif defined (USE_R3_INFO_OMP)
  r3_mype = omp_get_thread_num()
  r3_npes = omp_get_num_threads()
#else
  r3_mype = 0
  r3_npes = 1
#endif
  
  r3_root = 0
  
  ! Initialyze r3_info parameters
  
#if defined (USE_R3_INFO)
  r3_stak_n = 0
  
  r3_info_stak(:)      = 1
  r3_info_count(:,:)   = 0
  r3_info_start(:,:,:) = 0.0
  r3_info_sncds(:,:,:) = 0.0
  r3_info_name(:)      = ' '
  
  histx_start = -999
  histx_stop  = -999
  
  CALL getenv ('R3_HISTX_START', chvalue)
  IF (chvalue(1:2) .NE. '  ') READ (chvalue, *) histx_start
  
  CALL getenv ('R3_HISTX_STOP', chvalue)
  IF (chvalue(1:2) .NE. '  ') READ (chvalue, *) histx_stop
  
  r3_info_n        = 1
  r3_info_name(-6) = 'MKL-Libraries'
  r3_info_name(-5) = 'I/O-Input/Ouput'
  r3_info_name(-4) = 'BAR-Barrier'
  r3_info_name(-3) = 'COM-Communication'
  r3_info_name(-2) = 'OMP-Parallel'
  r3_info_name(-1) = 'CPU-Computation'
  r3_info_name(0)  = 'Other'
  r3_info_name(1)  = 'Total'
  
  r3_part_n          = 1
  r3_part_name(:,-4) = (/ '               ', ' TEXT          ' /)
  r3_part_name(:,-3) = (/ '               ', ' HEAP          ' /)
  r3_part_name(:,-2) = (/ '               ', ' STACK         ' /)
  r3_part_name(:,-1) = (/ '               ', ' SYM           ' /)
  r3_part_name(:,0)  = (/ '               ', ' START         ' /)
  r3_part_name(:,1)  = (/ 'Initialization ', ' INIT          ' /)
  r3_part_name(:,2)  = (/ 'First          ', ' FIRST         ' /)
  r3_part_name(:,3)  = (/ 'Solver         ', ' SOLV          ' /)
  r3_part_name(:,4)  = (/ 'Finish         ', 'FINISH         ' /)
  
  CALL r3_barrier (-2_4, ' ')
  CALL r3_info_print (-1_4, r3_info_index_0, 'START')
#endif
  
  RETURN
END SUBROUTINE r3_info_init_x
      
!
      
SUBROUTINE r3_info_print (mode, ind, text)
  
  ! Written : r. richter, sgi, oct 2006
  use globals, only : plocal_id  
  IMPLICIT NONE
  
  ! Input/output parameters
  
  INTEGER (4)        :: mode, ind
  CHARACTER (len=*)  :: text
  
  ! Functions
  
  LOGICAL (4)        :: omp_in_parallel
   
#if defined (USE_R3_INFO_OMP)
  IF (omp_in_parallel()) THEN
    CALL r3_info_print_x (mode, ind, text)
  ELSE
    
!$OMP PARALLEL  
    CALL r3_info_print_x (mode, ind, text)
!$OMP END PARALLEL
    
  END IF
#else
  CALL r3_info_print_x (mode, ind, text)
#endif
  
  RETURN
END SUBROUTINE r3_info_print
      
!
      
SUBROUTINE r3_info_print_x (mode, ind, text)
  
  ! Written : r. richter, sgi, oct 2006
  use globals, only : plocal_id, iter_run, nbiter
  use read_write_module    
  IMPLICIT NONE
      
#include "r3_info.h"
      
  ! Input/output parameters
      
  INTEGER (4)        :: mode, ind
  CHARACTER (len=*)  :: text
      
  ! Local parameters
      
  INTEGER (4)        :: i, j, jnd, ipid, ierror
  REAL (8)           :: sdiff
  CHARACTER (len=80) :: txt, name
  INTEGER            :: i_bb(8) ! Need to be kind independent
      
  ! Save parameter
      
  LOGICAL (4), SAVE  :: lfirst1 = .TRUE., lfirst2 = .TRUE.
  INTEGER (8), SAVE  :: flops(6) = 0
  REAL (8), SAVE     :: stamp(5)
      
#if defined (USE_R3_INFO_OMP)
!$OMP THREADPRIVATE (lfirst1, lfirst2, flops, stamp)
#endif
      
  ! Functions
      
  INTEGER (4)       :: getpid
  REAL (8)          :: r3_secondr
  character(LEN=20) :: file_name_r3
      
#if defined (USE_R3_INFO)
  CALL date_and_time (values = i_bb)
      
!#CP  
  if (plocal_id.eq.0) then
     !-> opening of the file if necessary
     if (.not.file_r3_open) then
        write(file_name_r3,'(A)')"gysela_r3.out"//char(0)
        open(uout_r3, file = file_name_r3, status = 'UNKNOWN', &
             position = 'APPEND', form = 'FORMATTED')
        file_r3_open = .true.
     end if
  end if
!#CP
      
  IF (lfirst1) THEN
    IF (r3_mype .EQ. r3_root) &
      WRITE (uout_r3, '(/,a,3(i2.2,a),i3.3,5a)') 'Info (', i_bb(5), ':', &
      i_bb(6), ':', i_bb(7), '.', i_bb(8), ') : ', &
      '                      OTHER     ', &
      'WALL 1     WALL 2     WALL 3    '
    lfirst1 = .FALSE.
  END IF
      
  IF (mode .GT. 1) THEN
    WRITE (txt, '(a,3(i2.2,a),i3.3,2a,i5)') 'Info (', i_bb(5), ':', &
      i_bb(6), ':', i_bb(7), '.', i_bb(8), ') : ', text, mode
  ELSE
    WRITE (txt, '(a,3(i2.2,a),i3.3,2a)') 'Info (', i_bb(5), ':', &
      i_bb(6), ':', i_bb(7), '.', i_bb(8), ') : ', text
  END IF
      
  IF (mode .NE. 0) THEN
    name = 'BAR_'
    DO i = 1, INDEX (text, ' ') - 1
      j = ICHAR (text(i:i))
      IF (j .GE. 65 .AND. j .LE. 90) THEN
        name(4+i:4+i) = CHAR (j + 32)
      ELSE
        name(4+i:4+i) = CHAR (j)
      END IF
    END DO
    IF (mode .NE. -1) CALL r3_barrier (ind, name)
  END IF
  
  stamp(5) = r3_secondr()
      
  jnd = r3_info_stak(1)
  sdiff = stamp(5) - r3_info_start(1,jnd,r3_part_n)
  r3_info_sncds(1,jnd,r3_part_n) = r3_info_sncds(1,jnd,r3_part_n) + sdiff
  r3_info_start(1,jnd,r3_part_n) = stamp(5)
      
  IF (mode .EQ. -1_4) THEN
    sdiff      = 0.0
    stamp(1:4) = stamp(5) - 1.0e-9
    flops(1:4) = flops(5)
    CALL r3_info_begin (r3_info_index_1, r3_part_name(1,1))
  ELSE IF (mode .LE. -2_4 .AND. mode .GE. -8_4) THEN
    CALL r3_info_end (MAX (r3_info_index_1, r3_info_index_2, &
      r3_info_index_3, r3_info_index_4, &
      r3_info_index_5, r3_info_index_6, &
      r3_info_index_7, r3_info_index_8))
    r3_info_sncds(:,1,r3_part_n) = stamp(5) - &
      r3_info_start(:,1,r3_part_n)
    r3_part_n = r3_part_n + 1
    r3_info_start(:,1,r3_part_n) = stamp(5)
    IF (r3_info_index_2 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_2, r3_part_name(1,2))
    ELSE IF (r3_info_index_3 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_3, r3_part_name(1,3))
    ELSE IF (r3_info_index_4 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_4, r3_part_name(1,4))
    ELSE IF (r3_info_index_5 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_5, r3_part_name(1,5))
    ELSE IF (r3_info_index_6 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_6, r3_part_name(1,6))
    ELSE IF (r3_info_index_7 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_7, r3_part_name(1,7))
    ELSE IF (r3_info_index_8 .EQ. -1_4) THEN
      CALL r3_info_begin (r3_info_index_8, r3_part_name(1,8))
    END IF
  ELSE IF (mode .EQ. -9_4) THEN
    CALL r3_info_end (MAX (r3_info_index_1, r3_info_index_2, &
      r3_info_index_3, r3_info_index_4, &
      r3_info_index_5, r3_info_index_6, &
      r3_info_index_7, r3_info_index_8))
    r3_info_sncds(:,1,r3_part_n) = r3_info_sncds(:,1,r3_part_n) + &
      stamp(5) - r3_info_start(:,1,r3_part_n)
  END IF
      
  IF (r3_mype .EQ. r3_root) &
    WRITE (uout_r3, '(a,2f11.4,2f11.3)') txt(1:38), sdiff, &
    stamp(5) - stamp(4), stamp(5) - stamp(3), &
    stamp(5) - stamp(1)
      
  stamp(4) = stamp(5)
  flops(4) = flops(5)
  IF (mode .EQ. -2) THEN
    stamp(3) = stamp(5)
    flops(3) = flops(5)
  END IF
      
  IF (mode .EQ. histx_start .OR. mode .EQ. histx_stop) THEN
    ipid = getpid()
#if defined (SGI_MIPS)
    CALL pxfkill (ipid, 12, ierror)
#endif
    IF (mode .EQ. histx_start .AND. r3_mype .EQ. r3_root) &
      WRITE (uout_r3, '(a,a)') txt(1:38), ' HISTX started'
    IF (mode .EQ. histx_stop .AND. r3_mype .EQ. r3_root) &
      WRITE (uout_r3, '(a,a)') txt(1:38), ' HISTX stopped'
  END IF
#endif
      
  !-> file closing
  if (iter_run.eq.nbiter) then
     close(uout_r3)
     file_r3_open = .false.
  end if
      
  RETURN
END SUBROUTINE r3_info_print_x
      
!
      
SUBROUTINE r3_info_begin (ind, name)
      
  ! Written : r. richter, sgi, oct 2006
  use read_write_module  
  IMPLICIT NONE
  
#include "r3_info.h"
  
  ! Input/output parameters
  
  INTEGER (4)        :: ind
  CHARACTER (len=*)  :: name
  
  ! Local parameters
  
  INTEGER (4)        :: i
  REAL (8)           :: stamp
  
  INTEGER (8), SAVE  :: flops = 0
  
#if defined (USE_R3_INFO_OMP)
!$OMP THREADPRIVATE (flops)
#endif
  
  ! Functions
      
  REAL (8)           :: r3_secondr      
      
#if defined (USE_R3_INFO)
      
#if defined (USE_R3_INFO_MIX)
  IF (r3_mytask .EQ. 0) THEN
#endif
      
    IF (ind .EQ. -1) THEN
      r3_info_n = r3_info_n + 1
      IF (r3_info_n .GT. r3_info_mx) THEN
        IF (r3_mype .EQ. r3_root) THEN
          DO i = 1, r3_info_n
            WRITE (uout_r3,*) i, TRIM (r3_info_name(i))
          END DO
          WRITE (uout_r3,*)
          WRITE (uout_r3,*) r3_info_mx + 1, TRIM (r3_info_name(i))
        END IF
        STOP 'Increase r3_info_mx'
      END IF
      ind = r3_info_n
      r3_info_name(ind) = TRIM (name)
    END IF
    
    r3_stak_n = r3_stak_n + 1
    IF (r3_stak_n .GT. r3_stak_mx) THEN
      IF (r3_mype .EQ. r3_root) THEN
        DO i = 1, r3_stak_mx
          WRITE (uout_r3,*) i, TRIM (r3_info_name(r3_info_stak(i)))
        END DO
        WRITE (uout_r3,*)
        WRITE (uout_r3,*) r3_stak_mx + 1, TRIM (r3_info_name(ind))
      END IF
      STOP 'Increase r3_stak_mx'
    END IF
    r3_info_stak(r3_stak_n) = ind
    
    r3_info_count(ind,r3_part_n) = r3_info_count(ind,r3_part_n) + 1
    
    stamp = r3_secondr()
    
    r3_info_start(:,ind,r3_part_n) = stamp
    
#if defined (USE_R3_INFO_MIX)
  END IF
#endif
      
#endif
      
  RETURN
END SUBROUTINE r3_info_begin
      
!
      
SUBROUTINE r3_info_end (ind)
      
  ! Written : r. richter, sgi, oct 2006
  use read_write_module    
  IMPLICIT NONE
      
#include "r3_info.h"
      
  ! Input/output parameters
      
  INTEGER (4)        :: ind
      
  ! Local parameters
      
  INTEGER (4)        :: i, jnd
  REAL (8)           :: stamp, sdiff(2)
  
  INTEGER (8), SAVE  :: flops = 0
      
#if defined (USE_R3_INFO_OMP)
!$OMP THREADPRIVATE (flops)
#endif
      
  ! Functions
      
  REAL (8)           :: r3_secondr
      
#if defined (USE_R3_INFO)
      
#if defined (USE_R3_INFO_MIX)
  IF (r3_mytask .EQ. 0) THEN
#endif
      
    stamp = r3_secondr()
    
    sdiff(:) = stamp - r3_info_start(:,ind,r3_part_n)
    r3_info_sncds(:,ind,r3_part_n) = &
      r3_info_sncds(:,ind,r3_part_n) + sdiff(:)
    
    IF (r3_info_stak(r3_stak_n) .NE. ind) &
      WRITE (uout_r3,*) 'Info ERROR : ', &
      TRIM (r3_info_name(r3_info_stak(r3_stak_n))), &
      ' Not equal to ', &
      TRIM (r3_info_name(r3_info_stak(ind)))
    
    r3_stak_n = r3_stak_n - 1
    DO i = 1, r3_stak_n
      jnd = r3_info_stak(i)
      r3_info_start(1,jnd,r3_part_n) = &
        r3_info_start(1,jnd,r3_part_n) + sdiff(1)
    END DO
    
#if defined (USE_R3_INFO_MIX)
  END IF
#endif
  
#endif
      
  RETURN
END SUBROUTINE r3_info_end
      
!
      
SUBROUTINE r3_info_summary ()
      
  ! Written : r. richter, sgi, oct 2006
      
  IMPLICIT NONE
      
#include "r3_info.h"
      
  ! Local parameters
      
  INTEGER (4)        :: omp_num(0:r3_npes-1), &
    omp_chars(80,-6:r3_info_mx,0:r3_npes-1)
  INTEGER (8)        :: omp_flops(2,-6:r3_info_mx,0:r3_npes-1), &
    omp_cnts(-6:r3_info_mx,0:r3_npes-1)
  REAL (8)           :: omp_sncds(2,-6:r3_info_mx,0:r3_npes-1)
  
  ! Functions
      
  LOGICAL (4)        :: omp_in_parallel
      
#if defined (USE_R3_INFO_OMP)
  IF (omp_in_parallel()) &
    STOP 'R3_INFO_SUMMARY in parrallel region is not allowed'
      
!$OMP PARALLEL DEFAULT (none) &
!$OMP  SHARED (omp_num, omp_chars, omp_flops, omp_cnts, omp_sncds)
  
  CALL r3_info_summary_x (omp_num, omp_chars, omp_flops, omp_cnts, &
    omp_sncds)
!$OMP END PARALLEL
      
#else
  CALL r3_info_summary_x (omp_num, omp_chars, omp_flops, omp_cnts, &
    omp_sncds)
#endif
      
  RETURN
END SUBROUTINE r3_info_summary
      
!
      
SUBROUTINE r3_info_summary_x (omp_num, omp_chars, omp_flops, &
  omp_cnts, omp_sncds)
      
  ! Written : r. richter, sgi, oct 2006
  use read_write_module    
  IMPLICIT NONE
      
#include <mpif.h>
#include "r3_info.h"
      
  ! Input/output parameters
  
  INTEGER (4)        :: omp_num(0:r3_npes-1), &
    omp_chars(80,-6:r3_info_mx,0:r3_npes-1)
  INTEGER (8)        :: omp_flops(2,-6:r3_info_mx,0:r3_npes-1), &
    omp_cnts(-6:r3_info_mx,0:r3_npes-1)
  REAL (8)           :: omp_sncds(2,-6:r3_info_mx,0:r3_npes-1)
  
  ! Local scalar parameters
  
  INTEGER (4)        :: i, j, k, j1, j2, ipart, num, nsub, info
  REAL (8)           :: dom_cum
  CHARACTER (len=16) :: txt, ttt
  LOGICAL (4)        :: found
  
  ! Local arrays parameters
  
  INTEGER (4)        :: ind(-6:r3_info_mx), request(4), &
    status(mpi_status_size,4)
  INTEGER (8)        :: flp_tot(2,-6:r3_info_mx), &
    flops(2,-6:r3_info_mx), &
    cnt_max(-6:r3_info_mx), cnts(-6:r3_info_mx)
  REAL (8)           :: cpu_min(-6:r3_info_mx), &
    cpu_max(-6:r3_info_mx), &
    cpu_tot(2,-6:r3_info_mx), &
    sncds(2,-6:r3_info_mx)
  CHARACTER (len=80) :: names(-6:r3_info_mx), &
    sub_name(-6:r3_info_mx)
      
#if defined (USE_R3_INFO)
  DO i = r3_stak_n, 2, -1
    CALL r3_info_end (r3_info_stak(i))
  END DO
  
  CALL r3_info_print (-9_4, -2, 'END')
  
  DO i = 2, r3_info_n
    
    IF (r3_info_name(i)(:3) .EQ. 'CPU') THEN
      j = -1
    ELSE IF (r3_info_name(i)(:3) .EQ. 'OMP') THEN
      j = -2
    ELSE IF (r3_info_name(i)(:3) .EQ. 'COM') THEN
      j = -3
    ELSE IF (r3_info_name(i)(:3) .EQ. 'BAR') THEN
      j = -4
    ELSE IF (r3_info_name(i)(:3) .EQ. 'I/O') THEN
      j = -5
    ELSE IF (r3_info_name(i)(:3) .EQ. 'MKL') THEN
      j = -6
    ELSE
      j = 0
    END IF
    
    DO k = 1, r3_part_n
      r3_info_sncds(1,j,k) = r3_info_sncds(1,j,k) + &
        r3_info_sncds(1,i,k) 
      r3_info_sncds(2,j,k) = r3_info_sncds(2,j,k) + &
        r3_info_sncds(2,i,k)
      r3_info_count(j,k)   = r3_info_count(j,k) + &
        r3_info_count(i,k)
    END DO
    
  END DO
      
  r3_part_n = r3_part_n + 1
  r3_part_name(:,r3_part_n) = (/ 'Global         ', &
    'GLOBAL         ' /)
  DO i = -6, r3_info_n
    r3_info_sncds(1,i,r3_part_n) = &
      SUM (r3_info_sncds(1,i,1:r3_part_n-1))
    r3_info_sncds(2,i,r3_part_n) = &
      SUM (r3_info_sncds(2,i,1:r3_part_n-1))
    r3_info_count(i,r3_part_n)   = &
      SUM (r3_info_count(i,1:r3_part_n-1))
  END DO
  
  CALL r3_barrier (-2_4, ' ')
  
  DO ipart = 1, r3_part_n
    
#if defined (USE_R3_INFO_MIX) || defined (USE_R3_INFO_MPI)
    CALL mpi_isend (r3_info_n, 1_4, INT (mpi_integer, 4), r3_root, &
      11_4, r3_comm, request(1), info)
    CALL mpi_isend (r3_info_sncds(1,-6,ipart), &
      INT (2 * (r3_info_n + 7), 4), INT (mpi_real8, &
      4), r3_root, 12_4, r3_comm, request(2), info)
    CALL mpi_isend (r3_info_count(-6,ipart), INT (r3_info_n + 7, &
      4), INT (mpi_integer8, 4), r3_root, 14_4, &
      r3_comm, request(3), info)
    CALL mpi_isend (r3_info_name(-6), INT (80 * (r3_info_n + 7), &
      4), INT (mpi_byte, 4), r3_root, 15_4, r3_comm, &
      request(4), info)
#elif defined (USE_R3_INFO_OMP)
    omp_num(r3_mype)       = r3_info_n
    omp_sncds(:,:,r3_mype) = r3_info_sncds(:,:,ipart)
    omp_cnts(:,r3_mype)    = r3_info_count(:,ipart)
    DO k = -6, r3_info_n
      DO j = 1, 80
        omp_chars(j,k,r3_mype) = ICHAR (r3_info_name(k)(j:j))
      END DO
    END DO
#endif
      
    IF (r3_mype .EQ. r3_root) THEN
      
      cpu_min(:)   = 1.0d+9
      cpu_max(:)   = -1.0d+9
      cpu_tot(:,:) = 0.0
      flp_tot(:,:) = 0
      cnt_max(:)   = 0
      
      nsub = r3_info_n
      sub_name(:) = r3_info_name(:)
      
      DO i = 0, r3_npes - 1
        
#if defined (USE_R3_INFO_MIX) || defined (USE_R3_INFO_MPI)
        CALL mpi_recv (num, 1_4, INT (mpi_integer, 4), &
          INT (r3_root + i, 4), 11_4, r3_comm, &
          status, info)
        CALL mpi_recv (sncds, INT (2 * (num + 7), 4), &
          INT (mpi_real8, 4), INT (r3_root + i, 4), &
          12_4, r3_comm, status, info)
        CALL mpi_recv (cnts, INT (num + 7, 4), INT (mpi_integer8, &
          4), INT (r3_root + i, 4), 14_4, r3_comm, &
          status, info)
        CALL mpi_recv (names, INT (80 * (num + 7), 4), &
          INT (mpi_byte, 4), INT (r3_root + i, 4), &
          15_4, r3_comm, status, info)
#elif defined (USE_R3_INFO_SMA)
        CALL shmem_get4 (num, r3_info_n, 1_4, INT (r3_root + i, &
          4))
        CALL shmem_get8 (sncds, r3_info_sncds(1,-6,ipart), &
          INT (2 * (num + 7), 4), &
          INT (r3_root + i, 4))
        CALL shmem_get8 (cnts, r3_info_count(-6,ipart), &
          INT (num + 7, 4), INT (r3_root + i, 4))
        CALL shmem_character_get (names, r3_info_name, &
          INT (80 * (num + 7), 4), &
          INT (r3_root + i, 4))
#elif defined (USE_R3_INFO_OMP)
        num = omp_num(i)
        sncds(:,:) = omp_sncds(:,:,i)
        flops(:,:) = omp_flops(:,:,i)
        cnts(:)    = omp_cnts(:,i)
        DO k = -6, num
          DO j = 1, 80
            names(k)(j:j) = CHAR (omp_chars(j,k,i))
          END DO
        END DO
#else
        num = r3_info_n
        sncds(:,:) = r3_info_sncds(:,:,ipart)
        cnts(:)    = r3_info_count(:,ipart)
        names(:)   = r3_info_name(:)
#endif
        ind(:)     = -100
        DO j = -6, num
          IF (names(j) .EQ. sub_name(j)) ind(j) = j
          IF (ind(j) .EQ. -100) THEN
            DO k = 2, nsub
              IF (names(j) .EQ. sub_name(k)) ind(j) = k
            END DO
          END IF
          IF (ind(j) .EQ. -100) THEN
            nsub = nsub + 1
            ind(j) = nsub
            sub_name(nsub) = names(j)
          END IF
        END DO
        
        DO j = -6, num
          k = ind(j)
          cpu_tot(:,k) = cpu_tot(:,k) + sncds(:,j)
          flp_tot(:,k) = flp_tot(:,k) + flops(:,j)
          cpu_min(k)   = MIN (cpu_min(k), sncds(1,j))
          IF (cnts(j) .EQ. 0) cpu_min(k) = 0.0d0
          cpu_max(k) = MAX (cpu_max(k), sncds(1,j))
          cnt_max(k) = MAX (cnt_max(k), cnts(j))
        END DO
        
        DO j = 2, nsub
          found = .FALSE.
          IF (j .LE. num .AND. sub_name(j) .EQ. names(j)) &
            found = .TRUE.
          IF (.NOT. found) THEN
            DO k = 2, num
              IF (sub_name(j) .EQ. names(k)) found = .TRUE.
            END DO
          END IF
          IF (.NOT. found) cpu_min(j) = 0.0d0
        END DO
        
      END DO
      
      DO j = -6, 0
        cpu_tot(:,j) = cpu_tot(1,j) / r3_npes
      END DO
      
      DO j = 1, nsub
        cpu_tot(:,j) = cpu_tot(:,j) / r3_npes
      END DO
      
      DO j = -6, nsub
        ind(j) = j
      END DO
      
      DO j = -6, 0
        DO i = -6, 0
          IF (ABS (cpu_tot(1,ind(i))) .LT. &
            ABS (cpu_tot(1,ind(j)))) THEN
            k      = ind(j)
            ind(j) = ind(i)
            ind(i) = k
          END IF
        END DO
      END DO
      
      DO j = 2, nsub
        DO i = 2, nsub
          IF (ABS (cpu_tot(1,ind(i))) .LT. &
            ABS (cpu_tot(1,ind(j)))) THEN
            k      = ind(j)
            ind(j) = ind(i)
            ind(i) = k
          END IF
        END DO
      END DO
      
      ttt = 'part : '//r3_part_name(2,ipart)(:6)
      txt = 'Info ('//ttt(8:13)//') : '
      
      WRITE (uout_r3, '(a,/,a,/,3a,/,3a,/,a,/,4a,/,a)') txt, txt, txt, &
        '                            ', '  Timing summary '// &
        ttt(1:13), txt, '                            ', &
        '  ****************************', txt, txt, &
        'NAME                      ', &
        'MIN        MAX        MEAN       INCL     ', &
        'MEAN%    CUMUL    INCL%    COUNT', txt
      
      WRITE (uout_r3, '(a,a20,3f11.3,11x,/a)') txt, sub_name(1), &
        cpu_tot(1,1), cpu_tot(1,1), cpu_tot(1,1), txt
      
      DO k = 1, 2
        IF (k .EQ. 1) THEN
          j1 = -6
          j2 = 0
        ELSE
          j1 = 2
          j2 = nsub
        END IF
        dom_cum = 0.0
        DO j = j1, j2
          i       = ind(j)
          dom_cum = dom_cum + ABS (cpu_tot(1,i))
          IF (cpu_tot(1,i) .GT. 1.0d-5 .AND. &
            cnt_max(i) .GT. 0) THEN
            IF (cpu_tot(1,i) .EQ. cpu_tot(2,i)) THEN
              WRITE (uout_r3, '(a,a20,3f11.3,11x,2f9.2,9x,i9)') txt, &
                sub_name(i), cpu_min(i), cpu_max(i), &
                cpu_tot(1,i), 100 * cpu_tot(1,i) / &
                cpu_tot(1,1), 100 * dom_cum / &
                cpu_tot(1,1), cnt_max(i)
            ELSE
              WRITE (uout_r3, '(a,a20,4f11.3,3f9.2,i9)') txt, &
                sub_name(i), cpu_min(i), cpu_max(i), &
                cpu_tot(1,i), cpu_tot(2,i), &
                100 * cpu_tot(1,i) / cpu_tot(1,1), &
                100 * dom_cum / cpu_tot(1,1), &
                100 * cpu_tot(2,i) / cpu_tot(1,1), &
                cnt_max(i)
            END IF
          END IF
        END DO
        WRITE (uout_r3, '(a)') txt
      END DO
    END IF
        
#if defined (USE_R3_INFO_MIX) || defined (USE_R3_INFO_MPI)
    CALL mpi_waitall (4_4, request, status, info)
#else
    CALL r3_barrier (-2_4, ' ')
#endif
        
  END DO
#endif
      
  RETURN
END SUBROUTINE r3_info_summary_x
    
! 
      
SUBROUTINE r3_barrier (ind, name)
      
  ! Written : r. richter, sgi, oct 2006
  
  IMPLICIT NONE
      
#include <mpif.h>
#include "r3_info.h"
      
  ! Input/output parameters
  
  INTEGER (4)        :: ind
  CHARACTER (len=*)  :: name
  
  ! Local parameters
      
  INTEGER (4)        :: info
      
  ! Statements
  
  IF (ind .GT. -2) CALL r3_info_begin (ind, name)
  
  IF (ind .GT. -3) THEN
    
#if defined (USE_R3_INFO_MIX) || defined (USE_R3_INFO_MPI)
    CALL mpi_barrier (r3_comm, info)
#elif defined (USE_R3_INFO_SMA)
    CALL shmem_barrier_all ()
#elif defined (USE_R3_INFO_OMP)
!$OMP BARRIER
#endif
      
  END IF
  
  IF (ind .GT. -2) CALL r3_info_end (ind)
  
  RETURN
END SUBROUTINE r3_barrier
