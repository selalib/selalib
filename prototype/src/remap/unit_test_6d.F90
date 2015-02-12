program remap_test_6d
  use sll_collective, only: sll_boot_collective, sll_halt_collective
  use sll_remapper
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none

#define RANK_TO_PRINT 0

  ! Test of the 6D remapper takes a 6D array whose global size is
  ! N1*N2*N3*N4*N5*N6 and is distributed among NPi*NPj*NPk*NPl*NPm*NPn
  ! processors.
  sll_real64, dimension(:,:,:,:,:,:), allocatable :: local_array1
  sll_real64, dimension(:,:,:,:,:,:), allocatable :: local_array2
!!$  sll_int32, dimension(:,:,:,:,:,:), allocatable :: local_array1
!!$  sll_int32, dimension(:,:,:,:,:,:), allocatable :: local_array2
  sll_real64, dimension(:,:,:,:,:,:), allocatable :: arrays_diff
!!$  sll_int32, dimension(:,:,:,:,:,:), allocatable :: arrays_diff
  ! Dimensions of the 6d array (global dimensions)
  integer, parameter                       :: ni = 32
  integer, parameter                       :: nj = 16
  integer, parameter                       :: nk = 16   ! change
  integer, parameter                       :: nl = 16   ! change
  integer, parameter                       :: nm = 16   ! change
  integer, parameter                       :: nn = 16   ! change
!!$  integer, parameter                       :: ni = 4
!!$  integer, parameter                       :: nj = 4
!!$  integer, parameter                       :: nk = 4   ! change
!!$  integer, parameter                       :: nl = 4   ! change
!!$  integer, parameter                       :: nm = 4   ! change
!!$  integer, parameter                       :: nn = 4   ! change
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_k_init
  integer                                   :: loc_sz_l_init
  integer                                   :: loc_sz_m_init
  integer                                   :: loc_sz_n_init

  integer                                   :: loc_sz_i_final
  integer                                   :: loc_sz_j_final
  integer                                   :: loc_sz_k_final
  integer                                   :: loc_sz_l_final
  integer                                   :: loc_sz_m_final
  integer                                   :: loc_sz_n_final
  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  integer                                   :: npl
  integer                                   :: npm
  integer                                   :: npn
  sll_int32                                 :: gi, gj, gk, gl, gm, gn
  integer                                   :: ierr
  integer                                   :: myrank
  sll_int64                                 :: colsz        ! collective size
  ! Remap stuff
  type(layout_6D), pointer                  :: layout1
  type(layout_6D), pointer                  :: layout2
  type(remap_plan_6D_real64), pointer       :: rmp6

  sll_real64                                :: rand_real
  integer, parameter                        :: nbtest = 25
  integer                                   :: i_test
  integer                                   :: i, j, k, l, m, n
  sll_int32, dimension(6)                   :: global_indices, g
  sll_real32   , dimension(1)               :: prod4test
  integer                                   :: ok
  sll_int32                                 :: lin_index
  sll_int32, dimension(6)                   :: theo_index
  sll_int32, dimension(6)                   :: tmp_array

  ! Boot parallel environment
  print *, 'Booting parallel environment...'
  call sll_boot_collective()

  !  end_result = .true.
  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  if( myrank .eq. RANK_TO_PRINT ) then
     print *, ' '
     print *, '--------------- REMAP 6D test ---------------------'
     print *, ' '
     print *, 'Running a test on ', colsz, 'processes'
     call flush(6)
  end if

  if (.not. is_power_of_two(colsz)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if

  ok = 1
  do, i_test=1, nbtest
     call flush(6)
     if( myrank .eq. 0 ) then
        print *, 'Iteration ', i_test, ' of ', nbtest
     end if
     layout1  => new_layout_6D( sll_world_collective )        
     call factorize_in_random_2powers_6d(colsz, npi, npj, npk, npl, npm, npn)
!!$     npi = 2
!!$     npj = 4
!!$     npk = 1
!!$     npl = 1
     if( myrank .eq. RANK_TO_PRINT ) then
        print *, 'source configuration: ', npi, npj, npk, npl, npm, npn
     end if

     call initialize_layout_with_distributed_6D_array( &
          ni, &
          nj, &
          nk, &
          nl, &
          nm, &
          nn, &
          npi, &
          npj, &
          npk, &
          npl, &
          npm, &
          npn, &
          layout1 )

!!$     if( myrank == RANK_TO_PRINT ) then
!!$        print *, '----------------------------------------'
!!$        print *, 'process: ', myrank, 'viewing layout1: '
!!$        call sll_view_lims_6D(layout1)
!!$        print *, '----------------------------------------'
!!$        call flush(6)
!!$     end if

     call compute_local_sizes_6d( &
          layout1, &
          loc_sz_i_init, &
          loc_sz_j_init, &
          loc_sz_k_init, &
          loc_sz_l_init, &
          loc_sz_m_init, &
          loc_sz_n_init )      
  
!!$     if( myrank == RANK_TO_PRINT) then
!!$        print *, 'local sizes: ', &
!!$             loc_sz_i_init, &
!!$             loc_sz_j_init, &
!!$             loc_sz_k_init, &
!!$             loc_sz_l_init, &
!!$             loc_sz_m_init, &
!!$             loc_sz_n_init 
!!$     end if

     SLL_ALLOCATE(local_array1(loc_sz_i_init, loc_sz_j_init, loc_sz_k_init, loc_sz_l_init, loc_sz_m_init, loc_sz_n_init), ierr)

     ! initialize the local data  
     do n=1, loc_sz_n_init
        do m=1, loc_sz_m_init  
           do l=1,loc_sz_l_init
              do k=1,loc_sz_k_init
                 do j=1,loc_sz_j_init 
                    do i=1,loc_sz_i_init
                       tmp_array(:) = (/i, j, k, l, m, n/)
                       global_indices = local_to_global_6D(layout1, tmp_array)
                       gi = global_indices(1)
                       gj = global_indices(2)
                       gk = global_indices(3)
                       gl = global_indices(4)
                       gm = global_indices(5)
                       gn = global_indices(6)
                       lin_index = gi + (gj-1)*ni + (gk-1)*ni*nj + &
                            (gl-1)*ni*nj*nk + (gm-1)*ni*nj*nk*nl + &
                            (gn-1)*ni*nj*nk*nl*nm
                       local_array1(i,j,k,l,m,n) = lin_index
                       ! Check if the theoretical index is what we expect:
                       theo_index = theoretical_global_6D_indices( &
                            lin_index, ni, nj, nk, nl, nm )

                       if( theo_index(1) .ne. gi ) then
                          print *, 'discrepancy found: gi = ', gi, &
                               'theoretical = ', theo_index(1)
                          print *, 'Rank: ', myrank
                          print *, '(i,j,k,l,m,n) = ', i,j,k,l,m,n
                          print *, '(gi,gj,gk,gl,gm,gn) = ',gi,gj,gk,gl,gm,gn
                          print *, 'linear index = ', lin_index
                          print *, 'theoretical = ', theo_index(:)
                          stop
                       end if

                       if( theo_index(2) .ne. gj ) then
                          print *, 'discrepancy found: gj = ', gj, &
                               'theoretical = ', theo_index(2)
                          print *, 'Rank: ', myrank
                          print *, '(i,j,k,l,m,n) = ', i,j,k,l,m,n
                          print *, '(gi,gj,gk,gl,gm,gn) = ',gi,gj,gk,gl,gm,gn
                          print *, 'linear index = ', lin_index
                          print *, 'theoretical = ', theo_index(:)
                          stop
                       end if

                       if( theo_index(3) .ne. gk ) then
                          print *, 'discrepancy found: gk = ', gk, &
                               'theoretical = ', theo_index(3)
                          print *, 'Rank: ', myrank
                          print *, '(i,j,k,l,m,n) = ', i,j,k,l,m,n
                          print *, '(gi,gj,gk,gl,gm,gn) = ',gi,gj,gk,gl,gm,gn
                          print *, 'linear index = ', lin_index
                          print *, 'theoretical = ', theo_index(:)
                          stop
                       end if

                       if( theo_index(4) .ne. gl ) then
                          print *, 'discrepancy found: gl = ', gl, &
                               'theoretical = ', theo_index(4)
                          print *, 'Rank: ', myrank
                          print *, '(i,j,k,l,m,n) = ', i,j,k,l,m,n
                          print *, '(gi,gj,gk,gl,gm,gn) = ',gi,gj,gk,gl,gm,gn
                          print *, 'linear index = ', lin_index
                          print *, 'theoretical = ', theo_index(:)
                          stop
                       end if

                       if( theo_index(5) .ne. gm ) then
                          print *, 'discrepancy found: gm = ', gm, &
                               'theoretical = ', theo_index(5)
                          print *, 'Rank: ', myrank
                          print *, '(i,j,k,l,m,n) = ', i,j,k,l,m,n
                          print *, '(gi,gj,gk,gl,gm,gn) = ',gi,gj,gk,gl,gm,gn
                          print *, 'linear index = ', lin_index
                          print *, 'theoretical = ', theo_index(:)
                          stop
                       end if

                       if( theo_index(6) .ne. gn ) then
                          print *, 'discrepancy found: gn = ', gn, &
                               'theoretical = ', theo_index(6)
                          print *, 'Rank: ', myrank
                          print *, '(i,j,k,l,m,n) = ', i,j,k,l,m,n
                          print *, '(gi,gj,gk,gl,gm,gn) = ',gi,gj,gk,gl,gm,gn
                          print *, 'linear index = ', lin_index
                          print *, 'theoretical = ', theo_index(:)
                          stop
                       end if

                    end do
                 end do
              end do
           end do
        end do
     end do
     
!!$     if( myrank == RANK_TO_PRINT ) then
!!$        print *, '*********************************************************'
!!$        print *, 'rank = ', myrank, 'initialized array:'
!!$        print *, local_array1(:,:,:,:,:,:)
!!$        call flush(6)
!!$        print *, '*********************************************************'
!!$     end if

     layout2  => new_layout_6D( sll_world_collective )
     call factorize_in_random_2powers_6d(colsz, npi, npj, npk, npl, npm, npn)
!!$     npi = 4
!!$     npj = 2
!!$     npk = 1
!!$     npl = 1

     if( myrank .eq. RANK_TO_PRINT ) then
        print *, 'target configuration: ', npi, npj, npk, npl, npm, npn
     end if

     call initialize_layout_with_distributed_6D_array( &
          ni, &
          nj, &
          nk, &
          nl, &
          nm, &
          nn, &
          npi, &
          npj, &
          npk, &
          npl, &
          npm, &
          npn, &
          layout2 )

     call reorganize_randomly_6d(layout2)
     
     call compute_local_sizes_6d( &
          layout2, &
          loc_sz_i_final, &
          loc_sz_j_final, &
          loc_sz_k_final, &
          loc_sz_l_final, &
          loc_sz_m_final, &
          loc_sz_n_final )

     SLL_ALLOCATE(local_array2(loc_sz_i_final, loc_sz_j_final, loc_sz_k_final, loc_sz_l_final, loc_sz_m_final, loc_sz_n_final), ierr )
    
     rmp6 => new_remap_plan( layout1, layout2, local_array1 )     

!!$     if( myrank .eq. RANK_TO_PRINT ) then
!!$        print *, myrank, 'just before remap, local_array1:'
!!$        print *, local_array1(:,:,:,:,:,:)
!!$        print *, myrank, 'just before remap, local_array2:'
!!$        print *, local_array2(:,:,:,:,:,:)
!!$     end if

     call apply_remap_6D( rmp6, local_array1, local_array2 ) 

!!$     if( myrank .eq. RANK_TO_PRINT ) then
!!$        print *, myrank, 'just after remap'
!!$        print *, local_array2(:,:,:,:,:,:)
!!$     end if

!#if 0
     SLL_ALLOCATE(arrays_diff(loc_sz_i_final, loc_sz_j_final, loc_sz_k_final, loc_sz_l_final, loc_sz_m_final, loc_sz_n_final),ierr )
 
#if 0
     if( myrank .eq. 0 ) then
        print *, i_test, myrank, 'Printing layout1: '
        call sll_view_lims_6D( layout1 )
        print *, i_test, myrank, 'Printing layout2: '
        call sll_view_lims_6D( layout2 )
     end if
#endif
     ! compare results with expected data
     do n=1, loc_sz_n_final
        do m=1, loc_sz_m_final
           do l=1, loc_sz_l_final
              do k=1,loc_sz_k_final
                 do j=1,loc_sz_j_final 
                    do i=1,loc_sz_i_final
                       tmp_array(:) = (/i,j,k,l,m,n/)
                       global_indices = local_to_global_6D( layout2, tmp_array )
                       gi = global_indices(1)
                       gj = global_indices(2)
                       gk = global_indices(3)
                       gl = global_indices(4)
                       gm = global_indices(5)
                       gn = global_indices(6)
                       arrays_diff(i,j,k,l,m,n) = local_array2(i,j,k,l,m,n) - &
                            (gi + (gj-1)*ni + (gk-1)*ni*nj + (gl-1)*ni*nj*nk + &
                            (gm-1)*ni*nj*nk*nl + (gn-1)*ni*nj*nk*nl*nm)
                       if (arrays_diff(i,j,k,l,m,n) /=0 ) then
                          ok = 0
                          print*, i_test, myrank, '"remap" unit test: FAIL'
                          print *, i_test, myrank, 'local indices: ', &
                               '(', i, j, k, l, m, n, ')'
                          print *, 'global indices wrt target layout'
                          print*, 'rank: ', myrank, &
                               'in global indices: (',gi,',', gj,',', gk,',', &
                               gl,',', gm, ',', gn, ')'
                          print *, 'rank: ', myrank, &
                               'local array1(', gi,',', gj, ',', gk, ',', gl, &
                               ',', gm, ',', gn, ')=', &
                               local_array1(i,j,k,l,m,n)  
                          print*, myrank, &
                               'local array2(',gi, ',', gj, ',', gk, ',', gl, &
                               ',', gm, ',', gn, ')=', &
                               local_array2(i,j,k,l,m,n)
                          g = theoretical_global_6D_indices(& 
                               int(local_array2(i,j,k,l,m,n),i32),&
                               ni, &
                               nj, &
                               nk, &
                               nl, &
                               nm ) 
                          print*, 'Theoretical indices: (',g(1), ',', g(2), &
                               ',',g(3),',',g(4),',', g(5), ',', g(6), ')'
                          !     if(myrank .eq. 1) then
                          print *, 'local_array2 = '
                          print *, local_array2(:,:,:,:,:,:)
                          !    end if
                    
                          print *, i_test, myrank, 'Printing layout1: '
                          call sll_view_lims_6D( layout1 )
                          print *, i_test, myrank, 'Printing layout2: '
                          call sll_view_lims_6D( layout2 )
                    
                          print*, 'program stopped by failure'
                          stop
                       end if
                       call flush(6)
                    end do
                 end do
              end do
           end do
        end do
     end do
     ! As the difference described above is null in each point, the sum of the 
     ! corresponding absolute values must be null. Each processor computes a 
     ! local 
     ! sum and all local sums are finally added and the result is sent to 
     ! processor 0 which will check if equal 0 to validate the test. (*)
     call sll_collective_reduce( &
          sll_world_collective, &
          (/ real(ok) /), &
          1, &
          MPI_PROD, &
          0, &
          prod4test )

     if( myrank .eq. 0) then
        print *, ' '
        print *, '-------------------------------------------'
        print *, ' '
        call flush(6)
     end if
     call flush(6) 
       
     call sll_collective_barrier(sll_world_collective)
  
     call sll_delete( layout1 )
     call sll_delete( layout2 )
     call sll_delete( rmp6 )
     SLL_DEALLOCATE_ARRAY(local_array1, ierr)
     SLL_DEALLOCATE_ARRAY(local_array2, ierr)
     SLL_DEALLOCATE_ARRAY(arrays_diff, ierr)
!#endif
  enddo

  if (myrank==0) then
     if (prod4test(1)==1.) then
        print*, 'TEST PASSED'
     endif
  endif



  call sll_halt_collective()
  
contains

  subroutine reorganize_randomly_6d(layout)
    implicit none
    type(layout_6D), pointer   :: layout
    integer                    :: i, colsz, proc_n, proc_p
    real                       :: rand_real
    colsz = sll_get_num_nodes(layout)
    do i=0, colsz-1
       call random_number(rand_real)
       proc_n = int(rand_real*(colsz-1))
       call random_number(rand_real)
       proc_p = int(rand_real*(colsz-1))
       call swap_box_6D(proc_n, proc_p, layout)
    enddo    
  end subroutine reorganize_randomly_6d

  subroutine swap_box_6D(proc_n, proc_p, layout)
    implicit none
    integer                  :: proc_n, proc_p
    type(layout_6D), pointer :: layout
    integer                  :: i_min_n, i_max_n
    integer                  :: j_min_n, j_max_n
    integer                  :: k_min_n, k_max_n 
    integer                  :: l_min_n, l_max_n
    integer                  :: m_min_n, m_max_n
    integer                  :: n_min_n, n_max_n

    integer                  :: i_min_p, i_max_p 
    integer                  :: j_min_p, j_max_p 
    integer                  :: k_min_p, k_max_p
    integer                  :: l_min_p, l_max_p    
    integer                  :: m_min_p, m_max_p    
    integer                  :: n_min_p, n_max_p    

    ! Get proc_n contents from layout
    i_min_n = get_layout_6D_i_min( layout, proc_n )
    i_max_n = get_layout_6D_i_max( layout, proc_n )
    j_min_n = get_layout_6D_j_min( layout, proc_n )
    j_max_n = get_layout_6D_j_max( layout, proc_n )
    k_min_n = get_layout_6D_k_min( layout, proc_n )
    k_max_n = get_layout_6D_k_max( layout, proc_n )
    l_min_n = get_layout_6D_l_min( layout, proc_n )
    l_max_n = get_layout_6D_l_max( layout, proc_n )
    m_min_n = get_layout_6D_m_min( layout, proc_n )
    m_max_n = get_layout_6D_m_max( layout, proc_n )
    n_min_n = get_layout_6D_n_min( layout, proc_n )
    n_max_n = get_layout_6D_n_max( layout, proc_n )

    ! Get proc_p contents from layout
    i_min_p = get_layout_6D_i_min( layout, proc_p )
    i_max_p = get_layout_6D_i_max( layout, proc_p )
    j_min_p = get_layout_6D_j_min( layout, proc_p )
    j_max_p = get_layout_6D_j_max( layout, proc_p )
    k_min_p = get_layout_6D_k_min( layout, proc_p )
    k_max_p = get_layout_6D_k_max( layout, proc_p )
    l_min_p = get_layout_6D_l_min( layout, proc_p )
    l_max_p = get_layout_6D_l_max( layout, proc_p )
    m_min_p = get_layout_6D_m_min( layout, proc_p )
    m_max_p = get_layout_6D_m_max( layout, proc_p )
    n_min_p = get_layout_6D_n_min( layout, proc_p )
    n_max_p = get_layout_6D_n_max( layout, proc_p )

    ! Set proc_n contents in layout
    call set_layout_6D_i_min( layout, proc_n, i_min_p )
    call set_layout_6D_i_max( layout, proc_n, i_max_p)
    call set_layout_6D_j_min( layout, proc_n, j_min_p )
    call set_layout_6D_j_max( layout, proc_n, j_max_p )
    call set_layout_6D_k_min( layout, proc_n, k_min_p )
    call set_layout_6D_k_max( layout, proc_n, k_max_p )
    call set_layout_6D_l_min( layout, proc_n, l_min_p )
    call set_layout_6D_l_max( layout, proc_n, l_max_p )
    call set_layout_6D_m_min( layout, proc_n, m_min_p )
    call set_layout_6D_m_max( layout, proc_n, m_max_p )
    call set_layout_6D_n_min( layout, proc_n, n_min_p )
    call set_layout_6D_n_max( layout, proc_n, n_max_p )

    ! Set proc_p contents in layout
    call set_layout_6D_i_min( layout, proc_p, i_min_n )
    call set_layout_6D_i_max( layout, proc_p, i_max_n )
    call set_layout_6D_j_min( layout, proc_p, j_min_n )
    call set_layout_6D_j_max( layout, proc_p, j_max_n )
    call set_layout_6D_k_min( layout, proc_p, k_min_n )
    call set_layout_6D_k_max( layout, proc_p, k_max_n )   
    call set_layout_6D_l_min( layout, proc_p, l_min_n )
    call set_layout_6D_l_max( layout, proc_p, l_max_n )
    call set_layout_6D_m_min( layout, proc_p, m_min_n )
    call set_layout_6D_m_max( layout, proc_p, m_max_n )
    call set_layout_6D_n_min( layout, proc_p, n_min_n )
    call set_layout_6D_n_max( layout, proc_p, n_max_n )
  end subroutine swap_box_6D

  function theoretical_global_6D_indices(d, ni, nj, nk, nl, nm)
    integer, dimension(1:6) :: theoretical_global_6D_indices
    integer, intent(in)     :: d, ni, nj, nk, nl, nm ! global array dims
    integer                 :: i, j, k, l, m, n
    integer                 :: q
    sll_real64              :: tmp
    ! We know that the linear index relates with the other
    ! indices by:
    !
    !          linear = i + ni*(j-1) + ni*nj*(k-1) + ni*nj*nk*(l-1) +
    !
    !                   ni*nj*nk*nl*(m-1) + ni*nj*nk*nl*nm*(n-1) 
    !
    ! where i,j,k,l,m,n are the indices sought. We start by working backwards.
    ! Note: the use of double precision is critical here, else the rounding
    ! will be wrong given the sizes of the numbers involved.
    tmp = real(d,f64)/real(ni*nj*nk*nl*nm,f64)
    n   = ceiling(tmp)
    ! reduce the size of the number we are working with:
    q   =  d - (n-1)*ni*nj*nk*nl*nm
    tmp = real(q,f64)/real(ni*nj*nk*nl,f64)
    m   = ceiling(tmp)
    ! reduce again
    q   = q - (m-1)*ni*nj*nk*nl
    tmp = real(q,f64)/real(ni*nj*nk,f64)
    l   = ceiling(tmp)
    ! reduce again
    q   = q - (l-1)*ni*nj*nk
    tmp = real(q,f64)/real(ni*nj,f64)
    k   = ceiling(tmp)
    ! reduce again
    q   = q - (k-1)*ni*nj
    tmp = real(q,f64)/real(ni,f64)
    j   = ceiling(tmp)
    ! final reduction
    q   = q - (j-1)*ni
    i   = q
    theoretical_global_6D_indices(1:6) = (/i,j,k,l,m,n/) 
  end function theoretical_global_6D_indices
  
  subroutine factorize_in_random_2powers_6d(n, n1, n2, n3, n4, n5, n6)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2, n3, n4, n5, n6
    integer   :: expo, expo1, expo2, expo3, expo4, expo5, expo6
    if (.not.is_power_of_two(n)) then   
       print*, 'factorize_in_random_2powers_6d(): ', &
            'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = int(rand_real*(expo-expo1))
    call random_number(rand_real)
    expo3 = int(rand_real*(expo-(expo1+expo2)))
    call random_number(rand_real)
    expo4 = int(rand_real*(expo-(expo1+expo2+expo3)))
    call random_number(rand_real)
    expo5 = int(rand_real*(expo-(expo1+expo2+expo3+expo4)))
    expo6 = expo - (expo1 + expo2 + expo3 + expo4 + expo5)
    n1 = 2**expo1
    n2 = 2**expo2
    n3 = 2**expo3
    n4 = 2**expo4
    n5 = 2**expo5
    n6 = 2**expo6
  end subroutine factorize_in_random_2powers_6d

      
end program remap_test_6d
    
