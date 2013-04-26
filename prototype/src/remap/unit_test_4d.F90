program remap_test_4d
  use sll_collective, only: sll_boot_collective, sll_halt_collective
  use sll_remapper
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none

  ! Test of the 4D remapper takes a 4D array whose global size N1*N2*N3*N4,
  ! distributed among NPi*NPj*NPk*NPl processors.
  sll_real64, dimension(:,:,:,:), allocatable :: local_array1
  sll_real64, dimension(:,:,:,:), allocatable :: local_array2
  sll_real64, dimension(:,:,:,:), allocatable :: arrays_diff
  ! Dimensions of the 4d array (global dimensions)
  integer, parameter                       :: ni = 64
  integer, parameter                       :: nj = 32
  integer, parameter                       :: nk = 16   ! change
  integer, parameter                       :: nl = 16   ! change
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_k_init
  integer                                   :: loc_sz_l_init
  integer                                   :: loc_sz_i_final
  integer                                   :: loc_sz_j_final
  integer                                   :: loc_sz_k_final
  integer                                   :: loc_sz_l_final
  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  integer                                   :: npl
  sll_int32                                 :: gi, gj, gk, gl
  integer                                   :: ierr
  integer                                   :: myrank
  sll_int64                                 :: colsz        ! collective size
  ! Remap stuff
  type(layout_4D), pointer                  :: layout1
  type(layout_4D), pointer                  :: layout2
  type(remap_plan_4D_real64), pointer       :: rmp4

  sll_real64                                :: rand_real
  integer, parameter                        :: nbtest = 25
  integer                                   :: i_test
  integer                                   :: i, j, k, l
  sll_int32, dimension(4)                   :: global_indices, g
  sll_real32   , dimension(1)               :: prod4test
  integer                                   :: ok
  sll_int32, dimension(4)                   :: tmpa


  ! Boot parallel environment
  call sll_boot_collective()
  !  end_result = .true.
  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  if( myrank .eq. 0) then
     print *, ' '
     print *, '--------------- REMAP 4D test ---------------------'
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
     layout1  => new_layout_4D( sll_world_collective )        
     call factorize_in_random_2powers_4d(colsz, npi, npj, npk, npl)
!!$     npi = 2
!!$     npj = 4
!!$     npk = 1
!!$     npl = 1
     if( myrank .eq. 0 ) then
        print *, 'source configuration: ', npi, npj, npk, npl
     end if

     call initialize_layout_with_distributed_4D_array( &
          ni, &
          nj, &
          nk, &
          nl, &
          npi, &
          npj, &
          npk, &
          npl, &
          layout1 )

     call compute_local_sizes_4d( &
          layout1, &
          loc_sz_i_init, &
          loc_sz_j_init, &
          loc_sz_k_init, &
          loc_sz_l_init )        

     SLL_ALLOCATE(local_array1(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init,loc_sz_l_init),ierr)
 
     ! initialize the local data    
     do l=1,loc_sz_l_init
        do k=1,loc_sz_k_init
           do j=1,loc_sz_j_init 
              do i=1,loc_sz_i_init
                 tmpa(:) = (/i, j, k, l/)
                 global_indices =  local_to_global_4D( layout1, tmpa )
                 gi = global_indices(1)
                 gj = global_indices(2)
                 gk = global_indices(3)
                 gl = global_indices(4)
                 local_array1(i,j,k,l) = &
                      gi + (gj-1)*ni + (gk-1)*ni*nj + (gl-1)*ni*nj*nk
              end do
           end do
        end do
     end do
     
     layout2  => new_layout_4D( sll_world_collective )
     call factorize_in_random_2powers_4d(colsz, npi, npj, npk, npl)
!!$     npi = 4
!!$     npj = 2
!!$     npk = 1
!!$     npl = 1

     if( myrank .eq. 0 ) then
        print *, 'target configuration: ', npi, npj, npk, npl
     end if

     call initialize_layout_with_distributed_4D_array( &
          ni, &
          nj, &
          nk, &
          nl, &
          npi, &
          npj, &
          npk, &
          npl, &
          layout2 )

     call reorganize_randomly_4d(layout2)
     
     call compute_local_sizes_4d( &
          layout2, &
          loc_sz_i_final, &
          loc_sz_j_final, &
          loc_sz_k_final, &
          loc_sz_l_final )

     SLL_ALLOCATE(local_array2(loc_sz_i_final,loc_sz_j_final,loc_sz_k_final,loc_sz_l_final), ierr )
    
     rmp4 => new_remap_plan( layout1, layout2, local_array1)     

     call apply_remap_4D( rmp4, local_array1, local_array2 ) 

     SLL_ALLOCATE(arrays_diff(loc_sz_i_final,loc_sz_j_final,loc_sz_k_final,loc_sz_l_final),ierr )
 
#if 0
     if( myrank .eq. 0 ) then
        print *, i_test, myrank, 'Printing layout1: '
        call sll_view_lims_4D( layout1 )
        print *, i_test, myrank, 'Printing layout2: '
        call sll_view_lims_4D( layout2 )
     end if
#endif
     ! compare results with expected data
     do l=1, loc_sz_l_final
        do k=1,loc_sz_k_final
           do j=1,loc_sz_j_final 
              do i=1,loc_sz_i_final
                 tmpa(:) = (/i,j,k,l/)
                 global_indices =  local_to_global_4D( layout2, tmpa )
                 gi = global_indices(1)
                 gj = global_indices(2)
                 gk = global_indices(3)
                 gl = global_indices(4)
                 arrays_diff(i,j,k,l) = local_array2(i,j,k,l) - &
                   (gi + (gj-1)*ni + (gk-1)*ni*nj + (gl-1)*ni*nj*nk)
                 if (arrays_diff(i,j,k,l) /=0 ) then
                    ok = 0
                    print*, i_test, myrank, '"remap" unit test: FAIL'
                    print *, i_test, myrank, 'local indices: ', &
                         '(', i, j, k, l, ')'
                    print *, 'global indices wrt target layout'
                    print*, myrank, &
                         'in global indices: (',gi,',', gj,',', gk,',',gl,')'
                    print*, myrank, &
                         'local array1(',gi,',', gj, ',', gk, ',', gl,')=', &
                         local_array1(i,j,k,l)  
                    print*, myrank, &
                         'local array2(',gi, ',', gj, ',', gk, ',', gl, ')=', &
                         local_array2(i,j,k,l)
                    g = theoretical_global_4D_indices(& 
                         int(local_array2(i,j,k,l),i32),&
                         ni, &
                         nj, &
                         nk ) ! por aqui
                    print*, 'Theoretical indices: (',g(1), ',', g(2),',',g(3),',',g(4), ')'
                    !     if(myrank .eq. 1) then
                    print *, local_array2(:,:,:,:)
                    !    end if
                    
                    print *, i_test, myrank, 'Printing layout1: '
                    call sll_view_lims_4D( layout1 )
                    print *, i_test, myrank, 'Printing layout2: '
                    call sll_view_lims_4D( layout2 )
                    
                    print*, 'program stopped by failure'
                    stop
                 end if
                 call flush(6)
              end do
           end do
        end do
     end do

     ! As the difference described above is null in each point, the sum of the 
     ! corresponding absolute values must be null. Each processor compute a local 
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
  
     call delete_layout_4D( layout1 )
     call delete_layout_4D( layout2 )
     SLL_DEALLOCATE_ARRAY(local_array1, ierr)
     SLL_DEALLOCATE_ARRAY(local_array2, ierr)
     SLL_DEALLOCATE_ARRAY(arrays_diff, ierr)
  enddo

  if (myrank==0) then
     if (prod4test(1)==1.) then
        print*, 'TEST PASSED'
     endif
  endif
  
  call sll_halt_collective()
  
contains

  subroutine reorganize_randomly_4d(layout)
    implicit none
    type(layout_4D), pointer   :: layout
    integer                    :: i, colsz, proc_n, proc_p
    real                       :: rand_real
    colsz = sll_get_num_nodes(layout)
    do i=0, colsz-1
       call random_number(rand_real)
       proc_n = int(rand_real*(colsz-1))
       call random_number(rand_real)
       proc_p = int(rand_real*(colsz-1))
       call swap_box_4D(proc_n, proc_p, layout)
    enddo    
  end subroutine reorganize_randomly_4d

  subroutine swap_box_4D(proc_n, proc_p, layout)
    implicit none
    integer                  :: proc_n, proc_p
    type(layout_4D), pointer :: layout
    integer                  :: i_min_n, i_max_n, j_min_n, j_max_n, &
    k_min_n, k_max_n, l_min_n, l_max_n, i_min_p, i_max_p, j_min_p, j_max_p, k_min_p, k_max_p, l_min_p, l_max_p    
    ! Get proc_n contents from layout
    i_min_n = get_layout_4D_i_min( layout, proc_n )
    i_max_n = get_layout_4D_i_max( layout, proc_n )
    j_min_n = get_layout_4D_j_min( layout, proc_n )
    j_max_n = get_layout_4D_j_max( layout, proc_n )
    k_min_n = get_layout_4D_k_min( layout, proc_n )
    k_max_n = get_layout_4D_k_max( layout, proc_n )
    l_min_n = get_layout_4D_l_min( layout, proc_n )
    l_max_n = get_layout_4D_l_max( layout, proc_n )

    ! Get proc_p contents from layout
    i_min_p = get_layout_4D_i_min( layout, proc_p )
    i_max_p = get_layout_4D_i_max( layout, proc_p )
    j_min_p = get_layout_4D_j_min( layout, proc_p )
    j_max_p = get_layout_4D_j_max( layout, proc_p )
    k_min_p = get_layout_4D_k_min( layout, proc_p )
    k_max_p = get_layout_4D_k_max( layout, proc_p )
    l_min_p = get_layout_4D_l_min( layout, proc_p )
    l_max_p = get_layout_4D_l_max( layout, proc_p )

    ! Set proc_n contents in layout
    call set_layout_4D_i_min( layout, proc_n, i_min_p )
    call set_layout_4D_i_max( layout, proc_n, i_max_p)
    call set_layout_4D_j_min( layout, proc_n, j_min_p )
    call set_layout_4D_j_max( layout, proc_n, j_max_p )
    call set_layout_4D_k_min( layout, proc_n, k_min_p )
    call set_layout_4D_k_max( layout, proc_n, k_max_p )
    call set_layout_4D_l_min( layout, proc_n, l_min_p )
    call set_layout_4D_l_max( layout, proc_n, l_max_p )

    ! Set proc_p contents in layout
    call set_layout_4D_i_min( layout, proc_p, i_min_n )
    call set_layout_4D_i_max( layout, proc_p, i_max_n )
    call set_layout_4D_j_min( layout, proc_p, j_min_n )
    call set_layout_4D_j_max( layout, proc_p, j_max_n )
    call set_layout_4D_k_min( layout, proc_p, k_min_n )
    call set_layout_4D_k_max( layout, proc_p, k_max_n )   
    call set_layout_4D_l_min( layout, proc_p, l_min_n )
    call set_layout_4D_l_max( layout, proc_p, l_max_n )
  end subroutine swap_box_4D

  function theoretical_global_4D_indices(d, ni, nj, nk)
    integer, dimension(1:4) :: theoretical_global_4D_indices
    integer, intent(in)     :: d, ni, nj, nk ! global array dims
    integer                 :: i, j, k, l
    integer                 :: q
    sll_real64              :: tmp
    ! We know that the linear index relates with the other
    ! indices by:
    !
    !          linear = i + ni*(j-1) + ni*nj*(k-1) + ni*nj*nk*(l-1)
    !
    ! where i,j,k,l are the indices sought. We start by working backwards.
    tmp = real(d,f64)/real(ni*nj*nk,f64)
    l   = ceiling(tmp)
    ! reduce the size of the number we are working with
    q   = d - (l-1)*ni*nj*nk
    tmp = real(q,f64)/real(ni*nj,f64)
    k   = ceiling(tmp)
    ! reduce again
    q   = q - (k-1)*ni*nj
    tmp = real(q,f64)/real(ni,f64)
    j   = ceiling(tmp)
    ! final reduction
    q   = q - (j-1)*ni
    i   = q
    theoretical_global_4D_indices(1:4) = (/i,j,k,l/) 
  end function theoretical_global_4D_indices
  
  subroutine factorize_in_random_2powers_4d(n, n1, n2, n3, n4)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2, n3, n4
    integer   :: expo, expo1, expo2, expo3, expo4
    if (.not.is_power_of_two(n)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = int(rand_real*(expo-expo1))
    call random_number(rand_real)
    expo3 = int(rand_real*(expo-(expo1+expo2)))
    expo4 = expo - (expo1 + expo2 + expo3)
    n1 = 2**expo1
    n2 = 2**expo2
    n3 = 2**expo3
    n4 = 2**expo4
  end subroutine factorize_in_random_2powers_4d

      
end program remap_test_4d
    
