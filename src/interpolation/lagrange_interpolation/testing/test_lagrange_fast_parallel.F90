! ---
! Program to demonstrate the simultaneous use of
! the domain decomposition module in concert with
! the fast Lagrange interpolation module.
! Interpolation of a distributed function is
! performed, a halo cell exchange is performed.
!
! To see what is going on,
! run the program with 2 or 4 MPI processes and
! plot the output files (intp_values_*.dat).
!
! (c) 2015, Klaus Reuter, khr@mpcdf.mpg.de
! ---


program test_lagrange_fast_parallel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_m_constants, only : sll_p_pi

  use sll_m_lagrange_fast, only : &
       sll_s_interpolate_array_disp_lagrange_fixed_halo_cells

  use sll_m_collective, only : &
       sll_s_boot_collective, &
       sll_f_get_collective_rank, &
       sll_f_get_collective_size, &
       sll_s_halt_collective, &
       sll_v_world_collective

  use sll_m_decomposition, only : &
       sll_o_apply_halo_exchange, &
       sll_t_cartesian_topology_6d, &
       sll_t_decomposition_6d, &
       sll_o_new_cartesian_domain_decomposition, &
       sll_o_new_cartesian_topology

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: nd = 6
  sll_int32 :: procs_per_dimension(nd)
  logical :: periodic(nd)
  type(sll_t_cartesian_topology_6d), pointer :: topology
  type(sll_t_decomposition_6d), pointer :: decomposition
  sll_int32 :: ierr
  sll_int32 :: world_size, my_rank
  sll_int32 :: global_grid_points_per_dimension(nd)
  sll_int32 :: halo_width_per_dimension(nd)
  sll_real64, dimension(:,:,:,:,:,:), pointer :: xi, fi, xp, fp
  sll_int32 :: i, stencil
  sll_real64 :: diff, alpha, diff_new
  character(len=128) :: filename
  character(len=6) :: rank_str

  ! --- initialize and perform 6d domain decomposition and interpolation

  call sll_s_boot_collective()
  world_size = sll_f_get_collective_size(sll_v_world_collective)
  my_rank = sll_f_get_collective_rank(sll_v_world_collective)
  write (rank_str,'(I6.6)') my_rank

  procs_per_dimension(1) = world_size
  procs_per_dimension(2:6) = 1

  periodic(1) = .true.
  periodic(2:6) = .false.
  topology => &
    sll_o_new_cartesian_topology(sll_v_world_collective, procs_per_dimension, periodic)


  global_grid_points_per_dimension(1) = 128
  global_grid_points_per_dimension(2:6) = 1
  halo_width_per_dimension(1) = 2
  halo_width_per_dimension(2:6) = 0
  decomposition => &
    sll_o_new_cartesian_domain_decomposition(topology, global_grid_points_per_dimension, halo_width_per_dimension)

  !write(*,*) decomposition%local%lo, decomposition%local%hi

  allocate( xi(decomposition%local%lo(1):decomposition%local%hi(1), &
    decomposition%local%lo(2):decomposition%local%hi(2), &
    decomposition%local%lo(3):decomposition%local%hi(3), &
    decomposition%local%lo(4):decomposition%local%hi(4), &
    decomposition%local%lo(5):decomposition%local%hi(5), &
    decomposition%local%lo(6):decomposition%local%hi(6)),&
    stat=ierr )
  SLL_ASSERT( ierr == 0 )
  allocate( fi(decomposition%local%lo(1):decomposition%local%hi(1), &
    decomposition%local%lo(2):decomposition%local%hi(2), &
    decomposition%local%lo(3):decomposition%local%hi(3), &
    decomposition%local%lo(4):decomposition%local%hi(4), &
    decomposition%local%lo(5):decomposition%local%hi(5), &
    decomposition%local%lo(6):decomposition%local%hi(6)),&
    stat=ierr )
  SLL_ASSERT( ierr == 0 )
  allocate( xp(decomposition%local%lo(1):decomposition%local%hi(1), &
    decomposition%local%lo(2):decomposition%local%hi(2), &
    decomposition%local%lo(3):decomposition%local%hi(3), &
    decomposition%local%lo(4):decomposition%local%hi(4), &
    decomposition%local%lo(5):decomposition%local%hi(5), &
    decomposition%local%lo(6):decomposition%local%hi(6)),&
    stat=ierr )
  SLL_ASSERT( ierr == 0 )
  allocate( fp(decomposition%local%lo(1):decomposition%local%hi(1), &
    decomposition%local%lo(2):decomposition%local%hi(2), &
    decomposition%local%lo(3):decomposition%local%hi(3), &
    decomposition%local%lo(4):decomposition%local%hi(4), &
    decomposition%local%lo(5):decomposition%local%hi(5), &
    decomposition%local%lo(6):decomposition%local%hi(6)),&
    stat=ierr )
  SLL_ASSERT( ierr == 0 )


  ! --- set up interpolation arrays

  xi = 0.0_f64
  fi = 0.0_f64
  xp = 0.0_f64
  fp = 0.0_f64

  ! --- offset from grid points
  alpha = 0.21_f64

  do i = decomposition%local%mn(1), decomposition%local%mx(1)
    xi(i,1,1,1,1,1) = real(i - 1, f64)
    fi(i,1,1,1,1,1) = f(xi(i,1,1,1,1,1), global_grid_points_per_dimension(1))
    xp(i,1,1,1,1,1) = xi(i,1,1,1,1,1) + alpha
  end do

  ! --- exchange boundaries
  call sll_o_apply_halo_exchange(topology, decomposition, fi)


  ! --- apply interpolation
  stencil = 2 * halo_width_per_dimension(1) + 1
  call sll_s_interpolate_array_disp_lagrange_fixed_halo_cells(fi(:,1,1,1,1,1), fp(:,1,1,1,1,1), alpha, stencil)


  ! --- write out original values
  filename = "grid_values_" // rank_str // ".dat"
  call dump(xi(decomposition%local%mn(1):decomposition%local%mx(1),1,1,1,1,1), &
    fi(decomposition%local%mn(1):decomposition%local%mx(1),1,1,1,1,1), &
    filename)

  ! --- write out interpolated values
  filename = "intp_values_" // rank_str // ".dat"
  call dump(xp(decomposition%local%mn(1):decomposition%local%mx(1),1,1,1,1,1), &
    fp(decomposition%local%mn(1):decomposition%local%mx(1),1,1,1,1,1), &
    filename)


  diff = 0.0_f64
  do i=decomposition%local%mn(1), decomposition%local%mx(1)
    diff_new = abs(f(xp(i,1,1,1,1,1), global_grid_points_per_dimension(1)) - fp(i,1,1,1,1,1))
    if (diff_new > diff) then
      ! print *, "### diff=", diff_new, ", i=", i
      diff = diff_new
    endif
  end do

  ! optional: do an MPI min reduction on diff

  deallocate(xi, fi, xp, fp)

  if (diff < 1.e-6) then
    print *, ""
    print *, "Fast Lagrange parallel interpolation unit test: PASSED"
    print *, "MPI rank=", my_rank, ", error =", diff
  else
    print *, ""
    print *, "Fast Lagrange parallel interpolation unit test: FAILED"
    print *, "MPI rank=", my_rank, ", error =", diff
  end if

  call sll_s_halt_collective()
  stop

contains

  function f(x, num_points)
    sll_int32, intent(in) :: num_points
    sll_real64, intent(in) :: x
    sll_real64 :: f
    f = cos(2*sll_p_pi*x/num_points)
  end function

  subroutine dump(x, y, filename)
    sll_real64, intent(in) :: x(:), y(:)
    character(len=*), intent(in) :: filename
    integer, parameter :: fd = 42
    integer :: n
    n = size(x)
    SLL_ASSERT( n == size(y) )
    open(unit=fd, file=trim(filename), action="write", status="replace")
    do i=1,n
      write (fd,*) x(i), y(i)
    enddo
    close(fd)
  end subroutine

end program
