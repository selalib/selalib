program poisson3d_with__fftw3



#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_remap.h"

  use sll_collective

  implicit none
  
  sll_int32 :: i, j, k, mx, my, mz
  sll_int32 :: nx, ny, nz
  sll_int32 :: ncount, t0, tn, istep, nstep
  sll_int32 :: nthreads = 4
  sll_int32 :: forward, backward
  sll_int64 :: col_size, myrank
  sll_int32 :: local_sz_x, local_sz_y, local_sz_z 
  sll_int32 :: npx, npy, npz ! numbers of procs in the directions
  sll_int32 :: e, e1, e2, e3
  sll_int32 :: ierr, gi, gj, gk

  sll_real64 :: x, y, z, tbegin, tend, err
  sll_real64 :: dx, dy, dz, time, dt, pi 
  sll_real64 :: cx, cy, cz, vpx, vpy, vpz 
  
  sll_real64, dimension(:,:,:), allocatable :: b, c, d, b_global
  sll_real64, dimension(:,:,:), allocatable :: u, f, g

  sll_int32, dimension(1:3) :: global

  type(layout_3d_t), pointer     :: layout1, layout2
  type(remap_plan_3d_t), pointer :: rmp3

  col_size = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  e = int(log(real(col_size))/log(2.))
  e1 = e/3
  e2 = (e-e1)/2
  e3 = e - (e1+e2)
  npx = 2**e1
  npy = 2**e2
  npz = 2**e3

  if ( col_size > min(nx-2,ny-2,nz-2) ) then
     if (myrank==0) then
        print*, 'the number of processors must be <=', min(nx,ny,nz)
     endif
     stop
  endif
  if ( .not.(is_power_of_two(col_size)) .or. (.not.is_power_of_two(int(nx,i64)-2)) .or. & 
       (.not.is_power_of_two(int(ny,i64)-2) ) .or. (.not.is_power_of_two(int(nz,i64)-2)) ) then
     if (myrank==0) then
        print*, 'this need to run with 2 power numbers of processors and of points.'
     endif
     stop
  endif  
  
  nx = 130; ny = 130; nz = 130
  dx = 1d0 / nx ; dy = 1d0 / ny ; dz = 1d0 / nz
  dt = 0.1; nstep = 10
  pi = 4d0 * datan(1d0)
  
  SLL_ALLOCATE(u(nx,ny,nz), ierr); u = 0.
  SLL_ALLOCATE(f(nx,ny,nz), ierr); f = 0.
  SLL_ALLOCATE(g(nx,ny,nz), ierr); g = 0.
  
  call init_poisson_fftw(nx,ny,nz)
  call system_clock(count=t0, count_rate=ncount)
  call cpu_time(tbegin)
  
  time = 0.0
  do istep = 1, nstep !loop over time
     
     do k = 1, nz
        z = (k-1)*dz
        do j = 1, ny
           y = (j-1)*dy
           do i = 1, nx
              x = (i-1)*dx
              f(i,j,k) = 12.*pi*pi*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)*cos(2*pi*time)
              g(i,j,k) = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)*cos(2*pi*time)
           enddo
        enddo
     enddo
     
     call solve_poisson_fftw(u, f, nx, ny, nz )
     
     write(*,'("istep = ",i4," time = ",g15.5," nx,ny,nz =",3i4)') &
          istep, time, nx, ny, nz
     
     err = maxval(u-g)
     print  '(3x,"norm(|ux - u_a|)               : ",2x,1pe10.3,//)', err
     
     time = time+dt
     
  end do !next time step
  
  call finalize_poisson_fftw()
  
  call system_clock(count=tn, count_rate=ncount)
  write(*,"(' elapsed time ', g15.5, ' s')") (tn - t0)/float(ncount)
  call cpu_time(tend)
  write(*,"(' cpu time ', g15.5, ' s')") tend-tbegin

  SLL_DEALLOCATE_ARRAY( u, ierr )
  SLL_DEALLOCATE_ARRAY( f, ierr )
  SLL_DEALLOCATE_ARRAY( g, ierr )
  
contains
  
  subroutine init_poisson_fftw(nx,ny,nz)
#include "fftw3.f"
    sll_int32, intent(in) :: nx,ny,nz
    
    mx = nx-2; my = ny-2; mz = nz-2
    
    cx = 1.0 / (dx*dx)
    cy = 1.0 / (dy*dy)
    cz = 1.0 / (dz*dz)

    layout1  => new_layout_3d( sll_world_collective )  
    call initialize_layout_with_distributed_3d_array( mx, my, mz, npx, npy, npz, layout1 )  
    call compute_local_sizes( layout1, local_sz_x, local_sz_y, local_sz_z )   
    
    !rhs
    SLL_ALLOCATE(b(local_sz_x, local_sz_y, local_sz_z ), ierr)
    SLL_ALLOCATE(d(local_sz_x, local_sz_y, local_sz_z ), ierr)
    SLL_ALLOCATE(c(local_sz_x, local_sz_y, local_sz_z ), ierr)
    
    !compute eigen values, build the matrix
    do k=1,local_sz_z
       do j=1,local_sz_y
          do i=1,local_sz_x
             global = local_to_global_3d( layout1, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             vpx=1.0-cos(float(gi)*pi/float(mx+1))
             vpy=1.0-cos(float(gj)*pi/float(my+1))
             vpz=1.0-cos(float(gk)*pi/float(mz+1))
             d(i,j,k)= 2.0 * (cz*vpz + cx*vpx + cy*vpy) &
                  * (8*(mx+1)*(my+1)*(mz+1))
          end do
       end do
    end do
    
    !initialize ffts
    call dfftw_init_threads(ierr)
    if (ierr == 0) stop 'fftw can''t use threads'
    call dfftw_plan_with_nthreads(nthreads)
    call dfftw_plan_r2r_3d(forward,mx,my,mz,c,b, &
         fftw_rodft00, fftw_rodft00,fftw_rodft00,&
         fftw_patient )
    call dfftw_plan_r2r_3d(backward,mx,my,mz,b,c, &
         fftw_rodft00, fftw_rodft00,fftw_rodft00,&
         fftw_patient )

    SLL_DEALLOCATE_ARRAY( b, ierr)
    SLL_DEALLOCATE_ARRAY( d, ierr)
    SLL_DEALLOCATE_ARRAY( c, ierr) 
    
  end subroutine init_poisson_fftw
  
  subroutine solve_poisson_fftw(psi, rhs, nx, ny, nz)
#include "fftw3.f"
    sll_int32, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: psi(nx,ny,nz)
    real(8), intent(in)  :: rhs(nx,ny,nz)
    
    !boundary conditions
    do k=1,local_sz_z
       do j=1,local_sz_y
          do i=1,local_sz_x
             global = local_to_global_3d( layout1, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             c(i,j,k) = rhs(gi+1,gj+1,gk+1)
             if ( gi == 1) then
                c(i,j,k) = c(i,j,k) + cx * psi( 1,gj+1,gk+1)
             else if (i == mx) then
                c(i,j,k) = c(i,j,k) + cx * psi(nx,gj+1,gk+1)
             else if ( j == 1) then
                c(i,j,k) = c(i,j,k) + cy * psi( gi+1,1,gk+1)
             else if (gj == my) then
                c(i,j,k) = c(i,j,k) + cy * psi(gi+1,ny,gk+1)
             else if ( k == 1) then
                c(i,j,k) = c(i,j,k) + cz * psi(gi+1,gj+1,1)
             else if (gk == mz) then
                c(i,j,k) = c(i,j,k) + cz * psi(gi+1,gj+1,nz)
             end if
          end do
       end do
    end do
!forward fft on the right hand side term
    call dfftw_execute_r2r(forward,c,b)
    
    b = b / d 

    layout2  => new_layout_3d( sll_world_collective )
    call initialize_layout_with_distributed_3d_array( mx, my, mz, 1, 1, 1, layout2 )
    call compute_local_sizes( layout2, local_sz_x, local_sz_y, local_sz_z )  
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, b)
    SLL_ALLOCATE( b_global(mx, my, mz), ierr )
    call apply_remap_3d( rmp3, b, b_global )
    
    !backward fft on the solution
    call dfftw_execute_r2r(backward,b_global,psi(2:nx-1,2:ny-1,2:nz-1))

    SLL_DEALLOCATE_ARRAY( b_global, ierr ) 
    
    return
    
  end subroutine solve_poisson_fftw
  
  subroutine finalize_poisson_fftw()
#include "fftw3.f"
    
    call dfftw_destroy_plan(forward)
    call dfftw_destroy_plan(backward)
    call dfftw_cleanup_threads(ierr)
    
  end subroutine finalize_poisson_fftw
 
  subroutine compute_local_sizes( layout, loc_sz_i, loc_sz_j, loc_sz_k )
    type(layout_3d_t), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_3d'
       print *, 'exiting...'
       stop
    end if
    my_rank = sll_get_collective_rank(get_layout_3d_collective(layout))
    i_min = get_layout_3d_i_min( layout, my_rank )
    i_max = get_layout_3d_i_max( layout, my_rank )
    j_min = get_layout_3d_j_min( layout, my_rank )
    j_max = get_layout_3d_j_max( layout, my_rank )
    k_min = get_layout_3d_k_min( layout, my_rank )
    k_max = get_layout_3d_k_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
  end subroutine compute_local_sizes
  
end program poisson3d_with__fftw3
