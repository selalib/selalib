
!*******************************************************************
!
! Selalib      
! Module: poisson_3d_periodic.F90
!
!> @brief 
!> 3D poisson solver with fftw3
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Pierre NAVARO (navaro@math.unistra.fr)
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!*******************************************************************

program poisson_3d

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_remap.h"
  
  use sll_collective
  use sll_fft

  implicit none
  
  sll_int32                                 :: i, j, k, mx, my, mz
  sll_int32                                 :: nx, ny, nz
  sll_int32                                 :: ncount, t0, tn, istep, nstep
  sll_int64                                 :: col_size, myrank
  sll_int32                                 :: npx, npy, npz ! numbers of procs in the directions
  sll_int32                                 :: e, e1, e2, e3
  sll_int32                                 :: ierr, gi, gj, gk
  sll_real64                                :: x, y, z, tbegin, tend, err
  sll_real64                                :: dx, dy, dz, time, dt, pi 
  sll_real64                                :: cx, cy, cz, vpx, vpy, vpz   
  sll_comp64, dimension(:,:,:), allocatable :: b
  sll_real64, dimension(:,:,:), allocatable :: c, d
  sll_comp64, dimension(:,:,:), allocatable :: tmp1, tmp2
  sll_real64, dimension(:,:,:), allocatable :: u, f, g
  sll_int32, dimension(1:3)                 :: global
  type(layout_3d_t), pointer                :: layout1, layout2
  type(remap_plan_3d_t), pointer            :: rmp3
  type(sll_fft_plan), pointer               :: p

  call sll_boot_collective()

  col_size = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  
  nx = 130; ny = 130; nz = 130
  dx = 1d0 / nx ; dy = 1d0 / ny ; dz = 1d0 / nz
  dt = 0.1; nstep = 10
  pi = 4d0 * datan(1d0)

  if ( col_size > min(nx-2,ny-2,nz-2) ) then
     if (myrank==0) then
        print*, 'the number of processors must be <=', min(nx-2,ny-2,nz-2)
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
  
  SLL_ALLOCATE(u(nx,ny,nz), ierr); u = 0.
  SLL_ALLOCATE(f(nx,ny,nz), ierr); f = 0.
  SLL_ALLOCATE(g(nx,ny,nz), ierr); g = 0.

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
     
     call solve_poisson_3d(u, f, nx, ny, nz )
     write(*,'("istep = ",i4," time = ",g15.5," nx,ny,nz =",3i4)') &
          istep, time, nx, ny, nz
     
     err = maxval(u-g)
     print  '(3x,"norm(|ux - u_a|)               : ",2x,1pe10.3,//)', err
     
     time = time+dt
     
  end do !next time step
  
  !call finalize_poisson()
  
  call system_clock(count=tn, count_rate=ncount)
  write(*,"(' elapsed time ', g15.5, ' s')") (tn - t0)/float(ncount)
  call cpu_time(tend)
  write(*,"(' cpu time ', g15.5, ' s')") tend-tbegin

  call sll_halt_collective()
  
contains
  
  subroutine solve_poisson_3d(psi, rhs, nx, ny, nz)

    sll_int32, intent(in)                                  :: nx, ny, nz
    sll_real64, intent(out), dimension(:,:,:), allocatable :: psi
    sll_real64, intent(in)                                 :: rhs(nx,ny,nz)

    mx = nx-2; my = ny-2; mz = nz-2
    
    cx = 1.0 / (dx*dx)
    cy = 1.0 / (dy*dy)
    cz = 1.0 / (dz*dz)

    e = int(log(real(col_size))/log(2.))
    e2 = e/2
    e3 = e - e2
    npx = 1
    npy = 2**e2
    npz = 2**e3

    layout1  => new_layout_3d( sll_world_collective )  
    call initialize_layout_with_distributed_3d_array( mx, my, mz, npx, npy, npz, layout1 )     

    ! Boundary conditions
    SLL_ALLOCATE(c(mx/npx, my/npy, mz/npz), ierr)
    do k=1,mz/npz
       do j=1,my/npy
          do i=1,mx/npx
             global = local_to_global_3d( layout1, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             c(i,j,k) = rhs(gi+1,gj+1,gk+1)
             if ( gi == 1) then
                c(i,j,k) = c(i,j,k) + cx * psi(1,gj+1,gk+1)
             else if (gi == mx) then
                c(i,j,k) = c(i,j,k) + cx * psi(nx,gj+1,gk+1)
             else if ( gj == 1) then
                c(i,j,k) = c(i,j,k) + cy * psi(gi+1,1,gk+1)
             else if (gj == my) then
                c(i,j,k) = c(i,j,k) + cy * psi(gi+1,ny,gk+1)
             else if ( gk == 1) then
                c(i,j,k) = c(i,j,k) + cz * psi(gi+1,gj+1,1)
             else if (gk == mz) then
                c(i,j,k) = c(i,j,k) + cz * psi(gi+1,gj+1,nz)
             end if
          end do
       end do
    end do

    ! FFTs in x-direction
    SLL_ALLOCATE(tmp1(mx/npx, my/npy, mz/npz), ierr)
    tmp1 = cmplx(c, 0_f64, kind=f64)
    p => new_plan_c2c_1d( mx, tmp1(:,1,1), tmp1(:,1,1), FFT_FORWARD )
    do k=1,mz/npz
       do j=1,my/npy
          call apply_fft_c2c_1d( p, tmp1(:,j,k), tmp1(:,j,k) )
       enddo
    enddo
    call delete(p)

    ! FFTs in y-direction
    e1 = e/2
    e3 = e - e1
    npx = 2**e1
    npy = 1
    npz = 2**e3
    layout2  => new_layout_3d( sll_world_collective )
    call initialize_layout_with_distributed_3d_array( mx, my, mz, npx, npy, npz, layout2 )
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, tmp1)
    SLL_ALLOCATE(tmp2(mx/npx, my/npy, mz/npz), ierr)
    call apply_remap_3d( rmp3, tmp1, tmp2 )
    p => new_plan_c2c_1d( my, tmp2(1,:,1), tmp2(1,:,1), FFT_FORWARD )
    do k=1,mz/npz
       do i=1,mx/npx
          call apply_fft_c2c_1d( p, tmp2(i,:,k), tmp2(i,:,k) )
       enddo
    enddo
    call delete(p)
    SLL_DEALLOCATE_ARRAY(tmp1, ierr)

    ! FFTs in z-direction
    e1 = e/2
    e2 = e - e1
    npx = 2**e1
    npy = 2**e2
    npz = 1
    call initialize_layout_with_distributed_3d_array( mx, my, mz, npx, npy, npz, layout1 )
    rmp3 => NEW_REMAPPER_PLAN_3D( layout2, layout1, tmp2)
    SLL_ALLOCATE(b(mx/npx, my/npy, mz/npz), ierr)
    call apply_remap_3d( rmp3, tmp2, b )
    p => new_plan_c2c_1d( mz, b(1,1,:), b(1,1,:), FFT_FORWARD )
    do j=1,my/npy
       do i=1,mx/npx
          call apply_fft_c2c_1d( p, b(i,j,:), b(i,j,:) )
       enddo
    enddo
    call delete(p)
    SLL_DEALLOCATE_ARRAY(tmp2, ierr)

    !Compute eigen values, build the matrix   
    SLL_ALLOCATE(d(mx/npx, my/npy, mz/npz), ierr) 
    do k=1,mz/npz
       do j=1,my/npz
          do i=1,mx/npx
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
    
    b = b / d
    
    ! Inverse FFTs in z-direction
    p => new_plan_c2c_1d( mz, b(1,1,:), b(1,1,:), FFT_INVERSE )
    do j=1,my/npy
       do i=1,mx/npx
          call apply_fft_c2c_1d( p, b(i,j,:), b(i,j,:) )
       enddo
    enddo
    call delete(p)

    ! Inverse FFTs in y-direction
    e1 = e/2
    e3 = e - e1
    npx = 2**e1
    npy = 1
    npz = 2**e3
    layout1  => new_layout_3d( sll_world_collective )
    call initialize_layout_with_distributed_3d_array( mx, my, mz, npx, npy, npz, layout2 )
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, b)
    SLL_ALLOCATE(tmp2(mx/npx, my/npy, mz/npz), ierr)
    call apply_remap_3d( rmp3, b, tmp2 )
    p => new_plan_c2c_1d( my, tmp2(1,:,1), tmp2(1,:,1), FFT_INVERSE )
    do k=1,mz/npz
       do i=1,mx/npx
          call apply_fft_c2c_1d( p, tmp2(i,:,k), tmp2(i,:,k) )
       enddo
    enddo
    call delete(p)
    SLL_DEALLOCATE_ARRAY(tmp1, ierr)

    ! Inverse FFTs in x-direction
    e2 = e/2
    e3 = e - e2
    npx = 1
    npy = 2**e2
    npz = 2**e3
    layout1  => new_layout_3d( sll_world_collective )
    call initialize_layout_with_distributed_3d_array( mx, my, mz, npx, npy, npz, layout1 )
    rmp3 => NEW_REMAPPER_PLAN_3D( layout2, layout1, tmp2)
    call apply_remap_3d( rmp3, tmp2, b )
    p => new_plan_c2c_1d( mx, b(:,1,1), b(:,1,1), FFT_INVERSE )
    do k=1,mz/npz
       do i=1,mx/npx
          call apply_fft_c2c_1d( p, b(:,j,k), b(:,j,k) )
       enddo
    enddo
    call delete(p)
    SLL_DEALLOCATE_ARRAY(tmp1, ierr)

    SLL_ALLOCATE(psi(mx/npx, my/npy, mz/npz), ierr)
    psi = real(b, f64)
    
    return
    
  end subroutine solve_poisson_3d
 
end program poisson_3d
