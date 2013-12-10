program poisson3d_with__fftw3_seq
  
  use fftw3
  implicit none
  
  integer :: i, j, k, mx, my, mz, ierr
  integer :: nx, ny, nz
  integer :: ncount, t0, tn, istep, nstep
  integer :: nthreads = 4
  
  integer(8) :: forward, backward
  
  real(8) :: x, y, z, tbegin, tend, err
  real(8) :: dx, dy, dz, time, dt, pi 
  real(8) :: cx, cy, cz, vpx, vpy, vpz
  
  real(8), dimension(:,:,:), allocatable :: b, c, d
  real(8), allocatable :: u(:,:,:)
  real(8), allocatable :: f(:,:,:)
  real(8), allocatable :: g(:,:,:)
  
  
  nx = 130; ny = 130; nz = 130
  dx = 1d0 / nx ; dy = 1d0 / ny ; dz = 1d0 / nz
  dt = 0.1; nstep = 10
  pi = 4d0 * datan(1d0)
  
  allocate(u(nx,ny,nz)); u = 0.
  allocate(f(nx,ny,nz)); f = 0.
  allocate(g(nx,ny,nz)); g = 0.
  
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
  
contains
  
  subroutine init_poisson_fftw(nx,ny,nz)

    integer, intent(in) :: nx,ny,nz
    
    mx = nx-2; my = ny-2; mz = nz-2
    
    cx = 1.0 / (dx*dx)
    cy = 1.0 / (dy*dy)
    cz = 1.0 / (dz*dz)
    
    !rhs
    allocate(b(mx,my,mz), d(mx,my,mz),c(mx,my,mz))
    
    !compute eigen values, build the matrix
    do k=1,mz
       do j=1,my
          do i=1,mx
             vpx=1.0-cos(float(i)*pi/float(mx+1))
             vpy=1.0-cos(float(j)*pi/float(my+1))
             vpz=1.0-cos(float(k)*pi/float(mz+1))
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
    
  end subroutine init_poisson_fftw
  
  subroutine solve_poisson_fftw(psi, rhs, nx, ny, nz)

    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: psi(nx,ny,nz)
    real(8), intent(in)  :: rhs(nx,ny,nz)
    
    !boundary conditions
    do i=1,mx
       do j=1,my
          do k = 1, mz
             c(i,j,k) = rhs(i+1,j+1,k+1)
             if ( i == 1) then
                c(i,j,k) = c(i,j,k) + cx * psi( 1,j+1,k+1)
             else if (i == mx) then
                c(i,j,k) = c(i,j,k) + cx * psi(nx,j+1,k+1)
             else if ( j == 1) then
                c(i,j,k) = c(i,j,k) + cy * psi( i+1,1,k+1)
             else if (j == my) then
                c(i,j,k) = c(i,j,k) + cy * psi(i+1,ny,k+1)
             else if ( k == 1) then
                c(i,j,k) = c(i,j,k) + cz * psi(i+1,j+1,1)
             else if (k == mz) then
                c(i,j,k) = c(i,j,k) + cz * psi(i+1,j+1,nz)
             end if
          end do
       end do
    end do
!forward fft on the right hand side term
    call dfftw_execute_r2r(forward,c,b)
    
    b = b / d 
    
    !backward fft on the solution
    call dfftw_execute_r2r(backward,b,psi(2:nx-1,2:ny-1,2:nz-1))
    
    
    return
    
  end subroutine solve_poisson_fftw
  
  subroutine finalize_poisson_fftw()
    
    call dfftw_destroy_plan(forward)
    call dfftw_destroy_plan(backward)
    call dfftw_cleanup_threads(ierr)
    
  end subroutine finalize_poisson_fftw
  
end program poisson3d_with__fftw3_seq
