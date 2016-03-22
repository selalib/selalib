!************************************************************************
!
! Selalib 2012     
! File : test_poisson_3d_par.F90
!
!> @brief 
!> Selalib poisson solvers (3D) unit test
!> Last modification: April 10, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!> Pierre NAVARO (navaro@math.unistra.fr)
!> Katharina Kormann
!                                  
!************************************************************************

program test_poisson_3d_periodic_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_o_collective_reduce, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_3d_periodic_par, only: &
    sll_s_poisson_3d_periodic_par_compute_e_from_phi, &
    sll_s_poisson_3d_periodic_par_free, &
    sll_s_poisson_3d_periodic_par_init, &
    sll_t_poisson_3d_periodic_par, &
    sll_s_poisson_3d_periodic_par_solve

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_3d, &
    sll_o_local_to_global, &
    sll_f_new_layout_3d

  use sll_mpi, only: &
    mpi_prod

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32                                    :: nx, ny, nz
  sll_int32                                    :: nx_loc, ny_loc, nz_loc
  sll_int32                                    :: ierr
  sll_real64                                   :: Lx, Ly, Lz
  sll_real64                                   :: dx, dy, dz
  sll_real64, dimension(:,:,:), allocatable    :: x, y, z
  sll_real64, dimension(:,:,:), allocatable    :: rho
  sll_real64, dimension(:,:,:), allocatable    :: phi_an
  sll_real64, dimension(:,:,:), allocatable    :: phi
  sll_real64, dimension(:,:,:), allocatable    :: ex
  sll_real64, dimension(:,:,:), allocatable    :: ey
  sll_real64, dimension(:,:,:), allocatable    :: ez
  sll_real64, dimension(:,:,:), allocatable    :: ex_an
  sll_real64, dimension(:,:,:), allocatable    :: ey_an
  sll_real64, dimension(:,:,:), allocatable    :: ez_an
  sll_int32                                    :: i, j, k
  type (sll_t_poisson_3d_periodic_par)         :: plan
  sll_real64                                   :: average_err
  sll_int32, dimension(1:3)                    :: global
  sll_int32                                    :: gi, gj, gk
  sll_int32                                    :: myrank
  type(sll_t_layout_3d), pointer               :: layout
  sll_int64                                    :: colsz ! collective size
  sll_int32                                    :: i_test
  sll_int32                                    :: npx, npy, npz
  sll_int32                                    :: e
  sll_real32                                   :: ok = 1.0
  sll_real32   , dimension(1)                  :: prod4test

  !Boot parallel environment
  call sll_s_boot_collective()

  nx = 64
  ny = 64
  nz = 64
  Lx = 2*sll_p_pi
  Ly = 2*sll_p_pi
  Lz = 2*sll_p_pi

  colsz  = int(sll_f_get_collective_size(sll_v_world_collective), i64)
  myrank = sll_f_get_collective_rank(sll_v_world_collective)


  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz

  e = int(log(real(colsz))/log(2.))

  ! Layout and local sizes for FFTs in x-direction
  layout => sll_f_new_layout_3d( sll_v_world_collective )
  npx = 1
  npy = 2**(e/2)
  npz = int(colsz)/npy
  call sll_o_initialize_layout_with_distributed_array( nx, ny, &
                                  nz, npx, npy, npz, layout )

  call sll_s_poisson_3d_periodic_par_init(layout, nx, ny, &
       nz, Lx, Ly, Lz, plan)


  call sll_o_compute_local_sizes( layout, nx_loc, ny_loc, nz_loc )

  SLL_ALLOCATE(rho(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(ex(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(ey(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(ez(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(ex_an(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(ey_an(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(ez_an(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(x(nx_loc,ny_loc,nz_loc),ierr)
  SLL_ALLOCATE(y(nx_loc,ny_loc,nz_loc),ierr)
  SLL_ALLOCATE(z(nx_loc,ny_loc,nz_loc),ierr)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc,nz_loc),ierr)

  do k=1,nz_loc
     do j=1,ny_loc
        do i=1,nx_loc
           global = sll_o_local_to_global( layout, (/i, j, k/))
           gi = global(1)
           gj = global(2)
           gk = global(3)
           x(i,j,k) = (gi-1)*dx
           y(i,j,k) = (gj-1)*dy
           z(i,j,k) = (gk-1)*dz
        end do
     end do
  end do

  do i_test = 1, 2

     if (i_test==1) then
        phi_an = cos(x)*sin(y)*cos(z)
        ex_an =  sin(x)*sin(y)*cos(z)
        ey_an =  -cos(x)*cos(y)*cos(z)
        ez_an =  cos(x)*sin(y)*sin(z)
     else if (i_test == 2) then
        phi_an = (4.0_f64/(sll_p_pi * sqrt(sll_p_pi)*Lx*Ly*Lz)) &
             * exp(-.5*(x-Lx/2)**2)                   &
             * exp(-.5*(y-Ly/2)**2) * sin(z)
        ex_an = phi_an * (x - Lx/2)
        ey_an = phi_an * (y - Ly/2)
        ez_an =  -(4.0_f64/(sll_p_pi * sqrt(sll_p_pi)*Lx*Ly*Lz)) &
             * exp(-.5*(x-Lx/2)**2)                   &
             * exp(-.5*(y-Ly/2)**2) * cos(z)
     end if

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              if (i_test == 1) then
                 rho(i,j,k) = 3*phi_an(i,j,k)
              else if(i_test == 2) then
                 rho(i,j,k) = phi_an(i,j,k)   &
                      *(3-((x(i,j,k)-Lx/2)**2 &
                      +(y(i,j,k)-Ly/2)**2))
              end if
           enddo
        enddo
     enddo

     SLL_ALLOCATE(phi(nx_loc,ny_loc,nz_loc), ierr)
     call sll_s_poisson_3d_periodic_par_solve(plan, rho, phi)

     average_err  = 0._f64

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              average_err  = average_err + abs( phi_an(i,j,k) &
                             - phi(i,j,k) )
           enddo
        enddo
     enddo

     average_err  = average_err  / (nx_loc*ny_loc*nz_loc)

     flush( output_unit ); print*, ' ------------------'
     flush( output_unit ); print*, ' myrank ', myrank
     flush( output_unit ); print*, 'local average error:', average_err
     flush( output_unit ); print*, 'dx*dy*dz =', dx*dy*dz
     flush( output_unit ); print*, ' ------------------'

     if (average_err> dx*dy*dz ) then
        print*, 'Test stopped by "sll_m_poisson_3d_periodic_par" failure'
        stop
     endif

     call sll_s_poisson_3d_periodic_par_compute_e_from_phi(plan, phi, ex, ey, ez)

     average_err  = 0._f64

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              average_err  = average_err + abs( ex_an(i,j,k) &
                             - ex(i,j,k) )
           enddo
        enddo
     enddo

     average_err  = average_err  / (nx_loc*ny_loc*nz_loc)

     flush( output_unit ); print*, ' ------------------'
     flush( output_unit ); print*, ' myrank ', myrank
     flush( output_unit ); print*, 'local average error in ex:', average_err
     flush( output_unit ); print*, ' ------------------'

     if (average_err> dx*dy*dz ) then
        print*, 'Test stopped by "sll_m_poisson_3d_periodic_par" failure'
        stop
     endif

     average_err  = 0._f64

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              average_err  = average_err + abs( ey_an(i,j,k) &
                             - ey(i,j,k) )
           enddo
        enddo
     enddo

     average_err  = average_err  / (nx_loc*ny_loc*nz_loc)

     flush( output_unit ); print*, ' ------------------'
     flush( output_unit ); print*, ' myrank ', myrank
     flush( output_unit ); print*, 'local average error in ey:', average_err
     flush( output_unit ); print*, ' ------------------'

     if (average_err> dx*dy*dz ) then
        print*, 'Test stopped by "sll_m_poisson_3d_periodic_par" failure'
        stop
     endif

     average_err  = 0._f64

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              average_err  = average_err + abs( ez_an(i,j,k) &
                             - ez(i,j,k) )
           enddo
        enddo
     enddo

     average_err  = average_err  / (nx_loc*ny_loc*nz_loc)

     flush( output_unit ); print*, ' ------------------'
     flush( output_unit ); print*, ' myrank ', myrank
     flush( output_unit ); print*, 'local average error in ez:', average_err
     flush( output_unit ); print*, ' ------------------'

     if (average_err> dx*dy*dz ) then
        print*, 'Test stopped by "sll_m_poisson_3d_periodic_par" failure'
        stop
     endif

     SLL_DEALLOCATE_ARRAY(phi, ierr)

  end do

     call sll_o_collective_reduce(sll_v_world_collective, (/ ok /), &
          1, MPI_PROD, 0, prod4test )
     if (myrank==0) then

        if (prod4test(1)==1.) then
           flush( output_unit )
           print*, ' '
           flush( output_unit )
           print*, '"sll_m_poisson_3d_periodic_par" test: PASSED'
           flush( output_unit )
           print*, ' '
        endif
     endif           

  call sll_s_poisson_3d_periodic_par_free(plan)

  SLL_DEALLOCATE_ARRAY(phi_an, ierr)
  SLL_DEALLOCATE_ARRAY(rho, ierr)
  call sll_s_halt_collective()

end program test_poisson_3d_periodic_par
