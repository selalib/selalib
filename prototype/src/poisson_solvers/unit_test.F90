
!***************************************************************************
!
! Selalib 2012     
! Module: unit_test.F90
!
!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Last modification: March 23, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!> Pierre NAVARO (navaro@math.unistra.fr)
!                                  
!***************************************************************************

program test_poisson_solvers

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_collective
  use sll_poisson_1d_periodic
  use sll_poisson_2d_periodic
  use sll_poisson_3d_periodic_seq
  use sll_poisson_3d_periodic_par

  implicit none

  sll_int32  :: nx, ny, nz
  ! nx, ny, nz are the numbers of points - 1 in directions x, y, z
  sll_real64 :: Lx, Ly, Lz
  sll_int64  :: colsz

  !Boot parallel environment
  call sll_boot_collective()

  nx = 64
  ny = 128
  nz = 256
  Lx = 2*sll_pi
  Ly = 2*sll_pi
  Lz = 2*sll_pi

  colsz  = sll_get_collective_size(sll_world_collective)

  if (colsz==1) then
     print*, ' '
     print*, 'Testing poisson_1d ...'
     print*, ' '
     call test_poisson_1d()
     print*, ' '
     print*, 'Testing poisson_2d ...'
     print*, ' '
     call test_poisson_2d()
     print*, ' '
  endif

  call test_sll_poisson_3d_periodic(nx, ny, nz, Lx, Ly, Lz)

  call sll_halt_collective()

contains

  subroutine test_poisson_1d()

    type (mesh_descriptor_1D), pointer :: geometry

    type (field_1D_vec1), pointer      :: ex
    type (field_1D_vec1), pointer      :: ex_exact
    type (field_1D_vec1), pointer      :: rho
    type (poisson_1d_periodic)         :: poisson

    sll_int32   :: nc_eta1
    sll_real64  :: eta1_min, eta1_max
    sll_real64  :: delta_eta1
    sll_int32   :: error
    sll_int32   :: mode
    sll_int32   :: i

    eta1_min = 0.0; eta1_max = 2*sll_pi;
    nc_eta1 = 128

    geometry => new_mesh_descriptor_1D( eta1_min, eta1_max, nc_eta1, PERIODIC )
    rho      => new_field_1D_vec1( geometry )
    ex       => new_field_1D_vec1( geometry )
    ex_exact => new_field_1D_vec1( geometry )

    call new(poisson, nc_eta1, error) 

    mode = 7
    delta_eta1 = geometry%delta_eta1
    do i=1,nc_eta1+1
       rho%data(i)      =  mode**2*sin(mode*(i-1)*delta_eta1)
       ex_exact%data(i) = -mode*cos(mode*(i-1)*delta_eta1)
    end do

    call solve(poisson, ex, rho)

    print*,'mode=',mode,'   error=',maxval(abs(FIELD_DATA(ex)-FIELD_DATA(ex_exact)))

  end subroutine test_poisson_1d

  subroutine test_poisson_2d()

    sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
    sll_int32   :: nc_eta1, nc_eta2
    sll_int32   :: error

    type(poisson_2d_periodic), pointer :: poisson_e_fields
    type(poisson_2d_periodic), pointer :: poisson_potential
    type(geometry_2D),         pointer :: geom
    type(mesh_descriptor_2D),  pointer :: mesh
    type(field_2D_vec2),       pointer :: exy, exy_exact
    type(field_2D_vec1),       pointer :: rho, phi, phi_exact
    sll_real64                         :: x1, x2
    sll_int32                          :: mode
    sll_int32                          :: i, j

    eta1_min = .0_f64; eta1_max = 2.0_f64*sll_pi
    eta2_min = .0_f64; eta2_max = 2.0_f64*sll_pi

    geom => new_geometry_2D ('cartesian')

    nc_eta1 = 127; nc_eta2 = 127

    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
         PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

    call write_mesh_2D(mesh)

    rho       => new_field_2D_vec1(mesh)
    phi       => new_field_2D_vec1(mesh)
    phi_exact => new_field_2D_vec1(mesh)
    exy       => new_field_2D_vec2(mesh)
    exy_exact => new_field_2D_vec2(mesh)

    poisson_e_fields  => new_poisson_2d_periodic(exy)
    poisson_potential => new_poisson_2d_periodic(phi)

    mode = 2
    do i = 1, nc_eta1+1
       do j = 1, nc_eta2+1
          x1 = eta1_min+(i-1)*mesh%delta_eta1
          x2 = eta2_min+(j-1)*mesh%delta_eta2
          phi%data(i,j) = mode * sin(mode*x1) * cos(mode*x2)
          rho%data(i,j) = -2_f64 * mode**3 * sin(mode*x1)*cos(mode*x2)
          exy%data(i,j)%v1 =  mode**2*cos(mode*x1)*cos(mode*x2)
          exy%data(i,j)%v2 = -mode**2*sin(mode*x1)*sin(mode*x2)
       end do
    end do

    call write_vec1d(phi%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"phi0","mesh",0)
    call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho0","mesh",0)
    call write_vec2d(exy%data%v1,exy%data%v2,mesh%nc_eta1+1,mesh%nc_eta2+1,"exy0","mesh",0)

    FIELD_DATA(exy_exact) = FIELD_DATA(exy)
    FIELD_DATA(phi_exact) = FIELD_DATA(phi)

    call solve_poisson_2d_periodic(poisson_potential,phi%data,rho%data,error)
    call solve_poisson_2d_periodic(poisson_e_fields,exy,rho,error)

    call write_vec1d(phi%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"phi1","mesh",0)
    call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho1","mesh",0)
    call write_vec2d(exy%data%v1,exy%data%v2,mesh%nc_eta1+1,mesh%nc_eta2+1,"exy1","mesh",0)

    write(*,*) " Ex Error = " , maxval(abs(exy_exact%data%v1-exy%data%v1))
    write(*,*) " Ey Error = " , maxval(abs(exy_exact%data%v2-exy%data%v2))
    write(*,*) " Po Error = " , maxval(abs(phi_exact%data+phi%data))

    call delete_poisson_2d_periodic(poisson_e_fields)
    call delete_poisson_2d_periodic(poisson_potential)

    call delete_field_2D_vec1( rho )
    call delete_field_2D_vec1( phi )
    call delete_field_2D_vec1( phi_exact )
    call delete_field_2D_vec2( exy )
    call delete_field_2D_vec2( exy_exact )

  end subroutine test_poisson_2d

  subroutine test_sll_poisson_3d_periodic(nx, ny, nz, Lx, Ly, Lz)

    sll_int32                                    :: nx, ny, nz
    ! nx, ny, nz are the numbers of points - 1 in directions x, y, z
    sll_int32                                    :: nx_loc, ny_loc, nz_loc
    sll_int32                                    :: ierr
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: dx, dy, dz
    sll_real64                                   :: x, y, z
    sll_real64, dimension(nx,ny,nz)              :: rho_seq_1, rho_seq_2
    sll_real64, dimension(:,:,:), allocatable    :: rho_par_1, rho_par_2
    sll_real64, dimension(nx,ny,nz)              :: phi_an_1, phi_an_2
    sll_real64, dimension(nx,ny,nz)              :: phi_seq_1, phi_seq_2
    sll_real64, dimension(:,:,:), allocatable    :: phi_par_1, phi_par_2
    sll_int32                                    :: i, j, k
    type (poisson_3d_periodic_plan_seq), pointer :: plan_seq
    type (poisson_3d_periodic_plan_par), pointer :: plan_par
    sll_real64                                   :: average_err_1
    sll_real64                                   :: average_err_2
    sll_real64                                   :: seq_par_diff_1
    sll_real64                                   :: seq_par_diff_2
    sll_int32, dimension(1:3)                    :: global
    sll_int32                                    :: gi, gj, gk
    sll_int32                                    :: myrank
    sll_real32                                   :: ok = 1.d0
    sll_real32, dimension(1)                     :: prod4test
    type(layout_3D_t), pointer                   :: layout
    sll_int64                                    :: colsz ! collective size
    sll_int32                                    :: npx, npy, npz
    ! npx, npy, npz are the numbers of processors in directions x, y, z
    sll_int32                                    :: e

    dx = Lx/nx
    dy = Ly/ny
    dz = Lz/nz

    do k=1,nz
       z = (k-1)*dz
       do j=1,ny
          y = (j-1)*dy
          do i=1,nx
             x = (i-1)*dx
             phi_an_1(i,j,k)  = cos(x)*sin(y)*cos(z)
             phi_an_2(i,j,k)  = (4/(sll_pi*sqrt(sll_pi)*Lx*Ly*Lz)) *  exp(-.5 &
                                * (x-Lx/2)**2) * exp(-.5*(y-Ly/2)**2) * sin(z)
             rho_seq_2(i,j,k) = phi_an_2(i,j,k) * ( 3 - ( (x-Lx/2)**2 + &
                                (y-Ly/2)**2 ) )
          enddo
       enddo
    enddo
    rho_seq_1 = 3*phi_an_1

    colsz  = sll_get_collective_size(sll_world_collective)
    myrank = sll_get_collective_rank(sll_world_collective)

    ! Test sequential periodic 3D poisson solver
    if (myrank==0) then
       call flush()
       print*, ' '
       print*, 'Testing poisson_3d (2 equations here ) in sequential ...'
    endif

    plan_seq => new_poisson_3d_periodic_plan_seq(nx, ny, nz, Lx, Ly, Lz)
    call solve_poisson_3d_periodic_seq(plan_seq, rho_seq_1, phi_seq_1)
    call solve_poisson_3d_periodic_seq(plan_seq, rho_seq_2, phi_seq_2)

    average_err_1 = sum( abs(phi_an_1-phi_seq_1) ) / (nx*ny*nz)
    average_err_2 = sum( abs(phi_an_2-phi_seq_2) ) / (nx*ny*nz)

    if (myrank==0) then
       call flush()
       print*, ' '
       call flush()
       print*, 'Average error for equation 1:', average_err_1
       call flush()
       print*, 'Average error for equation 2:', average_err_2
       call flush()
       print*, 'dx*dy*dz =', dx*dy*dz
    endif

    if ( max(average_err_1, average_err_2) <= dx*dx*dy) then
       if (myrank==0) then
          call flush()
          print*, ' '
          print*, '"sll_poisson_3d_periodic_seq" test: PASS'
       endif
    else
       call flush()
       print*, ' '
       print*, 'Test stopped by "sll_poisson_3d_periodic_seq" failure'
       print*, ' '
       stop
    endif

    ! Test parallel periodic 3D poisson solver

    if (myrank==0) then
       call flush()
       print*, ' '
       call flush()
       print*, 'Testing poisson_3d (2 equations here ) in parallel ...'
    endif

    colsz  = sll_get_collective_size(sll_world_collective)
    e = int(log(real(colsz))/log(2.))

    ! Layout and local sizes for FFTs in x-direction
    layout => new_layout_3D( sll_world_collective )
    npx = 1
    npy = 2**(e/2)
    npz = int(colsz)/npy
    call initialize_layout_with_distributed_3D_array( nx, ny, &
                                    nz, npx, npy, npz, layout )

    plan_par => new_poisson_3d_periodic_plan_par(layout, nx, &
                                           ny, nz, Lx, Ly, Lz)

    call compute_local_sizes( layout, nx_loc, ny_loc, nz_loc )
    SLL_ALLOCATE(rho_par_1(nx_loc,ny_loc,nz_loc), ierr)
    SLL_ALLOCATE(rho_par_2(nx_loc,ny_loc,nz_loc), ierr)

    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = local_to_global_3D( layout, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             rho_par_1(i,j,k) = rho_seq_1(gi,gj,gk)
             rho_par_2(i,j,k) = rho_seq_2(gi,gj,gk)
          enddo
       enddo
    enddo

    SLL_ALLOCATE(phi_par_1(nx_loc,ny_loc,nz_loc), ierr)
    SLL_ALLOCATE(phi_par_2(nx_loc,ny_loc,nz_loc), ierr)

    call solve_poisson_3d_periodic_par(plan_par, rho_par_1, phi_par_1)
    call solve_poisson_3d_periodic_par(plan_par, rho_par_2, phi_par_2)

    average_err_1  = 0.d0
    seq_par_diff_1 = 0.d0
    average_err_2  = 0.d0
    seq_par_diff_2 = 0.d0

    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = local_to_global_3D( layout, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             average_err_1  = average_err_1  + abs( phi_an_1 (gi,gj,gk) &
                              - phi_par_1(i,j,k) )
             seq_par_diff_1 = seq_par_diff_1 + abs( phi_seq_1(gi,gj,gk) &
                              - phi_par_1(i,j,k) )
             average_err_2  = average_err_2  + abs( phi_an_2 (gi,gj,gk) &
                              - phi_par_2(i,j,k) )
             seq_par_diff_2 = seq_par_diff_2 + abs( phi_seq_2(gi,gj,gk) &
                              - phi_par_2(i,j,k) )
          enddo
       enddo
    enddo

    average_err_1  = average_err_1  / (nx_loc*ny_loc*nz_loc)
    seq_par_diff_1 = seq_par_diff_1 / (nx_loc*ny_loc*nz_loc)
    average_err_2  = average_err_2  / (nx_loc*ny_loc*nz_loc)
    seq_par_diff_2 = seq_par_diff_2 / (nx_loc*ny_loc*nz_loc)

    call flush()
    print*, ' '
    call flush()
    print*, 'Average error for equation 1, in proc', myrank, &
            ':', average_err_1
    call flush()
    print*, 'Average error for equation 2, in proc', myrank, &
            ':', average_err_2
    call flush()
    print*, 'dx*dy*dz =', dx*dy*dz
    call flush()
    print*, 'Average diff between seq sol and par sol for ', &
            'equation 1, in proc', myrank, ':', seq_par_diff_1
    call flush()
    print*, 'Average diff between seq sol and par sol for', &
            'equation 2, in proc', myrank, ':', seq_par_diff_2


    if ( max(average_err_1, average_err_2) > dx*dx*dy) then
       ok = 1.d0
       call flush()
       print*, ' '
       call flush()
       print*, 'Test stopped by "sll_poisson_3d_periodic_par" failure'
       call flush()
       print*, 'myrank=', myrank
       call flush()
       print*, ' '
       stop
    endif

    call sll_collective_reduce(sll_world_collective, (/ ok /), &
                                     1, MPI_PROD, 0, prod4test )

    if (myrank==0) then
       if (prod4test(1)==1.d0) then
          call flush()
          print*, ' '
          call flush()
          print*, '"sll_poisson_3d_periodic_par" test: PASS'
          call flush()
          print*, ' '
       endif
    endif

    call delete_poisson_3d_periodic_plan_par(plan_par)
    call delete_poisson_3d_periodic_plan_seq(plan_seq)

    SLL_DEALLOCATE_ARRAY(phi_par_1, ierr)
    SLL_DEALLOCATE_ARRAY(phi_par_2, ierr)
    SLL_DEALLOCATE_ARRAY(rho_par_1, ierr)
    SLL_DEALLOCATE_ARRAY(rho_par_2, ierr)

  end subroutine test_sll_poisson_3d_periodic

end program test_poisson_solvers

