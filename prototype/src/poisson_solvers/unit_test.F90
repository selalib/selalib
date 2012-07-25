
!***************************************************************************
!
! Selalib 2012     
! Module: unit_test.F90
!
!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Last modification: April 10, 2012
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
  use sll_poisson_1d_periodic
  use sll_poisson_2d_periodic
  use sll_poisson_3d_periodic_seq

  implicit none


  print*, ' '
  print*, 'Testing poisson 1D ...'
  print*, ' '
  call test_poisson_1d()
  print*, ' '
  print*, 'Testing poisson 2D ...'
  print*, ' '
  call test_poisson_2d()
  print*, ' '
  print*, 'Testing poisson 3D ...'
  print*, ' '
  call test_poisson_3d()

contains

  subroutine test_poisson_1d()

    type (mesh_descriptor_1D), pointer :: geometry

    !type (field_1D_vec1), pointer      :: ex
    !type (field_1D_vec1), pointer      :: ex_exact
    !type (field_1D_vec1), pointer      :: rho
    sll_real64, dimension(:), allocatable :: ex
    sll_real64, dimension(:), allocatable :: ex_exact
    sll_real64, dimension(:), allocatable :: rho
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

    call new(poisson, eta1_min, eta1_max, nc_eta1, error) 

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

  subroutine test_poisson_3d()

    sll_int32                                    :: nc_eta1
    sll_int32                                    :: nc_eta2
    sll_int32                                    :: nc_eta3
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: delta_eta1
    sll_real64                                   :: delta_eta2
    sll_real64                                   :: delta_eta3
    sll_real64, dimension(:,:,:), allocatable    :: eta1
    sll_real64, dimension(:,:,:), allocatable    :: eta2
    sll_real64, dimension(:,:,:), allocatable    :: eta3
    sll_real64, dimension(:,:,:), allocatable    :: rho
    sll_real64, dimension(:,:,:), allocatable    :: phi_exact
    sll_real64, dimension(:,:,:), allocatable    :: phi
    sll_int32                                    :: i, j, k
    type (poisson_3d_periodic_plan_seq), pointer :: plan
    sll_real64                                   :: average_err
    sll_real64                                   :: vcell_volume
    sll_int32                                    :: error
    sll_int32                                    :: icase

    nc_eta1 = 64
    nc_eta2 = 128
    nc_eta3 = 256

    SLL_ALLOCATE(eta1(nc_eta1,nc_eta2,nc_eta3),      error)
    SLL_ALLOCATE(eta2(nc_eta1,nc_eta2,nc_eta3),      error)
    SLL_ALLOCATE(eta3(nc_eta1,nc_eta2,nc_eta3),      error)
    SLL_ALLOCATE(rho(nc_eta1,nc_eta2,nc_eta3),       error)
    SLL_ALLOCATE(phi(nc_eta1,nc_eta2,nc_eta3),       error)
    SLL_ALLOCATE(phi_exact(nc_eta1,nc_eta2,nc_eta3), error)

    Lx = 2*sll_pi
    Ly = 2*sll_pi
    Lz = 2*sll_pi

    delta_eta1 = Lx/nc_eta1
    delta_eta2 = Ly/nc_eta2
    delta_eta3 = Lz/nc_eta3

    do k=1,nc_eta3
       do j=1,nc_eta2
          do i=1,nc_eta1
             eta1(i,j,k) = (i-1)*delta_eta1
             eta2(i,j,k) = (j-1)*delta_eta2
             eta3(i,j,k) = (k-1)*delta_eta3
          enddo
       enddo
    enddo

    plan => new_poisson_3d_periodic_plan_seq(nc_eta1,nc_eta2,nc_eta3,Lx,Ly,Lz)

    do icase = 1, 2

       if (icase == 1) then
          phi_exact = cos(eta1)*sin(eta2)*cos(eta3)
          rho = 3*phi_exact
       else
          phi_exact = (4.0_f64/(sll_pi*sqrt(sll_pi)*Lx*Ly*Lz)) &
             * exp(-.5*(eta1-Lx/2)**2) * exp(-.5*(eta2-Ly/2)**2)*sin(eta3)
          rho = phi_exact * ( 3 - ( (eta1-Lx/2)**2 + (eta2-Ly/2)**2 ) )
       end if
   
       call solve_poisson_3d_periodic_seq(plan, rho, phi)
   
       average_err = sum( abs(phi-phi_exact) ) / (nc_eta1*nc_eta2*nc_eta3)

       vcell_volume = delta_eta1*delta_eta2*delta_eta3
       print*, ' '
       print*, 'Average error, vcell volume :', average_err, vcell_volume
   
       if ( average_err > vcell_volume) then
          print*, 'Test stopped by "sll_poisson_3d_periodic_seq" failure'
          stop
       endif

    end do

    print*, ' '
    print*, '"sll_poisson_3d_periodic_seq" test: PASS'
    print*, ' '

    SLL_DEALLOCATE_ARRAY(eta1,error)
    SLL_DEALLOCATE_ARRAY(eta2,error)
    SLL_DEALLOCATE_ARRAY(eta3,error)
    SLL_DEALLOCATE_ARRAY(phi,error)
    SLL_DEALLOCATE_ARRAY(phi_exact,error)
    SLL_DEALLOCATE_ARRAY(rho,error)
    call delete_poisson_3d_periodic_plan_seq(plan)

  end subroutine test_poisson_3d

end program test_poisson_solvers

