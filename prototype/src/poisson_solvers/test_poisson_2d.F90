program test_poisson_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"

   !-------------------------------------------------------------------
   !  test 2D Poisson solver based on FFT
   !-------------------------------------------------------------------


   use numeric_constants
   use sll_poisson_2D_periodic

   implicit none

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
         exy%data(i,j)%v1 =  1_f64*mode**2*cos(mode*x1)*cos(mode*x2)
         exy%data(i,j)%v2 = -1_f64*mode**2*sin(mode*x1)*sin(mode*x2)
      end do
   end do

   call write_vec1d(phi%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"phi0","mesh",0)
   call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho0","mesh",0)
   call write_vec2d(exy%data%v1,exy%data%v2,mesh%nc_eta1+1,mesh%nc_eta2+1,"exy0","mesh",0)

   FIELD_DATA(exy_exact) = FIELD_DATA(exy)
   FIELD_DATA(phi_exact) = FIELD_DATA(phi)
   
   call solve_poisson_2d_periodic(poisson_potential,phi%data,rho%data,error)
   call solve_poisson_2d_periodic(poisson_e_fields,exy,rho,error)
   
   write(*,*) " Ex Error = " , maxval(abs(exy_exact%data%v1-exy%data%v1))
   write(*,*) " Ey Error = " , maxval(abs(exy_exact%data%v2-exy%data%v2))
   write(*,*) " Po Error = " , maxval(abs(phi_exact%data+phi%data))

   do j = 1, nc_eta2
   do i = 1, nc_eta1
      write(13,*) i, j, exy_exact%data(i,j)%v1 - exy%data(i,j)%v1
   end do
   write(13,*)
   end do
   close(13)
   
   call write_vec1d(phi%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"phi1","mesh",0)
   call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho1","mesh",0)
   call write_vec2d(exy%data%v1,exy%data%v2,mesh%nc_eta1+1,mesh%nc_eta2+1,"exy1","mesh",0)
   
   call delete_poisson_2d_periodic(poisson_e_fields)
   call delete_poisson_2d_periodic(poisson_potential)
   
   call delete_field_2D_vec1( rho )
   call delete_field_2D_vec1( phi )
   call delete_field_2D_vec1( phi_exact )
   call delete_field_2D_vec2( exy )
   call delete_field_2D_vec2( exy_exact )

end program test_poisson_2d
