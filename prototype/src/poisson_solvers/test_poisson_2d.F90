program test_poisson_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_2d.h"
#include "sll_poisson_solvers.h"

   !-------------------------------------------------------------------
   !  test 2D Poisson solver based on FFT
   !-------------------------------------------------------------------

   use numeric_constants
   use sll_poisson_2D_periodic
   use geometry_functions
   use sll_module_mapped_meshes_2d_base
   use sll_module_mapped_meshes_2d_cartesian
   use sll_scalar_field_2d
   use sll_scalar_field_initializers_base

   implicit none

   sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
   sll_int32   :: nc_eta1, nc_eta2
   sll_int32   :: error

   type (sll_mapped_mesh_2d_cartesian), target :: mesh
   class(sll_mapped_mesh_2d_base), pointer     :: m
   type (scalar_field_2d)                      :: ex
   type (scalar_field_2d)                      :: ey
   type (scalar_field_2d)                      :: ex_exact
   type (scalar_field_2d)                      :: ey_exact
   type (scalar_field_2d)                      :: rho
   type (scalar_field_2d)                      :: phi
   type (scalar_field_2d)                      :: phi_exact
   type(poisson_2d_periodic)                   :: poisson

   sll_real64                         :: x1, x2
   sll_int32                          :: mode
   sll_int32                          :: i, j

   eta1_min = .0_f64; eta1_max = 2.0_f64*sll_pi
   eta2_min = .0_f64; eta2_max = 2.0_f64*sll_pi

   nc_eta1 = 127; nc_eta2 = 127

   call mesh%initialize("mesh", nc_eta1+1, nc_eta2+1)

   call mesh%write_to_file()
   m => mesh

   call initialize_scalar_field_2d( rho, "rho", &
                                    m, NODE_CENTERED_FIELD)

   call initialize_scalar_field_2d( phi, "phi", &
                                    m, NODE_CENTERED_FIELD)

   call initialize_scalar_field_2d( phi_exact, "phi_exact", &
                                    m, NODE_CENTERED_FIELD)

   call initialize_scalar_field_2d( ex, "ex", &
                                    m, NODE_CENTERED_FIELD)

   call initialize_scalar_field_2d( ex_exact, "ex_exact", &
                                    m, NODE_CENTERED_FIELD)

   call initialize_scalar_field_2d( ey, "ey", &
                                    m, NODE_CENTERED_FIELD)

   call initialize_scalar_field_2d( ey_exact, "ey_exact", &
                                    m, NODE_CENTERED_FIELD)

   write(*,*) " eta1_min, eta1_max, nc_eta1 ", eta1_min, eta1_max, nc_eta1
   write(*,*) " eta2_min, eta2_max, nc_eta2 ", eta2_min, eta2_max, nc_eta2

   mode = 2
   do i = 1, nc_eta1+1
      do j = 1, nc_eta2+1
         x1 = mesh%x1_at_node(i,j)*(eta1_max-eta1_min)
         x2 = mesh%x2_at_node(i,j)*(eta2_max-eta2_min)
         phi_exact%data(i,j) = mode * sin(mode*x1) * cos(mode*x2)
         ex_exact%data(i,j)  =  1_f64*mode**2*cos(mode*x1)*cos(mode*x2)
         ey_exact%data(i,j)  = -1_f64*mode**2*sin(mode*x1)*sin(mode*x2)
         rho%data(i,j) = -2_f64 * mode**3 * sin(mode*x1)*cos(mode*x2)
         write(13,*) x1, x2, phi_exact%data(i,j)
      end do
      write(13,*)
   end do
   close(13)

   call poisson%initialize( eta1_min, eta1_max, nc_eta1, &
                            eta2_min, eta2_max, nc_eta2, error) 

   call poisson%solve( phi%data, rho%data)
   write(*,*) " Po Error = " , maxval(abs(phi_exact%data+phi%data))
   call poisson%solve( phi%data, rho%data)
   write(*,*) " Po Error = " , maxval(abs(phi_exact%data+phi%data))

   call poisson%solve( ex%data, ey%data, rho%data)
   write(*,*) " Ex Error = " , maxval(abs(ex_exact%data-ex%data))
   write(*,*) " Ey Error = " , maxval(abs(ey_exact%data-ey%data))

   call poisson%solve( ex%data, ey%data, rho%data)
   write(*,*) " Ex Error = " , maxval(abs(ex_exact%data-ex%data))
   write(*,*) " Ey Error = " , maxval(abs(ey_exact%data-ey%data))

end program test_poisson_2d
