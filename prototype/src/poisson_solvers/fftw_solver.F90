program fft_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
use numeric_constants

use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

type(C_PTR) :: fw, bw
sll_int32 :: nx, ny
sll_int32 :: error
sll_int32 :: nc_eta1, nc_eta2
sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
type(geometry_2D),         pointer :: geom
type(mesh_descriptor_2D),  pointer :: mesh
type(field_2D_vec1),       pointer :: phi, phi_exact
type(field_2D_vec1),       pointer :: rho
sll_real64                         :: x1, x2
sll_int32                          :: mode
sll_int32                          :: i, j
sll_real64                         :: dx,dy
sll_real64                         :: kx0, kx
sll_real64                         :: ky0, ky
sll_int32                          :: ik, jk

sll_real64, dimension(:,:), allocatable :: kmod
sll_comp64, dimension(:,:), allocatable :: rhot
!complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: rhot
!type(C_PTR) :: p

eta1_min = .0_f64; eta1_max = 2.0_f64*sll_pi
eta2_min = .0_f64; eta2_max = 2.0_f64*sll_pi

geom => new_geometry_2D ('cartesian')

nc_eta1 = 256; nc_eta2 = 256

mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
        PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

call write_mesh_2D(mesh)

rho       => new_field_2D_vec1(mesh)
phi       => new_field_2D_vec1(mesh)
phi_exact => new_field_2D_vec1(mesh)

mode = 2
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      x1 = eta1_min+(i-1)*mesh%delta_eta1
      x2 = eta2_min+(j-1)*mesh%delta_eta2
      phi%data(i,j) = mode * sin(mode*x1) * cos(mode*x2)
      rho%data(i,j) = 2_f64 * mode**3 * sin(mode*x1)*cos(mode*x2)
   end do
end do

call write_vec1d(phi%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"phi0","mesh",0)
call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho0","mesh",0)

FIELD_DATA(phi_exact) = FIELD_DATA(phi)
FIELD_DATA(phi) = 0.0

nx = GET_FIELD_NC_ETA1(phi)
ny = GET_FIELD_NC_ETA2(phi)

dx = GET_FIELD_DELTA_ETA1(phi)
dy = GET_FIELD_DELTA_ETA2(phi)
write(*,*) "dx, dy = ", dx, dy

!p = fftw_alloc_complex(int((nx/2+1)*ny,C_SIZE_T))
SLL_ALLOCATE(rhot(nx/2+1,ny), error)
!call c_f_pointer(p, rhot, [nx/2+1,ny])

fw = fftw_plan_dft_r2c_2d(ny, nx, rho%data(1:nx,1:ny), rhot, FFTW_ESTIMATE);
bw = fftw_plan_dft_c2r_2d(ny, nx, rhot, phi%data(1:nx,1:ny), FFTW_ESTIMATE)

call fftw_execute_dft_r2c(fw, rho%data(1:nx,1:ny), rhot)

kx0=2._f64*sll_pi/(eta1_max-eta1_min)
ky0=2._f64*sll_pi/(eta2_max-eta2_min)

SLL_ALLOCATE(kmod(nx/2+1,ny), error)

do ik=1,nx/2+1
   kx  = (ik-1)*kx0
   do jk = 1, ny/2
      ky  = (jk-1)*ky0
      kmod(ik,jk) = kx*kx+ky*ky
   end do
   do jk = ny/2+1,ny     
      ky= (jk-1-ny)*ky0
      kmod(ik,jk) = kx*kx+ky*ky
   end do
end do
kmod(1,1) = 1.0_f64

rhot = rhot / kmod

call fftw_execute_dft_c2r(bw,rhot,phi%data(1:nx,1:ny))

!Normalize
phi%data = phi%data / (nx*ny)

phi%data(nx+1,:) = phi%data(1,:)
phi%data(:,ny+1) = phi%data(:,1)

call write_vec1d(phi%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"phi1","mesh",0)
call write_vec1d(rho%data,mesh%nc_eta1+1,mesh%nc_eta2+1,"rho1","mesh",0)

write(*,*) " E = ", maxval(phi_exact%data-phi%data)/maxval(phi_exact%data)

call delete_field_2D_vec1( rho )
call delete_field_2D_vec1( phi )
call delete_field_2D_vec1( phi_exact )

call fftw_destroy_plan(fw)
call fftw_destroy_plan(bw)
!call fftw_free(p)

end program fft_solver
