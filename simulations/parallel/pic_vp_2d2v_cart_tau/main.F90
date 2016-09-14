!NKS-two-scale solver for 4d VP
!2nd order EWI scheme
!two-scale E
program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_errors.h"

use m_zone
use m_particules, only: particle, plasma
use sll_m_fft
use sll_m_poisson_2d_base
use sll_m_poisson_2d_periodic
use sll_m_constants
use m_particules_m6, only: calcul_rho_m6, interpol_eb_m6
use mpi_module

implicit none

type(tm_mesh_fields)    :: f
type(particle)          :: p

type(sll_t_fft)         :: fw
type(sll_t_fft)         :: bw
type(sll_t_fft)         :: fft2d
sll_comp64, allocatable :: phi(:,:)
real(8)                 :: xxt(2)
real(8)                 :: epsq
real(8)                 :: dtau

real(8)    , allocatable :: tau   (:)
real(8)    , allocatable :: ltau  (:)
complex(8) , allocatable :: iltau (:)
complex(8) , allocatable :: eiltau(:)
complex(8) , allocatable :: pl    (:)
complex(8) , allocatable :: ql    (:)
complex(8) , allocatable :: gp1   (:,:)
complex(8) , allocatable :: gp2   (:,:)
complex(8) , allocatable :: gm1   (:,:)
complex(8) , allocatable :: gm2   (:,:)
complex(8) , allocatable :: wp1   (:)
complex(8) , allocatable :: wp2   (:)
complex(8) , allocatable :: wm1   (:)
complex(8) , allocatable :: wm2   (:)
complex(8) , allocatable :: xt1   (:,:)
complex(8) , allocatable :: xt2   (:,:)
complex(8) , allocatable :: Et1   (:,:)
complex(8) , allocatable :: Et2   (:,:)
complex(8) , allocatable :: temp1 (:)
complex(8) , allocatable :: temp2 (:)
complex(8) , allocatable :: z     (:)
complex(8) , allocatable :: fex   (:,:,:)
complex(8) , allocatable :: fey   (:,:,:)

complex(8) , allocatable :: up    (:,:,:)
complex(8) , allocatable :: um    (:,:,:)
complex(8) , allocatable :: up0   (:,:,:)
complex(8) , allocatable :: um0   (:,:,:)

complex(8) :: utmp
complex(8) :: vtmp

sll_real64            :: time
sll_real64            :: xmin
sll_real64            :: xmax
sll_real64            :: ymin
sll_real64            :: ymax
sll_int32             :: istep=1
sll_int32             :: iargc
sll_int32             :: n,m
sll_int32             :: i
sll_int32             :: j
sll_int32             :: error

sll_real64            :: aux1, aux2
sll_real64            :: s, dum

real       :: start_time, stop_time
integer    :: dat_file_id, ref_file_id
logical    :: file_exists
complex(8) :: cost, sint

character(len=272)    :: argv
class(sll_c_poisson_2d_base), pointer :: poisson

logical :: master = .false.
integer :: prank, psize, iproc
integer :: l1, l2, ll
complex(8), allocatable :: et1_loc(:,:), et2_loc(:,:)
complex(8), allocatable :: fex_loc(:,:,:), fey_loc(:,:,:)

call init_mpi(prank, psize)
if (prank == 0 ) master = .true.

if (master) then

  n = iargc()
  if (n == 0) stop 'Usage: ./bin/test_pic2d fichier-de-donnees.nml'
  do i = 1, n
    call getarg( i, argv)
    write(*,'(i2, 1x, a)') i, argv
  end do

  call readin( trim(argv) )
  call cpu_time(start_time)

end if

call mpi_global_master() !Brodcast global values

SLL_ALLOCATE(tau(0:ntau-1),           error)
SLL_ALLOCATE(ltau(0:ntau-1),          error)
SLL_ALLOCATE(iltau(0:ntau-1),         error)
SLL_ALLOCATE(eiltau(0:ntau-1),        error)
SLL_ALLOCATE(pl(0:ntau-1),            error)
SLL_ALLOCATE(ql(0:ntau-1),            error)
SLL_ALLOCATE(wp1(nbpart),             error)
SLL_ALLOCATE(wp2(nbpart),             error)
SLL_ALLOCATE(wm1(nbpart),             error)
SLL_ALLOCATE(wm2(nbpart),             error)
SLL_ALLOCATE(Et1(nbpart,0:ntau-1),    error)
SLL_ALLOCATE(Et2(nbpart,0:ntau-1),    error)

SLL_ALLOCATE(gp1(0:ntau-1,nbpart),    error)
SLL_ALLOCATE(gp2(0:ntau-1,nbpart),    error)
SLL_ALLOCATE(gm1(0:ntau-1,nbpart),    error)
SLL_ALLOCATE(gm2(0:ntau-1,nbpart),    error)
SLL_ALLOCATE(up(0:ntau-1,nbpart,2),   error)
SLL_ALLOCATE(um(0:ntau-1,nbpart,2),   error)
SLL_ALLOCATE(up0(0:ntau-1,nbpart,2),  error)
SLL_ALLOCATE(um0(0:ntau-1,nbpart,2),  error)
SLL_ALLOCATE(xt1(0:ntau-1,2),         error)
SLL_ALLOCATE(xt2(0:ntau-1,2),         error)

SLL_ALLOCATE(z(2),                    error)
SLL_ALLOCATE(fex(0:nx,0:ny,0:ntau-1), error)
SLL_ALLOCATE(fey(0:nx,0:ny,0:ntau-1), error)

SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny),   error)
SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny),   error)
SLL_CLEAR_ALLOCATE(f%bz(0:nx,0:ny),   error)
SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny),   error)

time   = 0._f64

epsq   = ep * ep
dtau   = 2._f64*sll_p_pi/ntau

m      = ntau/2
ltau   =(/ (n, n=0,m-1), (n, n=-m,-1 )/)
iltau  = 0.5_f64 * sll_p_i1 * cmplx(ltau,0.0_f64) / epsq
eiltau = exp(-iltau*dt) /ntau

pl(0)=dt
ql(0)=dt**2/2._f64
do i=1,ntau-1
  pl(i)=2._f64*epsq*sll_p_i1*(exp(-iltau(i)*dt)-1.0_f64)/ltau(i)
  ql(i)=2._f64*epsq*(2.0_f64*epsq*(1.0_f64-exp(-iltau(i)*dt)) &
                    -2.0_f64*epsq*iltau(i)*dt)/ltau(i)**2
enddo
do i=0,ntau-1
  tau(i) =i*dtau
enddo

do i=0,nx
  aux1 = alpha/kx * sin(kx*i*dx)
  aux2 = alpha * cos(kx*i*dx)
  do j=0,ny
    f%ex(i,j) = aux1
    f%r0(i,j) = aux2+sin(ky*j*dy)
    f%ey(i,j) = -cos(ky*j*dy)/ky!0._f64!
  enddo
enddo
      
xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p )

if (master) then
  print"('ep = ', g15.3)", ep
  print"('nbpart = ', g15.3)", nbpart
  print"('dt = ', g15.3)", dt
  print"('kx,ky = ', 2g15.3)", kx, ky
  print"('alpha = ', g15.3)", alpha
  print"('ntau = ', g15.3)", ntau
endif

poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)

wp1(:) = cmplx(  2._f64*((p%dpx+p%idx)*dx+ep*p%vpy),0.0,f64)
wp2(:) = cmplx(  2._f64*((p%dpy+p%idy)*dy-ep*p%vpx),0.0,f64)
wm1(:) = cmplx(- 2._f64*ep*p%vpy,0.0,f64)
wm2(:) = cmplx(  2._f64*ep*p%vpx,0.0,f64)

if (mod(ntau,psize) == 0) then

  l1 = prank*ntau/psize
  l2 = (prank+1)*ntau/psize-1
  ll = l2-l1+1

else
  
  if (master) then
    print*, "Remainder of the division of ntau by Nprocs must be zero." 
  end if

  call finish_mpi()
  stop

end if

print*, prank, psize, l1, l2, ll

allocate(et1_loc(nbpart,l1:l2))
allocate(et2_loc(nbpart,l1:l2))
allocate(fex_loc(0:nx,0:ny,l1:l2))
allocate(fey_loc(0:nx,0:ny,l1:l2))

do n = l1, l2

  cost = cmplx(cos(tau(n)), 0.0, f64)
  sint = cmplx(sin(tau(n)), 0.0, f64)
  do m=1,nbpart
    utmp   = 0.5_f64*(cost*wp1(m)-sint*wp2(m)+cost*wm1(m)+sint*wm2(m))
    vtmp   = 0.5_f64*(sint*wp1(m)+cost*wp2(m)-sint*wm1(m)+cost*wm2(m))
    xxt(1) = real( cost*utmp+sint*vtmp)
    xxt(2) = real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  call interpol_eb_m6( f, p )
  et1_loc(:,n) = cmplx(p%epx,0.0,f64) !g(0,tau,w(0))
  et2_loc(:,n) = cmplx(p%epy,0.0,f64)

enddo

call MPI_ALLGATHER(et1_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et1,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)
call MPI_ALLGATHER(et2_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et2,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)



SLL_ALLOCATE(temp1(0:ntau-1),         error)
SLL_ALLOCATE(temp2(0:ntau-1),         error)
call sll_s_fft_init_c2c_1d(fw,ntau,temp1,temp1,sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(bw,ntau,temp1,temp1,sll_p_fft_backward)

SLL_ALLOCATE(phi(nx/2+1,ny),          error)
call sll_s_fft_init_r2c_2d(fft2d,nx,ny,f%r0(0:nx-1,0:ny-1),phi)

do m=1,nbpart

  temp1= 2._f64*Et2(m,:)
  temp2=-2._f64*Et1(m,:)!g_+
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)

  do n=1,ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/ntau
  enddo

  temp1(0)= sll_p_i0
  temp2(0)= sll_p_i0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF+
  up(:,m,1)=wp1(m)+2._f64*epsq*(temp1-temp1(0))
  up(:,m,2)=wp2(m)+2._f64*epsq*(temp2-temp2(0))!1st ini data of U_+
  !---
  do n=0,ntau-1
    cost = cmplx(cos(2.0_f64*tau(n)), 0.0, f64) ! not the same dtau
    sint = cmplx(sin(2.0_f64*tau(n)), 0.0, f64)
    temp1(n)=-2._f64*( sint*Et1(m,n)+cost*Et2(m,n))
    temp2(n)= 2._f64*(-sint*Et2(m,n)+cost*Et1(m,n))!g_-
  enddo
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  do n=1,ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/ntau
  enddo
  temp1(0)=sll_p_i0
  temp2(0)=sll_p_i0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF-
  um(:,m,1)=wm1(m)+2._f64*epsq*(temp1-temp1(0))
  um(:,m,2)=wm2(m)+2._f64*epsq*(temp2-temp2(0))!1st ini data of U_-
enddo

!--corrected more initial data

do n=l1,l2

  cost = cmplx(cos(tau(n)), 0.0, f64)
  sint = cmplx(sin(tau(n)), 0.0, f64)

  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                   +cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                   -sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo

  call calcul_rho_m6( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
  fex_loc(:,:,n)= cmplx(f%ex,0.0,f64)
  fey_loc(:,:,n)= cmplx(f%ey,0.0,f64)!E_1st(0,x)

enddo

call MPI_ALLGATHER(fex_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fex,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)
call MPI_ALLGATHER(fey_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fey,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)

do n=l1,l2
  cost = cmplx(cos(tau(n)), 0.0, f64)
  sint = cmplx(sin(tau(n)), 0.0, f64)
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                   +cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                   -sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  f%ex= real(fex(:,:,n))
  f%ey= real(fey(:,:,n))
  call interpol_eb_m6( f, p )
  Et1_loc(:,n)= cmplx(p%epx,0.0,f64) !g_1st(0,tau,U_1st(0))
  Et2_loc(:,n)= cmplx(p%epy,0.0,f64)
enddo

call MPI_ALLGATHER(et1_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et1,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)
call MPI_ALLGATHER(et2_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et2,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)

do m=1,nbpart
  temp1= 2._f64*Et2(m,:)
  temp2=-2._f64*Et1(m,:)!g_+
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  up0(0,m,1)=temp1(0)/ntau!Pi g_+
  up0(0,m,2)=temp2(0)/ntau!Pi g_+
  do n=1,ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/ntau
  enddo
  temp1(0)=sll_p_i0
  temp2(0)=sll_p_i0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF+
  up(:,m,1)=wp1(m)+2._f64*epsq*(temp1-temp1(0))
  up(:,m,2)=wp2(m)+2._f64*epsq*(temp2-temp2(0))!3rd ini data of U_+
  !---
  do n=0,ntau-1
    cost = cmplx(cos(2_f64*tau(n)), 0.0, f64)
    sint = cmplx(sin(2_f64*tau(n)), 0.0, f64)
    temp1(n)=-2._f64*(sint*Et1(m,n)+cost*Et2(m,n))
    temp2(n)=2._f64*(-sint*Et2(m,n)+cost*Et1(m,n))!g_-
  enddo
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  um0(0,m,1)=temp1(0)/ntau!Pi g_-
  um0(0,m,2)=temp2(0)/ntau!Pi g_-
  do n=1,ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/ntau
  enddo
  temp1(0)=sll_p_i0
  temp2(0)=sll_p_i0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF-
  um(:,m,1)=wm1(m)+2._f64*epsq*(temp1-temp1(0))
  um(:,m,2)=wm2(m)+2._f64*epsq*(temp2-temp2(0))!3rd ini data of U_-

enddo

do n=l1,l2
  cost = cmplx(cos(tau(n)),0.0, f64)
  sint = cmplx(sin(tau(n)),0.0, f64)
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                   +cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                   -sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  call calcul_rho_m6( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
  fex_loc(:,:,n)= cmplx(f%ex,0.0,f64)
  fey_loc(:,:,n)= cmplx(f%ey,0.0,f64)!E_4(0,x)
enddo

call MPI_ALLGATHER(fex_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fex,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)
call MPI_ALLGATHER(fey_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fey,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)

!--time iteration---
time=dt
do n=l1,l2
  cost = cmplx(cos(tau(n)), 0.0, f64)
  sint = cmplx(sin(tau(n)), 0.0, f64)
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                   +cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                   -sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  f%ex= real(fex(:,:,n))
  f%ey= real(fey(:,:,n))
  call interpol_eb_m6( f, p )
  Et1_loc(:,n)= cmplx(p%epx,0.0,f64) !g_3rd(0,tau,U_3rd(0))
  Et2_loc(:,n)= cmplx(p%epy,0.0,f64)
enddo

call MPI_ALLGATHER(et1_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et1,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)
call MPI_ALLGATHER(et2_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et2,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)

do m=1,nbpart
  temp1= 2._f64*Et2(m,:)
  temp2=-2._f64*Et1(m,:)
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  xt1(:,1)=temp1/ntau!g_+tilde(t=0)
  xt1(:,2)=temp2/ntau!g_+tilde(t=0)
  !---
  do n=0,ntau-1
    cost = cmplx(cos(2_f64*tau(n)),0.0, f64)
    sint = cmplx(sin(2_f64*tau(n)),0.0, f64)
    temp1(n)=-2._f64*(sint*Et1(m,n)+cost*Et2(m,n))
    temp2(n)=2._f64*(-sint*Et2(m,n)+cost*Et1(m,n))
  enddo
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  xt2(:,1)=temp1/ntau!g_-tilde(t=0)
  xt2(:,2)=temp2/ntau!g_-tilde(t=0)

  call sll_s_fft_exec_c2c_1d(fw, up(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,2), temp2)
  do n=0,ntau-1
    temp1(n)=eiltau(n)*temp1(n)+pl(n)*xt1(n,1)!utilde_+^1,predict
    temp2(n)=eiltau(n)*temp2(n)+pl(n)*xt1(n,2)!utilde_+^1,predict
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, up0(:,m,1))!u_+(t1),predict
  call sll_s_fft_exec_c2c_1d(bw, temp2, up0(:,m,2))
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,2), temp2)
  do n=0,ntau-1
    temp1(n)=eiltau(n)*temp1(n)+pl(n)*xt2(n,1)!utilde_-^1,predict
    temp2(n)=eiltau(n)*temp2(n)+pl(n)*xt2(n,2)!utilde_-^1,predict
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, um0(:,m,1))!u_-(t1),predict
  call sll_s_fft_exec_c2c_1d(bw, temp2, um0(:,m,2))

  gp1(:,m) = xt1(:,1)
  gp2(:,m) = xt1(:,2)
  gm1(:,m) = xt2(:,1)
  gm2(:,m) = xt2(:,2)

enddo


do n= l1,l2
  cost = cmplx(cos(tau(n)),0.0,f64)
  sint = cmplx(sin(tau(n)),0.0,f64)
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up0(n,m,1)-sint*up0(n,m,2) &
                   +cost*um0(n,m,1)+sint*um0(n,m,2))
    vtmp = 0.5_f64*(sint*up0(n,m,1)+cost*up0(n,m,2) &
                   -sint*um0(n,m,1)+cost*um0(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  call calcul_rho_m6( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
  fex_loc(:,:,n)=cmplx(f%ex,0.0,f64)
  fey_loc(:,:,n)=cmplx(f%ey,0.0,f64)!prediction
enddo

call MPI_ALLGATHER(fex_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fex,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)
call MPI_ALLGATHER(fey_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fey,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)
!--correction--
do n=l1,l2
  cost = cmplx(cos(tau(n)),0.0,f64)
  sint = cmplx(sin(tau(n)),0.0,f64)
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up0(n,m,1)-sint*up0(n,m,2) &
                   +cost*um0(n,m,1)+sint*um0(n,m,2))
    vtmp = 0.5_f64*(sint*up0(n,m,1)+cost*up0(n,m,2) &
                   -sint*um0(n,m,1)+cost*um0(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  f%ex= real(fex(:,:,n))
  f%ey= real(fey(:,:,n))
  call interpol_eb_m6( f, p )
  Et1_loc(:,n)= cmplx(p%epx,0.0,f64) !g(t1,tau,U(t1))
  Et2_loc(:,n)= cmplx(p%epy,0.0,f64)
enddo

call MPI_ALLGATHER(et1_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et1,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)
call MPI_ALLGATHER(et2_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                   Et2,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)

do m=1,nbpart

  temp1 =  2._f64*Et2(m,:)
  temp2 = -2._f64*Et1(m,:)
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  xt1(:,1)=temp1/ntau!g_+tilde(t1) predict
  xt1(:,2)=temp2/ntau!g_+tilde(t1) predict
  !---
  do n=0,ntau-1
    cost = cmplx(cos(2_f64*tau(n)),0.0,f64)
    sint = cmplx(sin(2_f64*tau(n)),0.0,f64)
    temp1(n) = - 2._f64*( sint*Et1(m,n)+cost*Et2(m,n))
    temp2(n) =   2._f64*(-sint*Et2(m,n)+cost*Et1(m,n))
  enddo

  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)

  xt2(:,1)=temp1/ntau!g_-tilde(t1) predict
  xt2(:,2)=temp2/ntau!g_-tilde(t1) predict

  call sll_s_fft_exec_c2c_1d(fw, up(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,2), temp2)

  do n=0,ntau-1
    temp1(n)=eiltau(n)*temp1(n)+pl(n)*xt1(n,1)+ql(n)*(xt1(n,1)-gp1(n,m))/dt
    temp2(n)=eiltau(n)*temp2(n)+pl(n)*xt1(n,2)+ql(n)*(xt1(n,2)-gp2(n,m))/dt
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, up(:,m,1))!u_+(t1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, up(:,m,2))
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,2), temp2)
  do n=0,ntau-1
    temp1(n)=eiltau(n)*temp1(n)+pl(n)*xt2(n,1)+ql(n)*(xt2(n,1)-gm1(n,m))/dt
    temp2(n)=eiltau(n)*temp2(n)+pl(n)*xt2(n,2)+ql(n)*(xt2(n,2)-gm2(n,m))/dt
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, um(:,m,1))!u_-(t1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, um(:,m,2))

enddo

do n=l1,l2
  cost = cmplx(cos(tau(n)),0.0,f64)
  sint = cmplx(sin(tau(n)),0.0,f64)
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                   +cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                   -sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  call calcul_rho_m6( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
  fex_loc(:,:,n)= cmplx(f%ex,0.0,f64)
  fey_loc(:,:,n)= cmplx(f%ey,0.0,f64)
enddo

call MPI_ALLGATHER(fex_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fex,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)
call MPI_ALLGATHER(fey_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                   fey,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                   MPI_COMM_WORLD,code)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Loop over time ***
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call energy_use()

do istep = 2, nstep

  do n=l1,l2
    cost = cmplx(cos(tau(n)),0.0,f64)
    sint = cmplx(sin(tau(n)),0.0,f64)
    do m=1,nbpart
      utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                     +cost*um(n,m,1)+sint*um(n,m,2))
      vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                     -sint*um(n,m,1)+cost*um(n,m,2))
      xxt(1)=real( cost*utmp+sint*vtmp)
      xxt(2)=real(-sint*utmp+cost*vtmp)
      call apply_bc()
      p%idx(m) = floor(xxt(1)/dimx*nx)
      p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
      p%idy(m) = floor(xxt(2)/dimy*ny)
      p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    f%ex = real(fex(:,:,n))
    f%ey = real(fey(:,:,n))
    call interpol_eb_m6( f, p )
    Et1_loc(:,n)= cmplx(p%epx, 0.0, f64) 
    Et2_loc(:,n)= cmplx(p%epy, 0.0, f64)
  enddo

  call MPI_ALLGATHER(et1_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                     Et1,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)
  call MPI_ALLGATHER(et2_loc,nbpart*ll,MPI_DOUBLE_COMPLEX, &
                     Et2,nbpart*ll,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,code)

  do m=1,nbpart
    temp1=  2._f64*Et2(m,:)
    temp2= -2._f64*Et1(m,:)
    call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
    call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
    xt1(:,1)=temp1/ntau
    xt1(:,2)=temp2/ntau
    !---
    do n=0,ntau-1
      cost = cmplx(cos(2_f64*tau(n)),0.0,f64)
      sint = cmplx(sin(2_f64*tau(n)),0.0,f64)
      temp1(n) = -2._f64*( sint*Et1(m,n)+cost*Et2(m,n))
      temp2(n) =  2._f64*(-sint*Et2(m,n)+cost*Et1(m,n))
    enddo
    call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
    call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
    xt2(:,1)=temp1/ntau
    xt2(:,2)=temp2/ntau

    call sll_s_fft_exec_c2c_1d(fw, up(:,m,1), temp1)
    call sll_s_fft_exec_c2c_1d(fw, up(:,m,2), temp2)
    do n=0,ntau-1
      temp1(n)= eiltau(n)*temp1(n)+pl(n)*xt1(n,1) &
                +ql(n)*(xt1(n,1)-gp1(n,m))/dt
      temp2(n)= eiltau(n)*temp2(n)+pl(n)*xt1(n,2) &
                +ql(n)*(xt1(n,2)-gp2(n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(bw, temp1, up(:,m,1))
    call sll_s_fft_exec_c2c_1d(bw, temp2, up(:,m,2))
    call sll_s_fft_exec_c2c_1d(fw, um(:,m,1), temp1)
    call sll_s_fft_exec_c2c_1d(fw, um(:,m,2), temp2)
    do n=0,ntau-1
      temp1(n)=eiltau(n)*temp1(n)+pl(n)*xt2(n,1) &
              +ql(n)*(xt2(n,1)-gm1(n,m))/dt
      temp2(n)=eiltau(n)*temp2(n)+pl(n)*xt2(n,2) &
              +ql(n)*(xt2(n,2)-gm2(n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(bw, temp1, um(:,m,1))
    call sll_s_fft_exec_c2c_1d(bw, temp2, um(:,m,2))

    gp1(:,m)=xt1(:,1)
    gp2(:,m)=xt1(:,2)
    gm1(:,m)=xt2(:,1)
    gm2(:,m)=xt2(:,2)

  enddo

  !--update E--
  time = time + dt


  do n=l1,l2
    cost = cmplx(cos(tau(n)),0.0,f64)
    sint = cmplx(sin(tau(n)),0.0,f64)
    do m=1,nbpart
      utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2) &
                     +cost*um(n,m,1)+sint*um(n,m,2))
      vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2) &
                     -sint*um(n,m,1)+cost*um(n,m,2))
      xxt(1)=real( cost*utmp+sint*vtmp)
      xxt(2)=real(-sint*utmp+cost*vtmp)
      call apply_bc()
      p%idx(m) = floor(xxt(1)/dimx*nx)
      p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
      p%idy(m) = floor(xxt(2)/dimy*ny)
      p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    call calcul_rho_m6( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
    fex_loc(:,:,n) = cmplx(f%ex,0.0,f64)
    fey_loc(:,:,n) = cmplx(f%ey,0.0,f64)
  enddo

  
  call MPI_ALLGATHER(fex_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                     fex,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                     MPI_COMM_WORLD,code)
  call MPI_ALLGATHER(fey_loc,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX, &
                     fey,(nx+1)*(ny+1)*ll,MPI_DOUBLE_COMPLEX,     &
                     MPI_COMM_WORLD,code)


  do iproc=0, psize-1
    if (iproc == prank) then
      print *, ' Rank ', iproc, istep, time
    end if
    call MPI_Barrier(MPI_COMM_WORLD, code)
  enddo

  if (master) call energy_use()

enddo

call MPI_BARRIER(MPI_COMM_WORLD,code)

if (master) then

  call cpu_time(stop_time)
  
  s = 0.0_f64
  inquire( file= 'fh64.ref', exist=file_exists )
  if (file_exists) then
    open(newunit=ref_file_id,file='fh64.ref')
    do i=1,nx
      do j=1,ny
        read(ref_file_id,*) dum
        s = s + abs(dum-f%r0(i,j))
      enddo
      read(ref_file_id,*)
    enddo
    print *, "CPU time:", stop_time - start_time, "seconds, error = ", s
  else
    open(newunit=dat_file_id,file='fh64.dat')
    do i=1,nx
      do j=1,ny
        write(dat_file_id,*)f%r0(i,j)
      enddo
      write(dat_file_id,*)
    enddo
    close(851)
    print *, "CPU time:", stop_time - start_time, "seconds"
  endif
  
end if

deallocate(fex_loc)
deallocate(fey_loc)
deallocate(et1_loc)
deallocate(et2_loc)
call sll_s_fft_free(fw)
call sll_s_fft_free(bw)

call finish_mpi()

contains

subroutine apply_bc()
  do while ( xxt(1) > xmax )
    xxt(1) = xxt(1) - dimx
  enddo
  do while ( xxt(1) < xmin )
    xxt(1)= xxt(1) + dimx
  enddo
  do while ( xxt(2) > ymax )
    xxt(2)  = xxt(2)  - dimy
  enddo
  do while ( xxt(2)  < ymin )
    xxt(2) = xxt(2)  + dimy
  enddo
end subroutine apply_bc

subroutine energy_use()
sll_comp64 :: vp(2), vm(2), z(2), wp(2), wm(2)

cost = cmplx(cos(0.5_f64*time/epsq),0.0,f64)
sint = cmplx(sin(0.5_f64*time/epsq),0.0,f64)

do m=1,nbpart

  call sll_s_fft_exec_c2c_1d(fw, up(:,m,1),temp1)
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,2),temp2)
  wp = sll_p_i0
  do n=0,ntau-1
    wp(1)=wp(1)+temp1(n)/ntau*exp(iltau(n)*time)
    wp(2)=wp(2)+temp2(n)/ntau*exp(iltau(n)*time)
  enddo
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,1),temp1)
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,2),temp2)
  wm = sll_p_i0
  do n=0,ntau-1
    wm(1)=wm(1)+temp1(n)/ntau*exp(iltau(n)*time)
    wm(2)=wm(2)+temp2(n)/ntau*exp(iltau(n)*time)
  enddo
  vp(1)=cost*wp(1)-sint*wp(2)
  vp(2)=cost*wp(2)+sint*wp(1)
  vm(1)=cost*wm(1)+sint*wm(2)
  vm(2)=cost*wm(2)-sint*wm(1)
  z  =  0.5_f64*(vp+vm)
  xxt(1)=real(cost*z(1)+sint*z(2))
  xxt(2)=real(cost*z(2)-sint*z(1))
  call apply_bc()
  p%idx(m) = floor(xxt(1)/dimx*nx)
  p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
  p%idy(m) = floor(xxt(2)/dimy*ny)
  p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  p%vpx(m) =  (sint*vm(1)+cost*vm(2))/2.0d0/ep
  p%vpy(m) =-(-sint*vm(2)+cost*vm(1))/2.0d0/ep

enddo

call calcul_rho_m6( p, f )
call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)

open(10, file='energy.dat', position='append')
if (istep==1) rewind(10)
write(10,"(3g15.7)") time, 0.5*log(sum(f%ex*f%ex)), 0.5*log(sum(f%ey*f%ey))
!write(10,"(3g15.7)") time, energy_fourier_mode_xy(1,1,f%ex,f%ey), &
!                     energy_fourier_mode_xy(1,1,real(fex(:,:,0)),real(fey(:,:,0))) 

close(10)


!energy_p = 0.5_f64*sum(p%p*(p%vpx*p%vpx+p%vpy*p%vpy))
!energy_e = 0.5_f64*sum(f%ex*f%ex+f%ey*f%ey)
!call sll_s_fft_exec_r2c_2d(fft2d,f%r0(0:nx-1,0:ny-1),phi)
!xi       = real(sum(phi*conjg(phi))*dx*dy)

end subroutine energy_use

!---------------------------------------------------------
!> @brief Computes a Fourier mode of the electric field.
!> @details
!> @param[in] mode_x the mode to be computed.
!> @param[in] mode_y the mode to be computed.
!> @param[in] ex the electric field on the x-axis.
!> @param[in] ey the electric field on the y-axis.
!> @param[in] d the direction in which to compute the Fourier mode.
!> @return    the Fourier mode of the electric field.
function fourier_mode_xy(mode_x, mode_y, ex, ey, d)

sll_int32,  intent(in) :: mode_x
sll_int32,  intent(in) :: mode_y
sll_real64, intent(in) :: ex(:,:)
sll_real64, intent(in) :: ey(:,:)
sll_real64, intent(in) :: d(1:2)
sll_real64             :: fourier_mode_xy

sll_int32  :: i, j
sll_real64 :: term1, term2
sll_real64 :: sqr_norm_d

sqr_norm_d = sqrt(real(mode_x*mode_x+mode_y*mode_y))
if (sqr_norm_d < epsilon(sqr_norm_d)) then
  SLL_ERROR("pic_vp_2d2v_cart_tau", 'd is not a proper direction.')
end if

term1 = 0._f64
do j = 1, ny
   do i = 1, nx
      term1 = term1 + (d(1) * ex(i, j) + d(2) * ey(i, j)) *         &
        cos(real(mode_x * (i-1),f64) * sll_p_twopi / real(nx,f64) + &
            real(mode_y * (j-1),f64) * sll_p_twopi / real(ny,f64))
   enddo
enddo
term1 = term1**2
term2 = 0._f64
do j = 1, ny
   do i = 1, nx
      term2 = term2 + (d(1) * ex(i, j) + d(2) * ey(i, j)) *         &
        sin(real(mode_x * (i-1),f64) * sll_p_twopi / real(nx,f64) + &
            real(mode_y * (j-1),f64) * sll_p_twopi / real(ny,f64))
   enddo
enddo
term2 = term2**2
fourier_mode_xy = sqrt(2._f64 * dx * dy * (term1 + term2) / &
    (sqr_norm_d * real(nx * ny, f64)))

end function fourier_mode_xy

!---------------------------------------------------------
!> @brief Computes the energy of the Fourier mode of the electric field.
!> @details
!> @param[in] mode_x the mode to be computed.
!> @param[in] mode_y the mode to be computed.
!> @param[in] ex the electric field on the x-axis.
!> @param[in] ey the electric field on the y-axis.
!> @return    the energy of the Fourier mode of the electric field.
function energy_fourier_mode_xy( mode_x, mode_y, ex, ey)

sll_int32,  intent(in) :: mode_x
sll_int32,  intent(in) :: mode_y
sll_real64, intent(in) :: ex(:,:)
sll_real64, intent(in) :: ey(:,:)

sll_real64 :: energy_fourier_mode_xy

energy_fourier_mode_xy = sqrt( &
  fourier_mode_xy(mode_x, mode_y, ex, ey, (/ 1._f64, 0._f64 /)**2) + &
  fourier_mode_xy(mode_x, mode_y, ex, ey, (/ 0._f64, 1._f64 /)**2))

end function energy_fourier_mode_xy

end program test_pic2d

