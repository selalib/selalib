! conservation de la masse
!plot 'diag.dat' u 1:3 w lp title 'BSL', 'diag.dat' u 1:7 w lp title 'BSL NC', 'diag.dat' u 1:11 w lp title 'FSL', 'diag.dat' u 1:15 w lp title 'FSL NC'
! convergence en espace
!plot 'Conv_collela_rot_f3.dat' u 1:2 w lp title 'BSL', 'Conv_collela_rot_f3.dat' u 1:3 w lp title 'BSL NC', 'Conv_collela_rot_f3.dat' u 1:4 w lp title 'FSL', 'Conv_collela_rot_f3.dat' u 1:5 w lp title 'FSL NC'

program test_deposit_cubic_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_splines, only: &
    sll_s_compute_cubic_spline_2d, &
    sll_s_deposit_value_2d, &
    sll_f_interpolate_value_2d, &
    sll_s_init_cubic_spline_2d, &
    sll_t_cubic_spline_2d

  use sll_m_fft, only: &
    sll_s_fft_exec_c2r_1d, &
    sll_s_fft_exec_r2c_1d, &
    sll_s_fft_free, &
    sll_s_fft_init_c2r_1d, &
    sll_s_fft_init_r2c_1d, &
    sll_t_fft

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_t_cubic_spline_2d), pointer :: spl_bsl
type(sll_t_cubic_spline_2d), pointer :: spl_fsl

sll_int32  :: N,Neta1,Neta2,mesh_case,test_case,step,nb_step,visu_step,field_case
sll_int32  :: i,j,bc1_type,bc2_type,err
sll_int32  :: i1,i2,i3
sll_real64 :: eta1,delta_eta1,eta1_min,eta1_max,eta2,delta_eta2,eta2_min,eta2_max
sll_real64 :: x1,x2,x1_min,x2_min,x1_max,x2_max,x1c,x2c,dt
sll_real64 :: T
sll_real64 :: val,val_bsl,val_fsl,tmp1
sll_real64 :: val_spe 
sll_real64 :: a1,a2,eta1c,eta2c
sll_real64,dimension(:,:), pointer :: f
sll_real64,dimension(:,:), pointer :: fh_bsl
sll_real64,dimension(:,:), pointer :: fh_fsl
sll_real64,dimension(:,:), pointer :: fh_spe
sll_real64,dimension(:,:), pointer :: x1_array,x2_array,eta1feet,eta2feet,eta1tot,eta2tot,x1tot,x2tot
sll_real64,dimension(:,:), pointer :: diag
character(len=20) :: conv_name, mass_name
character(len=3)  :: mesh_name, field_name, time_name
  
sll_int32                         :: nc_eta1, nc_eta2
sll_real64, dimension(:), pointer :: d_dx1, d_dx2
sll_real64, dimension(:), pointer :: kx1, kx2
type(sll_t_fft)                   :: fwx1, fwx2
type(sll_t_fft)                   :: bwx1, bwx2
sll_comp64, dimension(:), pointer :: fk1, fk2

sll_int32     :: error
sll_real64    :: kx10
sll_real64    :: kx20

N=128
Neta1 = N
Neta2 = N

! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic
eta1_min = 0._f64
eta1_max = 2._f64*sll_p_pi
eta2_min = 0._f64
eta2_max = 2._f64*sll_p_pi


nc_eta1    = Neta1
nc_eta2    = Neta2

SLL_CLEAR_ALLOCATE(d_dx1(1:nc_eta1), error)
SLL_CLEAR_ALLOCATE(d_dx2(1:nc_eta2), error)
SLL_ALLOCATE(fk1(1:nc_eta1/2+1), error)
fk1 = cmplx(0.0,0.0, kind=f64)
SLL_ALLOCATE(fk2(1:nc_eta2/2+1), error)
fk2 = cmplx(0.0,0.0, kind=f64)

call sll_s_fft_init_r2c_1d(fwx1, nc_eta1, d_dx1, fk1)
call sll_s_fft_init_c2r_1d(bwx1, nc_eta1,   fk1, d_dx1)
call sll_s_fft_init_r2c_1d(fwx2, nc_eta2, d_dx2, fk2)
call sll_s_fft_init_c2r_1d(bwx2, nc_eta2,   fk2, d_dx2)

SLL_CLEAR_ALLOCATE(kx1(1:nc_eta1/2+1), error)
SLL_CLEAR_ALLOCATE(kx2(1:nc_eta2/2+1), error)
 
kx10 = 2._f64*sll_p_pi/(eta1_max-eta1_min)
kx20 = 2._f64*sll_p_pi/(eta2_max-eta2_min)

kx1(1) = 1.0_f64
do i=2,nc_eta1/2+1
   kx1(i) = (i-1)*kx10
end do
kx2(1) = 1.0_f64
do j=2,nc_eta2/2+1
   kx2(j) = (j-1)*kx20
end do

! ---- * Parameters * ----

! --- Space and time parameters ---
! For the python script curvilinear-exe.py
!namelist /param/ N,T,mesh_case,test_case,field_case
!read(*,NML=param);open(unit=900,file="gyrof.param");write(900,NML=param);close(900)


! Final time
T = .1_f64

! -- mesh type --
! 1 : cartesian
mesh_case = 1

! -- distribution function --
! 4 : centered sll_m_gaussian in eta1 and eta2
test_case = 4

! -- advecton field --
! 1 : translation of vector (a1,a2)
! 2 : rotation
! 3 : non homogeneous rotation
! 4 : divergence free complex symmetric field (polar mesh only)
field_case = 1

a1=4._f64
a2=1._f64

! -- visualization parameters --
visu_step = 1
  
! ---- * Construction of the mesh * ----
  
bc1_type = sll_p_periodic
bc2_type = sll_p_periodic
  
eta1c = 0.5_f64*(eta1_max+eta1_min)
eta2c = 0.5_f64*(eta2_max+eta2_min)

! ---- * Time and space steps * ----

! space steps
delta_eta1 = (eta1_max-eta1_min)/real(Neta1,f64)
delta_eta2 = (eta2_max-eta2_min)/real(Neta2,f64)

! time step and number of steps
dt = T*delta_eta1 
nb_step = floor(T/dt)

! ---- * Messages * ----
  
print *,'# N=',N
print *,'# T=',real(nb_step,f64)*dt
print *,'# mesh_case=',mesh_case
print *,'# test_case=',test_case
print *,'# field_case=',field_case
    
! ---- * Allocation and creation of the splines * ----
	
! allocations of the arrays
SLL_ALLOCATE(f(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(fh_bsl(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(fh_fsl(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(fh_spe(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(x1_array(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(x2_array(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta1feet(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta2feet(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta1tot(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta2tot(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(x1tot(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(x2tot(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(diag(10,0:nb_step), err)
	
! creation of the splines
call sll_s_init_cubic_spline_2d(spl_bsl, Neta1+1, Neta2+1, &
  eta1_min, eta1_max, &
  0._f64, 2._f64*sll_p_pi, &
  bc1_type, bc2_type)

call sll_s_init_cubic_spline_2d(spl_fsl, Neta1+1, Neta2+1, &
  eta1_min, eta1_max, &
  0._f64, 2._f64*sll_p_pi, &
  bc1_type, bc2_type)
  
! ---- * Initializations * ----
  
! Analytic distribution function and data for the mesh
open(unit=900,file='f0.dat')
do i=1,Neta1+1
  eta1 = eta1_min + (i-1)*delta_eta1
  do j=1,Neta2+1
    eta2 = eta2_min + (j-1)*delta_eta2
    
    x1_min = eta1_min
    x2_min = eta2_min
    x1_max = eta1_max
    x2_max = eta2_max
    
    x1c = 0.5_f64*(x1_min+x1_max)
    x2c = 0.5_f64*(x2_min+x2_max)
    
    x1_array(i,j) = eta1
    x2_array(i,j) = eta2
    
    f(i,j) = exp(-2_f64*(eta1-eta1c)**2)*exp(-2_f64*(eta2-eta2c)**2)

    write(900,*) eta1, eta2, f(i,j)
  enddo
  write(900,*)
enddo
close(900)

! Distribution functions for the four methods
fh_bsl    = f
fh_fsl    = f
fh_spe    = f
  
! Diagnostic with the weighted mass (temp1) and the "classical" mass (temp2) 
diag = 0._f64

tmp1 = sum(f)*delta_eta1*delta_eta2
  
diag = tmp1  ! analytic solution

! total displacement
eta1tot = x1_array
eta2tot = x2_array

do step=1,nb_step ! ---- * Evolution in time * ----
  
  fh_bsl(:,Neta2+1)    = fh_bsl(:,1)
  fh_fsl(:,Neta2+1)    = fh_fsl(:,1)
      
  call sll_s_compute_cubic_spline_2d(fh_bsl,spl_bsl)
  call sll_s_compute_cubic_spline_2d(fh_fsl,spl_fsl)
    
  do i=1,Neta1+1
    do j=1,Neta2+1
              
      ! ------------ Analytic part -----------------
      
      ! --- Total displacement init ---
      
      eta1 = eta1tot(i,j)
      eta2 = eta2tot(i,j)
      
      x1tot(i,j) = eta1
      x2tot(i,j) = eta2
    
      ! translation
      if (field_case==1) then
        x1 = x1tot(i,j) + 0.01_f64*a1*dt
        x2 = x2tot(i,j) + 0.01_f64*a2*dt
      else if (field_case==2) then                   
        ! rotation
        x1 = cos(dt)*(x1tot(i,j)-x1c)  + sin(dt)*(x2tot(i,j)-x2c) + x1c
        x2 = -sin(dt)*(x1tot(i,j)-x1c) + cos(dt)*(x2tot(i,j)-x2c) + x2c
      else if (field_case==3) then
        ! non homogeneous equation
        x1 = (x1tot(i,j)-x1c)*cos(sqrt(a1*a2)*dt) - (x2tot(i,j)-x2c)*sqrt(a2/a1)*sin(sqrt(a1*a2)*dt) + x1c
        x2 = (x1tot(i,j)-x1c)*sqrt(a1/a2)*sin(sqrt(a1*a2)*dt) + (x2tot(i,j)-x2c)*cos(sqrt(a1*a2)*dt) + x2c
      end if
        
      eta1 = x1 
      eta2 = x2

      ! --- Evaluation of f ---
      f(i,j) = exp(-2_f64*(eta1-eta1c)**2)*exp(-2_f64*(eta2-eta2c)**2)
      
      ! --- Total displacement update ---
      eta1tot(i,j) = eta1
      eta2tot(i,j) = eta2
       
      ! ------------ BSL part -----------------
       
      eta1 = eta1_min + real(i-1,f64)*delta_eta1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
      
      call displacement(dt, eta1, eta2)
      
      call apply_bc()
      
      ! --- Interpolation ---
      fh_bsl(i,j)    = sll_f_interpolate_value_2d(eta1,eta2,spl_bsl)
        
      ! ------------ FSL part -----------------
        
      eta1 = eta1_min + real(i-1,f64)*delta_eta1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
        
      ! --- Displacement ---
      call displacement(-dt, eta1, eta2)
        
      call apply_bc()
        
      eta1feet(i,j) = eta1
      eta2feet(i,j) = eta2
				
    enddo
  enddo

  ! --- Deposition FSL ---
    
  call sll_s_deposit_value_2d(eta1feet,eta2feet,spl_fsl,fh_fsl)

  ! --- Spectral ---

  do j = 1, nc_eta2+1

    d_dx1 = fh_spe(1:nc_eta1,j)
    call sll_s_fft_exec_r2c_1d(fwx1, d_dx1, fk1)
    do i = 2, nc_eta1/2+1
      fk1(i) = fk1(i)*cmplx(cos(kx1(i)*0.01*a1*dt),sin(kx1(i)*0.01*a1*dt),kind=f64)
    end do
    call sll_s_fft_exec_c2r_1d(bwx1, fk1, d_dx1)
  
    fh_spe(1:nc_eta1,j) = d_dx1 / nc_eta1

  end do

  fh_spe(nc_eta1+1,:) = fh_spe(1,:)
                    
  do i = 1, nc_eta1+1

    d_dx2 = fh_spe(i,1:nc_eta2)
    call sll_s_fft_exec_r2c_1d(fwx2, d_dx2, fk2)
    do j = 2, nc_eta2/2+1
      fk2(j) = fk2(j)*cmplx(cos(kx2(j)*0.01*a2*dt),sin(kx2(j)*0.01*a2*dt),kind=f64)
    end do
    call sll_s_fft_exec_c2r_1d(bwx2, fk2, d_dx2)
  
    fh_spe(i,1:nc_eta2) = d_dx2 / nc_eta2

  end do

  fh_spe(:,nc_eta2+1) = fh_spe(:,1)
				
  ! ------------ Diagnostics -----------------
    
  diag(1,step) = sum(f)*delta_eta1*delta_eta2
  diag(2,step) = sum(fh_bsl)*delta_eta1*delta_eta2
  diag(3,step) = sum(fh_fsl)*delta_eta1*delta_eta2
  diag(4,step) = sum(fh_spe)*delta_eta1*delta_eta2

enddo
  
! --- diagnostics ---

! File name

mesh_name = "crt"

SELECT CASE (field_case)
  CASE (1)
    field_name = "trs"
  CASE (2)
    field_name = "rot"
  CASE (3)
    field_name = "rnh"
END SELECT

i1 = int(T)/100
i2 =(int(T)-100*i1)/10
i3 = int(T)-100*i1-10*i2
time_name = char(i1+48)//char(i2+48)//char(i3+48)

conv_name = 'Conv_'//mesh_name//'_'//field_name//'_'//time_name//'.dat'
mass_name = 'Mass_'//mesh_name//'_'//field_name//'_'//time_name//'.dat'
  
val_bsl    = maxval(abs(f-fh_bsl))
val_fsl    = maxval(abs(f-fh_fsl))
val_spe    = maxval(abs(f-fh_spe))
val        = maxval(abs(fh_fsl-fh_bsl))

write(*,*) N,'bsl:',val_bsl,'fsl:',val_fsl,'spe:',val_spe

open(unit=850,file='fh.dat')  
do i=1,Neta1+1
  eta1 = eta1_min+real(i-1,f64)*delta_eta1
  do j=1,Neta2+1
    eta2 = eta2_min + real(j-1,f64)*delta_eta2
    write(850,*) eta1,eta2,x1_array(i,j),x2_array(i,j),f(i,j),fh_bsl(i,j),fh_fsl(i,j), fh_spe(i,j)
  enddo
  write(850,*) ' '
enddo
close(850)

call sll_s_fft_free(fwx1)
call sll_s_fft_free(bwx1)
call sll_s_fft_free(fwx2)
call sll_s_fft_free(bwx2)

contains

subroutine apply_bc()

  ! --- Corrections on the BC ---
  if (bc1_type.eq.sll_p_hermite) eta1 = min(max(eta1,eta1_min),eta1_max)
  if (bc2_type.eq.sll_p_hermite) eta2 = min(max(eta2,eta2_min),eta2_max)

  if (bc1_type==sll_p_periodic) then
    do while (eta1>eta1_max)
      eta1 = eta1-(eta1_max-eta1_min)
    enddo
    do while (eta1<eta1_min)
      eta1 = eta1+(eta1_max-eta1_min)
    enddo
  endif

  if (bc2_type==sll_p_periodic) then
    do while (eta2>eta2_max)
      eta2 = eta2-(eta2_max-eta2_min)
    enddo
    do while (eta2<eta2_min)
      eta2 = eta2+(eta2_max-eta2_min)
    enddo
  endif

end subroutine apply_bc

subroutine displacement(dt, x1, x2)
sll_real64, intent(in)  :: dt
sll_real64, intent(out) :: x1
sll_real64, intent(out) :: x2

  ! --- Displacement ---
  
  ! translation
  if (field_case==1) then
    x1 = x1_array(i,j) + 0.01_f64*a1*dt
    x2 = x2_array(i,j) + 0.01_f64*a2*dt
  else if (field_case==2) then                   
    x1 = cos(dt)*(x1_array(i,j)-x1c)  + sin(dt)*(x2_array(i,j)-x2c) + x1c
    x2 = -sin(dt)*(x1_array(i,j)-x1c) + cos(dt)*(x2_array(i,j)-x2c) + x2c
  else if (field_case==3) then
    x1 = (x1_array(i,j)-x1c)*cos(sqrt(a1*a2)*dt) - (x2_array(i,j)-x2c)*sqrt(a2/a1)*sin(sqrt(a1*a2)*dt) + x1c
    x2 = (x1_array(i,j)-x1c)*sqrt(a1/a2)*sin(sqrt(a1*a2)*dt) + (x2_array(i,j)-x2c)*cos(sqrt(a1*a2)*dt) + x2c
  end if

end subroutine displacement


end program

