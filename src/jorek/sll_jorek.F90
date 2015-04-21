module jorek_model

use typedef
use spm_def
use indices_def
use jorek_param_def
use mesh_def
use blackbox_def
use greenbox_def
use quadratures_def 
use linear_solver_def
use basis_def
use febasis
use febasis2d_bezier
use space_def
use field_def
use matrix_def
use fem_def

!  use febasis2d_splines
!  use febasis2d_boxsplines

implicit none

! ... quad mesh
type(def_fem_quad_2d),    target  :: fem_model 
type(def_space_quad_2d),  target  :: space_trial
class(def_space_quad_2d), pointer :: ptr_space_trial
type(def_space_quad_2d),  target  :: space_test
class(def_space_quad_2d), pointer :: ptr_space_test
! ...

!  ! ... triangle mesh
!  type(def_fem_triangle_2d)  , target :: fem_model 
!
!  type(def_space_triangle_2d), target :: space_trial
!  class(def_space_triangle_2d), pointer :: ptr_space_trial
!
!  type(def_space_triangle_2d), target :: space_test
!  class(def_space_triangle_2d), pointer :: ptr_space_test
!  ! ...

type(def_mesh_2d),     target  :: fem_mesh
type(def_0_form_2d),   target  :: field_u
class(def_field_2d),   pointer :: ptr_field
type(def_matrix_2d) ,  target  :: matrix_stiffnes
class(def_matrix_2d),  pointer :: ptr_system

type(def_linear_solver_2d) :: solver
integer                    :: n_var_sys 
integer                    :: n_var_unknown

real(kind=rk), dimension(:,:), pointer  :: diagnostics
integer(kind=jorek_ints_kind)           :: nstep_max 
integer, parameter, private             :: dtllevel_base=0
integer, parameter                      :: nvar = 1
integer, parameter                      :: matrix_a_id = 0
integer                                 :: mode_m1
integer                                 :: mode_n1
real(kind=jorek_coef_kind)              :: a
real(kind=jorek_coef_kind)              :: acenter
real(kind=jorek_coef_kind)              :: r0
real(kind=jorek_coef_kind)              :: z0

end module jorek_model

module sll_jorek

use typedef
use jorek_model
use field
use mesh
use spm
use field
use linear_solver
use matrix
use space
use fem

implicit none

type, public :: sll_jorek_solver

end type sll_jorek_solver

interface sll_create
module procedure initialize_jorek
end interface sll_create
interface sll_solve
module procedure solve_jorek
end interface sll_solve

contains

subroutine initialize_jorek(jorek)

type(sll_jorek_solver) :: jorek

integer          :: nb_args
integer          :: ierr
integer          :: tracelogdetail 
integer          :: tracelogoutput

character(len = 1024)         :: exec
character(len = 1024)         :: rootname
character(len = *), parameter :: mod_name = "main"

character(len=1024) :: filename_parameter

integer, parameter  :: n_dim = 2
integer             :: err
character(len=1024) :: dirname
integer             :: myrank
character(len=1024) :: argname


argname = "--parameters"
call jorek_get_arguments(argname, filename_parameter, ierr)
myrank = 0
#ifdef MPI_ENABLED
call mpi_comm_rank ( mpi_comm_world, myrank, err)
#endif
if (myrank==0) print *, "Parameters text file ", trim(filename_parameter)
call initialize_jorek_parameters(filename_parameter)

#ifdef DEBUG_TRACE   
! ... only used in 3d case
call jorek_param_getint(int_trace_detail_id,tracelogdetail,ierr)
call jorek_param_getint(int_trace_output_id,tracelogoutput,ierr)
ll_stdoutput = (tracelogoutput==1)
call opentracelog(al_stdoutput = ll_stdoutput, ai_dtllevel = tracelogdetail)

call printlog("=========================================", ai_dtllevel = 0)
#ifdef PETSC_ENABLED
call printlog("== petsc enabled                       ==", ai_dtllevel = 0)
#else
call printlog("== petsc disabled                      ==", ai_dtllevel = 0)
#endif
#ifdef DEBUG_ELEMENT
call printlog("== debug element assembly mode         ==", ai_dtllevel = 0)
#endif 
call printlog("=========================================", ai_dtllevel = 0)
#endif
 
call getarg(0, exec)
rootname = "output"

nb_args = iargc()

print *, "runname : ", trim(exec)

! ... initialize spm, mpi 
call spm_initialize(ierr)

! ... define model parameters 
! ..................................................

current_model       = 1
n_var_unknown       = 1
n_var_sys           = 1
i_vp_rho            = 1
nmatrices           = 1
i_vu_rho            = 1 

call jorek_param_getint(int_modes_m1_id, mode_m1, err)
call jorek_param_getint(int_modes_n1_id, mode_n1, err)

call jorek_param_getreal(real_rgeo_id,r0,err)
call jorek_param_getreal(real_zgeo_id,z0,err)
call jorek_param_getreal(real_amin_id,a,err) 
call jorek_param_getreal(real_acenter_id,acenter,err)
call jorek_param_getint(int_nstep_max_id,nstep_max,err)

! ...
call space_create(space_trial, fem_mesh)
ptr_space_trial => space_trial
call space_create(space_test, fem_mesh)
ptr_space_test => space_test
! ...

! ...
argname = "--geometry"
call jorek_get_arguments(argname, dirname, err)
call create_mesh(fem_mesh, dirname)
! ...

call model_create(fem_model,       &
                  fem_mesh,        &
                  ptr_space_trial, &
                  ptr_space_test,  &
                  n_var_unknown,   &
                  n_var_sys)

call field_create(field_u, fem_model%space_trial, field_name="density")
ptr_field => field_u
call model_append_field(fem_model, ptr_field)

call matrix_create(matrix_stiffnes)
ptr_system => matrix_stiffnes

call matrix_append_unknown_field(ptr_system, ptr_field)
call model_append_matrix(fem_model, ptr_system)

call model_initialize(fem_model)    

! output file
open(unit=li_file_stream_norm, file='output_var_diag.dat', status='unknown')
open(unit=li_file_stream_visu, file='output_var_visu.dat', status='unknown')

end subroutine initialize_jorek

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solve_jorek(jorek)

type(sll_jorek_solver) :: jorek
real(kind=rk), dimension(:), allocatable :: y
real(kind=rk), dimension(1:1,1:2)        :: x
real(kind=rk), dimension(1:1,1:2)        :: v
integer                                  :: i   
real(kind=rk)       :: t_fin, t_deb
character(len=1024) :: argname
integer             :: err
integer(kind=spm_ints_kind)          :: nrows !number of rows
character(len=20)   :: stamp_default
character(len=20)   :: msg
integer             :: myrank
character(len=1024) :: filename
integer             :: ierr

! ... define weak formulation 
ptr_system%ptr_matrix_contribution => matrix_for_vi_vj
ptr_system%ptr_rhs_contribution    => rhs_for_vi
! ...

! ... loop over elements 
call cpu_time(t_deb)
call model_assembly(fem_model, ptr_system)
call cpu_time(t_fin)
! ...

print *, "rhs ", ptr_system%opr_global_rhs(1:6)

! ... example of solver calls
argname = "--solver"
call jorek_get_arguments(argname, filename, err)
call linear_solver_create(solver, ptr_system, filename)
call linear_solver_solve(solver,                     &
                         ptr_system%opr_global_rhs,  &
                         ptr_system%opr_global_unknown) 

call linear_solver_free(solver)

! ... update related fields to the linear system
call update_fields(fem_model, ptr_system)

! ... evaluate field
x = 1.0
v = 0.0
call field_evaluate(field_u, 1, x, v)
print *, "field value at ", x, " is ", v
! ...

call get_nr_matrix(ptr_system, nrows)

write(msg,*) myrank
stamp_default = "-proc_" // trim ( adjustl ( msg ) )
stamp_default = trim ( adjustl ( adjustr ( stamp_default ) ) )

call getarg(0, filename)

allocate(y(nrows))
y=0.d0 
call spm_matmult(matrix_a_id, ptr_system%opr_global_unknown, y, err)

print *, maxval(y-ptr_system%opr_global_rhs), &
         minval(y-ptr_system%opr_global_rhs) 

! ... evaluate unknowns on vertecies and compte model norms 
call model_diagnostics(fem_model,             &
                       analytical_model,      &
                       assembly_diagnostics,  &
                       plot_diagnostics,      &
                       ptr_system)

call model_save(fem_model, analytical_model, filename, 0)

open(unit=12,file=trim("rhs"//adjustl(adjustr(stamp_default)))//".txt"&
        & , action="write", status="replace")
do i=1,nrows
   write(12,*) ptr_system%opr_global_rhs(i) 
end do
close(12)
open(unit=13,file=trim ("unknown"//adjustl(adjustr(stamp_default)))//".txt"&
        & , action="write", status="replace")
do i=1,nrows
   write(13,*) ptr_system%opr_global_unknown(i) 
end do
close(13)
open(unit=14, file=trim ("var"//adjustl(adjustr(stamp_default)))//".txt"&
        & , action="write", status="replace")
do i=1, fem_model%ptr_mesh%oi_n_nodes
   write(14,*) fem_model%opr_global_var(:,i) 
end do
close(14)
! ...

#ifdef DEBUG_TRACE
call concatmsg(" min of var ", ai_dtllevel = 0)
call concatmsg(minval(fem_model%opr_global_var(1,:)), ai_dtllevel = 0)
call concatmsg(" max of var ", ai_dtllevel = 0)
call concatmsg(maxval(fem_model%opr_global_var(1,:)), ai_dtllevel = 0)
call printmsg(ai_dtllevel = 0)
#endif

call model_free(fem_model)

close(li_file_stream_norm)
close(li_file_stream_visu)

call spm_cleanall(ierr)
call spm_finalize(ierr)

#ifdef DEBUG_TRACE   
call closetracelog()
#endif

end subroutine solve_jorek


subroutine  rhs_for_vi(ptr_matrix, bbox2di, gbox2d)
class(def_matrix_2d), pointer :: ptr_matrix
type(def_blackbox_2d) :: bbox2di
type(def_greenbox_2d)                :: gbox2d
! local
real(kind=rk) :: f_rhs
real(kind=rk) :: contribution
integer       :: ijg
real(kind=rk) :: wvol
real(kind=rk) :: vi_0
real(kind=rk) :: vi_r
real(kind=rk) :: vi_z
real(kind=rk) :: vi_rr
real(kind=rk) :: vi_rz
real(kind=rk) :: vi_zz
real(kind=rk), dimension(:), pointer :: rhs_contribution

rhs_contribution => gbox2d%rhs_contribution

ijg   = bbox2di%ijg
wvol  = bbox2di%wvol(ijg)
vi_0  = bbox2di%b_0(ijg)
vi_r  = bbox2di%b_x1(ijg)
vi_z  = bbox2di%b_x2(ijg)
vi_rr = bbox2di%b_x1x1(ijg)
vi_rz = bbox2di%b_x1x2(ijg)
vi_zz = bbox2di%b_x2x2(ijg)

! ... add l2 contribution
f_rhs       = analytical_rhs(bbox2di)

contribution                = vi_0 *wvol*f_rhs
rhs_contribution(i_vu_rho)  =  contribution 
! ...

end subroutine rhs_for_vi
! ............................................................... 

! ............................................................... 
subroutine matrix_for_vi_vj(ptr_matrix, bbox2di, bbox2dj, gbox2d)

class(def_matrix_2d), pointer :: ptr_matrix
type(def_blackbox_2d)         :: bbox2di
type(def_blackbox_2d)         :: bbox2dj
type(def_greenbox_2d)         :: gbox2d
! local
real(kind=rk) :: contribution
integer       :: ijg
real(kind=rk) :: wvol
real(kind=rk) :: vi_0
real(kind=rk) :: vi_r
real(kind=rk) :: vi_z
real(kind=rk) :: vi_rr
real(kind=rk) :: vi_rz
real(kind=rk) :: vi_zz
real(kind=rk) :: vj_0
real(kind=rk) :: vj_r
real(kind=rk) :: vj_z
real(kind=rk) :: vj_rr
real(kind=rk) :: vj_rz
real(kind=rk) :: vj_zz
real(kind=rk), dimension(:,:), pointer :: matrix_contribution

matrix_contribution => gbox2d%matrix_contribution

ijg   = bbox2di%ijg
wvol  = bbox2di%wvol(ijg)
vi_0  = bbox2di%b_0(ijg)    ; vj_0  = bbox2dj%b_0(ijg)
vi_r  = bbox2di%b_x1(ijg)   ; vj_r  = bbox2dj%b_x1(ijg)
vi_z  = bbox2di%b_x2(ijg)   ; vj_z  = bbox2dj%b_x2(ijg)
vi_rr = bbox2di%b_x1x1(ijg) ; vj_rr = bbox2dj%b_x1x1(ijg)
vi_rz = bbox2di%b_x1x2(ijg) ; vj_rz = bbox2dj%b_x1x2(ijg)
vi_zz = bbox2di%b_x2x2(ijg) ; vj_zz = bbox2dj%b_x2x2(ijg)

! ... add stiffness contribution
contribution = ( vi_r*vj_r + vi_z*vj_z ) * wvol

matrix_contribution(i_vu_rho, i_vu_rho) =  contribution 
! ...

end subroutine matrix_for_vi_vj

subroutine assembly_diagnostics(bbox2d, gbox2d,nstep)

type(def_blackbox_2d) :: bbox2d
type(def_greenbox_2d)                :: gbox2d
integer       :: ijg
integer       :: nstep
real(kind=rk) :: wvol

if(nstep .ge. 0) then

  ijg   = bbox2d%ijg
  wvol  = bbox2d%wvol(ijg)

  fem_model%opr_diagnostics(1,nstep+1) = fem_model%opr_diagnostics(1,nstep+1) + &
 	  & gbox2d%varn_0(1, ijg) * &
 	  & wvol

  fem_model%opr_diagnostics(2,nstep+1) = fem_model%opr_diagnostics(2,nstep+1) + &
 	  & gbox2d%varn_0(1, ijg) * &
 	  & gbox2d%varn_0(1, ijg) * &
 	  & wvol

  fem_model%opr_diagnostics(3,nstep+1) = fem_model%opr_diagnostics(3,nstep+1) + &
             & gbox2d%varn_x1(1, ijg) * &
 	  & gbox2d%varn_x1(1, ijg) * &
 	  & wvol  + &
 	  & gbox2d%varn_x2(1, ijg) * &
 	  & gbox2d%varn_x2(1, ijg) * &
 	  & wvol 

  !... diff norms
  fem_model%opr_diagnostics(4,nstep+1) = fem_model%opr_diagnostics(4,nstep+1) + &
             & ( gbox2d%varn_0(1, ijg) - gbox2d%sol_analytical(ijg, 1, 1) ) * &
 	  & ( gbox2d%varn_0(1, ijg) - gbox2d%sol_analytical(ijg, 1, 1) ) * &
 	  & wvol 

  fem_model%opr_diagnostics(5,nstep+1) = fem_model%opr_diagnostics(5,nstep+1) + &
 	  & ( gbox2d%varn_x1(1, ijg) - gbox2d%sol_analytical(ijg, 1, 2) ) * &
 	  & ( gbox2d%varn_x1(1, ijg) - gbox2d%sol_analytical(ijg, 1, 2) ) * &
 	  & wvol + &
 	  & ( gbox2d%varn_x2(1, ijg) - gbox2d%sol_analytical(ijg, 1, 3) ) * &
 	  & ( gbox2d%varn_x2(1, ijg) - gbox2d%sol_analytical(ijg, 1, 3) ) * &
 	  & wvol 

end if

end subroutine assembly_diagnostics


subroutine plot_diagnostics(nstep)
implicit none
integer ::coordinates_poloidal,ierr,nstep
real(kind=rk) :: dt, time

call jorek_param_getint(int_coordinates_poloidal_id,coordinates_poloidal,ierr)
call jorek_param_getreal(real_dt_id,dt,ierr)

if(nstep .ge. 0) then

  fem_model%opr_diagnostics(2,nstep+1) = sqrt(fem_model%opr_diagnostics(2,nstep+1))
  fem_model%opr_diagnostics(3,nstep+1) = sqrt(fem_model%opr_diagnostics(3,nstep+1))
 
  fem_model%opr_diagnostics(4,nstep+1) = sqrt(fem_model%opr_diagnostics(4,nstep+1))
  fem_model%opr_diagnostics(5,nstep+1) = sqrt(fem_model%opr_diagnostics(5,nstep+1))
  
  print *, "======      masse ======"
  print *,'masse :',fem_model%opr_diagnostics(1,nstep+1)
  print *, "======      norms ======"
  print *,'norm l2 :',fem_model%opr_diagnostics(2,nstep+1)
  print *,'norm h1 :',fem_model%opr_diagnostics(3,nstep+1)
  print *, "====== diff-norms ======"
  print *,'error l2 :',fem_model%opr_diagnostics(4,nstep+1)/fem_model%opr_diagnostics(2,nstep+1)
  print *,'error h1 :',fem_model%opr_diagnostics(5,nstep+1)/fem_model%opr_diagnostics(3,nstep+1)
  print *, "========================"
  
  if(nstep .eq. 0) then
     write(li_file_stream_norm,*) '# ,time, var_id, masse,  norm_l2, semi_norm_h1, diff_norm_l2, diff_semi_norm_h1'
  end if
  time=nstep*dt
  write(li_file_stream_norm,*) time,                         &
                          1,                                 &
                       fem_model%opr_diagnostics(1,nstep+1), &
                       fem_model%opr_diagnostics(2,nstep+1), &
                       fem_model%opr_diagnostics(3,nstep+1), &
                       fem_model%opr_diagnostics(4,nstep+1), &
                       fem_model%opr_diagnostics(5,nstep+1)

end if

return 
end subroutine plot_diagnostics
! .........................................


!! test case 1 :
!!             - mesh: square, collela or square-periodic
!!             - r0=z0=1 (lenght square)
!!             - cylindrical or cartesian coordindate
!!             - dirichet or periodic boundary condition
!!             - f(x,y)= 8pi**2 * sin(2pi*r)*sin(2pi*z)

!! test case 2 : validate this case

!! test case 3 : validate this case

!! test case 4 :
!!             - mesh: square-periodic
!!             - r0=z0=1 (lenght square)
!!             - cylindrical or cartesian coordindate
!!             - periodic boundary condition
!!             - f(x,y)= 8pi**2 * cos(2pi*r)*cos(2pi*z)


function analytical_rhs(bbox2d)

real(kind=rk)         :: analytical_rhs 
type(def_blackbox_2d) :: bbox2d
real(kind=rk)         :: r, z

real(kind=rk)         :: k1
real(kind=rk)         :: k2
integer               :: testcase
integer               :: ierr
real(kind=rk)         :: r0
real(kind=rk)         :: a
real(kind=rk)         :: acenter
integer               :: ijg

ijg = bbox2d%ijg

r   = bbox2d%xp_0(1,ijg)
z   = bbox2d%xp_0(2,ijg)

k1 = 2.0 * pi * float(mode_m1)
k2 = 2.0 * pi * float(mode_n1)

call jorek_param_getint(int_testcase_id, testcase, ierr)
if(testcase .eq. 1) then
   analytical_rhs = (k1**2 + k2**2) * sin(k1*r)*sin(k2*z)
endif    
if(testcase .eq. 2) then
   analytical_rhs = 4.0 * ( r**2 + z**2 ) * sin(1.0-r**2-z**2) &
 	& + 4.0 * cos ( 1.0 - r**2 - z**2 ) 
endif
if(testcase .eq. 3) then
   analytical_rhs =8*z**2/a**2 - 4 + 2*(2*r - 2*r0)**2/a**2 &
 	&+ 4*(z**2 + (r - r0)**2)/a**2 &
 	& + 4*(z**2 - acenter**2 + (r - r0)**2)/a**2 
endif

! periodic but not dirichet homogeneous
 if(testcase .eq. 4) then
   analytical_rhs = (k1**2 + k2**2) * cos(k1*r)*cos(k2*z)
endif 

end function analytical_rhs

subroutine analytical_model(x, v,info,api_info,n_variable,n_dimension)

integer                                          :: n_dimension
integer                                          :: n_variable
real(kind=8), dimension(n_dimension), intent(in) :: x
real(kind=8), dimension(n_variable, n_dimension+1), intent(out) :: v
real(kind=8), dimension(10), intent(in)          :: info
integer, dimension(10), intent(in)               :: api_info
! local
real(kind=rk) :: r
real(kind=rk) :: z
real(kind=rk) :: k1
real(kind=rk) :: k2
integer       :: testcase
integer       :: ierr

r = x(1)
z = x(2)

k1 = 2.0 * pi * float(mode_m1)
k2 = 2.0 * pi * float(mode_n1)

call jorek_param_getint(int_testcase_id, testcase, ierr)
if(testcase .eq. 1) then
  v(1,1) = sin(k1*r)*sin(k2*z)    ! ... u
  v(1,2) = k1*cos(k1*r)*sin(k2*z) ! ... u_r
  v(1,3) = k2*sin(k1*r)*cos(k2*z) ! ... u_z
endif    
if(testcase .eq. 2) then
  v(1,1) = sin ( 1.0 - r**2 - z**2 ) ! ... u
  v(1,2) = - 2.0 * r * cos ( 1.0 - r**2 - z**2) ! ... u_r
  v(1,3) = - 2.0 * z * cos ( 1.0 - r**2 - z**2) ! ... u_z
endif
if(testcase .eq. 3) then
  v(1,1) = (1-(z**2+(r-r0)**2)/a**2)*(z**2-acenter**2+(r-r0)**2) ! ... u
  v(1,2) = (1-(z**2+(r-r0)**2)/a**2)*(2*r-2*r0) &                ! ... u_r
     & - (2*r-2*r0)*(z**2-acenter**2+(r-r0)**2)/a**2
  v(1,3) = 2*z*(1-(z**2+(r-r0)**2)/a**2) &
     & - 2*z*(z**2 - acenter**2+(r-r0)**2)/a**2 ! ... u_z
endif

! periodic but not dirichet homogeneous
if(testcase .eq. 4) then
  v(1, 1) = cos(k1*r)*cos(k2*z)     ! ... u
  v(1, 2) = -k1*sin(k1*r)*cos(k2*z) ! ... u_r
  v(1, 3) = -k2*cos(k1*r)*sin(k2*z) ! ... u_z
endif    

end subroutine analytical_model

end module sll_jorek
