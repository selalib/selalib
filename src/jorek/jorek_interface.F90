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
use space_def
use field_def
use matrix_def
use fem_def

!  use febasis2d_splines
!  use febasis2d_boxsplines

implicit none

! ... quad mesh
type(def_fem_quad_2d),    target        :: fem_model 
type(def_space_quad_2d),  target        :: space_trial
class(def_space_quad_2d), pointer       :: ptr_space_trial
type(def_space_quad_2d),  target        :: space_test
class(def_space_quad_2d), pointer       :: ptr_space_test
type(def_mesh_2d),        target        :: fem_mesh
type(def_0_form_2d),      target        :: field_u
class(def_field_2d),      pointer       :: ptr_field
type(def_matrix_2d) ,     target        :: matrix_stiffnes
class(def_matrix_2d),     pointer       :: ptr_system

type(def_linear_solver_2d)              :: solver
integer                                 :: n_var_sys 
integer                                 :: n_var_unknown

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
integer, parameter                      :: one=1
real(8), dimension(:), pointer          :: rho
real(8), dimension(:), pointer          :: e_x
real(8), dimension(:), pointer          :: e_y

end module jorek_model

module jorek_interface

use typedef
use jorek_model
use linear_solver
use field
use mesh
use spm_def
use matrix
use space
use fem

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_jorek_model()

integer                     :: nb_args
integer(kind=spm_ints_kind) :: ierror
integer                     :: tracelogdetail 
integer                     :: tracelogoutput

character(len=1024) :: progname
character(len=1024) :: rootname

character(len=1024) :: filename_parameter

integer, parameter  :: n_dim = 2
integer             :: ierr
character(len=1024) :: dirname
integer             :: myrank
character(len=1024) :: argname

type(def_element), pointer :: elmt => null()

argname = "--parameters"
call jorek_get_arguments(argname, filename_parameter, ierr)
print *, "Parameters text file ", trim(filename_parameter)
call initialize_jorek_parameters(filename_parameter)

call getarg(0, progname)
rootname = "output"

nb_args = iargc()


! ... initialize spm, mpi 
call spm_initialize(ierror)
myrank = 0
#ifdef MPI_ENABLED
call mpi_comm_rank ( mpi_comm_world, myrank, ierr)
#endif
if (myrank == 0) print *, "runname : ", trim(progname)

! ... define model parameters 
! ..................................................

current_model       = 1
n_var_unknown       = 1
n_var_sys           = 1
i_vp_rho            = 1
nmatrices           = 1
i_vu_rho            = 1 

call jorek_param_getint(int_modes_m1_id, mode_m1, ierr)
call jorek_param_getint(int_modes_n1_id, mode_n1, ierr)

call jorek_param_getreal(real_rgeo_id,r0,ierr)
call jorek_param_getreal(real_zgeo_id,z0,ierr)
call jorek_param_getreal(real_amin_id,a,ierr) 
call jorek_param_getreal(real_acenter_id,acenter,ierr)
call jorek_param_getint(int_nstep_max_id,nstep_max,ierr)

call space_create(space_trial, fem_mesh)
ptr_space_trial => space_trial
call space_create(space_test, fem_mesh)
ptr_space_test => space_test

argname = "--geometry"
call jorek_get_arguments(argname, dirname, ierr)
call create_mesh(fem_mesh, dirname=dirname)

call model_create(fem_model,       &
                  fem_mesh,        &
                  ptr_space_trial, &
                  ptr_space_test)

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

end subroutine create_jorek_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_jorek_model()

real(kind=rk), dimension(1:1,1:2)        :: x
real(kind=rk), dimension(1:1,1:2)        :: v
integer                                  :: i   
character(len=1024) :: argname
integer             :: err
integer(kind=spm_ints_kind)          :: nrows !number of rows
character(len=20)   :: stamp_default
character(len=20)   :: msg
integer             :: myrank = 0
character(len=1024) :: filename
integer             :: ierr

! ... define weak formulation 
ptr_system%ptr_matrix_contribution => matrix_for_vi_vj
ptr_system%ptr_rhs_contribution    => rhs_for_vi

! ... loop over elements 
call model_assembly(fem_model, ptr_system)
! ...

! opr_global_var( & 
print*, 'order =', fem_model%space_trial%basis%oi_n_order
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
stamp_default = "-proc_" // trim(adjustl(msg))
stamp_default = trim(adjustl(adjustr(stamp_default)))

! ... evaluate unknowns on vertecies and compte model norms 
ptr_system%ptr_assembly_diagnostics => assembly_diagnostics
call model_diagnostics(fem_model,             &
                       analytical_model,      &
                       plot_diagnostics,      &
                       ptr_system, 0)

call model_save(fem_model, analytical_model, "jorek_plot", 0)

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

close(li_file_stream_norm)
close(li_file_stream_visu)


end subroutine run_jorek_model

subroutine delete_jorek_model()

integer              :: ierr

call model_free(fem_model)
call spm_cleanall(ierr)
call spm_finalize(ierr)

end subroutine delete_jorek_model


subroutine rhs_for_vi(ptr_matrix, bbox2di, gbox2d)

class(def_matrix_2d), pointer :: ptr_matrix
type(def_blackbox_2d)         :: bbox2di
type(def_greenbox_2d)         :: gbox2d

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

rhs_contribution => ptr_matrix%rhs_contribution

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

matrix_contribution => ptr_matrix%matrix_contribution

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

subroutine assembly_diagnostics(ptr_matrix, bbox2d, gbox2d,nstep)

class(def_matrix_2d), pointer :: ptr_matrix
type(def_blackbox_2d)         :: bbox2d
type(def_greenbox_2d)         :: gbox2d

integer                       :: ijg
integer                       :: nstep
real(kind=rk)                 :: wvol

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
else if(testcase .eq. 2) then
  v(1,1) = sin(1.0-r**2-z**2 ) ! ... u
  v(1,2) = -2.0*r*cos(1.0-r**2-z**2) ! ... u_r
  v(1,3) = -2.0*z*cos(1.0-r**2-z**2) ! ... u_z
else if(testcase .eq. 3) then
  v(1,1) = (1-(z**2+(r-r0)**2)/a**2)*(z**2-acenter**2+(r-r0)**2) ! ... u
  v(1,2) = (1-(z**2+(r-r0)**2)/a**2)*(2*r-2*r0) &                ! ... u_r
     & - (2*r-2*r0)*(z**2-acenter**2+(r-r0)**2)/a**2
  v(1,3) = 2*z*(1-(z**2+(r-r0)**2)/a**2) &
     & - 2*z*(z**2 - acenter**2+(r-r0)**2)/a**2 ! ... u_z
else if(testcase .eq. 4) then ! periodic but not dirichet homogeneous
  v(1,1) = cos(k1*r)*cos(k2*z)     ! ... u
  v(1,2) = -k1*sin(k1*r)*cos(k2*z) ! ... u_r
  v(1,3) = -k2*cos(k1*r)*sin(k2*z) ! ... u_z
endif    

end subroutine analytical_model

subroutine compute_electric_fields(bbox2d, gbox2d)
integer              :: prank
integer              :: ierr
type(def_blackbox_2d) :: bbox2d
type(def_greenbox_2d) :: gbox2d
integer               :: ijg
real(kind=rk)         :: wvol

!call mpi_comm_rank(prank, ierr)
ijg   = bbox2d%ijg
wvol  = bbox2d%wvol(ijg)


!call evaluate_on_elmts(ptr_matrix,                      &
!&                      prank,                           &
!&                      fem_model%ptr_mesh,              &
!&                      fem_model%greenbox,              &
!&                      fem_model%space_trial%basis,     &
!&                      ptr_matrix%opr_global_unknown,   &
!&                      fem_model%opr_global_var,        &
!&                      func_analytical,                 &
!&                      Assembly_Diags,                  &
!&                      Plot_Diags,                      &
!&                      0                                )

end subroutine compute_electric_fields

subroutine assembly_electric_fields(bbox2d, gbox2d)

type(def_blackbox_2d) :: bbox2d
type(def_greenbox_2d) :: gbox2d
integer               :: ijg
integer               :: nstep
real(kind=rk)         :: wvol

ijg   = bbox2d%ijg
wvol  = bbox2d%wvol(ijg)

e_x = e_x + gbox2d%varn_x1(1, ijg) * wvol
e_y = e_y + gbox2d%varn_x2(1, ijg) * wvol

end subroutine assembly_electric_fields


!function interpolate_rho_value(bbox_2d)
!real(kind=rk)         :: interpolate_rho_value 
!type(def_blackbox_2d) :: bbox2d
!real(kind=rk)         :: x, y
!integer               :: ijg
!integer, parameter    :: one = 1
!
!ijg = bbox2d%ijg
!is1 = fem_mesh%opo_elements(ijg)%opi_vertices(1)
!is2 = fem_mesh%opo_elements(ijg)%opi_vertices(1)
!is3 = fem_mesh%opo_elements(ijg)%opi_vertices(1)
!is4 = fem_mesh%opo_elements(ijg)%opi_vertices(1)
!
!xs1 = fem_mesh%thenodes%coor2d(1, one, is1 )
!xs2 = fem_mesh%thenodes%coor2d(1, one, is2 )
!xs3 = fem_mesh%thenodes%coor2d(1, one, is3 )
!xs4 = fem_mesh%thenodes%coor2d(1, one, is4 )
!
!ys1 = fem_mesh%thenodes%coor2d(2, one, is1 )
!ys2 = fem_mesh%thenodes%coor2d(2, one, is2 )
!ys3 = fem_mesh%thenodes%coor2d(2, one, is3 )
!ys4 = fem_mesh%thenodes%coor2d(2, one, is4 )
!
!x   = bbox2d%xp_0(1,ijg)
!y   = bbox2d%xp_0(2,ijg)
!
!!call  the interpolator here
!
!end function interpolate_rho_value

end module jorek_interface
