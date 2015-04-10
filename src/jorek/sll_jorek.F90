module jorek_data
#include "sll_jorek.h"

  type(def_mesh_2d)                          :: mesh2d
  type(def_greenbox_2d)                      :: gbox2d
  type(def_quadrature_square), target        :: quad
  type(def_linear_solver)                    :: solver
  integer                                    :: n_var_sys 
  integer                                    :: n_var_unknown
  real(kind=rk), dimension(:),   allocatable :: global_rhs 
  real(kind=rk), dimension(:),   allocatable :: global_unknown
  real(kind=rk), dimension(:,:), allocatable :: var
  real(kind=rk), dimension(:,:), pointer     :: diagnostics
  integer(kind=jorek_ints_kind)              :: nstep_max 
  integer, parameter                         :: dtllevel_base = 0
  integer, parameter                         :: nvar = 1
  integer, parameter                         :: matrix_a_id = 0
  integer                                    :: mode_m1
  integer                                    :: mode_n1
  real(kind=jorek_coef_kind)                 :: a
  real(kind=jorek_coef_kind)                 :: acenter
  real(kind=jorek_coef_kind)                 :: r0
  real(kind=jorek_coef_kind)                 :: z0
  integer,       dimension(:), allocatable   :: rvars
  integer,       dimension(:), allocatable   :: cvars

end module jorek_data

module sll_jorek

#include "sll_working_precision.h"

use jorek_data

implicit none

type, public :: sll_jorek_solver


end type sll_jorek_solver

interface sll_create
  module procedure :: initialize_jorek
end interface sll_create

interface sll_solve
  module procedure :: solve_jorek
end interface sll_solve

interface sll_delete
  module procedure :: delete_jorek
end interface sll_delete

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine solve_jorek(this)
  
  type(sll_jorek_solver)                   :: this
  real(kind=rk)                            :: t_fin
  real(kind=rk)                            :: t_deb
  character(len = 1024)                    :: filename
  integer(kind=spm_ints_kind)              :: nrows 
  integer                                  :: ierr
  integer                                  :: myrank
  character(len=20)                        :: stamp_default
  character(len=20)                        :: msg
  integer                                  :: i   

  call mpi_comm_rank ( mpi_comm_world, myrank, ierr)

  call spm_getnr(matrix_a_id, nrows, ierr)

  write(msg,*) myrank
  stamp_default = "-proc_"//trim(adjustl(msg))
  stamp_default = trim(adjustl(adjustr(stamp_default)))

  call getarg(0, filename)

  ! ... loop over elements 
  call cpu_time(t_deb)
  call loop_on_elmts( matrix_a_id,      &
                      myrank,           &
                      mesh2d,           & 
                      gbox2d,           &
                      rhs_for_vi,       &
                      matrix_for_vi_vj, &
                      global_unknown,   &
                      var,              &
                      global_rhs,       &
                      rvars,            &
                      cvars)

  call cpu_time(t_fin)

  print *, "rhs ", global_rhs(1:6)

  ! ... example of solver calls
  call linear_solver_new_with_matrix_id(solver, matrix_a_id)
  call linear_solver_solve(solver, global_rhs, global_unknown) 
  call linear_solver_free(solver)

  print *, "unknown ", global_unknown(1:6)

  ! ... evaluate unknowns on verticies and compute model norms 
  call evaluate_on_elmts(matrix_a_id,          &
                         myrank,               &
                         mesh2d,               &
                         gbox2d,               &
                         global_unknown,       &
                         var,                  &
                         analytical_model,     &
                         assembly_diagnostics, &
                         plot_diagnostics,     &
                         rvars,                &
                         cvars,0)

  filename = trim(filename) 
  call savemeshes(mesh2d, gbox2d, var, analytical_model, filename,0)

  open(unit=12, file=trim("rhs"//adjustl(adjustr(stamp_default)))//".txt"&
            & , action="write", status="replace")
  do i=1,nrows
     write(12,*) global_rhs(i) 
  end do
  close(12)
  open(unit=13, file=trim("unknown"//adjustl(adjustr(stamp_default)))//".txt"&
            & , action="write", status="replace")
  do i=1,nrows
     write(13,*) global_unknown(i) 
  end do
  close(13)
  open(unit=14, file=trim("var"//adjustl(adjustr(stamp_default)))//".txt"&
            & , action="write", status="replace")
  do i=1, mesh2d%oi_n_nodes
     write(14,*) var(:,i) 
  end do
  close(14)

end subroutine solve_jorek

subroutine delete_jorek(this)

  type(sll_jorek_solver) :: this
  sll_int32 :: ierr

  deallocate(global_rhs)
  deallocate(global_unknown)
  deallocate(var)
  deallocate(diagnostics)

  call free_greenbox(gbox2d)
  call spm_cleanall(ierr)
  call spm_finalize(ierr)

end subroutine delete_jorek

subroutine initialize_jorek(this)

  type(sll_jorek_solver)             :: this
  character(len=1024)                :: filename_parameter

  sll_int32                          :: imesh=1
  sll_int32                          :: ierr
  integer, parameter                 :: n_dim = 2

  character(len=*), parameter        :: mod_name = "poisson_2d"
  sll_int32                          :: n_gauss_rp
  sll_int32                          :: n_gauss_zp
  character(len=1024)                :: flag 
  integer(kind=spm_ints_kind)        :: nrows 
  integer(kind=spm_ints_kind)        :: ncols
  integer                            :: i
  integer                            :: n_nodes_global

  flag = "--parameters"
  call jorek_get_arguments(flag, filename_parameter, ierr)

  call initialize_jorek_parameters(filename_parameter)

  ! ... define model parameters 

  current_model = 1
  n_var_unknown = 1
  n_var_sys     = 1
  i_vp_rho      = 1
  nmatrices     = 1
  i_vu_rho      = 1 

  allocate(namesvaru(n_var_unknown)) 
  namesvaru(i_vu_rho)  = "density"
  allocate(namesvarp(n_var_unknown)) 
  namesvarp(i_vp_rho)  = "density"

  call jorek_param_getint(int_modes_m1_id, mode_m1, ierr)
  call jorek_param_getint(int_modes_n1_id, mode_n1, ierr)

  call jorek_param_getreal(real_rgeo_id,r0,ierr)
  call jorek_param_getreal(real_zgeo_id,z0,ierr)
  call jorek_param_getreal(real_amin_id,a,ierr) 
  call jorek_param_getreal(real_acenter_id,acenter,ierr)
  call jorek_param_getint(int_nstep_max_id,nstep_max,ierr)

  call create_quadrature(quad, n_dim, 0, 3, 3)

  mesh2d%ptr_quad => quad

  ! ... define basis functions 
  call initbasis( mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  n_gauss_rp,                    &
                  n_gauss_zp,                    &
                  mesh2d%ptr_quad%oi_n_points)

  call jorek_param_getint(int_typemesh_id,imesh,ierr)

  call initgrid(mesh2d,                          &
                imesh,                           &
                mesh2d%oi_n_nodes,               &
                mesh2d%oi_n_elmts,               &
                mesh2d%oi_n_max_vtex_per_elmt,   &
                mesh2d%oi_n_max_order,           &
                mesh2d%oi_n_max_vtex_per_elmt,   &
                mesh2d%oi_n_max_nbnet,           &
                n_var_unknown,                   &
                mesh2d%oi_n_max_order,           &
                mesh2d%oi_n_nzero_bloc,          &
                mesh2d%ptr_quad%oi_n_points      )

  call spm_initialize(nmatrices, ierr)

  allocate(rvars(0:n_var_sys))
  allocate(cvars(0:n_var_sys))

  rvars(0) = n_var_sys
  do i = 1, n_var_sys
     rvars(i) = i
  end do

  cvars(0) = n_var_sys
  do i = 1, n_var_sys
     cvars(i) = i
  end do

  call initialize_matrix(matrix_a_id, mesh2d, rvars, cvars)

  ! ... get the new size of the matrix 
  call spm_getnr(matrix_a_id, nrows, ierr)
  call spm_getnc(matrix_a_id, ncols, ierr)
  ! ...

  n_nodes_global = mesh2d % oi_n_nodes

  allocate(global_rhs(ncols))
  allocate(global_unknown(ncols))
  allocate(var(mesh2d%oi_n_max_order*n_var_unknown, n_nodes_global) )
  allocate(diagnostics(n_diag,nstep_max+1))

  global_rhs     = 0.0
  global_unknown = 0.0
  var            = 0.0
  diagnostics    = 0.0

  call create_greenbox(gbox2d,                      &
                       n_var_unknown,               &
                       n_var_sys,                   &
                       mesh2d%ptr_quad%oi_n_points)

  call coordinates_initialize(2, mesh2d%ptr_quad%oi_n_points, ierr)

end subroutine initialize_jorek



!! test case 1 :
!!   - mesh: square, collela or square-periodic
!!   - r0=z0=1 (lenght square)
!!   - cylindrical or cartesian coordindate
!!   - dirichet or periodic boundary condition
!!   - f(x,y)= 8pi**2 * sin(2pi*r)*sin(2pi*z)
!! test case 2 : validate this case
!! test case 3 : validate this case
!! test case 4 :
!!   - mesh: square-periodic
!!   - r0=z0=1 (lenght square)
!!   - cylindrical or cartesian coordindate
!!   - periodic boundary condition
!!   - f(x,y)= 8pi**2 * cos(2pi*r)*cos(2pi*z)
function analytical_rhs(bbox2d)

  real(kind=rk)         :: analytical_rhs 
  type(def_blackbox_2d) :: bbox2d

  real(kind=rk)         :: r
  real(kind=rk)         :: z
  integer               :: testcase
  integer               :: ierr
  integer               :: ijg
  real(kind=rk)         :: k1
  real(kind=rk)         :: k2

  ijg = bbox2d%ijg

  r   = bbox2d%xp_0(1,ijg)
  z   = bbox2d%xp_0(2,ijg)

  k1 = 2.0 * pi * float(mode_m1)
  k2 = 2.0 * pi * float(mode_n1)
 
  call jorek_param_getint(int_testcase_id, testcase, ierr)

  select case (testcase)

  case(1) 
     analytical_rhs = (k1**2+k2**2)*sin(k1*r)*sin(k2*z)
  case(2)
     analytical_rhs = 4.0*(r**2+z**2)*sin(1.0-r**2-z**2) &
   	            & + 4.0*cos(1.0-r**2-z**2) 
  case(3)
    analytical_rhs =8*z**2/a**2-4+2*(2*r-2*r0)**2/a**2 &
                 & +4*(z**2+(r-r0)**2)/a**2            &
                 & +4*(z**2-acenter**2+(r-r0)**2)/a**2
  case(4)
    analytical_rhs = (k1**2 + k2**2) * cos(k1*r)*cos(k2*z) ! periodic but not dirichet homogeneous
  end select
   
end function analytical_rhs

subroutine analytical_model( x,          &
                             v,          &
                             apr_info,   &
                             api_info,   &
                             n_variable, &
                             n_dimension )

  integer                                          :: n_variable
  integer                                          :: n_dimension
  real(kind=8), dimension(n_dimension), intent(in) :: x
  real(kind=8), dimension(n_variable, n_dimension+1), intent(out) :: v
  real(kind=8), dimension(10), intent(in) :: apr_info
  integer,      dimension(10), intent(in) :: api_info

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

  select case(testcase)
  case(1)
     ! ... u
     v(1,1) = sin(k1*r)*sin(k2*z)
     ! ... u_r
     v(1,2) = k1*cos(k1*r)*sin(k2*z)
     ! ... u_z
     v(1,3) = k2*sin(k1*r)*cos(k2*z)
     ! ...
  case(2)
     ! ... u
     v(1,1) = sin(1.0-r**2-z**2)
     ! ... u_r
     v(1,2) = -2.0*r*cos(1.0-r**2-z**2)
     ! ... u_z
     v(1,3) = -2.0*z*cos(1.0-r**2-z**2)
     ! ...
  case(3)
     ! ... u
     v(1,1) = (1-(z**2+(r-r0)**2)/a**2)*(z**2-acenter**2+(r-r0)**2)
     ! ... u_r
     v(1,2) = (1-(z**2+(r-r0)**2)/a**2)*(2*r-2*r0) &
   	    & - (2*r-2*r0)*(z**2-acenter**2+(r-r0)**2)/a**2
     ! ... u_z
     v(1,3) = 2*z*(1-(z**2+(r-r0)**2)/a**2) &
   	    & - 2*z*(z**2-acenter**2+(r-r0)**2)/a**2
     ! ...

  case(4) ! periodic but not dirichet homogeneous
     ! ... u
     v(1,1) = cos(k1*r)*cos(k2*z)
     ! ... u_r
     v(1,2) = -k1*sin(k1*r)*cos(k2*z)
     ! ... u_z
     v(1,3) = -k2*cos(k1*r)*sin(k2*z)
     ! ...
  end select
  
end subroutine analytical_model

subroutine assembly_diagnostics(bbox2d, gbox2d,nstep)

  type(def_blackbox_2d) :: bbox2d
  type(def_greenbox_2d) :: gbox2d
  integer               :: ijg
  integer               :: nstep
  real(kind=rk)         :: wvol
  
  if(nstep .ge. 0) then
  
    ijg   = bbox2d % ijg
    wvol  = bbox2d % wvol(ijg)
  
    diagnostics(1,nstep+1) = diagnostics(1,nstep+1) + &
   	  & gbox2d%varn_0(1, ijg) * &
   	  & wvol
  
    diagnostics(2,nstep+1) = diagnostics(2,nstep+1) + &
   	  & gbox2d%varn_0(1, ijg) * &
   	  & gbox2d%varn_0(1, ijg) * &
   	  & wvol
  
    diagnostics(3,nstep+1) = diagnostics(3,nstep+1) + &
               & gbox2d%varn_x1(1, ijg) * &
   	  & gbox2d%varn_x1(1, ijg) * &
   	  & wvol  + &
   	  & gbox2d%varn_x2(1, ijg) * &
   	  & gbox2d%varn_x2(1, ijg) * &
   	  & wvol 
  
    !... diff norms
    diagnostics(4,nstep+1) = diagnostics(4,nstep+1) + &
        & ( gbox2d%varn_0(1, ijg) - gbox2d%sol_analytical(ijg, 1, 1) ) * &
   	  & ( gbox2d%varn_0(1, ijg) - gbox2d%sol_analytical(ijg, 1, 1) ) * &
   	  & wvol 
  
    diagnostics(5,nstep+1) = diagnostics(5,nstep+1) + &
   	  & ( gbox2d%varn_x1(1, ijg) - gbox2d%sol_analytical(ijg, 1, 2) ) * &
   	  & ( gbox2d%varn_x1(1, ijg) - gbox2d%sol_analytical(ijg, 1, 2) ) * &
   	  & wvol + &
   	  & ( gbox2d%varn_x2(1, ijg) - gbox2d%sol_analytical(ijg, 1, 3) ) * &
   	  & ( gbox2d%varn_x2(1, ijg) - gbox2d%sol_analytical(ijg, 1, 3) ) * &
   	  & wvol 
  
  end if
  
  return 

end subroutine assembly_diagnostics

subroutine plot_diagnostics(nstep)

  integer ::coordinates_poloidal,ierr,nstep
  real(kind=rk) :: dt, time
  
  call jorek_param_getint(int_coordinates_poloidal_id,coordinates_poloidal,ierr)
  call jorek_param_getreal(real_dt_id,dt,ierr)

  if(nstep .ge. 0) then

    diagnostics(2,nstep+1) = sqrt(diagnostics(2,nstep+1))
    diagnostics(3,nstep+1) = sqrt(diagnostics(3,nstep+1))

    diagnostics(4,nstep+1) = sqrt(diagnostics(4,nstep+1))
    diagnostics(5,nstep+1) = sqrt(diagnostics(5,nstep+1))
    
    print *, "======      masse ======"
    print *,'masse :',diagnostics(1,nstep+1)
    print *, "======      norms ======"
    print *,'norm l2 :',diagnostics(2,nstep+1)
    print *,'norm h1 :',diagnostics(3,nstep+1)
    print *, "====== diff-norms ======"
    print *,'error l2 :',diagnostics(4,nstep+1)/diagnostics(2,nstep+1)
    print *,'error h1 :',diagnostics(5,nstep+1)/diagnostics(3,nstep+1)
    print *, "========================"
    
    if(nstep .eq. 0) then
      write(6,*) '# ,time,var_id,masse,norm_l2,semi_norm_h1,diff_norm_l2,diff_semi_norm_h1'
    end if
    time=nstep*dt
    write(6,*) time, 1, diagnostics(1:5,nstep+1)

  end if

end subroutine plot_diagnostics

subroutine rhs_for_vi(bbox2di, gbox2d)

  type(def_blackbox_2d) :: bbox2di
  type(def_greenbox_2d) :: gbox2d
  
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

  rhs_contribution => gbox2d % rhs_contribution

  ijg   = bbox2di%ijg
  wvol  = bbox2di%wvol(ijg)
  vi_0  = bbox2di%b_0(ijg)
  vi_r  = bbox2di%b_x1(ijg)
  vi_z  = bbox2di%b_x2(ijg)
  vi_rr = bbox2di%b_x1x1(ijg)
  vi_rz = bbox2di%b_x1x2(ijg)
  vi_zz = bbox2di%b_x2x2(ijg)

  ! ... add l2 contribution
  f_rhs = analytical_rhs(bbox2di)

  contribution                = vi_0 *wvol*f_rhs
  rhs_contribution(i_vu_rho)  =  contribution 

end subroutine rhs_for_vi

subroutine matrix_for_vi_vj(bbox2di, bbox2dj, gbox2d)

  type(def_blackbox_2d) :: bbox2di
  type(def_blackbox_2d) :: bbox2dj
  type(def_greenbox_2d) :: gbox2d

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

  matrix_contribution => gbox2d % matrix_contribution

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

end subroutine matrix_for_vi_vj

end module sll_jorek
