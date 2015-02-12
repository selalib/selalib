!! Copyright (C) 2009,2010,2011,2012  Marco Restelli
!!
!! This file is part of:
!!   FEMilaro -- Finite Element Method toolkit
!!
!! FEMilaro is free software; you can redistribute it and/or modify it
!! under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 3 of the License, or
!! (at your option) any later version.
!!
!! FEMilaro is distributed in the hope that it will be useful, but
!! WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with FEMilaro; If not, see <http://www.gnu.org/licenses/>.
!!
!! author: Marco Restelli                   <marco.restelli@gmail.com>

!>\brief
!!
!! Simplified interface to UMFPACK.
!!
!! \n
!!
!! This module provides a simplified interface to the MUMPS solver,
!! according to the general layout given in \c mod_linsolver_base.
!!
!! The structure of this module is very similar to \c mod_mumpsintf,
!! with the simplification that UMFpack is serial, so this module does
!! not include any MPI related object.
!<----------------------------------------------------------------------
module mod_umfintf

!-----------------------------------------------------------------------

 use mod_utils, only: &
   t_realtime, my_second

 use mod_messages, only: &
   mod_messages_initialized, &
   error,   &
   warning, &
   info

 use mod_kinds, only: &
   mod_kinds_initialized, &
   wp

 use mod_sparse, only: &
   mod_sparse_initialized, &
   ! sparse types
   t_col,       &
   transpose,   &
   matmul,      &
   diag,        &
   spdiag,      &
   clear

 use mod_umfpack, only: &
   mod_umfpack_constructor, &
   mod_umfpack_destructor,  &
   mod_umfpack_initialized, &
   umfpack_control, &
   umfpack_prl,     &
   umfpack_info,    &
   umf_int,         &
   umf_dp,          &
   umf_void,        &
   umfpack_a,       &
   umfpack_at,      &
   umfpack_aat,     &
   umfpack_pt_l,    &
   umfpack_l,       &
   umfpack_lt_p,    &
   umfpack_lat_p,   &
   umfpack_lt,      &
   umfpack_lat,     &
   umfpack_u_qt,    &
   umfpack_u,       &
   umfpack_q_ut,    &
   umfpack_q_uat,   &
   umfpack_ut,      &
   umfpack_uat,     &
   umf4def,         &
   umf4pcon,        &
   umf4pinf,        &
   umf4sym,         &
   umf4num,         &
   umf4sol,         &
   umf4solr,        &
   umf4fsym,        &
   umf4fnum

 use mod_state_vars, only: &
   mod_state_vars_initialized, &
   c_stv

 use mod_linsolver_base, only: &
   mod_linsolver_base_initialized, &
   c_linpb

!-----------------------------------------------------------------------
 
 implicit none

!-----------------------------------------------------------------------

! Module interface

 public :: &
   mod_umfintf_constructor, &
   mod_umfintf_destructor,  &
   mod_umfintf_initialized, &
   c_umfpackpb

 private

!-----------------------------------------------------------------------

 !> Linear UMFpack solver problem
 !!
 !! This type describes a linear system from the UMFpack viewpoint.
 type, extends(c_linpb), abstract :: c_umfpackpb
  !> Set the <tt>[UMFPACK PRL]</tt> parameter to control the amount of
  !! output; see the UMFpack documentation for details.
  integer :: print_level = 1
  !> Control parameters
  real(umf_dp) :: umf_contr(umfpack_control)
  !> UMFpack diagnostics and error codes
  real(umf_dp) :: umf_info ( umfpack_info  )
  !> System matrix
  type(t_col), pointer :: m
  !> solve \f$A^Tx=b\f$
  logical :: transposed_mat = .false.
  !> right-hand side
  real(wp), pointer :: rhs(:)
  !> Pointers to UMFpack objects
  integer(umf_void), private :: symbolic, numeric
  logical, private :: an_set = .false. !< internal consistency check
  logical, private :: fc_set = .false. !< internal consistency check
 contains
  procedure, pass(s) :: factor => umfpack_factor
  procedure, pass(s) :: solve  => umfpack_solve
  procedure, pass(s) :: clean  => umfpack_clean
  procedure(i_xassign), deferred, pass(s) :: xassign
 end type c_umfpackpb

 !> Convert the UMFpack solution into a \c c_stv object
 abstract interface
  subroutine i_xassign(x,s,umfpack_x)
   import :: wp, c_stv, c_umfpackpb
   implicit none
   real(wp),           intent(in) :: umfpack_x(:)
   class(c_umfpackpb), intent(inout) :: s
   class(c_stv),       intent(inout) :: x
  end subroutine i_xassign
 end interface

 real(t_realtime) :: t0, t1
 logical, protected :: &
   mod_umfintf_initialized = .false.
 character(len=*), parameter :: &
   this_mod_name = 'mod_umfintf'

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

 subroutine mod_umfintf_constructor()
  character(len=*), parameter :: &
    this_sub_name = 'constructor'

   !Consistency checks ---------------------------
   if( (mod_messages_initialized.eqv..false.) .or. &
          (mod_kinds_initialized.eqv..false.) .or. &
         (mod_sparse_initialized.eqv..false.) ) then
     call error(this_sub_name,this_mod_name, &
                'Not all the required modules are initialized.')
   endif
   if(mod_umfintf_initialized.eqv..true.) then
     call warning(this_sub_name,this_mod_name, &
                  'Module is already initialized.')
   endif
   !----------------------------------------------

   call mod_umfpack_constructor()

   mod_umfintf_initialized = .true.
 end subroutine mod_umfintf_constructor

!-----------------------------------------------------------------------
 
 subroutine mod_umfintf_destructor()
  character(len=*), parameter :: &
    this_sub_name = 'destructor'
   
   !Consistency checks ---------------------------
   if(mod_umfintf_initialized.eqv..false.) then
     call error(this_sub_name,this_mod_name, &
                'This module is not initialized.')
   endif
   !----------------------------------------------

   call mod_umfpack_destructor()

   mod_umfintf_initialized = .false.
 end subroutine mod_umfintf_destructor

!-----------------------------------------------------------------------

 subroutine umfpack_factor(s,phase)
  class(c_umfpackpb), intent(inout) :: s
  character(len=*), intent(in), optional :: phase
 
  logical :: do_analysis, do_factorization
  character(len=*), parameter :: &
    this_sub_name = 'umfpack_factor'

   ! Check whether only one phase is required
   if(present(phase)) then
     if(trim(phase).eq.'analysis') then
       do_analysis      = .true.
       do_factorization = .false.
     elseif(trim(phase).eq.'factorization') then
       do_analysis      = .false.
       do_factorization = .true.
     else
       call error(this_sub_name,this_mod_name,    &
              'Unknown phase "'//trim(phase)//'".')
     endif
   else
     do_analysis      = .true.
     do_factorization = .true.
   endif

   ! Inizialize umf_contr with the default UMFPACK parameters
   call umf4def(s%umf_contr)
   ! set the verbose mode
   s%umf_contr(umfpack_prl) = real( s%print_level , umf_dp )
   call umf4pcon(s%umf_contr)

   ! Analysis
   if(do_analysis) then
     if(s%an_set) call error(this_sub_name,this_mod_name, &
         'The symbolic analysis is already defined')
     call umf4sym( int(s%m%m,umf_int) , int(s%m%n,umf_int) ,          &
       int(s%m%ap,umf_int), int(s%m%ai,umf_int), real(s%m%ax,umf_dp), &
       s%symbolic , s%umf_contr , s%umf_info )
     call umf4pinf(s%umf_contr,s%umf_info) ! print output
     s%an_set = .true.
   endif

   ! Factorization
   if(do_factorization) then
     if(.not.s%an_set) call error(this_sub_name,this_mod_name, &
         'The symbolic analysis must preceed the factorization')
     if(s%fc_set) then
       call umf4fnum(s%numeric)
       s%fc_set = .false.
     endif
     call umf4num( int(s%m%ap,umf_int) , int(s%m%ai,umf_int) ,        &
                                                 real(s%m%ax,umf_dp), &
       s%symbolic,s%numeric , s%umf_contr , s%umf_info )
     call umf4pinf(s%umf_contr,s%umf_info) ! output
     s%fc_set = .true.
   endif

 end subroutine umfpack_factor
 
!-----------------------------------------------------------------------
 
 subroutine umfpack_solve(x,s)
  class(c_umfpackpb), intent(inout) :: s
  class(c_stv),       intent(inout) :: x
 
  integer(umf_int) :: system_type
  real(umf_dp), allocatable :: x_dp(:)
  character(len=*), parameter :: &
    this_sub_name = 'umfpack_solve'

   if(.not.s%fc_set) call error(this_sub_name,this_mod_name, &
     'A system must be factored beforeit can be solved.')

   system_type = umfpack_a
   if(s%transposed_mat) system_type = umfpack_aat

   ! solve
   allocate(x_dp(s%m%n))
   call umf4solr(system_type,                       &
     int(s%m%ap,umf_int), int(s%m%ai,umf_int),      &
     real(s%m%ax,umf_dp), x_dp, real(s%rhs,umf_dp), &
     s%numeric, s%umf_contr, s%umf_info)
   call umf4pinf(s%umf_contr,s%umf_info) ! print output
   call s%xassign(x,real(x_dp,wp))
   deallocate(x_dp)

 end subroutine umfpack_solve
 
!-----------------------------------------------------------------------
 
 subroutine umfpack_clean(s)
  class(c_umfpackpb), intent(inout) :: s

  character(len=*), parameter :: &
    this_sub_name = 'umfpack_clean'
 
   if(s%an_set) then
     call umf4fsym(s%symbolic); s%an_set = .false.
   endif
   if(s%fc_set) then
     call umf4fnum(s%numeric ); s%fc_set = .false.
   endif

 end subroutine umfpack_clean
 
!-----------------------------------------------------------------------

end module mod_umfintf

