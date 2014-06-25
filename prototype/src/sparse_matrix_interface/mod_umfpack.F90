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
!! UMFPACK interface.
!!
!! \n
!!
!! This module is the fortran interface to <em>umf4_f77wrapper.c</em>.
!! Essentially, we provide here an explicit interface to the UMFPACK C
!! functions. There are two issues in interfacing fortran to C:
!! <ul>
!!  <li> the mane-mangling
!!  <li> the type correspondence.
!! </ul>
!! The name-mangling is handles in <em>umf4_f77wrapper.c</em>, so we
!! don't consider it here. Instead, we address the problem of type and
!! kind compatibility. There are three types involved in the UMFPACK
!! interface, handled as follows.
!! <ul>
!!  <li> integer type: as discussed in the UMFPACK manual, sections
!!  5.1 and 7, UMFPACK can use either \c int or <tt>UF_long</tt>
!!  integers, where <tt>UF_long</tt> is defined in <em>UFconfig.h</em>
!!  and is always a 64-bit integer in a 64-bit code. The file
!!  <em>umf4_f77wrapper.c</em> defines the integer type used in the
!!  fortran interface as
!!  <ul>
!!   <li> \c Int \f$\mapsto\f$ \c UF_long if compiled with
!!   <tt>-DDLONG</tt>
!!   <li> \c Int \f$\mapsto\f$ \c int otherwise,
!!  </ul>
!!  and adjusts the C calls consistently. The present module thus must
!!  match the integer type defined in <em>umf4_f77wrapper.c</em>, and
!!  this is done by the fortran kind <tt>integer(umf_int)<tt>.
!!  Summarizing, one has to set at compile time one of the following
!!  combinations:
!!  <ul>
!!   <li> <em>small matrices:</em>
!!    <ul>
!!     <li> do NOT use <tt>-DDLONG</tt> in the compilation of
!!     <em>umf4_f77wrapper.c</em> and
!!     <li> define in the present module \c umf_int as
!!     <tt>integer*4</tt>;
!!    </ul>
!!   <li> <em>large matrices:</em>
!!    <ul>
!!     <li> use <tt>-DDLONG</tt> in the compilation of
!!     <em>umf4_f77wrapper.c</em> and
!!     <li> define in the present module \c umf_int as
!!     <tt>integer*8</tt>.
!!    </ul>
!!  </ul>
!!  <li> real type: this is the simplest one, since we use double
!!  precision and the fortran kind <tt>real(umf_dp)</tt>, for a
!!  convenient definition of \c umf_dp.
!!  <li> pointer to void: as fortran doesn't have a corresponding
!!  type, this is probably the trickiest part. As suggested in the
!!  UMFPACK manual, section 7, we use the kind <tt>umf_void</tt>
!!  defined as <tt>integer*8</tt>.
!! </ul>
!! \warning This module has been tested with UMFPACK-5.4.0; different
!! versions may require different settings, as explicitly stated also
!! in <em>umf4_f77wrapper.c</em>.
!! \note One should consider rewriting this module by making use of
!! the standardized fortran/C interoperability.
!<
module mod_umfpack
#include "sll_working_precision.h"

!-----------------------------------------------------------------------

!!$ use mod_messages, only: &
!!$   mod_messages_initialized, &
!!$   error,   &
!!$   warning, &
!!$   info

!-----------------------------------------------------------------------
 
 implicit none

!-----------------------------------------------------------------------

! Module interface

 public :: &
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

 private

!-----------------------------------------------------------------------

! Module types and parameters

 ! Choose here the integer type and kind (make sure that what you
 ! select here is consistent with the compile options of
 ! umf4_f77wrapper.c, see the Makefile).
 !integer*4, parameter ::           & ! small matrices
 sll_int32, parameter ::           & 
   umf_int_example = 1
 integer, parameter ::             &
   umf_int = kind(umf_int_example)

 double precision, parameter ::    &
   umf_real_example = 1.0d0
 integer, parameter ::             &
   & umf_dp = kind(umf_real_example)!, &
  ! & umf_int = kind(umf_real_example), &
  ! & umf_void = kind(umf_real_example)

 sll_int64, parameter ::           &
   umf_void_example = 1
 integer, parameter ::             &
   umf_void = kind(umf_void_example)

 ! public members
 integer(umf_int), parameter :: &
   umfpack_a     = 0,  &! Ax = b
   umfpack_at    = 1,  &! AHx = b
   umfpack_aat   = 2,  &! ATx = b
   umfpack_pt_l  = 3,  &! PTLx = b
   umfpack_l     = 4,  &! Lx = b
   umfpack_lt_p  = 5,  &! LHPx = b
   umfpack_lat_p = 6,  &! LTPx = b
   umfpack_lt    = 7,  &! LHx = b
   umfpack_lat   = 8,  &! LTx = b
   umfpack_u_qt  = 9,  &! UQTx = b
   umfpack_u     = 10, &! Ux = b
   umfpack_q_ut  = 11, &! QUHx = b
   umfpack_q_uat = 12, &! QUTx = b
   umfpack_ut    = 13, &! UHx = b
   umfpack_uat   = 14   ! UTx = b

 ! These parameters are taken from umfpack.h. Notice that different
 ! versions of UMFPACK may use different values for these parameteres,
 ! so one shoud check that this module is consistent with the system
 ! umfpack.h.
 integer(umf_int), parameter ::      &
   umfpack_control = 20, &
   umfpack_info = 90,    &
   umfpack_prl = 0 + 1     ! add 1 when going from C to FORTRAN indexes

! Module variables

 ! public members
 logical, protected ::               &
   mod_umfpack_initialized = .false.
 ! private members

 interface 
   subroutine umf4def(control)
    import :: umfpack_control, umf_dp
    implicit none
     real(umf_dp), intent(out) :: control(umfpack_control)
   end subroutine umf4def
 end interface

 interface 
   subroutine umf4fsym(symbolic)
    import :: umf_void
    implicit none
     integer(umf_void), intent(inout) :: symbolic
   end subroutine umf4fsym
 end interface

 interface 
   subroutine umf4fnum(numeric)
    import :: umf_void
    implicit none
     integer(umf_void), intent(inout) :: numeric
   end subroutine umf4fnum
 end interface

 interface 
   subroutine umf4pcon(control)
    import :: umfpack_control, umf_dp
    implicit none
     real(umf_dp), intent(in) :: control(umfpack_control)
   end subroutine umf4pcon
 end interface

 interface 
   subroutine umf4sym(n_row,n_col,ap,ai,ax,symbolic,control,info)
    import :: umfpack_control, umfpack_info, umf_int, umf_dp, umf_void
    implicit none
     integer(umf_int), intent(in) :: n_row
     integer(umf_int), intent(in) :: n_col
     integer(umf_int), intent(in) :: ap(*)
     integer(umf_int), intent(in) :: ai(*)
     real(umf_dp), intent(in) :: ax(*)
     integer(umf_void), intent(out) :: symbolic
     real(umf_dp), intent(in) :: control(umfpack_control)
     real(umf_dp), intent(out) :: info(umfpack_info)
   end subroutine umf4sym
 end interface

 interface 
   subroutine umf4num(ap,ai,ax,symbolic,numeric,control,info)
    import :: umfpack_control, umfpack_info, umf_int, umf_dp, umf_void
    implicit none
     integer(umf_int), intent(in) :: ap(*)
     integer(umf_int), intent(in) :: ai(*)
     real(umf_dp), intent(in) :: ax(*)
     integer(umf_void), intent(in) :: symbolic
     integer(umf_void), intent(out) :: numeric
     real(umf_dp), intent(in) :: control(umfpack_control)
     real(umf_dp), intent(out) :: info(umfpack_info)
   end subroutine umf4num
 end interface

 interface 
   subroutine umf4sol(sys,x,b,numeric,control,info)
    import :: umfpack_control, umfpack_info, umf_int, umf_dp, umf_void
    implicit none
     integer(umf_int), intent(in) :: sys
     real(umf_dp), intent(out) :: x(*)
     real(umf_dp), intent(in) :: b(*)
     integer(umf_void), intent(in) :: numeric
     real(umf_dp), intent(inout) :: control(umfpack_control)
     real(umf_dp), intent(out) :: info(umfpack_info)
   end subroutine umf4sol
 end interface
     
 interface 
   subroutine umf4solr(sys,ap,ai,ax,x,b,numeric,control,info)
    import :: umfpack_control, umfpack_info, umf_int, umf_dp, umf_void
    implicit none
     integer(umf_int), intent(in) :: sys
     integer(umf_int), intent(in) :: ap(*)
     integer(umf_int), intent(in) :: ai(*)
     real(umf_dp), intent(in) :: ax(*)
     real(umf_dp), intent(out) :: x(*)
     real(umf_dp), intent(in) :: b(*)
     integer(umf_void), intent(in) :: numeric
     real(umf_dp), intent(in) :: control(umfpack_control)
     real(umf_dp), intent(out) :: info(umfpack_info)
   end subroutine umf4solr
 end interface

 interface
   subroutine umf4pinf(control,info)
    import :: umfpack_control, umfpack_info, umf_dp
    implicit none
     real(umf_dp), intent(in) :: control(umfpack_control)
     real(umf_dp), intent(in) :: info(umfpack_info)
   end subroutine umf4pinf
 end interface

 character(len=*), parameter :: &
   this_mod_name = 'mod_umfpack'

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

 subroutine mod_umfpack_constructor()

  character(len=100) message(6)
  character(len=*), parameter :: &
    this_sub_name = 'constructor'

   !Consistency checks ---------------------------
!!$   if( (mod_messages_initialized.eqv..false.) ) then
!!$     call error(this_sub_name,this_mod_name, &
!!$                'Not all the required modules are initialized.')
!!$   endif
!!$   if(mod_umfpack_initialized.eqv..true.) then
!!$     call warning(this_sub_name,this_mod_name, &
!!$                  'Module is already initialized.')
!!$   endif
   !----------------------------------------------

   write(message(1),'(A,I3)')  'UMFPACK integers: umf_int = ', umf_int
   write(message(2),'(A,I20)') 'huge  = ',  huge(1_umf_int)
   write(message(3),'(A,I20)') 'range = ', range(1_umf_int)
   write(message(4),'(A,I3)')  'UMFPACK reals:     umf_dp = ', umf_dp
   write(message(5),'(A,I5)') 'precision = ', precision(1.0_umf_dp)
   write(message(6),'(A,I5)') 'range     = ',     range(1.0_umf_dp)
!!$   call info(this_sub_name,this_mod_name,message)

   mod_umfpack_initialized = .true.
 end subroutine mod_umfpack_constructor

!-----------------------------------------------------------------------
 
 subroutine mod_umfpack_destructor()
  character(len=*), parameter :: &
    this_sub_name = 'destructor'
   
   !Consistency checks ---------------------------
!!$   if(mod_umfpack_initialized.eqv..false.) then
!!$     call error(this_sub_name,this_mod_name, &
!!$                'This module is not initialized.')
!!$   endif
   !----------------------------------------------

   mod_umfpack_initialized = .false.
 end subroutine mod_umfpack_destructor

!-----------------------------------------------------------------------

end module mod_umfpack

