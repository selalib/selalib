!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

module sll_poisson_polar_parallel
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

  use sll_fft
  use sll_tridiagonal
  use sll_collective
  use sll_remapper
  use sll_boundary_condition_descriptors

  implicit none
  !>type sll_poisson_polar
  !>type for the Poisson solver in polar coordinate
  type sll_poisson_polar
   sll_real64                          :: dr, rmin, rmax
   sll_int32                           :: nr, nt
   sll_int32                           :: bc(2)
   type(sll_fft_plan), pointer         :: pfwd,pinv
   sll_comp64, dimension(:),   pointer :: fk,phik
   sll_real64, dimension(:),   pointer :: a,cts
   sll_int32,  dimension(:),   pointer :: ipiv
   type(layout_2D),  pointer           :: layout_r !< layout sequential in x
   type(layout_2D),  pointer           :: layout_t !< layout sequential in y
   type(remap_plan_2D_real64), pointer :: rmp_rt   !< remap r->theta 
   type(remap_plan_2D_real64), pointer :: rmp_tr   !< remap theta->r
   sll_real64, dimension(:,:), pointer :: f_r      !< array sequential in x
   sll_real64, dimension(:,:), pointer :: f_t      !< array sequential in y
  end type sll_poisson_polar

  interface initialize
     module procedure initialize_poisson_polar
  end interface initialize

  interface solve
     module procedure solve_poisson_polar
  end interface solve

  integer, private :: nr_loc
  integer, private :: nt_loc

contains

  !> Initialize the Poisson solver in polar coordinates
  subroutine initialize_poisson_polar(this, &
             layout_r, layout_t,        &
             rmin,rmax,nr,nt,bc_rmin,bc_rmax)

    implicit none
    type(sll_poisson_polar)  :: this
    type(layout_2D), pointer :: layout_r !< sequential in r direction
    type(layout_2D), pointer :: layout_t !< sequential in theta direction

    sll_real64               :: rmin    !< rmin
    sll_real64               :: rmax    !< rmax
    sll_int32                :: nr      !< number of cells radial
    sll_int32                :: nt  !< number of cells angular
    sll_int32, optional      :: bc_rmin !< radial boundary conditions
    sll_int32, optional      :: bc_rmax !< radial boundary conditions
    sll_int32                :: error
    sll_real64, dimension(:), allocatable :: buf

    SLL_ALLOCATE(this%fk(nr+1),error)
    SLL_ALLOCATE(this%phik(nr+1),error)
    SLL_ALLOCATE(this%a(3*(nr-1)),error)
    SLL_ALLOCATE(this%cts(7*(nr-1)),error)
    SLL_ALLOCATE(this%ipiv(nr-1),error)

    this%rmin = rmin
    this%rmax = rmax
    this%dr   = (rmax-rmin)/nr
    this%nr   = nr
    this%nt   = nt

    if (present(bc_rmin) .and. present(bc_rmax)) then
      this%bc(1)=bc_rmin
      this%bc(2)=bc_rmax
    else
      this%bc(1)=-1
      this%bc(2)=-1
    end if

    SLL_ALLOCATE(buf(nt),error)
    this%pfwd => fft_new_plan(nt,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(nt,buf,buf,FFT_INVERSE)
    SLL_DEALLOCATE_ARRAY(buf,error)

    ! Layout and local sizes for FFTs in r-direction
    this%layout_r => layout_r
    call compute_local_sizes_2d(this%layout_r,nr_loc,nt_loc)
    SLL_CLEAR_ALLOCATE(this%f_r(1:nr_loc,1:nt_loc),error)

    ! Layout and local sizes for FFTs in theta-direction
    this%layout_t => layout_t
    call compute_local_sizes_2d(this%layout_t,nr_loc,nt_loc)
    SLL_CLEAR_ALLOCATE(this%f_t(1:nr_loc,1:nt_loc),error)

    this%rmp_rt => new_remap_plan(this%layout_r, this%layout_t, this%f_r)
    this%rmp_tr => new_remap_plan(this%layout_t, this%layout_r, this%f_t)

  end subroutine initialize_poisson_polar

  !>delete a sll_poisson_polar object
  subroutine delete_poisson_polar(this)

    implicit none

    type(sll_poisson_polar), pointer :: this
    sll_int32 :: err
    if (associated(this)) then
       call fft_delete_plan(this%pfwd)
       call fft_delete_plan(this%pinv)
       SLL_DEALLOCATE_ARRAY(this%fk,err)
       SLL_DEALLOCATE_ARRAY(this%phik,err)
       SLL_DEALLOCATE_ARRAY(this%a,err)
       SLL_DEALLOCATE_ARRAY(this%cts,err)
       SLL_DEALLOCATE_ARRAY(this%ipiv,err)
       SLL_DEALLOCATE(this,err)
    end if

  end subroutine delete_poisson_polar

!===================
!  Poisson solver
!===================

  !>poisson solver for polar system : -\Delta (phi)=rhs
  !>initialization must be done outside the solver
  subroutine solve_poisson_polar(this,rhs,phi)

    implicit none

    type(sll_poisson_polar) :: this !< contains data for the solver
    !>rhs field, size (nr+1)*(nt+1), input
    sll_real64, dimension(this%nr+1,this%nt+1), intent(in)  :: rhs
    !>unknown field, size (nr+1)*(nt+1), output
    sll_real64, dimension(this%nr+1,this%nt+1), intent(out) :: phi

    sll_real64 :: rmin,dr
    sll_int32  :: nr, nt,bc(2)

    sll_real64 :: r
    sll_int32  :: i, k
    sll_real64 :: kval

    nr   = this%nr
    nt   = this%nt
    rmin = this%rmin
    dr   = this%dr
    bc   = this%bc

    call verify_argument_sizes_par(this%layout_t, rhs, phi)

    this%f_t = rhs
    
    call compute_local_sizes_2d( this%layout_t, nr_loc, nt_loc )
    
    do i=1,nr_loc+1
      call fft_apply_plan(this%pfwd,this%f_t(i,1:nt),this%f_t(i,1:nt))
    end do

    call apply_remap_2D( this%rmp_tr, this%f_t, this%f_r )
    call compute_local_sizes_2d( this%layout_r, nr_loc, nt_loc )

    do k = 0,nt_loc/2

      kval=real(k,f64)

      do i=2,nr
        r=rmin+(i-1)*dr
        this%a(3*(i-1)  )=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
        this%a(3*(i-1)-1)= 2.0_f64/dr**2+(kval/r)**2
        this%a(3*(i-1)-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
        this%fk(i)=fft_get_mode(this%pfwd,this%f_r(i,1:nt),k)
      enddo

      this%phik=0.0_f64

      if(bc(1)==SLL_DIRICHLET) then
        this%a(1)=0.0_f64
      else if(bc(1)==SLL_NEUMANN) then
        this%a(2)=this%a(2)+this%a(1) 
        this%a(1)=0._f64
      else if(bc(1)==SLL_NEUMANN_MODE0)then 
        if(k==0) then
          this%a(2)=this%a(2)+this%a(1)
          this%a(1)=0._f64
        else 
          this%a(1)=0._f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==SLL_DIRICHLET) then
        this%a(3*(nr-1))=0.0_f64
      else if(bc(2)==SLL_NEUMANN)then
        this%a(3*(nr-1)-1)=this%a(3*(nr-1)-1)+this%a(3*(nr-1)) 
        this%a(3*(nr-1))=0.0_f64
      else if(bc(2)==SLL_NEUMANN_MODE0)then 
        if(k==0)then
          this%a(3*(nr-1)-1)=this%a(3*(nr-1)-1)+this%a(3*(nr-1))
          this%a(3*(nr-1))=0.0_f64
        else 
          this%a(3*(nr-1))=0.0_f64
        endif
      endif

      call setup_cyclic_tridiag(this%a,nr-1,this%cts,this%ipiv)
      call solve_cyclic_tridiag(this%cts,this%ipiv,this%fk(2:nr),nr-1,this%phik(2:nr))

      !boundary condition at rmin
      if(bc(1)==1)then
        this%phik(1)=0.0_f64
      else if(bc(1)==2)then
        this%phik(1)=this%phik(2) 
      else if(bc(1)==3)then 
        if(k==0)then
          this%phik(1)=this%phik(2)
        else 
          this%phik(1)=0.0_f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==1)then
        this%phik(nr+1)=0.0_f64
      else if(bc(2)==2)then
        this%phik(nr+1)=this%phik(nr)
      else if(bc(2)==3)then 
        if(k==0)then
          this%phik(nr+1)=this%phik(nr)
        else 
          this%phik(nr+1)=0.0_f64
        endif
      endif

      do i=1,nr+1
        call fft_set_mode(this%pinv,phi(i,1:nt),this%phik(i),k)
      end do

    end do

    call apply_remap_2D( this%rmp_rt, this%f_r, this%f_t )
    call compute_local_sizes_2d( this%layout_t, nr_loc, nt_loc )
    do i=1,nr_loc+1
      call fft_apply_plan(this%pinv,phi(i,1:nt),phi(i,1:nt))
      phi(i,nt+1)=phi(i,1)
    end do


  end subroutine solve_poisson_polar


end module sll_poisson_polar_parallel
