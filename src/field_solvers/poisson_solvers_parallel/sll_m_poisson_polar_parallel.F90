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

!> @ingroup poisson_solvers
!> Module to sll_o_solve Poisson equation on polar mesh using FFT transform
module sll_m_poisson_polar_parallel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_neumann_mode_0

  use sll_m_collective, only: &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_fft, only: &
    sll_s_fft_apply_plan_c2c_1d, &
    sll_p_fft_backward, &
    sll_s_fft_delete_plan, &
    sll_p_fft_forward, &
    sll_f_fft_new_plan_c2c_1d, &
    sll_t_fft_plan

  use sll_m_remapper, only: &
    sll_o_apply_remap_2d, &
    sll_o_compute_local_sizes, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_2d_comp64

  use sll_m_tridiagonal, only: &
    sll_s_setup_cyclic_tridiag, &
    sll_o_solve_cyclic_tridiag

  implicit none

  public :: &
    sll_o_initialize, &
    sll_f_new_poisson_polar, &
    sll_t_poisson_polar, &
    sll_o_solve, &
    sll_s_solve_poisson_polar

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>Class for the Poisson solver in polar coordinate
  type sll_t_poisson_polar
   sll_real64                          :: rmin     !< left corner of r dimension
   sll_real64                          :: rmax     !< right corner of r dimension
   sll_real64                          :: dr       !< step size along r
   sll_int32                           :: nr       !< number of nodes along r
   sll_int32                           :: nt       !< number of nodes along theta
   sll_int32                           :: bc(2)    !< boundary conditions options
   type(sll_t_fft_plan), pointer         :: fw       !< Forward FFT plan
   type(sll_t_fft_plan), pointer         :: bw       !< Inverse FFT plan
   sll_comp64, dimension(:),   pointer :: fk       !< RHSf fft
   sll_comp64, dimension(:),   pointer :: phik     !< Potential fft
   sll_real64, dimension(:),   pointer :: mat      !< Matrix
   sll_real64, dimension(:),   pointer :: cts      !< Lapack coefficient
   sll_int32,  dimension(:),   pointer :: ipiv     !< Lapack pivot indices
   type(sll_t_layout_2d),  pointer           :: layout_r !< layout sequential in r
   type(sll_t_layout_2d),  pointer           :: layout_a !< layout sequential in theta
   type(sll_t_remap_plan_2d_comp64), pointer :: rmp_ra   !< remap r->theta 
   type(sll_t_remap_plan_2d_comp64), pointer :: rmp_ar   !< remap theta->r
   sll_comp64, dimension(:,:), pointer :: f_r      !< array sequential in r
   sll_comp64, dimension(:,:), pointer :: f_a      !< array sequential in theta
   sll_real64, dimension(:), pointer :: dlog_density,inv_Te !<for quasi neutral solver

  end type sll_t_poisson_polar

  !> sll_o_initialize the Poisson solver on polar mesh
  interface sll_o_initialize
     module procedure initialize_poisson_polar
  end interface sll_o_initialize

  !> Compute the potential solving the Poisson equation on polar mesh
  interface sll_o_solve
     module procedure sll_s_solve_poisson_polar
  end interface sll_o_solve

contains
  !> new: wrapper of sll_o_initialize
  function sll_f_new_poisson_polar( &
    layout_r, &
    layout_a, &
    rmin, &
    rmax, &
    nr, &
    ntheta, &
    bc_rmin, &
    bc_rmax, &
    dlog_density, &
    inv_Te) &
    result(this)
    implicit none
    type(sll_t_poisson_polar), pointer  :: this     !< Poisson solver class
    type(sll_t_layout_2d), pointer :: layout_r !< sequential in r direction
    type(sll_t_layout_2d), pointer :: layout_a !< sequential in theta direction


    sll_real64               :: rmin    !< rmin
    sll_real64               :: rmax    !< rmax
    sll_int32                :: nr      !< number of cells radial
    sll_int32                :: ntheta  !< number of cells angular
    sll_int32, optional      :: bc_rmin !< radial boundary conditions
    sll_int32, optional      :: bc_rmax !< radial boundary conditions
    sll_real64,dimension(:),optional :: dlog_density !< for quasi neutral solver
    sll_real64,dimension(:),optional :: inv_Te       !< for quasi neutral solver
    
    !local variables
    sll_int32 :: ierr
    
    SLL_ALLOCATE(this,ierr)
   
    call initialize_poisson_polar( &
      this, &
      layout_r, &
      layout_a, &
      rmin, &
      rmax, &
      nr, &
      ntheta, &
      bc_rmin, &
      bc_rmax, &
      dlog_density, &
      inv_Te)

  end function sll_f_new_poisson_polar

  !> sll_o_initialize the Poisson solver in polar coordinates
  subroutine initialize_poisson_polar(this, layout_r, layout_a, &
             rmin,rmax,nr,ntheta,bc_rmin,bc_rmax,dlog_density,inv_Te)

    implicit none
    type(sll_t_poisson_polar)  :: this     !< Poisson solver class
    type(sll_t_layout_2d), pointer :: layout_r !< sequential in r direction
    type(sll_t_layout_2d), pointer :: layout_a !< sequential in theta direction


    sll_real64               :: rmin    !< rmin
    sll_real64               :: rmax    !< rmax
    sll_int32                :: nr      !< number of cells radial
    sll_int32                :: ntheta  !< number of cells angular
    sll_int32, optional      :: bc_rmin !< radial boundary conditions
    sll_int32, optional      :: bc_rmax !< radial boundary conditions
    sll_real64,dimension(:),optional :: dlog_density !< for quasi neutral solver
    sll_real64,dimension(:),optional :: inv_Te       !< for quasi neutral solver

    sll_int32                :: error
    sll_comp64, dimension(:), allocatable :: buf
    sll_int32 :: nr_loc
    sll_int32 :: na_loc
    sll_int32 :: psize


    SLL_ALLOCATE(this%f_r(nr+1,ntheta+1),error)
    SLL_ALLOCATE(this%fk(nr+1),error)
    SLL_ALLOCATE(this%phik(nr+1),error)
    SLL_ALLOCATE(this%mat(3*(nr-1)),error)
    SLL_ALLOCATE(this%cts(7*(nr-1)),error)
    SLL_ALLOCATE(this%ipiv(nr-1),error)
    SLL_ALLOCATE(this%dlog_density(nr+1),error)
    SLL_ALLOCATE(this%inv_Te(nr+1),error)

    this%dlog_density = 0._f64
    this%inv_Te = 0._f64
    
    if(present(dlog_density))then
      this%dlog_density = dlog_density
    endif
    if(present(inv_Te))then
      this%inv_Te = inv_Te
    endif



    this%rmin=rmin
    this%rmax=rmax
    this%dr=(rmax-rmin)/nr
    this%nr=nr
    this%nt=ntheta

    if (present(bc_rmin) .and. present(bc_rmax)) then
      this%bc(1)=bc_rmin
      this%bc(2)=bc_rmax
    else
      this%bc(1)=-1
      this%bc(2)=-1
    end if

    SLL_ALLOCATE(buf(ntheta),error)
    this%fw => sll_f_fft_new_plan_c2c_1d(ntheta,buf,buf,sll_p_fft_forward,normalized=.true.)
    this%bw => sll_f_fft_new_plan_c2c_1d(ntheta,buf,buf,sll_p_fft_backward)
    SLL_DEALLOCATE_ARRAY(buf,error)

    psize = sll_f_get_collective_size(sll_v_world_collective)

    this%layout_r => layout_r
    call sll_o_compute_local_sizes(layout_r,nr_loc,na_loc)
    SLL_ALLOCATE(this%f_r(1:nr_loc,1:na_loc),error)
    this%f_r = (0.0_f64,0.0_f64)

    ! Layout and local sizes for FFTs in theta-direction
    this%layout_a => layout_a
    call sll_o_compute_local_sizes(layout_a,nr_loc,na_loc)
    SLL_ALLOCATE(this%f_a(1:nr_loc,1:na_loc),error)
    this%f_a = (0.0_f64,0.0_f64)

    this%rmp_ra => sll_o_new_remap_plan(this%layout_r, this%layout_a, this%f_r)
    this%rmp_ar => sll_o_new_remap_plan(this%layout_a, this%layout_r, this%f_a)



  end subroutine initialize_poisson_polar

  !>delete a sll_t_poisson_polar object
  subroutine delete_poisson_polar(this)

    implicit none

    type(sll_t_poisson_polar), pointer :: this !< Poisson solver object
    sll_int32 :: err
    if (associated(this)) then
       call sll_s_fft_delete_plan(this%fw)
       deallocate(this%fw)
       call sll_s_fft_delete_plan(this%bw)
       deallocate(this%bw)
       SLL_DEALLOCATE_ARRAY(this%fk,err)
       SLL_DEALLOCATE_ARRAY(this%phik,err)
       SLL_DEALLOCATE_ARRAY(this%mat,err)
       SLL_DEALLOCATE_ARRAY(this%cts,err)
       SLL_DEALLOCATE_ARRAY(this%ipiv,err)
       SLL_DEALLOCATE(this,err)
    end if

  end subroutine delete_poisson_polar


  !> sll_o_solve the Poisson equation and get the potential
  subroutine sll_s_solve_poisson_polar(this,rhs,phi)

    implicit none

    type(sll_t_poisson_polar) :: this !< Poisson solver object
    sll_real64, dimension(:,:), intent(in)  :: rhs !< Charge density
    sll_real64, dimension(:,:), intent(out) :: phi !< Potential

    sll_real64 :: rmin,dr
    sll_int32  :: nr, ntheta,bc(2)

    sll_real64 :: r
    sll_int32  :: i, j, k
    sll_int32  :: nr_loc
    sll_int32  :: na_loc
    sll_int32  :: global(2)

    nr     = this%nr
    ntheta = this%nt
    rmin   = this%rmin
    dr     = this%dr
    bc     = this%bc

    call verify_argument_sizes_par(this%layout_a, rhs)
    this%f_a = cmplx(rhs, 0_f64, kind=f64)

    call sll_o_compute_local_sizes( this%layout_a, nr_loc, na_loc )

    do i=1,nr_loc
      call sll_s_fft_apply_plan_c2c_1d(this%fw,this%f_a(i,1:ntheta),this%f_a(i,1:ntheta))
    end do

    
    call sll_o_apply_remap_2d( this%rmp_ar, this%f_a, this%f_r )
    
    
    call sll_o_compute_local_sizes( this%layout_r, nr_loc, na_loc )
    !global = sll_o_local_to_global( this%layout_r, (/1, 1/))
    !k = global(2)/2 - 1
    
    !print *,'#na_loc=',sll_f_get_collective_rank(sll_v_world_collective),na_loc
    !do i=1,na_loc
    !  print *,sll_f_get_collective_rank(sll_v_world_collective),i,this%f_r(nr_loc/2,i)
    !enddo
    
    !do j = 1, na_loc-1, 2
    do j = 1, na_loc  
      !k = k + 1
      global = sll_o_local_to_global( this%layout_r, (/1, j/))
      k = global(2)
      if (k<=ntheta/2) then
        k = k-1
      else
        k = ntheta-(k-1)
      endif
      
      !print *,'#k=',k,sll_f_get_collective_rank(sll_v_world_collective)

      do i=2,nr

        r = rmin + (i-1)*dr
        this%mat(3*(i-1)  ) = -1.0_f64/dr**2-1.0_f64/(2._f64*dr*r) &
          -this%dlog_density(i)/(2._f64*dr)
        this%mat(3*(i-1)-1) =  2.0_f64/dr**2+(k/r)**2 &
          +this%inv_Te(i)
        this%mat(3*(i-1)-2) = -1.0_f64/dr**2+1.0_f64/(2*dr*r) &
          +this%dlog_density(i)/(2._f64*dr)

        this%fk(i) = this%f_r(i,j)
        !if( j == 1 ) then
        !   this%fk(i) = cmplx(this%f_r(i,1),0.0_f64,kind=f64)
        !else if( j == ntheta ) then
        !   this%fk(i) = cmplx(this%f_r(i,2),0.0_f64,kind=f64)
        !else
        !   this%fk(i) = cmplx(this%f_r(i,j), &
        !                      this%f_r(i,j+1),kind=f64)
        !endif

      enddo

      this%phik=(0.0_f64,0.0_f64)

      !boundary condition at rmin
      if (bc(1)==sll_p_dirichlet) then !Dirichlet
        this%mat(1)=0.0_f64
      else if (bc(1)==sll_p_neumann) then
        this%mat(2)=this%mat(2)+this%mat(1) !Neumann
        this%mat(1)=0._f64
      else if(bc(1)==sll_p_neumann_mode_0)then 
        if(k==0)then!Neumann for mode zero
          this%mat(2)=this%mat(2)+this%mat(1)
          this%mat(1)=0._f64
        else !Dirichlet for other modes
          this%mat(1)=0._f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==sll_p_dirichlet)then !Dirichlet
        this%mat(3*(nr-1))=0.0_f64
      else if(bc(2)==sll_p_neumann)then
        this%mat(3*(nr-1)-1)=this%mat(3*(nr-1)-1)+this%mat(3*(nr-1)) !Neumann
        this%mat(3*(nr-1))=0.0_f64
      else if(bc(2)==sll_p_neumann_mode_0)then 
        if(k==0)then!Neumann for mode zero
          this%mat(3*(nr-1)-1)=this%mat(3*(nr-1)-1)+this%mat(3*(nr-1))
          this%mat(3*(nr-1))=0.0_f64
        else !Dirichlet for other modes
          this%mat(3*(nr-1))=0.0_f64
        endif
      endif

      call sll_s_setup_cyclic_tridiag(this%mat,nr-1,this%cts,this%ipiv)
      call sll_o_solve_cyclic_tridiag(this%cts,this%ipiv,this%fk(2:nr),nr-1,this%phik(2:nr))

      !boundary condition at rmin
      if(bc(1)==sll_p_dirichlet)then !Dirichlet
        this%phik(1)=(0.0_f64,0.0_f64)
      else if (bc(1)==sll_p_neumann) then
        this%phik(1)=this%phik(2) !Neumann
      else if (bc(1)==sll_p_neumann_mode_0) then 
        if (k==0) then!Neumann for mode zero
          this%phik(1)=this%phik(2)
        else !Dirichlet for other modes
          this%phik(1)=(0.0_f64,0.0_f64)
        endif
      endif

      !boundary condition at rmax
      if (bc(2)==sll_p_dirichlet) then !Dirichlet
        this%phik(nr+1)=(0.0_f64,0.0_f64)
      else if (bc(2)==sll_p_neumann) then
        this%phik(nr+1)=this%phik(nr) !Neumann
      else if (bc(2)==sll_p_neumann_mode_0) then 
        if(k==0)then!Neumann for mode zero
          this%phik(nr+1)=this%phik(nr)
        else !Dirichlet for other modes
          this%phik(nr+1)=(0.0_f64,0.0_f64)
        endif
      endif

      do i=1,nr+1
        this%f_r(i,j) = this%phik(i)   
         !if( j == 1 ) then
         !   this%f_r(i,1) = real(this%phik(i),kind=f64)
         !else if( j == ntheta ) then
         !   this%f_r(i,2) = real(this%phik(i),kind=f64)
         !else
         !   this%f_r(i,j) = real(this%phik(i),kind=f64)
         !   this%f_r(i,j+1) = aimag(this%phik(i))
         !endif

      end do

    end do

    call sll_o_apply_remap_2d( this%rmp_ra, this%f_r, this%f_a )
    call sll_o_compute_local_sizes( this%layout_a, nr_loc, na_loc )

    call verify_argument_sizes_par(this%layout_a, phi)
    
    do i=1,nr_loc
      !call fft_apply_plan(this%bw,this%f_a(i,1:ntheta),phi(i,1:ntheta))
      call sll_s_fft_apply_plan_c2c_1d(this%bw,this%f_a(i,1:ntheta),this%f_a(i,1:ntheta))
      phi(i,1:ntheta) = real(this%f_a(i,1:ntheta))
    end do

  end subroutine sll_s_solve_poisson_polar


  !> Check if array sizes are compatble with the layout 
  subroutine verify_argument_sizes_par(layout, array)

    type(sll_t_layout_2d), pointer       :: layout
    sll_real64, dimension(:,:)     :: array
    sll_int32,  dimension(2)       :: n ! nx_loc, ny_loc
    sll_int32                      :: i

    ! Note that this checks for strict sizes, not an array being bigger
    ! than a certain size, but exactly a desired size... This may be a bit
    ! too stringent.
    call sll_o_compute_local_sizes( layout, n(1), n(2) )

    do i=1,2
       if ( (n(i)/=size(array,i))) then
          print*, 'ERROR: solve_poisson_polar_parallel()', &
               'size of either rhs or phi does not match expected size. '
          if (i==1) then
             print*, 'solve_poisson_polar_parallel(): ', &
                  'mismatch in direction r'
          else if (i==2) then
             print*, 'solve_poisson_polar_parallel(): ', &
                  'mismatch in direction theta'
          endif
          print *, 'Exiting...'
          call sll_s_halt_collective()
          stop
       endif
    enddo

  end subroutine verify_argument_sizes_par

end module sll_m_poisson_polar_parallel
