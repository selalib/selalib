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

!> @author Eric MADAULE
!> @brief Poisson equation solver in polar coordinate
!> @details Solver for the Poisson equation 
!> \f[ \Delta \phi = f \f]
!> in polar coordinate
!> using a fft in direction \f$ \theta \f$ and final differencies 
!> in direction \f$ r \f$.
!> This way we solve a tridiagonal system with the solver from SELALIB.
!>
!>\section how How to use the Poisson polar solver?
!>
!>You must add \code sll_poisson_2d_polar \endcode to the list of linked libraries.
!>The Poisson solver uses the FFT, so you also need to link the FFT
!>
!>1. Declare a Poisson polar plan
!>\code type(sll_plan_poisson_polar), pointer :: plan \endcode
!>2. Initialize the plan
!>\code plan => new_plan_poisson_polar(dr,rmin,nr,ntheta,boundary_conditions) \endcode
!>nr and ntheta are the number of step in direction r and theta
!> 3. Execute the plan
!> \code call poisson_solve_polar(plan,in,out) \endcode
!> 4. Delete the plan
!> \code call delete_plan_poisson_polar(plan) \endcode
!>
!>\section bc Boundary conditions :
!>
!>The boundary conditions define the potential behaviour in r_{min} (BOT_) and r_{max} (TOP_)
!>They can take the value DIRICHLET or NEUMANN
!>
!>Summary
!>The boundary conditions define the potential behaviour in r_{min} (BOT_) and r_{max} (TOP_)
!>They can take the value DIRICHLET or NEUMANN
!>
!>Summary :
!>
!>The different boundary conditions are :
!> - BOT_DIRICHLET
!> - BOT_NEUMANN
!> - TOP_DIRICHLET
!> - TOP_NEUMANN
!>
!>You must combine BOT_ and TOP_ conditions with '+'.
!>
!>\section examples Example :
!>
!>Full code :
!>\code
!>integer :: nr, ntheta
!>real :: dr, rmin
!>integer :: bc
!>type(sll_plan_poisson_polar), pointer :: plan
!>real, dimension(:,:), allocatable :: in, out
!>
!>!define all parameters
!>nr     = 100
!>ntheta = 64 !ntheta must be a power of 2
!>rmin   = 1.0
!>dr     = 0.05
!>bc     = BOT_NEUMANN+TOP_DIRICHLET
!>
!>allocate(in(nr+1,ntheta+1))
!>allocate(out(nr+1,ntheta+1))
!>
!>in= !definition of in
!>
!>!initialization of plan
!>plan => new_plan_poisson_polar(dr,rmin,nr,ntheta,bc)
!>
!>!computation of Poisson
!>call poisson_solve_polar(plan,in,out)
!>
!>!deletion of plan
!>call delete_plan_poisson_polar(plan)
!>\endcode

module sll_poisson_2d_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use sll_tridiagonal
  use sll_constants
  use sll_boundary_condition_descriptors

  implicit none
  !>type sll_plan_poisson_polar
  !>type for the Poisson solver in polar coordinate
  type sll_plan_poisson_polar
     sll_real64                          :: rmin   !< r min
     sll_real64                          :: rmax   !< r max
     sll_real64                          :: dr     !< step size
     sll_int32                           :: nr     !< number of points in r
     sll_int32                           :: ntheta !< number of points in theta
     sll_int32                           :: bc(2)  !< boundary conditon type
     type(sll_fft_plan), pointer         :: pfwd   !< fft plan in theta
     type(sll_fft_plan), pointer         :: pinv   !< inverse fft plan in theta
     sll_real64, dimension(:,:), pointer :: f_fft  !< potential fft in theta
     sll_comp64, dimension(:),   pointer :: fk     !< \f$ f_k \f$
     sll_comp64, dimension(:),   pointer :: phik   !< \f$ phi_k \f$
     sll_real64, dimension(:), pointer   :: a      !< data for the tridiagonal solver
     sll_real64, dimension(:), pointer   :: cts    !< lapack array
     sll_int32, dimension(:),  pointer   :: ipiv   !< lapack pivot data
     sll_real64, dimension(:), pointer ::dlog_density,inv_Te !<for quasi neutral solver


  end type sll_plan_poisson_polar

  !> Initialize the polar poisson solver
  interface initialize
     module procedure initialize_poisson_polar
  end interface initialize

  !> Get potential from the polar poisson solver
  interface solve
     module procedure solve_poisson_polar
  end interface solve

contains

!> Creation of sll_plan_poisson_polar object for the 
!> Poisson solver in polar coordinate
  function new_plan_poisson_polar(dr,rmin,nr,ntheta,bc,dlog_density,inv_Te) result(this)

    implicit none

    sll_real64 :: dr             !< size of space in direction r
    sll_real64 :: rmin           !< interior radius
    sll_int32  :: nr             !< number of space in direction r
    sll_int32  :: ntheta         !< number of space in direction theta
    sll_int32, optional :: bc(2) !< Boundary conditions, can be combined with +
                                 !< bc is optionnal and default is Dirichlet condition in rmin and rmax
    type(sll_plan_poisson_polar), pointer :: this !< Poisson solver structure
    sll_real64,dimension(:),optional ::dlog_density,inv_Te !< for quasi neutral solver

    sll_int32 :: err
    sll_real64, dimension(:), allocatable :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fk(nr+1),err)
    SLL_ALLOCATE(this%phik(nr+1),err)
!!$    SLL_ALLOCATE(this%a(3*(nr+1)),err)
!!$    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
!!$    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%a(3*(nr-1)),err)
    SLL_ALLOCATE(this%cts(7*(nr-1)),err)
    SLL_ALLOCATE(this%ipiv(nr-1),err)
    
    SLL_ALLOCATE(this%dlog_density(nr+1),err)
    SLL_ALLOCATE(this%inv_Te(nr+1),err)
    
    this%dlog_density = 0._f64
    this%inv_Te = 0._f64
    
    if(present(dlog_density))then
      this%dlog_density = dlog_density
    endif
    if(present(inv_Te))then
      this%inv_Te = inv_Te
    endif
    
    
    
    this%dr=dr
    this%rmin=rmin
    this%nr=nr
    this%ntheta=ntheta

    if (present(bc)) then
      this%bc=bc
    else
      this%bc(1)=-1
      this%bc(2)=-1
    end if

 

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)
    
    
    
    SLL_DEALLOCATE_ARRAY(buf,err)
    
    

  end function new_plan_poisson_polar

  !> Initialize the Poisson solver in polar coordinates
  subroutine initialize_poisson_polar(this, rmin,rmax,nr,ntheta,bc_rmin,bc_rmax,dlog_density,inv_Te)

    implicit none
    type(sll_plan_poisson_polar) :: this !< Poisson solver structure

    sll_real64               :: rmin     !< rmin
    sll_real64               :: rmax     !< rmax
    sll_int32                :: nr       !< number of cells radial
    sll_int32                :: ntheta   !< number of cells angular
    sll_int32, optional      :: bc_rmin  !< radial boundary conditions
    sll_int32, optional      :: bc_rmax  !< radial boundary conditions
    sll_int32                :: error
    sll_real64, dimension(:), allocatable :: buf
    sll_real64,dimension(:),optional ::dlog_density,inv_Te

    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),error)
    SLL_ALLOCATE(this%fk(nr+1),error)
    SLL_ALLOCATE(this%phik(nr+1),error)
    SLL_ALLOCATE(this%a(3*(nr-1)),error)
    SLL_ALLOCATE(this%cts(7*(nr-1)),error)
    SLL_ALLOCATE(this%ipiv(nr-1),error)

    this%rmin=rmin
    this%rmax=rmax
    this%dr=(rmax-rmin)/nr
    this%nr=nr
    this%ntheta=ntheta

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



    if (present(bc_rmin) .and. present(bc_rmax)) then
      this%bc(1)=bc_rmin
      this%bc(2)=bc_rmax
    else
      this%bc(1)=-1
      this%bc(2)=-1
    end if

    SLL_ALLOCATE(buf(ntheta),error)
    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)
    SLL_DEALLOCATE_ARRAY(buf,error)

  end subroutine initialize_poisson_polar

!======================================
! deletion of sll_plan_poisson_polar
!======================================

  !>delete_plan_poisson_polar(plan)
  !>delete a sll_plan_poisson_polar object
  subroutine delete_plan_poisson_polar(this)

    implicit none

    type(sll_plan_poisson_polar), pointer :: this
    sll_int32 :: err
    if (associated(this)) then
       call fft_delete_plan(this%pfwd)
       call fft_delete_plan(this%pinv)
       SLL_DEALLOCATE_ARRAY(this%f_fft,err)
       SLL_DEALLOCATE_ARRAY(this%fk,err)
       SLL_DEALLOCATE_ARRAY(this%phik,err)
       SLL_DEALLOCATE_ARRAY(this%a,err)
       SLL_DEALLOCATE_ARRAY(this%cts,err)
       SLL_DEALLOCATE_ARRAY(this%ipiv,err)
       SLL_DEALLOCATE(this,err)
    end if

  end subroutine delete_plan_poisson_polar

!===================
!  Poisson solver
!===================

  !>subroutine solve_poisson_polar(plan,f,phi)
  !>poisson solver for polar system : -\Delta (phi)=f
  !>plan : sll_plan_poisson_polar, contains data for the solver
  !>f : distribution function, size (nr+1)*(ntheta+1), input
  !>phi : unknown field, size (nr+1)*(ntheta+1), output
  !>initialization must be done outside the solver
  subroutine solve_poisson_polar(plan,f,phi)

    implicit none

    type(sll_plan_poisson_polar) :: plan
    sll_real64, dimension(:,:), intent(in)  :: f
    sll_real64, dimension(:,:), intent(out) :: phi

    sll_real64 :: rmin,dr
    sll_int32  :: nr, ntheta,bc(2)

    sll_real64 :: r
    sll_int32  :: i, k, ind_k
    sll_real64 :: kval

    sll_comp64 :: err_loc
    !sll_int32  :: ierr_sup_1em12
    sll_real64 :: err

    nr     = plan%nr
    ntheta = plan%ntheta
    rmin   = plan%rmin
    dr     = plan%dr
    
    !print *,'#nr=',nr,ntheta
    !stop

    !do k=1,ntheta+1!/2
    !  print *,k,f(2,k)!,plan%f_fft(2,k)!fft_get_mode(plan%pfwd,plan%f_fft(2,1:ntheta),k)
    !enddo
    
    !stop



    bc         = plan%bc
    plan%f_fft = f

   

    do i=1,nr+1
      call fft_apply_plan(plan%pfwd,plan%f_fft(i,1:ntheta),plan%f_fft(i,1:ntheta))
    end do

    !do k=0,ntheta/2
    !  print *,k,fft_get_mode(plan%pfwd,plan%f_fft(2,1:ntheta),k)
    !enddo
    
    
    
    

    ! poisson solver
    do k = 0,ntheta/2

      ind_k=k

      kval=real(ind_k,f64)

      do i=2,nr
        r = rmin + (i-1)*dr
        plan%a(3*(i-1)  ) = -1.0_f64/dr**2-1.0_f64/(2._f64*dr*r)-plan%dlog_density(i)/(2._f64*dr)
        !plan%a(3*(i-1)  ) = -1.0_f64/dr**2-1.0_f64/(2*dr*r)
        plan%a(3*(i-1)-1) =  2.0_f64/dr**2+(kval/r)**2+plan%inv_Te(i)
        plan%a(3*(i-1)-2) = -1.0_f64/dr**2+1.0_f64/(2._f64*dr*r)+plan%dlog_density(i)/(2._f64*dr)
        
        !print *,'before1'
        
        plan%fk(i)=fft_get_mode(plan%pfwd,plan%f_fft(i,1:ntheta),k)

        !print *,'after1'

      enddo
      
      
      
      !print *,k,maxval(abs(plan%fk))
      
      plan%phik=0.0_f64

      !boundary condition at rmin
      if(bc(1)==SLL_DIRICHLET)then !Dirichlet
        plan%a(1)=0.0_f64
      endif
      if(bc(1)==SLL_NEUMANN)then
        plan%a(2)=plan%a(2)+plan%a(1) !Neumann
        plan%a(1)=0._f64
      endif
      if(bc(1)==SLL_NEUMANN_MODE_0)then 
        if(k==0)then!Neumann for mode zero
          plan%a(2)=plan%a(2)+plan%a(1)
          plan%a(1)=0._f64
        else !Dirichlet for other modes
          plan%a(1)=0._f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==SLL_DIRICHLET)then !Dirichlet
        plan%a(3*(nr-1))=0.0_f64
      endif
      if(bc(2)==SLL_NEUMANN)then
        plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1)) !Neumann
        plan%a(3*(nr-1))=0.0_f64
      endif
      if(bc(2)==SLL_NEUMANN_MODE_0)then 
        if(k==0)then!Neumann for mode zero
          plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1))
          plan%a(3*(nr-1))=0.0_f64
        else !Dirichlet for other modes
          plan%a(3*(nr-1))=0.0_f64
        endif
      endif

      call setup_cyclic_tridiag(plan%a,nr-1,plan%cts,plan%ipiv)
      call solve_cyclic_tridiag(plan%cts,plan%ipiv,plan%fk(2:nr),nr-1,plan%phik(2:nr))

      !boundary condition at rmin
      if(bc(1)==1)then !Dirichlet
        plan%phik(1)=0.0_f64
      endif
      if(bc(1)==2)then
        plan%phik(1)=plan%phik(2) !Neumann
      endif
      if(bc(1)==3)then 
        if(k==0)then!Neumann for mode zero
          plan%phik(1)=plan%phik(2)
        else !Dirichlet for other modes
          plan%phik(1)=0.0_f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==1)then !Dirichlet
        plan%phik(nr+1)=0.0_f64
      endif
      if(bc(2)==2)then
        plan%phik(nr+1)=plan%phik(nr) !Neumann
      endif
      if(bc(2)==3)then 
        if(k==0)then!Neumann for mode zero
          plan%phik(nr+1)=plan%phik(nr)
        else !Dirichlet for other modes
          plan%phik(nr+1)=0.0_f64
        endif
      endif

      do i=1,nr+1
        call fft_set_mode(plan%pinv,phi(i,1:ntheta),plan%phik(i),k)
      end do
    end do


      err = 0._f64
      do i=4,nr-4
        r=rmin+real(i-1,f64)*dr
        err_loc=(plan%phik(i+1)-2*plan%phik(i)+plan%phik(i-1))/dr**2
        err_loc=err_loc-plan%phik(i)*plan%inv_Te(i)
        err_loc=err_loc+(plan%phik(i+1)-plan%phik(i-1))/(2._f64*r*dr)
        err_loc=err_loc+(plan%phik(i+1)-plan%phik(i-1))/(2._f64*dr)*plan%dlog_density(i)
        err_loc=-err_loc+kval**2/r**2*plan%phik(i)
        err_loc=(err_loc-plan%fk(i))
        if(abs(err_loc)>err)then
          err=abs(err_loc)
        endif
      enddo
      
      if(err>1.e-12)then 
        print *,'#err for QNS=',err 
      endif
    ! FFT INVERSE
    do i=1,nr+1
      call fft_apply_plan(plan%pinv,phi(i,1:ntheta),phi(i,1:ntheta))
    end do

    phi(:,ntheta+1)=phi(:,1)
    
    
    
    

  end subroutine solve_poisson_polar

  !>subroutine poisson_solve_polar(plan,f,phi)
  !>poisson solver for polar system : -\Delta (phi)=f
  !>plan : sll_plan_poisson_polar, contains data for the solver
  !>f : distribution function, size (nr+1)*(ntheta+1), input
  !>phi : unknown field, size (nr+1)*(ntheta+1), output
  !>initialization must be done outside the solver
  subroutine poisson_solve_polar(plan,f,phi,ierr)

    implicit none

    type(sll_plan_poisson_polar), pointer :: plan
    sll_real64, dimension(plan%nr+1,plan%ntheta+1), intent(in)  :: f
    sll_real64, dimension(plan%nr+1,plan%ntheta+1), intent(out) :: phi
    sll_int32 ,                           optional, intent(out) :: ierr

    sll_real64 :: rmin,dr
    sll_int32  :: nr, ntheta,bc(2)

    sll_real64 :: r
    sll_int32  :: i, k, ind_k
    sll_real64 :: kval, err
    sll_comp64 :: err_loc
    sll_int32  :: ierr_sup_1em12

    nr     = plan%nr
    ntheta = plan%ntheta
    rmin   = plan%rmin
    dr     = plan%dr
    
    bc         = plan%bc
    plan%f_fft = f

    do i=1,nr+1
      call fft_apply_plan(plan%pfwd,plan%f_fft(i,1:ntheta),plan%f_fft(i,1:ntheta))
    end do

    ierr_sup_1em12 = 0
    ! poisson solver
    do k = 0,ntheta/2
      ind_k=k
      kval=real(ind_k,f64)

      do i=2,nr
        r=rmin+real(i-1,f64)*dr
        plan%a(3*(i-1))=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
        plan%a(3*(i-1)-1)=2.0_f64/dr**2+(kval/r)**2
        plan%a(3*(i-1)-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)

        plan%fk(i)=fft_get_mode(plan%pfwd,plan%f_fft(i,1:ntheta),k)!ind_k)          
      enddo

      plan%phik=0.0_f64

      !boundary condition at rmin
      if(bc(1)==SLL_DIRICHLET)then !Dirichlet
        plan%a(1)=0.0_f64
      endif
      if(bc(1)==SLL_NEUMANN)then
        plan%a(2)=plan%a(2)+plan%a(1) !Neumann
        plan%a(1)=0._f64
      endif
      if(bc(1)==SLL_NEUMANN_MODE_0)then 
        if(k==0)then!Neumann for mode zero
          plan%a(2)=plan%a(2)+plan%a(1)
          plan%a(1)=0._f64
        else !Dirichlet for other modes
          plan%a(1)=0._f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==SLL_DIRICHLET)then !Dirichlet
        plan%a(3*(nr-1))=0.0_f64
      endif
      if(bc(2)==SLL_NEUMANN)then
        plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1)) !Neumann
        plan%a(3*(nr-1))=0.0_f64
      endif
      if(bc(2)==SLL_NEUMANN_MODE_0)then 
        if(k==0)then!Neumann for mode zero
          plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1))
          plan%a(3*(nr-1))=0.0_f64
        else !Dirichlet for other modes
          plan%a(3*(nr-1))=0.0_f64
        endif
      endif

      call setup_cyclic_tridiag(plan%a,nr-1,plan%cts,plan%ipiv)
      call solve_cyclic_tridiag(plan%cts,plan%ipiv,plan%fk(2:nr),nr-1,plan%phik(2:nr))

      !boundary condition at rmin
      if(bc(1)==SLL_DIRICHLET)then !Dirichlet
        plan%phik(1)=0.0_f64
      endif
      if(bc(1)==SLL_NEUMANN)then
        plan%phik(1)=plan%phik(2) !Neumann
      endif
      if(bc(1)==SLL_NEUMANN_MODE_0)then 
        if(k==0)then!Neumann for mode zero
          plan%phik(1)=plan%phik(2)
        else !Dirichlet for other modes
          plan%phik(1)=0.0_f64
        endif
      endif

      !boundary condition at rmax
      if(bc(2)==SLL_DIRICHLET)then !Dirichlet
        plan%phik(nr+1)=0.0_f64
      endif
      if(bc(2)==SLL_NEUMANN)then
        plan%phik(nr+1)=plan%phik(nr) !Neumann
      endif
      if(bc(2)==SLL_NEUMANN_MODE_0)then 
        if(k==0)then!Neumann for mode zero
          plan%phik(nr+1)=plan%phik(nr)
        else !Dirichlet for other modes
          plan%phik(nr+1)=0.0_f64
        endif
      endif

      err = 0._f64
      do i=4,nr-4
        r=rmin+real(i-1,f64)*dr
        err_loc=(plan%phik(i+1)-2*plan%phik(i)+plan%phik(i-1))/dr**2
        err_loc=err_loc+(plan%phik(i+1)-plan%phik(i-1))/(2._f64*r*dr)
        err_loc=-err_loc+kval**2/r**2*plan%phik(i)
        err_loc=(err_loc-plan%fk(i))
        if(abs(err_loc)>err)then
          err=abs(err_loc)
        endif
      enddo

      if (err>1e-12) then
        ierr_sup_1em12 = ierr_sup_1em12 + 1
      endif

      if(err>1e-4)then
        do i=2,nr
          r=rmin+real(i-1,f64)*dr
          err_loc=(plan%phik(i+1)-2*plan%phik(i)+plan%phik(i-1))/dr**2
          err_loc=err_loc+(plan%phik(i+1)-plan%phik(i-1))/(2._f64*r*dr)
          err_loc=-err_loc+kval**2/r**2*plan%phik(i)
          print *,r,real(err_loc),aimag(err_loc),real(plan%fk(i)),aimag(plan%fk(i))
        enddo
        stop
      endif

      do i=1,nr+1
        call fft_set_mode(plan%pinv,phi(i,1:ntheta),plan%phik(i),k)!ind_k)
      end do
    end do

    ! FFT INVERSE
    do i=1,nr+1
      call fft_apply_plan(plan%pinv,phi(i,1:ntheta),phi(i,1:ntheta))
    end do

    phi(:,ntheta+1)=phi(:,1)

    if (ierr_sup_1em12.ne.0) then
      if (present(ierr)) &
        ierr = ierr_sup_1em12 
    end if
  end subroutine poisson_solve_polar



end module sll_poisson_2d_polar
