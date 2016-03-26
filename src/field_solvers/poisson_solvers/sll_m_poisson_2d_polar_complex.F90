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
!> @brief Poisson equation solver in polar coordinates
!> @details Solver for the Poisson equation 
!> \f[ \Delta \phi = f \f]
!> in polar coordinates
!> using a fft in direction \f$ \theta \f$ and finite differences 
!> in direction \f$ r \f$.
!> This way we solve a tridiagonal system with the cyclic reduction solver.
!>
!> <b>How to use the Poisson polar solver?</b>
!>
!>You must add \code sll_m_poisson_2d_polar \endcode to the list of linked libraries.
!>The Poisson solver uses the FFT, so you also need to link the FFT
!>
!>1. Declare a Poisson polar plan
!>\code type(sll_t_plan_poisson_polar), pointer :: plan \endcode
!>2. Initialize the plan
!>\code plan => new_poisson_polar(dr,rmin,nr,ntheta,boundary_conditions) \endcode
!> `nr` and `ntheta` are step numbers in direction r and theta
!> 3. Execute the plan
!> \code call sll_s_poisson_solve_polar(plan,in,out) \endcode
!> 4. Delete the plan
!> \code call sll_o_delete(plan) \endcode
!>
!>\section bc Boundary conditions :
!>
!>The boundary conditions define the potential behaviour in \f$ r_{min} \f$
!> (BOT_) and \f$ r_{max} \f$ (TOP_)
!>They can take the value DIRICHLET or NEUMANN
!>
!>Summary
!>The boundary conditions define the potential behaviour in \f$ r_{min} \f$
!> (BOT_) and \f$ r_{max} \f$ (TOP_)
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
!>type(sll_t_plan_poisson_polar), pointer :: plan
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
!>plan => new_poisson_polar(dr,rmin,nr,ntheta,bc)
!>
!>!computation of Poisson
!>call sll_s_poisson_solve_polar(plan,in,out)
!>
!>!deletion of plan
!>call sll_o_delete(plan)
!>\endcode

module sll_m_poisson_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
  sll_p_dirichlet,                              &
  sll_p_neumann,                                &
  sll_p_neumann_mode_0

use sll_m_fft, only:         &
  sll_s_fft_exec_c2c_1d,     &
  sll_p_fft_backward,        &
  sll_s_fft_free,            &
  sll_p_fft_forward,         &
  sll_s_fft_init_c2c_1d,     &
  sll_t_fft



use sll_m_tridiagonal, only:  &
  sll_s_setup_cyclic_tridiag, &
  sll_o_solve_cyclic_tridiag

use sll_m_poisson_2d_base, only: &
  sll_c_poisson_2d_base, &
  sll_i_function_of_position

implicit none

public ::                       &
  sll_f_new_plan_poisson_polar, &
  sll_s_poisson_solve_polar,    &
  sll_o_create,                 &
  sll_o_delete,                 &
  sll_t_plan_poisson_polar,     &
  sll_o_solve,                  &
  sll_s_solve_poisson_polar,    &
  sll_f_new_poisson_2d_polar,   &
  sll_p_poisson_drift_kinetic

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!>type for the Poisson solver in polar coordinate
type :: sll_t_plan_poisson_polar

  sll_real64           :: rmin            !< r min
  sll_real64           :: rmax            !< r max
  sll_real64           :: dr              !< step size
  sll_int32            :: nr              !< number of points in r
  sll_int32            :: ntheta          !< number of points in theta
  sll_int32            :: bc(2)           !< boundary conditon type
  type(sll_t_fft)      :: pfwd            !< fft plan in theta
  type(sll_t_fft)      :: pinv            !< inverse fft plan in theta
  sll_comp64,  pointer :: f_fft(:,:)      !< potential fft in theta
  sll_comp64,  pointer :: fk(:)           !< \f$ f_k \f$
  sll_comp64,  pointer :: phik(:)         !< \f$ phi_k \f$
  sll_real64,  pointer :: a(:)            !< data for the tridiagonal solver
  sll_real64,  pointer :: cts(:)          !< lapack array
  sll_int32,   pointer :: ipiv(:)         !< lapack pivot data
  sll_real64,  pointer :: dlog_density(:) !< for quasi neutral solver
  sll_real64,  pointer :: inv_Te(:)       !< for quasi neutral solver

end type sll_t_plan_poisson_polar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Classic Poisson solver
sll_int32, parameter :: SLL_POISSON_CLASSIC = 0
!> Poisson solver for drift kinetic simulation
sll_int32, parameter :: sll_p_poisson_drift_kinetic = 1

!> Poisson solver in polar coordinates
type, extends(sll_c_poisson_2d_base) :: poisson_2d_polar_solver     

  type(sll_t_plan_poisson_polar), pointer :: solver       !< workspace for polar solver
  sll_int32                               :: poisson_case !< classic or drift kinetic
  sll_real64,dimension(:), pointer        :: dlog_density !< QNS parameter
  sll_real64,dimension(:), pointer        :: inv_Te       !< QNS parameter

contains

  !> Initialize and allocate arrays.
  procedure, pass(poisson) :: initialize => initialize_poisson_2d_polar_solver
  !> Solves \f$ -\Delta phi(x,y) = rho(x,y) \f$
  procedure, pass(poisson) :: compute_phi_from_rho => compute_phi_from_rho_2d_polar
  !> Solves \f$ -\Delta phi(x,y) = rho(x,y) \f$ and \f$ E = \nabla  \phi \f$
  procedure, pass(poisson) :: compute_E_from_rho => compute_E_from_rho_2d_polar
  
  !> Compute the squared L_2 for given coefficients
  procedure :: &
       l2norm_squared => l2norm_squared_2d_polar
  !> Compute the right hand side from a given function
  procedure :: &
       compute_rhs_from_function => compute_rhs_from_function_2d_polar
  !> Destructor
  procedure :: &
       free => delete_2d_polar

end type poisson_2d_polar_solver

!> Initialize the polar poisson solver
interface sll_o_create
  module procedure initialize_poisson_polar
end interface sll_o_create

!> Get potential from the polar poisson solver
interface sll_o_solve
  module procedure sll_s_solve_poisson_polar
end interface sll_o_solve

!> Deallocate memory
interface sll_o_delete
  module procedure delete_plan_poisson_polar
end interface sll_o_delete

contains

!> Creation of sll_t_plan_poisson_polar object for the 
!> Poisson solver in polar coordinate
!> @returns a Poisson solver object for polar coordinates
function sll_f_new_plan_poisson_polar(dr,           &
                                      rmin,         &
                                      nr,           &
                                      ntheta,       &
                                      bc,           &
                                      dlog_density, &
                                      inv_Te) result(self)

  sll_real64 :: dr             !< size of space in direction r
  sll_real64 :: rmin           !< interior radius
  sll_int32  :: nr             !< number of space in direction r
  sll_int32  :: ntheta         !< number of space in direction theta
  sll_int32, optional :: bc(2) !< Boundary conditions, can be combined with +
                               !< optional and default is Dirichlet in rmin and rmax

  type(sll_t_plan_poisson_polar), pointer  :: self         !< Poisson solver structure
  sll_real64,dimension(:),        optional :: dlog_density !< for quasi neutral solver
  sll_real64,dimension(:),        optional :: inv_Te       !< for quasi neutral solver

  sll_int32 :: err
  sll_comp64, dimension(:), allocatable :: buf

  SLL_ALLOCATE(self,err)
  SLL_ALLOCATE(buf(ntheta),err)
  SLL_ALLOCATE(self%f_fft(nr+1,ntheta+1),err)
  SLL_ALLOCATE(self%fk(nr+1),err)
  SLL_ALLOCATE(self%phik(nr+1),err)
  SLL_ALLOCATE(self%a(3*(nr-1)),err)
  SLL_ALLOCATE(self%cts(7*(nr-1)),err)
  SLL_ALLOCATE(self%ipiv(nr-1),err)
  
  SLL_ALLOCATE(self%dlog_density(nr+1),err)
  SLL_ALLOCATE(self%inv_Te(nr+1),err)
  
  self%dlog_density = 0._f64
  self%inv_Te       = 0._f64
  
  if (present(dlog_density)) self%dlog_density = dlog_density
  if (present(inv_Te))       self%inv_Te = inv_Te
  
  self%dr     = dr
  self%rmin   = rmin
  self%nr     = nr
  self%ntheta = ntheta

  if (present(bc)) then
    self%bc=bc
  else
    self%bc(1)=-1
    self%bc(2)=-1
  end if

  call sll_s_fft_init_c2c_1d(&
       self%pfwd, ntheta,buf,buf,sll_p_fft_forward,normalized = .TRUE.)
  call sll_s_fft_init_c2c_1d( &
       self%pinv, ntheta,buf,buf,sll_p_fft_backward)
  
  SLL_DEALLOCATE_ARRAY(buf,err)

end function sll_f_new_plan_poisson_polar

!> Initialize the Poisson solver in polar coordinates
subroutine initialize_poisson_polar(self,         &
                                    rmin,         &
                                    rmax,         &
                                    nr,           &
                                    ntheta,       &
                                    bc_rmin,      &
                                    bc_rmax,      &
                                    dlog_density, &
                                    inv_Te)

  type(sll_t_plan_poisson_polar) :: self !< Poisson solver object

  sll_real64, intent(in) :: rmin            !< r min
  sll_real64, intent(in) :: rmax            !< r max
  sll_int32,  intent(in) :: nr              !< number of cells radial
  sll_int32,  intent(in) :: ntheta          !< number of cells angular
  sll_int32,  optional   :: bc_rmin         !< radial boundary conditions
  sll_int32,  optional   :: bc_rmax         !< radial boundary conditions
  sll_real64, optional   :: dlog_density(:) !< For quasi neutral solver
  sll_real64, optional   :: inv_Te(:)       !< For quasi neutral solver

  sll_int32               :: error
  sll_comp64, allocatable :: buf(:)

  SLL_ALLOCATE(self%f_fft(nr+1,ntheta+1),error)
  SLL_ALLOCATE(self%fk(nr+1),            error)
  SLL_ALLOCATE(self%phik(nr+1),          error)
  SLL_ALLOCATE(self%a(3*(nr-1)),         error)
  SLL_ALLOCATE(self%cts(7*(nr-1)),       error)
  SLL_ALLOCATE(self%ipiv(nr-1),          error)

  self%rmin   = rmin
  self%rmax   = rmax
  self%dr     = (rmax-rmin)/nr
  self%nr     = nr
  self%ntheta = ntheta

  SLL_ALLOCATE(self%dlog_density(nr+1),error)
  SLL_ALLOCATE(self%inv_Te(nr+1),error)
  
  self%dlog_density = 0._f64
  self%inv_Te = 0._f64
  
  if(present(dlog_density)) self%dlog_density = dlog_density
  if(present(inv_Te))       self%inv_Te       = inv_Te

  if (present(bc_rmin) .and. present(bc_rmax)) then
    self%bc(1) = bc_rmin
    self%bc(2) = bc_rmax
  else
    self%bc(1) = -1
    self%bc(2) = -1
  end if

  SLL_ALLOCATE(buf(ntheta),error)
  call sll_s_fft_init_c2c_1d(self%pfwd, ntheta, &
   buf,buf,sll_p_fft_forward,normalized = .TRUE.)
  call sll_s_fft_init_c2c_1d(self%pinv, ntheta,buf,buf,sll_p_fft_backward)
  SLL_DEALLOCATE_ARRAY(buf,error)

end subroutine initialize_poisson_polar

!=====================================
!deletion of sll_t_plan_poisson_polar
!=====================================

!>delete a sll_t_plan_poisson_polar object
subroutine delete_plan_poisson_polar(self)

  type(sll_t_plan_poisson_polar), pointer :: self
  sll_int32 :: err

  if (associated(self)) then
    call sll_s_fft_free(self%pfwd)
    call sll_s_fft_free(self%pinv)
    SLL_DEALLOCATE_ARRAY(self%f_fft,err)
    SLL_DEALLOCATE_ARRAY(self%fk,err)
    SLL_DEALLOCATE_ARRAY(self%phik,err)
    SLL_DEALLOCATE_ARRAY(self%a,err)
    SLL_DEALLOCATE_ARRAY(self%cts,err)
    SLL_DEALLOCATE_ARRAY(self%ipiv,err)
    SLL_DEALLOCATE(self,err)
  end if

end subroutine delete_plan_poisson_polar

!===================
!  Poisson solver
!===================

!>subroutine sll_s_solve_poisson_polar(plan,f,phi)
!>poisson solver for polar system : \f$ -\Delta (phi)=f \f$
!>@param plan : sll_t_plan_poisson_polar, contains data for the solver
!>@param f : distribution function, size (nr+1)*(ntheta+1), input
!>@param phi : unknown field, size (nr+1)*(ntheta+1), output
!>initialization must be done outside the solver
subroutine sll_s_solve_poisson_polar(plan,f,phi)

  type(sll_t_plan_poisson_polar) :: plan
  sll_real64, dimension(:,:), intent(in)  :: f
  sll_real64, dimension(:,:), intent(out) :: phi

  sll_real64 :: rmin,dr
  sll_int32  :: nr, ntheta,bc(2)

  sll_real64 :: r
  sll_int32  :: i, k, ind_k
  sll_real64 :: kval

  sll_comp64 :: err_loc
  sll_real64 :: err

  nr     = plan%nr
  ntheta = plan%ntheta
  rmin   = plan%rmin
  dr     = plan%dr
  
  bc         = plan%bc
  plan%f_fft = f

  do i=1,nr+1
    call sll_s_fft_exec_c2c_1d(plan%pfwd,plan%f_fft(i,1:ntheta), &
      plan%f_fft(i,1:ntheta))
  end do

  !do k = 0,ntheta/2
  do k = 1,ntheta

    ind_k=k
    if (ind_k<=ntheta/2) then
        ind_k = ind_k-1
    else
        ind_k = ntheta-(ind_k-1)
    endif



    kval=real(ind_k,f64)

    do i=2,nr
      r = rmin + (i-1)*dr
      plan%a(3*(i-1)  ) = -1.0_f64/dr**2-1.0_f64/(2._f64*dr*r) &
       -plan%dlog_density(i)/(2._f64*dr)
      plan%a(3*(i-1)-1) =  2.0_f64/dr**2+(kval/r)**2+plan%inv_Te(i)
      plan%a(3*(i-1)-2) = -1.0_f64/dr**2+1.0_f64/(2._f64*dr*r) &
       +plan%dlog_density(i)/(2._f64*dr)
      
      !plan%fk(i)=sll_f_fft_get_mode_r2c_1d(plan%pfwd,plan%f_fft(i,1:ntheta),k)
      plan%fk(i) = plan%f_fft(i,k)
    enddo
    
    plan%phik=(0.0_f64,0.0_f64)

    !boundary condition at rmin
    if(bc(1)==sll_p_dirichlet)then !Dirichlet
      plan%a(1)=0.0_f64
    endif
    if(bc(1)==sll_p_neumann)then
      plan%a(2)=plan%a(2)+plan%a(1) !Neumann
      plan%a(1)=0._f64
    endif
    if(bc(1)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%a(2)=plan%a(2)+plan%a(1)
        plan%a(1)=0._f64
      else !Dirichlet for other modes
        plan%a(1)=0._f64
      endif
    endif

    !boundary condition at rmax
    if(bc(2)==sll_p_dirichlet)then !Dirichlet
      plan%a(3*(nr-1))=0.0_f64
    endif
    if(bc(2)==sll_p_neumann)then
      plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1)) !Neumann
      plan%a(3*(nr-1))=0.0_f64
    endif
    if(bc(2)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1))
        plan%a(3*(nr-1))=0.0_f64
      else !Dirichlet for other modes
        plan%a(3*(nr-1))=0.0_f64
      endif
    endif

    call sll_s_setup_cyclic_tridiag(plan%a,nr-1,plan%cts,plan%ipiv)
    call sll_o_solve_cyclic_tridiag(plan%cts,plan%ipiv,plan%fk(2:nr), &
                   nr-1,plan%phik(2:nr))

    !boundary condition at rmin
    if(bc(1)==sll_p_dirichlet)then !Dirichlet
      plan%phik(1)=(0.0_f64,0.0_f64)
    endif
    if(bc(1)==sll_p_neumann)then
      plan%phik(1)=plan%phik(2) !Neumann
    endif
    if(bc(1)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%phik(1)=plan%phik(2)
      else !Dirichlet for other modes
        plan%phik(1)=(0.0_f64,0.0_f64)
      endif
    endif

    !boundary condition at rmax
    if(bc(2)==sll_p_dirichlet)then !Dirichlet
      plan%phik(nr+1)=(0.0_f64,0.0_f64)
    endif
    if(bc(2)==sll_p_neumann)then
      plan%phik(nr+1)=plan%phik(nr) !Neumann
    endif
    if(bc(2)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%phik(nr+1)=plan%phik(nr)
      else !Dirichlet for other modes
        plan%phik(nr+1)=(0.0_f64,0.0_f64)
      endif
    endif

    do i=1,nr+1
      plan%f_fft(i,k) = plan%phik(i)
    !  call sll_s_fft_set_mode_c2r_1d(plan%pinv,phi(i,1:ntheta),plan%phik(i),k)
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
    !call sll_s_fft_exec_r2r_1d(plan%pinv,phi(i,1:ntheta),phi(i,1:ntheta))
    call sll_s_fft_exec_c2c_1d(plan%pinv,plan%f_fft(i,1:ntheta),plan%f_fft(i,1:ntheta))
    phi(i,1:ntheta) = real(plan%f_fft(i,1:ntheta),f64)
  end do

  phi(:,ntheta+1)=phi(:,1)
  
end subroutine sll_s_solve_poisson_polar

!>subroutine sll_s_poisson_solve_polar(plan,f,phi)
!>poisson solver for polar system : \f$ -\Delta (phi)=fa\f$
!>@param plan : sll_t_plan_poisson_polar, contains data for the solver
!>@param f : distribution function, size (nr+1)*(ntheta+1), input
!>@param phi : unknown field, size (nr+1)*(ntheta+1), output
!>initialization must be done outside the solver
subroutine sll_s_poisson_solve_polar(plan,f,phi,ierr)

  type(sll_t_plan_poisson_polar), pointer :: plan
  sll_real64, dimension(plan%nr+1,plan%ntheta+1), intent(in)  :: f
  sll_real64, dimension(plan%nr+1,plan%ntheta+1), intent(out) :: phi
  sll_int32 ,                           optional              :: ierr !< error code

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
    call sll_s_fft_exec_c2c_1d(plan%pfwd,plan%f_fft(i,1:ntheta), &
     plan%f_fft(i,1:ntheta))
  end do

  ierr_sup_1em12 = 0
  ! poisson solver
  !do k = 0,ntheta/2
  do k=1,ntheta
    ind_k=k
    if (ind_k<=ntheta/2) then
        ind_k = ind_k-1
    else
        ind_k = ntheta-(ind_k-1)
    endif
    
    kval=real(ind_k,f64)

    do i=2,nr
      r=rmin+real(i-1,f64)*dr
      plan%a(3*(i-1))=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
      plan%a(3*(i-1)-1)=2.0_f64/dr**2+(kval/r)**2
      plan%a(3*(i-1)-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)

      !plan%fk(i)=sll_f_fft_get_mode_r2c_1d(plan%pfwd,plan%f_fft(i,1:ntheta),k)
      plan%fk(i) = plan%f_fft(i,k)
    enddo

    plan%phik=(0.0_f64,0.0_f64)

    !boundary condition at rmin
    if(bc(1)==sll_p_dirichlet)then !Dirichlet
      plan%a(1)=0.0_f64
    endif
    if(bc(1)==sll_p_neumann)then
      plan%a(2)=plan%a(2)+plan%a(1) !Neumann
      plan%a(1)=0._f64
    endif
    if(bc(1)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%a(2)=plan%a(2)+plan%a(1)
        plan%a(1)=0._f64
      else !Dirichlet for other modes
        plan%a(1)=0._f64
      endif
    endif

    !boundary condition at rmax
    if(bc(2)==sll_p_dirichlet)then !Dirichlet
      plan%a(3*(nr-1))=0.0_f64
    endif
    if(bc(2)==sll_p_neumann)then
      plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1)) !Neumann
      plan%a(3*(nr-1))=0.0_f64
    endif
    if(bc(2)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1))
        plan%a(3*(nr-1))=0.0_f64
      else !Dirichlet for other modes
        plan%a(3*(nr-1))=0.0_f64
      endif
    endif

    call sll_s_setup_cyclic_tridiag(plan%a,nr-1,plan%cts,plan%ipiv)
    call sll_o_solve_cyclic_tridiag(plan%cts,plan%ipiv, &
      plan%fk(2:nr),nr-1,plan%phik(2:nr))

    !boundary condition at rmin
    if(bc(1)==sll_p_dirichlet)then !Dirichlet
      plan%phik(1)=(0.0_f64,0.0_f64)
    endif
    if(bc(1)==sll_p_neumann)then
      plan%phik(1)=plan%phik(2) !Neumann
    endif
    if(bc(1)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%phik(1)=plan%phik(2)
      else !Dirichlet for other modes
        plan%phik(1)=(0.0_f64,0.0_f64)
      endif
    endif

    !boundary condition at rmax
    if(bc(2)==sll_p_dirichlet)then !Dirichlet
      plan%phik(nr+1)=(0.0_f64,0.0_f64)
    endif
    if(bc(2)==sll_p_neumann)then
      plan%phik(nr+1)=plan%phik(nr) !Neumann
    endif
    if(bc(2)==sll_p_neumann_mode_0)then 
      if(k==0)then!Neumann for mode zero
        plan%phik(nr+1)=plan%phik(nr)
      else !Dirichlet for other modes
        plan%phik(nr+1)=(0.0_f64,0.0_f64)
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
      !call sll_s_fft_set_mode_c2r_1d(plan%pinv,phi(i,1:ntheta),plan%phik(i),k)
      plan%f_fft(i,k) = plan%phik(i)
    end do
  end do

  ! FFT INVERSE
  do i=1,nr+1
    call sll_s_fft_exec_c2c_1d(plan%pinv,plan%f_fft(i,1:ntheta),plan%f_fft(i,1:ntheta))
    phi(i,1:ntheta) = real(plan%f_fft(i,1:ntheta),f64)
    !call sll_s_fft_exec_r2r_1d(plan%pinv,phi(i,1:ntheta),phi(i,1:ntheta))
  end do

  phi(:,ntheta+1)=phi(:,1)

  if (ierr_sup_1em12.ne.0) then
    if (present(ierr)) &
      ierr = ierr_sup_1em12 
  end if

end subroutine sll_s_poisson_solve_polar

!> Allocate a new Poisson solver in polar coordinates
!> @returns a pointer to the derived type
function sll_f_new_poisson_2d_polar( &
  eta1_min,                          &
  eta1_max,                          &
  nc_eta1,                           &
  nc_eta2,                           &
  bc,                                &
  dlog_density,                      &
  inv_Te,                            &
  poisson_case)                      &     
  result(poisson)
    
  type(poisson_2d_polar_solver), pointer :: poisson
  sll_real64, intent(in)                 :: eta1_min
  sll_real64, intent(in)                 :: eta1_max
  sll_int32,  intent(in)                 :: nc_eta1
  sll_int32,  intent(in)                 :: nc_eta2
  sll_int32,  intent(in)                 :: bc(2)
  sll_real64, intent(in), optional       :: dlog_density(:)
  sll_real64, intent(in), optional       :: inv_Te(:)
  sll_int32,  optional                   :: poisson_case
  sll_int32                              :: ierr
    
  SLL_ALLOCATE(poisson,ierr)
  call initialize_poisson_2d_polar_solver( &
    poisson,                               &
    eta1_min,                              &
    eta1_max,                              &
    nc_eta1,                               &
    nc_eta2,                               &
    bc,                                    &
    dlog_density,                          &
    inv_Te,                                &
    poisson_case)
  
end function sll_f_new_poisson_2d_polar

subroutine initialize_poisson_2d_polar_solver( &
  poisson,                                     &
  eta1_min,                                    &
  eta1_max,                                    &
  nc_eta1,                                     &
  nc_eta2,                                     &
  bc,                                          &
  dlog_density,                                &
  inv_Te,                                      &
  poisson_case)

  class(poisson_2d_polar_solver) :: poisson

  sll_real64, intent(in) :: eta1_min
  sll_real64, intent(in) :: eta1_max
  sll_int32,  intent(in) :: nc_eta1
  sll_int32,  intent(in) :: nc_eta2
  sll_int32,  intent(in) :: bc(2)
  sll_real64, optional   :: dlog_density(:)
  sll_real64, optional   :: inv_Te(:)
  sll_int32,  optional   :: poisson_case

  sll_int32  :: ierr
  sll_real64 :: delta_eta
  
  delta_eta = (eta1_max-eta1_min)/real(nc_eta1,f64)
  
  if(present(poisson_case)) then
    poisson%poisson_case = poisson_case  
  else   
    poisson%poisson_case = SLL_POISSON_CLASSIC
  endif
  
  select case(poisson%poisson_case)

    case (SLL_POISSON_CLASSIC)
      poisson%solver => sll_f_new_plan_poisson_polar( delta_eta,& 
                                                      eta1_min, &
                                                      nc_eta1,  &
                                                      nc_eta2,  &
                                                      bc)
   case (sll_p_poisson_drift_kinetic)    

     SLL_ALLOCATE(poisson%dlog_density(nc_eta1+1),ierr)
     SLL_ALLOCATE(poisson%inv_Te(nc_eta1+1),ierr)

     if(.not.(present(dlog_density)))then
       print *,'#dlog_density should be present in initialize_poisson_2d_polar_solver'
       stop
     endif

     if(size(dlog_density)<nc_eta1+1)then
       print *,'#Bad size for dlog_density',size(dlog_density)
       stop
     endif

     if(.not.(present(inv_Te)))then
       print *,'#dlog_density should be present in initialize_poisson_2d_polar_solver'
       stop
     endif

     if(size(inv_Te)<nc_eta1+1)then
       print *,'#Bad size for dlog_density',size(inv_Te)
       stop
     endif

     poisson%dlog_density(1:nc_eta1+1)=dlog_density(1:nc_eta1+1)
     poisson%inv_Te(1:nc_eta1+1)=inv_Te(1:nc_eta1+1)
     poisson%solver => sll_f_new_plan_poisson_polar( &
       delta_eta,                                    & 
       eta1_min,                                     &
       nc_eta1,                                      &
       nc_eta2,                                      &
       bc,                                           &
       poisson%dlog_density,                         &
       poisson%inv_Te)

   case default

      print *,'#bad value of poisson_case=', poisson%poisson_case
      print *,'#not implemented'
      print *,'#in initialize_poisson_2d_polar_solver'
      stop

   end select   

end subroutine initialize_poisson_2d_polar_solver

subroutine compute_phi_from_rho_2d_polar( poisson, phi, rho )

  class(poisson_2d_polar_solver), target      :: poisson
  sll_real64,dimension(:,:),      intent(in)  :: rho
  sll_real64,dimension(:,:),      intent(out) :: phi
  
  select case(poisson%poisson_case)
    case (SLL_POISSON_CLASSIC)
      call sll_s_poisson_solve_polar(poisson%solver,rho,phi)            
    case (sll_p_poisson_drift_kinetic)    
      call sll_s_solve_poisson_polar(poisson%solver,rho,phi)
    case default
      print *,'#bad value of poisson_case=', poisson%poisson_case
      print *,'#not implemented'
      print *,'in compute_phi_from_rho_2d_polar'
      stop
  end select   

end subroutine compute_phi_from_rho_2d_polar

!> Solves \f$ \vec{E} = -\nabla \phi \f$ with \f$ -\Delta \phi(x,y) = rho(x,y) \f$.
subroutine compute_E_from_rho_2d_polar( poisson, E1, E2, rho )

  class(poisson_2d_polar_solver)        :: poisson
  sll_real64,dimension(:,:),intent(in)  :: rho
  sll_real64,dimension(:,:),intent(out) :: E1
  sll_real64,dimension(:,:),intent(out) :: E2
    
  print *,'#compute_E_from_rho_2d_polar'      
  print *,'#not implemented for the moment'
    
  E1 = 0._f64
  E2 = 0._f64
  print *,maxval(rho)
    
  if(.not.(associated(poisson%solver)))then
    stop '#poisson%solver is not associated'
  endif

  stop
      
end subroutine compute_E_from_rho_2d_polar

subroutine compute_rhs_from_function_2d_polar(poisson, func, coefs_dofs)
  class( poisson_2d_polar_solver)                    :: poisson !< Maxwell solver object.
  procedure(sll_i_function_of_position)          :: func !< Function to be projected.
  sll_real64, intent(out)                        :: coefs_dofs(:) !< Coefficients of the projection.

  SLL_ERROR('compute_rhs_from_function_2d_polar', 'Procedure not implemented.')
  
end subroutine compute_rhs_from_function_2d_polar

function l2norm_squared_2d_polar( poisson, coefs_dofs) result(r)
  class( poisson_2d_polar_solver), intent(in) :: poisson !< Poisson solver object.
       sll_real64 , intent(in)                :: coefs_dofs(:,:) !< Values of the coefficient vectors for each DoF
       sll_real64                             :: r
  
  SLL_ERROR('l2norm_squared_2d_polar', 'Procedure not implemented.')
  r = 0.0_f64

end function l2norm_squared_2d_polar
  
subroutine delete_2d_polar( poisson )
  class( poisson_2d_polar_solver) :: poisson !< Poisson solver object.

end subroutine delete_2d_polar
  

end module sll_m_poisson_2d_polar
