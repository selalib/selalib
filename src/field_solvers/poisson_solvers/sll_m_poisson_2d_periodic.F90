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
!> @brief
!> Regular cartesian two dimensional mesh with periodic bounday conditions.
!> @details
!> Numerical method uses Fast Fourier Transform and periodic
!> boundary conditions.
!> @snippet poisson_solvers/test_poisson_2d_fft.F90 example
module sll_m_poisson_2d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_fftw.h"

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base, &
    sll_f_function_of_position

  use sll_m_constants, only: sll_p_pi

use iso_c_binding, only: &
    c_associated, &
    c_double, &
    c_double_complex, &
    c_f_pointer, &
    c_ptr, &
    c_size_t

#ifdef FFTW_F2003
  use sll_m_fftw3, only: &
    fftw_alloc_complex, &
    fftw_alloc_real, &
    fftw_destroy_plan, &
    fftw_estimate, &
    fftw_execute_dft_c2r, &
    fftw_execute_dft_r2c, &
    fftw_free, &
    fftw_plan_dft_c2r_2d, &
    fftw_plan_dft_r2c_2d
#else
  use sll_m_fftw3, only : &
       fftw_estimate
#endif

implicit none

private

public ::                            &
  sll_f_new_poisson_2d_periodic, &
  sll_f_new_poisson_2d_periodic_fftpack, &
  sll_t_poisson_2d_periodic,     &
  sll_t_poisson_2d_periodic_fftpack, &
#ifdef FFTW
  sll_t_poisson_2d_periodic_fftw,    &
#endif
  sll_o_initialize,                  &
  sll_o_solve

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> fft type use to do fft with fftpack library
type :: fftclass
   sll_real64, dimension(:), pointer :: coefc !< data for complex fft
   sll_real64, dimension(:), pointer :: work  !< work data for fft
   sll_real64, dimension(:), pointer :: workc !< workc complex
   sll_real64, dimension(:), pointer :: coefd !< data for double fft
   sll_real64, dimension(:), pointer :: workd !< work data
   sll_real64, dimension(:), pointer :: coefcd!< data for complex fft
   sll_int32  :: n  !< number of samples in each sequence
end type fftclass

!> Object with data to sll_o_solve Poisson equation on 2d domain with
!> periodic boundary conditions
type :: sll_t_poisson_2d_periodic_fftpack
  sll_int32   :: nc_x  !< number of cells direction x
  sll_int32   :: nc_y  !< number of cells direction y
  sll_real64  :: x_min !< left corner direction x
  sll_real64  :: x_max !< right corner direction x
  sll_real64  :: y_min !< left corner direction y
  sll_real64  :: y_max !< right corner direction y
  sll_real64  :: dx    !< step size direction x
  sll_real64  :: dy    !< step size direction y
  sll_comp64, dimension(:,:), pointer :: rhst !< rhs fft
  sll_comp64, dimension(:,:), pointer :: ext  !< x electric field fft
  sll_comp64, dimension(:,:), pointer :: eyt  !< y electric field fft
  sll_real64, dimension(:,:), pointer :: kx   !< wave number x
  sll_real64, dimension(:,:), pointer :: ky   !< wave number y
  sll_real64, dimension(:,:), pointer :: k2   !< \f$ k_x^2+k_y^2 \f$
  type(fftclass)                      :: fftx !< fft plan in direction x
  type(fftclass)                      :: ffty !< fft plan in direction y
end type sll_t_poisson_2d_periodic_fftpack

#ifdef FFTW
!> derived type to sll_o_solve the Poisson equation on 2d regular 
!> cartesian mesh with periodic boundary conditions on both sides
type :: sll_t_poisson_2d_periodic_fftw

   private
   sll_real64, dimension(:,:), pointer :: kx       !< wave number in x
   sll_real64, dimension(:,:), pointer :: ky       !< wave number in y
   sll_real64, dimension(:,:), pointer :: k2       !< \f[ k_x^2 + k_y^2 \f]
   sll_int32                           :: nc_x     !< cells number in x
   sll_int32                           :: nc_y     !< cells number in y
   sll_real64                          :: dx       !< x step size
   sll_real64                          :: dy       !< y step size
   fftw_plan                           :: fw       !< forward fftw plan
   fftw_plan                           :: bw       !< backward fftw plan
   fftw_comp                , pointer  :: rht(:,:) !< fft(rho)
   fftw_comp                , pointer  :: exy(:,:) !< fft(ex and ey)
   fftw_plan                           :: p_rho    !< C array pointer
   fftw_plan                           :: p_exy    !< C array pointer
   fftw_plan                           :: p_tmp    !< C array pointer
   fftw_real                , pointer  :: tmp(:,:)

end type sll_t_poisson_2d_periodic_fftw

#endif

type, extends(sll_c_poisson_2d_base) :: sll_t_poisson_2d_periodic

  type(sll_t_poisson_2d_periodic_fftpack), private, pointer :: solver

contains

  !> Create the Poisson solver
  procedure, public, pass(poisson) :: initialize => &
    initialize_poisson_2d_periodic
  !> Compute potential solving the Poisson equation
  procedure, public, pass(poisson) :: compute_phi_from_rho => &
    compute_phi_from_rho_2d_fft
  !> Compute electric fields solving the Poisson equation
  procedure, public, pass(poisson) :: compute_E_from_rho => &
    compute_E_from_rho_2d_fft

  !> Compute the squarred L_2 for given coefficients
  procedure :: &
       l2norm_squared => l2norm_squarred_2d_periodic
  !> Compute the right hand side from a given function
  procedure :: &
       compute_rhs_from_function => compute_rhs_from_function_2d_periodic
  !> Destructor
  procedure :: free => delete_2d_periodic
    
end type sll_t_poisson_2d_periodic

interface sll_o_initialize
  module procedure initialize_poisson_2d_periodic_fftpack
#ifdef FFTW
  module procedure initialize_poisson_2d_periodic_fftw
#endif
end interface

interface sll_o_solve
   module procedure solve_potential_poisson_2d_periodic_fftpack
   module procedure solve_e_fields_poisson_2d_periodic_fftpack
#ifdef FFTW
   module procedure solve_potential_poisson_2d_periodic_fftw
   module procedure solve_e_fields_poisson_2d_periodic_fftw
#endif
end interface

interface delete
   module procedure free_poisson_2d_periodic_fftpack
#ifdef FFTW
   module procedure delete_poisson_2d_periodic_fftw
#endif
end interface

contains

  
  function l2norm_squarred_2d_periodic(poisson, coefs_dofs) result(r)
    class( sll_t_poisson_2d_periodic) , intent(in)        :: poisson !< Poisson solver object.
    sll_real64   , intent(in)                                  :: coefs_dofs(:,:) !< Values of the coefficient vectors for each DoF
    sll_real64                                     :: r
    
    r = 0.0_f64
    print*, 'l2norm_squared not implemented for sll_t_poisson_2d_periodic.'
    
  end function l2norm_squarred_2d_periodic
  
  subroutine compute_rhs_from_function_2d_periodic(poisson, func, coefs_dofs)
    class( sll_t_poisson_2d_periodic)                    :: poisson !< Maxwell solver object.
    procedure(sll_f_function_of_position)          :: func !< Function to be projected.
    sll_real64, intent(out)                        :: coefs_dofs(:) !< Coefficients of the projection.
    
    print*, 'compute_rhs_from_function not implemented for sll_t_poisson_2d_periodic.'
    
  end subroutine compute_rhs_from_function_2d_periodic
  
  subroutine delete_2d_periodic(poisson)
    class( sll_t_poisson_2d_periodic)                    :: poisson !< Maxwell solver object.
  end subroutine delete_2d_periodic



  !> @returns a pointer to the derived type sll_t_poisson_2d_periodic.
  function sll_f_new_poisson_2d_periodic( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2) &     
    result(poisson)
      
    type(sll_t_poisson_2d_periodic),pointer :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_periodic( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    
  end function sll_f_new_poisson_2d_periodic
  
  subroutine initialize_poisson_2d_periodic( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    class(sll_t_poisson_2d_periodic) :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
    
    SLL_ALLOCATE(poisson%solver,ierr)
    
    call sll_o_initialize( &
      poisson%solver, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      ierr) 

  end subroutine initialize_poisson_2d_periodic
  
  !> solves \f$ -\Delta phi(x,y) = rho (x,y) \f$
  subroutine compute_phi_from_rho_2d_fft( poisson, phi, rho )
    class(sll_t_poisson_2d_periodic), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    call sll_o_solve( poisson%solver, phi, rho)
    
  end subroutine compute_phi_from_rho_2d_fft

  !> @brief
  !> sll_o_solve Poisson equation to compute electric fields
  !> @details
  !> solves 
  !> \f[ 
  !> E(x,y) = -\nabla \phi(x,y) \\
  !> -\Delta \phi(x,y) = \rho(x,y)
  !> \f]
  subroutine compute_E_from_rho_2d_fft( poisson, E1, E2, rho )
    class(sll_t_poisson_2d_periodic) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: E1
    sll_real64,dimension(:,:),intent(out) :: E2
      
    call sll_o_solve( poisson%solver, E1, E2, rho)
      
  end subroutine compute_E_from_rho_2d_fft

  !> Create a sll_o_new solver
  !> @return
  function sll_f_new_poisson_2d_periodic_fftpack(&
    x_min, &
    x_max, &
    nc_x, &
    y_min, &
    y_max, &
    nc_y, &
    error) &
    result(self)

   type(sll_t_poisson_2d_periodic_fftpack),pointer :: self   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code

   SLL_ALLOCATE(self, error)
   call initialize_poisson_2d_periodic_fftpack( &
           self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

  end function sll_f_new_poisson_2d_periodic_fftpack 


!> delete sll_t_poisson_2d_periodic_fftpack
subroutine free_poisson_2d_periodic_fftpack( self, error )
   type(sll_t_poisson_2d_periodic_fftpack) :: self   !< self object
   sll_int32,  intent(out)   :: error  !< error code

   SLL_DEALLOCATE(self%rhst, error)
   SLL_DEALLOCATE(self%ext, error)
   SLL_DEALLOCATE(self%eyt, error)
   SLL_DEALLOCATE(self%kx, error)
   SLL_DEALLOCATE(self%ky, error)
   SLL_DEALLOCATE(self%k2, error)
   
   SLL_DEALLOCATE(self%fftx%coefc, error)
   SLL_DEALLOCATE(self%fftx%work, error)
   SLL_DEALLOCATE(self%fftx%workc, error)
   SLL_DEALLOCATE(self%fftx%coefd, error)
   SLL_DEALLOCATE(self%fftx%workd, error)
   SLL_DEALLOCATE(self%fftx%coefcd, error)

   SLL_DEALLOCATE(self%ffty%coefc, error)
   SLL_DEALLOCATE(self%ffty%work, error)
   SLL_DEALLOCATE(self%ffty%workc, error)
   SLL_DEALLOCATE(self%ffty%coefd, error)
   SLL_DEALLOCATE(self%ffty%workd, error)
   SLL_DEALLOCATE(self%ffty%coefcd, error)

      
end subroutine free_poisson_2d_periodic_fftpack

!> Create an object to sll_o_solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
subroutine initialize_poisson_2d_periodic_fftpack( &
           self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

   type(sll_t_poisson_2d_periodic_fftpack) :: self   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code
   
   self%nc_x = nc_x
   self%nc_y = nc_y

   self%x_min = x_min
   self%x_max = x_max
   self%y_min = y_min
   self%y_max = y_max
   self%dx   = (x_max-x_min) / real(nc_x, f64)
   self%dy   = (y_max-y_min) / real(nc_y, f64)

   SLL_ALLOCATE(self%rhst(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%ext (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%eyt (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%kx  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%ky  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%k2  (nc_y,nc_x/2+1), error)

#ifdef DEBUG
   print*, " FFTPACK version of poisson 2d periodic solver "
#endif

   call initdfft(self%fftx, nc_x)
   call initcfft(self%ffty, nc_y)

   call wave_number_vectors(self)

end subroutine initialize_poisson_2d_periodic_fftpack

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fftpack(self,sol,rhs)

   type(sll_t_poisson_2d_periodic_fftpack)       :: self !< self object
   sll_real64, dimension(:,:), intent(in)  :: rhs  !< charge density
   sll_real64, dimension(:,:), intent(out) :: sol  !< electric potential
   sll_int32                               :: nc_x !< number of cells direction x
   sll_int32                               :: nc_y !< number of cells direction y
   sll_int32                               :: i, j

   nc_x = self%nc_x
   nc_y = self%nc_y

   sol(1:nc_x,1:nc_y) = rhs(1:nc_x,1:nc_y)
   do j=1,nc_y
      call dfftf(nc_x, sol(1:nc_x,j), self%fftx%coefd)
   end do

   call transpose_r2c(sol(1:nc_x,1:nc_y), self%rhst)

   do i=1,nc_x/2+1
      call zfftf( nc_y, self%rhst(:,i), self%ffty%coefcd)
   end do

   self%rhst = self%rhst / self%k2

   do i=1,nc_x/2+1
      call zfftb( nc_y, self%rhst(:,i), self%ffty%coefcd )
   end do

   call transpose_c2r(self%rhst, sol(1:nc_x,1:nc_y))

   do j=1,nc_y
      call dfftb( nc_x, sol(1:nc_x,j),  self%fftx%coefd )
   end do

   sol(1:nc_x,1:nc_y) = sol(1:nc_x,1:nc_y) / real(nc_x*nc_y, f64)     ! normalize FFTs

   if (size(sol,1) == nc_x+1) sol(nc_x+1,:) = sol(1,:)
   if (size(sol,2) == nc_y+1) sol(:,nc_y+1) = sol(:,1)

end subroutine solve_potential_poisson_2d_periodic_fftpack

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftpack(self,field_x,field_y,rhs,nrj)
! THIS routine changes the RHS despite its declaration as intent(in)
! It should be fixed !!!
   type(sll_t_poisson_2d_periodic_fftpack) :: self    !< self object
   sll_real64, dimension(:,:), intent(in)  :: rhs     !< charge density
   sll_real64, dimension(:,:), intent(out) :: field_x !< electric field direction x
   sll_real64, dimension(:,:), intent(out) :: field_y !< electric field direction y
   sll_int32                               :: nc_x    !< number of cells direction x
   sll_int32                               :: nc_y    !< number of cells direction y
   sll_int32                               :: i, j
   sll_real64, optional                    :: nrj     !< \f$ \sqrt{e_x^2+e_y^2} \f$

   nc_x = self%nc_x
   nc_y = self%nc_y

   self%rhst = cmplx(0.0_f64,0.0,kind=f64)
   self%ext  = cmplx(0.0_f64,0.0,kind=f64)
   self%eyt  = cmplx(0.0_f64,0.0,kind=f64)
   field_x   = 0.0_f64
   field_y   = 0.0_f64

   do j=1,nc_y
      call dfftf(nc_x, rhs(1:nc_x,j), self%fftx%coefd)
   end do

   call transpose_r2c(rhs(1:nc_x,1:nc_y), self%rhst)

   do i=1,nc_x/2+1
      call zfftf( nc_y, self%rhst(:,i), self%ffty%coefcd)
   end do

   self%ext(1,1) = (0.0_f64,0.0_f64)
   self%eyt(1,1) = (0.0_f64,0.0_f64)
   self%ext = -cmplx(0.0_f64,self%kx/self%k2,kind=f64)*self%rhst
   self%eyt = -cmplx(0.0_f64,self%ky/self%k2,kind=f64)*self%rhst

   do i=1,nc_x/2+1
      call zfftb( nc_y, self%ext(:,i), self%ffty%coefcd )
      call zfftb( nc_y, self%eyt(:,i), self%ffty%coefcd )
   end do

   call transpose_c2r(self%ext, field_x(1:nc_x,1:nc_y))
   call transpose_c2r(self%eyt, field_y(1:nc_x,1:nc_y))

   do j=1,nc_y
      call dfftb( nc_x, field_x(1:nc_x,j), self%fftx%coefd )
      call dfftb( nc_x, field_y(1:nc_x,j), self%fftx%coefd )
   end do

   field_x(1:nc_x,1:nc_y) = field_x(1:nc_x,1:nc_y) / (nc_x*nc_y)
   field_y(1:nc_x,1:nc_y) = field_y(1:nc_x,1:nc_y) / (nc_x*nc_y)

   if (size(field_x,1) == nc_x+1) field_x(nc_x+1,:) = field_x(1,:)
   if (size(field_x,2) == nc_y+1) field_x(:,nc_y+1) = field_x(:,1)
   if (size(field_y,1) == nc_x+1) field_y(nc_x+1,:) = field_y(1,:)
   if (size(field_y,2) == nc_y+1) field_y(:,nc_y+1) = field_y(:,1)

   if (present(nrj)) then 
      nrj=sum(field_x(1:nc_x,1:nc_y)*field_x(1:nc_x,1:nc_y) &
        +field_y(1:nc_x,1:nc_y)*field_y(1:nc_x,1:nc_y))*self%dx*self%dy
   end if

end subroutine solve_e_fields_poisson_2d_periodic_fftpack

subroutine wave_number_vectors(self)

   type(sll_t_poisson_2d_periodic_fftpack) :: self
   sll_int32  :: ik, jk
   sll_int32  :: nc_x, nc_y
   sll_real64 :: kx, ky, kx0, ky0
   
   nc_x = self%nc_x
   nc_y = self%nc_y
   
   kx0 = 2._f64*sll_p_pi/(self%x_max-self%x_min)
   ky0 = 2._f64*sll_p_pi/(self%y_max-self%y_min)
   
   do ik=1,nc_x/2+1
      kx = (ik-1)*kx0
      do jk = 1, nc_y/2
         ky = (jk-1)*ky0
         self%kx(jk,ik) = kx
         self%ky(jk,ik) = ky
      end do
      do jk = nc_y/2+1 , nc_y     
         ky = (jk-1-nc_y)*ky0
         self%kx(jk,ik) = kx
         self%ky(jk,ik) = ky
      end do
   end do
   self%kx(1,1) = 1.0_f64
   
   self%k2 = self%kx*self%kx+self%ky*self%ky

end subroutine wave_number_vectors

!> convert real array to complex and transpose
subroutine transpose_r2c(real_array, comp_array)

   sll_real64, dimension(:,:), intent(in)  :: real_array
   sll_comp64, dimension(:,:), intent(out) :: comp_array
   sll_int32 :: i, j, n1, n2

   n1 = size(real_array,1)
   n2 = size(real_array,2)

   SLL_ASSERT(size(comp_array,1)==n2)
   SLL_ASSERT(size(comp_array,2)==n1/2+1)

   do j=1,n2
      comp_array(j,1) = cmplx(real_array(1,j),0._f64,kind=f64)
      do i=2, n1/2
         comp_array(j,i) = cmplx(real_array(2*i-2,j),real_array(2*i-1,j),kind=f64)
      end do
      comp_array(j,n1/2+1) = cmplx(real_array(n1,j),0._f64,kind=f64)
   end do

end subroutine transpose_r2c

!> convert complex array to real and transpose
subroutine transpose_c2r(comp_array, real_array)

   sll_comp64, dimension(:,:), intent(in)  :: comp_array
   sll_real64, dimension(:,:), intent(out) :: real_array
   sll_int32 :: i, j, n1, n2

   n1 = size(real_array,1)
   n2 = size(real_array,2)

   SLL_ASSERT((n2==size(comp_array,1)))
   SLL_ASSERT((size(comp_array,2)==n1/2+1))

   do j=1,n2
      real_array(1,j) = real(comp_array(j,1),kind=f64)
      do i=2,n1/2
         real_array(2*i-2,j) = real(comp_array(j,i),kind=f64)
         real_array(2*i-1,j) = aimag(comp_array(j,i))
      end do
      real_array(n1,j) = real(comp_array(j,n1/2+1),kind=f64)
   end do

end subroutine transpose_c2r

subroutine initdfft(self,l)

   type(fftclass) :: self
   sll_int32 :: l 
   self%n = l 
   allocate(self%coefd(2*self%n+15))
   call dffti(self%n,self%coefd)

end subroutine initdfft

subroutine initcfft(self,l)

   type(fftclass) :: self
   sll_int32 :: l 
   self%n = l
   allocate(self%coefcd(4*self%n+15))
   call zffti(self%n,self%coefcd)

end subroutine initcfft

#ifdef FFTW

!> Create a sll_o_new solver
!> @return a pointer to the solver derived type
function new_poisson_2d_periodic_fftw( &
   x_min,                              &
   x_max,                              &
   nc_x,                               &
   y_min,                              &
   y_max,                              &
   nc_y,                               &
   error)                              &
   result(self)

  type(sll_t_poisson_2d_periodic_fftw),pointer :: self   !< self object
  sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
  sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
  sll_real64, intent(in)    :: x_min  !< left corner direction x
  sll_real64, intent(in)    :: x_max  !< right corner direction x
  sll_real64, intent(in)    :: y_min  !< left corner direction y
  sll_real64, intent(in)    :: y_max  !< right corner direction y
  sll_int32,  intent(out)   :: error  !< error code

  SLL_ALLOCATE(self, error)
  call initialize_poisson_2d_periodic_fftw( &
          self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

end function new_poisson_2d_periodic_fftw 

!> sll_o_initialize the Poisson solver
subroutine initialize_poisson_2d_periodic_fftw(self, &
  x_min,                                             &
  x_max,                                             &
  nc_x,                                              &
  y_min,                                             &
  y_max,                                             &
  nc_y,                                              &
  error )

  type(sll_t_poisson_2d_periodic_fftw) :: self   !< Self data object
  sll_real64, intent(in)         :: x_min  !< left corner direction x
  sll_real64, intent(in)         :: x_max  !< right corner direction x
  sll_real64, intent(in)         :: y_min  !< left corner direction y
  sll_real64, intent(in)         :: y_max  !< right corner direction y
  sll_int32                      :: error  !< error code
  sll_int32                      :: nc_x   !< number of cells direction x
  sll_int32                      :: nc_y   !< number of cells direction y
  sll_int32                      :: ik, jk
  sll_real64                     :: kx1, kx0, ky0

  fftw_int                       :: sz


  self%nc_x = nc_x
  self%nc_y = nc_y

  self%dx   = (x_max-x_min) / real( nc_x, f64)
  self%dy   = (y_max-y_min) / real( nc_y, f64)

#ifdef DEBUG
  print*, " FFTW version of poisson 2d periodic solver "
#endif

#ifdef FFTW_F2003

  sz = int((nc_x/2+1)*nc_y,C_SIZE_T)

  self%p_rho = fftw_alloc_complex(sz)
  call c_f_pointer(self%p_rho, self%rht, [nc_x/2+1,nc_y])

  self%p_exy = fftw_alloc_complex(sz)
  call c_f_pointer(self%p_exy, self%exy, [nc_x/2+1,nc_y])

  self%p_tmp = fftw_alloc_real(int(nc_x*nc_y,C_SIZE_T))
  call c_f_pointer(self%p_tmp, self%tmp, [nc_x,nc_y])

  self%fw = fftw_plan_dft_r2c_2d(nc_y,nc_x,self%tmp,self%exy,FFTW_ESTIMATE)
  self%bw = fftw_plan_dft_c2r_2d(nc_y,nc_x,self%exy,self%tmp,FFTW_ESTIMATE)

#else

  SLL_ALLOCATE(self%rht(1:nc_x/2+1,1:nc_y),error)
  SLL_ALLOCATE(self%exy(1:nc_x/2+1,1:nc_y),error)
  SLL_ALLOCATE(self%tmp(nc_x,nc_y),error)
  call dfftw_plan_dft_r2c_2d(self%fw,nc_x,nc_y,self%tmp,self%rht,FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_2d(self%bw,nc_x,nc_y,self%rht,self%tmp,FFTW_ESTIMATE)
   
#endif

  SLL_ALLOCATE(self%kx(nc_x/2+1,nc_y), error)
  SLL_ALLOCATE(self%ky(nc_x/2+1,nc_y), error)
  SLL_ALLOCATE(self%k2(nc_x/2+1,nc_y), error)

  kx0 = 2._f64*sll_p_pi/(x_max-x_min)
  ky0 = 2._f64*sll_p_pi/(y_max-y_min)
  
  do ik=1,nc_x/2+1
     kx1 = (ik-1)*kx0
     do jk = 1, nc_y/2
        self%kx(ik,jk) = kx1
        self%ky(ik,jk) = (jk-1)*ky0
     end do
     do jk = nc_y/2+1 , nc_y     
        self%kx(ik,jk) = kx1
        self%ky(ik,jk) = (jk-1-nc_y)*ky0
     end do
  end do

  self%kx(1,1) = 1.0_f64
  self%k2 = self%kx*self%kx+self%ky*self%ky
  self%kx = self%kx/self%k2
  self%ky = self%ky/self%k2

end subroutine initialize_poisson_2d_periodic_fftw

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fftw(self, phi, rho)

  type(sll_t_poisson_2d_periodic_fftw)          :: self !< self data object
  sll_real64, dimension(:,:), intent(in)  :: rho  !< charge density
  sll_real64, dimension(:,:), intent(out) :: phi  !< electric potential
  sll_int32                               :: nc_x !< number of cells direction x
  sll_int32                               :: nc_y !< number of cells direction y

  nc_x = self%nc_x
  nc_y = self%nc_y

  self%tmp = rho(1:nc_x,1:nc_y)
  call fftw_execute_dft_r2c(self%fw, self%tmp, self%rht)

  self%rht = self%rht / self%k2

  call fftw_execute_dft_c2r(self%bw, self%rht, self%tmp)

  phi(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)   
  
  !Node centered case
  if(size(phi,1) == nc_x+1) phi(nc_x+1,:) = phi(1,:)
  if(size(phi,2) == nc_y+1) phi(:,nc_y+1) = phi(:,1)

end subroutine solve_potential_poisson_2d_periodic_fftw

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftw(self,e_x,e_y,rho,nrj)

  type(sll_t_poisson_2d_periodic_fftw),intent(inout) :: self !< Self data object
  sll_real64, dimension(:,:),    intent(in)    :: rho  !< Charge density
  sll_real64, dimension(:,:),    intent(out)   :: e_x  !< Electric field x
  sll_real64, dimension(:,:),    intent(out)   :: e_y  !< Electric field y
  sll_real64, optional                         :: nrj  !< Energy 

  sll_int32  :: nc_x, nc_y
  sll_real64 :: dx, dy

  nc_x = self%nc_x
  nc_y = self%nc_y

  self%tmp = rho(1:nc_x,1:nc_y)
  call fftw_execute_dft_r2c(self%fw, self%tmp, self%rht)

  self%exy(1,1) = (0.0_f64,0.0_f64)
  self%exy = -cmplx(0.0_f64,self%kx,kind=f64)*self%rht
  call fftw_execute_dft_c2r(self%bw, self%exy, self%tmp)
  e_x(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)

  self%exy(1,1) = (0.0_f64,0.0_f64)
  self%exy = -cmplx(0.0_f64,self%ky,kind=f64)*self%rht
  call fftw_execute_dft_c2r(self%bw, self%exy, self%tmp)

  e_y(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)

  !Node centered case
  if (size(e_x,1) == nc_x+1) e_x(nc_x+1,:) = e_x(1,:)
  if (size(e_x,2) == nc_y+1) e_x(:,nc_y+1) = e_x(:,1)
  if (size(e_y,1) == nc_x+1) e_y(nc_x+1,:) = e_y(1,:)
  if (size(e_y,2) == nc_y+1) e_y(:,nc_y+1) = e_y(:,1)

  if (present(nrj)) then 
     dx = self%dx
     dy = self%dy
     nrj=sum(e_x(1:nc_x,1:nc_y)*e_x(1:nc_x,1:nc_y) &
       +e_y(1:nc_x,1:nc_y)*e_y(1:nc_x,1:nc_y))*dx*dy
  end if

end subroutine solve_e_fields_poisson_2d_periodic_fftw

!> Delete the Poisson object
subroutine delete_poisson_2d_periodic_fftw(self)

  type(sll_t_poisson_2d_periodic_fftw) :: self

#ifdef FFTW_F2003
  call fftw_free(self%p_rho)
  if (c_associated(self%p_exy)) call fftw_free(self%p_exy)
  if (c_associated(self%p_tmp)) call fftw_free(self%p_tmp)
  if (c_associated(self%p_rho)) call fftw_free(self%p_rho)
#endif

  call fftw_destroy_plan(self%fw)
  call fftw_destroy_plan(self%bw)

end subroutine delete_poisson_2d_periodic_fftw

#endif /* FFTW */
 
end module sll_m_poisson_2d_periodic
