!> \brief Implements a conservative and consistant algorithm based on an idea by Lin and Rood
!> for the solution in 2D of \f$\partial_t u + \partial_x (au) = 0\f$, with a divergence free.
!> This algorithm is based on an advective backward semi-Lagrangian method for predicting 
!> the fluxes in a conservative finite difference scheme which leads to conservativity
!> Consitancy (conservation of constant states) is achieved by prediction using 
!> the advective form and by using a reconstruction of the stream function associated
!> to a identical to the reconstruction used for the fluxes
!> It is described in Qiu - Sonnendrucker

module sll_linrood
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_2d.h"

  use sll_constants
  use weno_recon
  use weno_interp
  use sll_advection_field
  use distribution_function
  implicit none

  type linrood_plan
     type (weno_recon_1D), pointer   :: recon_eta1
     type (weno_recon_1D), pointer   :: recon_eta2
     type (weno_interp_1D), pointer  :: interp_eta1
     type (weno_interp_1D), pointer  :: interp_eta2
     sll_int32  :: nc_eta1, nc_eta2
     sll_int32  :: order
     sll_real64, dimension(:,:), pointer   :: dist_func_2d, f_temp
     sll_real64, dimension(:,:), pointer   :: advfield_1, advfield_2
     sll_real64, dimension(:,:), pointer   :: flux
  end type linrood_plan

contains
    
!> initialize opaque pointer that contains information used by the algorithm
!> \param[in] dist_func_2D 2D distribution function object
!> \return pointer to opaque data type 
!> Only periodic boundary conditions are implemented up to now
  subroutine new_linrood_plan(this, dist_func_2d)
    type (linrood_plan) :: this
    type (sll_distribution_function_2D)  :: dist_func_2D 
    !type (weno_recon_1d) :: recon_eta1, recon_eta2
    !type (weno_interp_1d) :: interp_eta1, interp_eta2
    type(sll_logical_mesh_2d), pointer :: mesh
    sll_int32  :: ierr
    sll_int32  :: i_weno
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    
    ! get dimensions
!#ifdef STDF95
!    nc_eta1    = GET_FIELD_NC_ETA1( dist_func_2D%extend_type ) 
!    nc_eta2    = GET_FIELD_NC_ETA2( dist_func_2D%extend_type ) 
!#else
    mesh => dist_func_2d%transf%get_logical_mesh()
    nc_eta1    = mesh%num_cells1
    nc_eta2    = mesh%num_cells2

    ! save dimensions in plan
    this%nc_eta1 = nc_eta1  
    this%nc_eta2 = nc_eta2
!#endif

    ! allocate arrays
    SLL_ALLOCATE(this%dist_func_2d(nc_eta1,nc_eta1),ierr)
    SLL_ALLOCATE(this%flux(nc_eta1,nc_eta1),ierr)
    SLL_ALLOCATE(this%advfield_1(nc_eta1,nc_eta1),ierr)
    SLL_ALLOCATE(this%advfield_2(nc_eta1,nc_eta1),ierr)
    SLL_ALLOCATE(this%f_temp(nc_eta1,nc_eta1),ierr)

    ! intialize order of interpolation and reconstruction
    this%order = 5
    
    i_weno = 0  ! WENO (0 for WENO)
    
    ! initialize interpolators and reconstructors
    this%recon_eta1 =>  new_WENO_recon_1D(nc_eta1, 0.0_f64, 1.0_f64, this%order, i_weno)
    this%recon_eta2 =>  new_WENO_recon_1D(nc_eta2, 0.0_f64, 1.0_f64, this%order, i_weno)
    this%interp_eta1 =>  new_WENO_interp_1D(nc_eta1, 0.0_f64, 1.0_f64, this%order, i_weno) 
    this%interp_eta2 =>  new_WENO_interp_1D(nc_eta2, 0.0_f64, 1.0_f64, this%order, i_weno) 
    

  end subroutine new_linrood_plan


 subroutine linrood_step(plan, dist_func_2d, advfield, t, deltat)
   type (linrood_plan)                  :: plan
   type (sll_distribution_function_2d)  :: dist_func_2D  
   ! advection field defined by its stream function
   type(scalar_field_2d), intent(in) :: advfield
   sll_real64  :: deltat
   sll_real64  :: t
   type(sll_logical_mesh_2d), pointer :: mesh
   ! local variables
   sll_int32 :: i1, i2
   sll_int32 :: nc_eta1, nc_eta2
   sll_real64 :: eta1, eta2
   sll_real64 :: delta_eta1, delta_eta2
   sll_real64, dimension(max(plan%nc_eta1, plan%nc_eta2)) :: adt
   sll_real64, dimension(max(plan%nc_eta1, plan%nc_eta2)) :: aux_in, aux_out
    
   ! get dimensions
!#ifdef STDF95
!   nc_eta1    = GET_FIELD_NC_ETA1( dist_func_2D%extend_type ) 
!   nc_eta2    = GET_FIELD_NC_ETA2( dist_func_2D%extend_type )
!#else
   mesh => dist_func_2D%transf%get_logical_mesh()
   nc_eta1    = mesh%num_cells1
   nc_eta2    = mesh%num_cells2
!#endif
   delta_eta1 = 1.0_8 / nc_eta1
   delta_eta2 = 1.0_8 / nc_eta2

   ! Compute the two components of advection field by reconstruction
   !----------------------------------------------------------------
   ! First component: a_1 = d H/ d eta_2
   do i1 = 1, nc_eta1
      do i2 = 1, nc_eta2
!#ifdef STDF95
!         aux_in(i2) =  FIELD_2D_AT_I( advfield%extend_type, i1, i2 )
!#else
         aux_in(i2) =  FIELD_2D_AT_I( advfield, i1, i2 )
!#endif
      end do
      call FD_WENO_recon_1D(plan%recon_eta2, nc_eta2, aux_in, aux_out)
      do i2 = 1, nc_eta2
         plan%advfield_1( i1, i2 ) = cos(t*sll_pi/1.5_f64)*aux_out(i2)
      end do
   end do
   ! Second component: a_2 = - d H/ d eta_1
   do i2 = 1, nc_eta2
      do i1 = 1, nc_eta1
!#ifdef STDF95
!         aux_in(i1) =  FIELD_2D_AT_I( advfield%extend_type, i1, i2 )
!#else
         aux_in(i1) =  FIELD_2D_AT_I( advfield, i1, i2 )
!#endif
      end do
      call FD_WENO_recon_1D(plan%recon_eta1, nc_eta1, aux_in, aux_out)
      do i1 = 1, nc_eta1
         plan%advfield_2( i1, i2 ) =  - cos(t*sll_pi/1.5_f64)*aux_out(i1)
      end do
   end do

   ! First step of algorithm
   ! Advance distribution function (not multiplied by Jacobian on delta/2 using 
   ! first order split backward semi-Lagrangian method
   ! For conservative method sqrt(g) f is stored at cell centers
   !---------------------------------------------------------------------------
   
   ! compute f by dividing by Jacobian
   eta2 = 0.5_8 * delta_eta2
   do i2 = 1, nc_eta2     
      eta1 = 0.5_8 * delta_eta1
      do i1 = 1, nc_eta1
         ! FIELD_2D_AT_I( dist_func_2d, i1, i2 ) = dist_func_2d%data(i1,i2)
         plan%f_temp(i1,i2) = FIELD_2D_AT_I( dist_func_2d, i1, i2 ) / &
              FIELD_2D_JACOBIAN_AT_I( dist_func_2d, eta1, eta2 )
!#endif
         eta1 = eta1 + delta_eta1
      end do
      eta2 = eta2 + delta_eta2
   end do
   
   ! BSL Advection in eta1
   do i2 = 1, nc_eta2
      ! extract 1d distribution function and displacement field
      do i1 = 1, nc_eta1
         aux_in(i1) = plan%f_temp(i1,i2)
         adt(i1) = plan%advfield_1(i1,i2) * 0.5_8 * deltat
      end do
      ! do 1d bsl advection
      call bsl_1d(aux_in, aux_out, adt, nc_eta1, plan%interp_eta1)
      ! copy line back to 2d distribution function
      do i1 = 1, nc_eta1
         plan%f_temp(i1,i2) = aux_out(i1) 
      end do
   end do

   ! BSL Advection in eta2
   do i1 = 1, nc_eta1
      ! extract 1d distribution function and displacement field
      do i2 = 1, nc_eta2
!         aux_in(i2) = plan%dist_func_2d(i1,i2)
         aux_in(i2) = plan%f_temp(i1,i2)
         adt(i2) = plan%advfield_2(i1,i2) * 0.5_8 * deltat
      end do
      call bsl_1d(aux_in, aux_out, adt, nc_eta2, plan%interp_eta2)
      do i2 = 1, nc_eta2
         plan%f_temp(i1,i2) = aux_out(i2) 
      end do
   end do
!!$do i1 = 1, nc_eta1
!!$do i2 = 1, nc_eta2
!!$FIELD_2D_AT_I( dist_func_2d, i1, i2 ) = plan%f_temp(i1,i2)
!!$end do
!!$end do
! Second step of algorithm
! conservative update from dimension 1

   do i1 = 1, nc_eta1
      do i2 = 1, nc_eta2
         aux_in(i2) = plan%f_temp(i1,i2)*plan%advfield_1(i1,i2)
      end do 
         call deriv_f_1d(aux_in, aux_out, nc_eta2, plan%recon_eta2)
      do i2 = 1, nc_eta2
         FIELD_2D_AT_I( dist_func_2d, i1, i2 ) = &
              FIELD_2D_AT_I( dist_func_2d, i1, i2 ) - deltat*aux_out(i2)
       enddo
  end do

! conservative update from dimension 2
 do i2 = 1, nc_eta2
   do i1 = 1, nc_eta1
         aux_in(i1) = plan%f_temp(i1,i2)*plan%advfield_2(i1,i2)
   end do 
         call deriv_f_1d(aux_in, aux_out, nc_eta1, plan%recon_eta1)
      do i1 = 1, nc_eta1
         FIELD_2D_AT_I( dist_func_2d, i1, i2 ) = &
              FIELD_2D_AT_I( dist_func_2d, i1, i2 ) - deltat*aux_out(i1)
       enddo
  end do

  end subroutine linrood_step

  !> Does a 1d first order semi-Lagrangian advection of f in the advection field a
  subroutine bsl_1d(f_in, f_out, adt, nc, interp)
    sll_real64, dimension(:) :: f_in, f_out  ! function to be advected
    sll_real64, dimension(:) :: adt  ! advection field times time-step
    sll_int32 :: nc       ! number of cells
    type (weno_interp_1D), pointer       :: interp

    ! local variables
    sll_real64, dimension(nc)   :: origin
    sll_real64  :: delta, eta
    sll_int32   :: i
    
    delta = 1.0_8 / nc
    eta = 0.5_8 * delta
    do i = 1, nc
       ! periodic domain of period 1
       origin(i) = modulo(eta - adt(i),1.0_f64)
      ! print*, i, eta - adt(i), origin(i)
       eta = eta + delta
    end do

   call interpolate_WENO_1D( interp, nc, f_in, origin, f_out)

  end subroutine bsl_1d
    
  !> Does a 1d conservative update of f: the derivative of flux functions
  subroutine deriv_f_1d(f_in, f_out, nc, recon)
    sll_real64, dimension(:) :: f_in, f_out  ! function to be differentiate 
    sll_int32 :: nc       ! number of cells
    type (WENO_recon_1D), pointer       :: recon
    call FD_WENO_recon_1D(recon, nc, f_in, f_out)
  end subroutine deriv_f_1d
    
    

end module sll_linrood
