!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
!
! MODULE: sll_m_pic_1d_particle_loading
!
! DESCRIPTION:
!> @author Jakob Ameres
!> @brief Loading utility for particle-in-cell method in 1d
!> @details This module provides tools as distribution sampling for the PIC loading mechanism.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_m_pic_1d_particle_loading
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_collective, only: &
    sll_s_collective_barrier, &
    sll_t_collective_t, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_charge, &
    sll_p_mass, &
    sll_p_epsilon_0, &
    sll_p_kb, &
    sll_p_kx, &
    sll_p_pi, &
    sll_p_proton_mass

  use sll_m_particle_1d_description, only: &
    sll_t_particle_1d_group

  use sll_m_prob, only: &
    sll_s_normal_cdf_inv

  use sll_m_sobol, only: &
    sll_s_i8_sobol_generate

  implicit none

  public :: &
    sll_f_control_variate_xv, &
    sll_v_enable_deltaf, &
    sll_s_load_particle_species, &
    sll_v_num_species, &
    sll_s_set_loading_parameters, &
    sll_s_initialize_intrinsic_mpi_random, &
    sll_s_pic1d_ensure_boundary_conditions, &
    sll_s_pic1d_ensure_boundary_conditions_species, &
    sll_f_pic1d_ensure_periodicity, &
    sll_p_pic1d_testcase_bumpontail, &
    sll_p_pic1d_testcase_ionbeam, &
    sll_p_pic1d_testcase_ionbeam_electrons, &
    sll_p_pic1d_testcase_landau, &
    sll_p_pic1d_testcase_quiet

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !Definitions for different loadings
    sll_int32, parameter :: sll_p_pic1d_testcase_landau=1
    sll_int32, parameter :: sll_p_pic1d_testcase_bumpontail=2
    sll_int32, parameter :: SLL_PIC1D_TESTCASE_TWOSTREAM=3
    sll_int32, parameter :: sll_p_pic1d_testcase_ionbeam=4
    sll_int32, parameter :: sll_p_pic1d_testcase_quiet=5
    sll_int32, parameter :: sll_p_pic1d_testcase_ionbeam_electrons=6

    !Definitions for different Loaders
    sll_int32, parameter :: SLL_PIC1D_SAMPLER_SOBOL= 1
    sll_int32, parameter :: SLL_PIC1D_SAMPLER_HAMMERSLY= 2
    sll_int32, parameter :: SLL_PIC1D_SAMPLER_SYSTEM= 3

    !
    sll_int32, parameter :: SLL_PIC1D_SAMPLING_DIST_UNIFORM= 1
    sll_int32, parameter :: SLL_PIC1D_SAMPLING_DIST_IMPORTANT= 2

    !Flags can be combined
    sll_int32, parameter :: SLL_PIC1D_LOADER_NONE       = 0
    sll_int32, parameter :: SLL_PIC1D_LOADER_LANDAU   = 2**1
    sll_int32, parameter :: SLL_PIC1D_LOADER_TWOSTREAM   = 2**2
    sll_int32, parameter :: SLL_PIC1D_LOADER_BUMPONTAIL   = 2**3
    sll_int32, parameter :: SLL_PIC1D_SAMPLING_INVERSE   = 2**4

    !    sll_real64, parameter :: sll_p_kb = 1.3806488D-23
    !    sll_real64, parameter :: PLASMA_FREQUENCY=sqrt(sll_p_kb*T/sll_p_mass)

    integer :: ierr
    sll_int32 :: coll_rank, coll_size

    !Parameters for different loading types
    !sll_int32 ::  pic1d_testcase = sll_p_pic1d_testcase_landau
    sll_int32 :: sll_v_num_species
    sll_real64 :: landau_alpha=0.01_f64
    sll_real64 :: landau_mode=0.4_f64

    !sll_real64 :: bumpontail_a=0.4_f64
    sll_real64 :: bumpontail_a=0.04_f64
    sll_real64 :: bumpontail_v0=4.0_f64
    sll_real64 :: bumpontail_sigma=0.5_f64
    sll_real64 :: plasma_size=0.25_f64 !Relative size of plasma
    !sll_real64 :: sll_pic_boundary_condition=real(sll_p_periodic,f64)

    sll_int32 :: numberof_streams=1
    logical  :: sll_v_enable_deltaf=.FALSE.

    sll_real64 :: interval_length
    sll_real64 :: interval_a, interval_b

    !Probability
    abstract interface
        function sll_pdf_1d(x) result(p)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(size(x)) :: p
        endfunction
    endinterface

    abstract interface
        function sll_pdf_1d1t(x,t) result(p)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(:),intent(in) :: t
            sll_real64, dimension(size(x),size(t)) :: p
        endfunction
    endinterface

    !<Ony mostly supposes that the spatial and velocity density are independent
    !<which means they factorize. So one does not or should not necessarily
    !<ever need this interface. We keep it here for sake of completeness
    abstract interface
        function sll_pdf_1d1v(x,v) result(p)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(:) ,intent(in):: v
            sll_real64, dimension(size(x)) :: p
        endfunction
    endinterface

    abstract interface
        function sll_pdf_1d1v1t(x,v,t) result(p)
            use sll_m_working_precision
            sll_real64, dimension(:),intent(in) :: x
            sll_real64, dimension(:),intent(in) :: v
            sll_real64, dimension(:),intent(in) :: t
            sll_real64, dimension(size(x),size(t)) :: p
        endfunction
    endinterface

    !Probability density
    !<Here we suppose that f(x,v)=f(x)*f(v)
    procedure (sll_pdf_1d), pointer :: sampling_dist_x
    procedure (sll_pdf_1d), pointer :: sampling_dist_v
    procedure (sll_pdf_1d), pointer :: initial_dist_x
    procedure (sll_pdf_1d), pointer :: initial_dist_v
    procedure (sll_pdf_1d), pointer :: control_variate_x
    procedure (sll_pdf_1d), pointer :: control_variate_v

contains

  subroutine sll_s_set_loading_parameters( landau_alpha_user, landau_mode_user, numberof_streams_user )
    sll_real64, intent(in) :: landau_alpha_user
    sll_real64, intent(in) :: landau_mode_user
    sll_int32 , intent(in) :: numberof_streams_user

    !logical, intent(in)::enable_deltaf_user
    landau_mode      = landau_mode_user
    landau_alpha     = landau_alpha_user
    numberof_streams = numberof_streams_user
    !sll_v_enable_deltaf=enable_deltaf_user

  end subroutine sll_s_set_loading_parameters

  function gaussianrnd( mu, sigma ) result(x)
    sll_real64, intent(in) :: mu    !< mean
    sll_real64, intent(in) :: sigma !< standard deviation
    sll_real64 :: x

    sll_real64 :: u1, u2, r, phi, y
#ifdef DEBUG
    sll_real64 :: dummy
    dummy = sigma+mu
#endif

    call random_number(u1)
    call random_number(u2)

    r   = sqrt(-2.0_f64*log(u1))
    phi = sll_p_kx*u2
    x   = r*cos(phi)
    y   = r*sin(phi)

    ! box muller wiener
!    sll_real :: r1, r2, x ,y
!
!
!    call random_number(r1)
!    call random_number(r2)
!
!    r1= 1.0-r1
!    r1 = -alog(real(r1))
!    r1 = sqrt(2.0*r1)
!    r2 = 2.0*sll_p_pi*r2
!
!    x  = r1*cos(r2)- mu
!    y  = r1*sin(r2) -mu
!
  end function gaussianrnd

!
!function box_mueller_antithetic( mu , sigma , uniform_random_numbers  ) RESULT(gaussian_random_numbers)
!#include "sll_m_working_precision.h"
!    use sll_m_constants
!    IMPLICIT NONE
!    sll_real , intent(in):: mu    !< mean
!    sll_real, intent(in) :: sigma !< standard deviation
!    sll_real64, dimension(:), intent(in) :: uniform_random_numbers
!    sll_real64, dimension(2*size(uniform_random_numbers)) :: gaussian_random_numbers
!
!    integer :: idx
!
!    !! SLL_ASSERT( mod(size(uniform_random_numbers),2)==0)
!
!    ! Box Muller Wiener
!    sll_real :: R1, R2, X ,Y
!    do idx=1,size(uniform_random_numbers)
!        R1=gaussian_random_numbers(idx)
!        R2=1-gaussian_random_numbers(idx)
!        R1= 1.0_f64-R1
!        R1 = -ALOG(real(R1))
!        R1 = SQRT(2.0*R1)
!        R2 = 2.0*sll_p_pi*R2
!        gaussian_random_numbers(2*idx-1)= R1*COS(R2)- mu
!        gaussian_random_numbers(2*idx) = R1*SIN(R2) -mu
!    enddo
!
!endfunction


  function gaussian_from_rnd( mu, sigma, uniform_random_numbers ) result(gaussian_random_numbers)
    sll_real64, intent(in) :: mu    !< mean
    sll_real64, intent(in) :: sigma !< standard deviation
    sll_real64, intent(in) :: uniform_random_numbers(:)

    sll_real64 :: gaussian_random_numbers(size(uniform_random_numbers))
    sll_int32  :: idx

    !! SLL_ASSERT( mod(size(uniform_random_numbers),2)==0)

    ! Box Muller Wiener
    sll_real64 :: r1, r2!, x ,y
#ifdef debug
    sll_real64 :: dummy
    dummy = sigma+mu
#endif

    do idx=1,size(uniform_random_numbers)/2
        r1=uniform_random_numbers(2*idx-1)
        r2=uniform_random_numbers(2*idx)
        r1= 1.0_f64-r1
        r1 = -log(real(r1,f64))
        r1 = sqrt(2.0*r1)
        r2 = 2.0*sll_p_pi*r2
        gaussian_random_numbers(2*idx-1)= r1*cos(r2)- mu
        gaussian_random_numbers(2*idx) = r1*sin(r2) -mu
    end do
  end function gaussian_from_rnd


  function sll_normal_prb_kernel( mu, sigma, x  ) result(fx)
    sll_real64, intent(in):: mu    !< mean
    sll_real64, intent(in) :: sigma !< standard deviation
    sll_real64, intent(in) :: x(:)

    sll_real64 :: fx(size(x))
    !integer :: idx

    fx = (1.0_f64/(sqrt(2.0_f64*sll_p_pi)*sigma)) * &
         exp( -(x-mu)**2/(2.0_f64*sigma**2) )
  end function sll_normal_prb_kernel

!function sll_pdf_1d_signature( x)  result(y)
!    sll_real64, dimension(:) ::x
!    sll_real64, dimension(size(x)) :: y
!endfunction
!
!function sll_pdf_1d1v_signature( x ,y)  result(y)
!    sll_real64, dimension(:) ::x
!    sll_real64, dimension(size(x)) :: y
!endfunction

!<Dummy function for probability density in
!function sll_pdf_1d1v( x,v)  result(y)
!    sll_real64, dimension(:) ::x
!    sll_real64, dimension(:) ::v
!    SLL_ASSERT(size(x)==size(v))
!    sll_real64, dimension(size(x)) :: y
!endfunction
!
!function sll_pdf_1d( pdf_1d1v ,v)  result(y)
!    sll_real64, dimension(:) ::x
!    procedure(sll_pdf_1d1v) :: pdf_1d1v
!    y=pdf_1d1v( )
!    sll_real64, dimension(size(x)) :: y
!endfunction
!
!function sll_pdf_1v( x,v)  result(y)
!    sll_real64, dimension(:) ::x
!    sll_real64, dimension(:) ::v
!    SLL_ASSERT(size(x)==size(v))
!    sll_real64, dimension(size(x)) :: y
!endfunction
!

  function sll_normal_rnd( mu, sigma, uniformrandom ) result( normal)
    sll_real64, intent(in)    :: mu
    sll_real64, intent(in)    :: sigma
    sll_real64, intent(inout) :: uniformrandom(:)

    sll_real64 :: normal(size(uniformrandom))
    sll_int32  :: idx, N

    N = size(uniformrandom)

    SLL_ASSERT(sigma>0.0_f64)
    do idx=1,N
        if (uniformrandom(idx)==0.0_f64) uniformrandom(idx)=0.00001_f64
        call sll_s_normal_cdf_inv( uniformrandom(idx), &
            0.0_f64 , 1.0_f64, normal(idx))
    end do

    normal = normal*sigma + mu

  end function sll_normal_rnd


  function sll_normal_landaudamp_prb_kernel( mu, sigma, alpha, k, x ) result(fx)
    sll_real64, intent(in) :: mu    !< mean
    sll_real64, intent(in) :: sigma !< standard deviation
    sll_real64, intent(in) :: alpha !< landau damping factor
    sll_real64, intent(in) :: k     !< Wave vector k=2pi/lambda
    sll_real64, intent(in) :: x(:)

    sll_real64 :: fx(size(x))
    !integer :: idx

    fx = (1.0_f64/(sqrt(2.0_f64*sll_p_pi)*sigma)) * &
      exp( (x-mu)/(2.0_f64*sigma**2))*(1.0_f64 + alpha*cos(k*x))

  end function sll_normal_landaudamp_prb_kernel


  function sll_maxwellboltzmann1d( m, T, velocity ) result(prob)
    sll_real64, intent(in) :: T           !< temperature in K
    sll_real64, intent(in) :: m           !< particle mass in kg
    sll_real64, intent(in) :: velocity(:) !< velocity in m/s

    sll_real64 :: prob(size(velocity))

    !> @param Boltzmann constant (def) J/K
    sll_int64  :: idx, ndx
    sll_real64 :: scale_par

    scale_par=sqrt(sll_p_kb*T/m)

    ndx = int(size(velocity),i64)
    forall (idx=1:ndx)
        prob(idx)=(sqrt(2*sll_p_pi )*scale_par) *  &
            4.0_f64*sll_p_pi* velocity(idx)**2 &
            *exp(-  velocity(idx)**2 /(2*scale_par**2))
    end forall
  end function sll_maxwellboltzmann1d


!< takes M random numbers
  function birdsall_normal_1d( uniform_random_numbers, M_user ) &
      result(normal_random_numbers)
    sll_real64, intent(in)           :: uniform_random_numbers(:)
    sll_int32 , intent(in), optional :: M_user

    sll_int32  :: M
    sll_real64 :: normal_random_numbers(size(uniform_random_numbers))
    sll_int32  :: idx

    if (present(M_user)) then
        M=M_user
    else
        M=12
    endif

    do idx=1, size(normal_random_numbers)
        normal_random_numbers(idx) = &
          sqrt(M/12.0_f64) * sum( uniform_random_numbers(1:idx) - M/2.0_f64 )
    enddo

  end function birdsall_normal_1d


  !<Loads the number velocity streams in at the given offset
  !<Default behavoiour is to load 2 streams at -1 and 1
  !< stream_offsets_user = (/ -5,5 /)
  subroutine sll_pic1d_load_stream( uniform_random_numbers, particlespeed, &
      numberofstreams, stream_offsets_user  )
    sll_real64, intent(inout)        :: uniform_random_numbers(:)
    sll_real64, intent(inout)        :: particlespeed(:)
    sll_int32 , intent(in)           :: numberofstreams
    sll_real64, intent(in), optional :: stream_offsets_user(:)
    !sll_real64, intent(in), optional :: stream_widths_user(:)

    sll_real64, allocatable :: stream_offsets(:)
    sll_real64, allocatable :: stream_widths(:)
    integer :: nparticles
    integer :: i, stream_idx

    nparticles=size(uniform_random_numbers)
    SLL_ASSERT(nparticles==size(particlespeed) )

    SLL_CLEAR_ALLOCATE(stream_offsets(1:numberofstreams), ierr)
    SLL_CLEAR_ALLOCATE(stream_widths(1:numberofstreams), ierr)

    if  (present(stream_offsets_user)) then
        SLL_ASSERT(numberofstreams==size(stream_offsets_user))
        stream_offsets=stream_offsets_user
    else
        if (numberofstreams>1) then
            do i=1, numberofstreams
                stream_offsets(i)= (i-1.0_f64)/(numberofstreams-1) -0.5_f64
            enddo
            !stream_offsets=stream_offsets*1.0_f64
            stream_widths=(1.0_f64)/(numberofstreams**2)
        else
            stream_offsets(1)=0.0_f64
            stream_widths(1)=1.0_f64
        endif
    endif

    do stream_idx=1,numberofstreams
        do i=(nparticles)*(stream_idx-1)/numberofstreams +1, (nparticles)*(stream_idx)/numberofstreams
            if (uniform_random_numbers(i)==0.0_f64) uniform_random_numbers(i)=0.00001_f64
            call sll_s_normal_cdf_inv( uniform_random_numbers(i) , 0.0_f64 ,&
                stream_widths(stream_idx) , particlespeed(i) )
            particlespeed(i) = particlespeed(i) +stream_offsets(stream_idx)
        enddo
    enddo

    SLL_DEALLOCATE_ARRAY(stream_offsets, ierr )
    SLL_DEALLOCATE_ARRAY(stream_widths, ierr )

  end subroutine sll_pic1d_load_stream


  !<Loads the velocity component of bump on tail as
  !<$ f_0(,v)=
  !< \frac{1}{1+a}\left(  \frac{1}{\sqrt{2\pi}} e^{-\frac{v^2}{2}} +
  !< \frac{a}{\sqrt{2\pi} \sigma}  e^{-\frac{ (v-v_0)^2}{2\sigma^2}}
  !< \right)
  subroutine sll_pic1d_load_bumpontail_velocity( uniform_random_numbers, &
      particlespeed, a, v0, sigma )
    sll_real64, intent(inout) :: uniform_random_numbers(:)
    sll_real64, intent(inout) :: particlespeed(:)
    sll_real64, intent(in)    :: a
    sll_real64, intent(in)    :: v0
    sll_real64, intent(in)    :: sigma

    sll_real64 :: percentage
    sll_int32  :: np !number of particles
    sll_int32  :: i

    np = size(particlespeed)

    SLL_ASSERT(size(uniform_random_numbers)==np)
    particlespeed=0.0_f64
    percentage=1.0_f64/(1.0_f64+a)
    do i=1,floor(percentage*np*1.0_f64  )
        if (uniform_random_numbers(i)==0.0_f64) uniform_random_numbers(i)=0.00001_f64
        call sll_s_normal_cdf_inv( uniform_random_numbers(i) , 0.0_f64 ,&
            1.0_f64 , particlespeed(i) )
    enddo

    do i=floor(percentage*np*1.0_f64  )+1,np
        if (uniform_random_numbers(i)==0.0_f64) uniform_random_numbers(i)=0.00001_f64
        call sll_s_normal_cdf_inv( uniform_random_numbers(i) , v0,&
            sigma , particlespeed(i) )
    enddo

  end subroutine sll_pic1d_load_bumpontail_velocity


  function ionbeam(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    y=0.0_f64
    where (x>=interval_a .and. x-interval_a<=(interval_b-interval_a)/3.0_f64)
      y = 3.0_f64
    end where
  end function ionbeam


!subroutine  load_particles (nparticles, interval_a_user, interval_b_user,steadyparticleposition, &
    !        particleposition, particlespeed, particleweight, particleweight_constant, particle_qm)
!#include "sll_m_working_precision.h"
!#include "sll_m_memory.h"
!    use sll_m_constants
!    implicit none
!    integer, intent(in) :: nparticles
!    sll_real64, intent(in) ::interval_a_user
!    sll_real64, intent(in) :: interval_b_user
!    sll_real64, DIMENSION(:), intent(inout):: particleposition
!    sll_real64, DIMENSION(:), intent(inout) :: particlespeed
!    sll_real64, DIMENSION(:), intent(inout):: particleweight
!    sll_real64, DIMENSION(:), intent(inout):: particleweight_constant
!    sll_real64, DIMENSION(:), intent(inout):: particle_qm
!    sll_real64 :: funlandau,x
!
!    !Ionbeam scenario
!    sll_real64 :: electron_ratio, hplus_ratio, hminus_ratio
!    sll_real64 :: electron_temp, hplus_temp, hminus_temp
!    sll_int32 :: idx_up, idx_low
!
!    sll_real64, DIMENSION(:), intent(inout) :: steadyparticleposition
!    integer :: ierr
!    sll_real64 :: mu, sigma
!    integer :: i
!    real ( kind = 8 ) :: tmp
!    sll_real64 :: maxwellian_a=1.0_f64
!    sll_real64, DIMENSION(:,:), allocatable :: phasespace
!
!
!    !!sll_real64 :: landau_damping=0.01_f64
!    !Plasma temperature in Kelvin 150000273.15
!    !maxwellian_a=sqrt(sll_p_kb*150000273.15_f64/sll_p_mass)
!    !funlandau(x)=cos(x)
!
!    interval_a=interval_a_user
!    interval_b=interval_b_user
!    interval_length=interval_b-interval_a
!
!    !print *, funlandau(particleposition(1))
!
!    mu=0.0_f64
!    sigma=1.0_f64
!
!    !!!!!!!!!REMOVE THIS FROM HERE!!!! LATER WE HAVE ONE PIC MODULE
!    coll_rank = sll_f_get_collective_rank( sll_v_world_collective )
!    coll_size = sll_f_get_collective_size( sll_v_world_collective )
!    call sll_s_collective_barrier(sll_v_world_collective)
!
!    selectcase (pic1d_testcase)
!        !######################################################################
!        case(sll_p_pic1d_testcase_landau)
!            !Generate random numbers
!            SLL_CLEAR_ALLOCATE(phasespace(1:3,1:nparticles),ierr)
!            call sll_s_i8_sobol_generate ( 3_f64, nparticles, coll_rank*nparticles, phasespace )
!                particle_qm=-1.0_f64 !Electrons
!
!            particleweight_constant=0
!            if (sll_v_enable_deltaf .eqv. .TRUE.) then
!                !Landau damping
!                !Load deltaf and set weights
!
!                call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, 1)
!
!                sampling_dist_v=>sll_pic1d_normalPDF
!                initial_dist_x=>sll_pic_1d_landaudamp_PDF
!                initial_dist_v=>sll_pic1d_normalPDF
!                control_variate_v=>initial_dist_v
!
!
!!                control_variate_x=>sll_pic1d_constantPDFx
!!                sampling_dist_x=>sll_pic1d_abscosPDFlandau
!!                particleposition=sll_pic1d_abscosiCDF(interval_length,landau_mode,phasespace(3,:))
!
!
!
!                !control_variate_x=>sll_pic1d_constantPDFx
!                control_variate_x=>sll_pic_1d_landaudamp_PDF
!                sampling_dist_x=>sll_pic_1d_landaudamp_PDF
!                call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, phasespace(3,:) ,particleposition)
!                particleposition=sll_f_pic1d_ensure_periodicity(particleposition, interval_a, interval_b)
!                particleweight_constant=initial_dist_x(particleposition)/sampling_dist_x(particleposition)
!
!                particleweight_constant=initial_dist_xv(particleposition,particlespeed)/sampling_dist_xv(particleposition,particlespeed)
!                particleweight=particleweight_constant                             &
    !                                -sll_f_control_variate_xv(particleposition,particlespeed) /&
    !                                  sampling_dist_xv(particleposition,particlespeed)
!
!
!                particleweight_constant=particleweight_constant/(coll_size*nparticles)
!                particleweight=particleweight/(coll_size*nparticles)
!
!
!
!                !Load disturbance proportional to the actual distortion for landau damping
!                !particleweight=(1.0_f64)*landau_alpha/(real(nparticles*coll_size, i64))
!
!
!                !call sll_pic1d_normalizeweights( particleweight_constant, 1.0_f64 )
!
!!                particleweight=(initial_dist_x(particleposition)*initial_dist_v(particlespeed)&
    !!                                -control_variate_x(particleposition)*control_variate_v(particlespeed))&
    !!                                /(sampling_dist_x(particleposition)*sampling_dist_v(particlespeed))
!
!
!                !call sll_pic1d_normalizeweights( particleweight_constant, 1.0_f64 )
!
!
!                !call sll_pic1d_normalizeweights( particleweight, 1.0_f64 )
!
!                !particleweight=particleweight*abs(cos(landau_mode*particleposition))/cos(landau_mode*particleposition)
!
!                !particleweight=particleweight
!                !normalize particleweight
!                !particleweight=particleweight*0.5_f64+0.5_f64
!
!                !particleweight_constant=interval_length/(1.0_f64 + cos(landau_mode*particleposition)&
    !                    !    /(real(nparticles*coll_size, i64)
!
!                !particleweight=(1.0_f64)*landau_alpha/(real(nparticles*coll_size, i64))
!
!                !        &
    !                    !                *abs(cos(landau_mode*particleposition))/cos(landau_mode*particleposition)
!                !       ! Normalize weights
!                !particleweight=(particleweight/sum(particleweight))/coll_size
!
!
!                !*&(landau_alpha*cos(landau_mode*particleposition))/(particleposition)
!
!                !where (particleposition==0) particleweight=1.0_f64
!                !particleposition=phasespace(3,:)*(interval_b-interval_a)+ interval_a
!
!                ! particleposition=sll_f_pic1d_ensure_periodicity(particleposition, interval_a, interval_b)
!
!
!                print *, "DELTA-F LOADING DONE"
!
!                !project from 0-1 on actual interval
!                steadyparticleposition=interval_a + phasespace(2,:)*(interval_b-interval_a)
!                steadyparticleposition=sll_f_pic1d_ensure_periodicity(steadyparticleposition,  interval_a, interval_b)
!            else
!                steadyparticleposition=interval_a + phasespace(2,:)*(interval_b-interval_a)
!                steadyparticleposition=sll_f_pic1d_ensure_periodicity(steadyparticleposition,  interval_a, interval_b)
!
!                if  (landau_alpha/=0) then
!                    call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, &
    !                        phasespace(3,:) ,particleposition)
!                else
!                    particleposition=phasespace(3,:)*(interval_b-interval_a) - interval_a
!                endif
!                particleweight=(1.0_f64/(nparticles*coll_size*1.0_f64))
!                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
!
!                call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, numberof_streams)
!            endif
!            SLL_DEALLOCATE_ARRAY(phasespace,ierr)
!            !######################################################################
!        case(sll_p_pic1d_testcase_bumpontail)
!                particle_qm=-1.0_f64 !Electrons
!
!            !Generate random numbers
!            SLL_CLEAR_ALLOCATE(phasespace(1:2,1:nparticles),ierr)
!            call sll_s_i8_sobol_generate ( 2_f64, nparticles, coll_rank*nparticles, phasespace )
!
!            if (sll_v_enable_deltaf .eqv. .TRUE.) then
!
!
!                print *, "DELTA-F NOT IMPLEMENTED YET"
!                stop
!            else
!                particlespeed=0
!                call sll_pic1d_load_bumpontail_velocity(phasespace(1,:) , particlespeed, &
    !                    bumpontail_a, bumpontail_v0  ,bumpontail_sigma)
!                if  (landau_alpha/=0) then
!                    call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, &
    !                        phasespace(2,:) ,particleposition)
!                else
!                    particleposition=phasespace(2,:)*(interval_b-interval_a) - interval_a
!                endif
!                steadyparticleposition=interval_a + phasespace(2,:)*(interval_b-interval_a)
!                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
!                steadyparticleposition=sll_f_pic1d_ensure_periodicity(steadyparticleposition,  interval_a, interval_b)
!                particleweight=(1.0_f64/(nparticles*coll_size*1.0_f64))
!            endif
!            SLL_DEALLOCATE_ARRAY(phasespace,ierr)
!
!        !######################################################################
!        case(sll_p_pic1d_testcase_ionbeam)
!            !Generate random numbers
!!            SLL_CLEAR_ALLOCATE(phasespace(1:4,1:nparticles/2),ierr)
!!            call sll_s_i8_sobol_generate ( 6_f64, nparticles/2, coll_rank*(nparticles/2), phasespace )
!
!            SLL_CLEAR_ALLOCATE(phasespace(1:6,1:nparticles/2),ierr)
!            call sll_s_i8_sobol_generate ( 6_f64, nparticles/2, coll_rank*(nparticles/2), phasespace )
!
!            if (sll_v_enable_deltaf .eqv. .TRUE.) then
!                print *, "DELTA-F NOT IMPLEMENTED YET"
!                stop
!            else
!
!                hplus_ratio=1.0_f64      !fixed
!
!                electron_ratio=0.99_f64
!                hminus_ratio=1.0_f64-electron_ratio
!
!                hplus_temp= 0.0148_f64
!                hminus_temp= 0.0148_f64
!                electron_temp=1.0_f64
!
!                particleweight=0
!
!                particleweight=(2.0_f64/(nparticles*coll_size*1.0_f64))
!
!                !Load all particles at beginning and set the Hminus to zero.
!
!                !Load Slow H+
!                particlespeed(nparticles/2+1:nparticles)=sll_normal_rnd(0.0_f64, &
    !                                        sqrt(sll_p_mass/sll_p_proton_mass ),phasespace(3,:))
!   
!                particleposition(nparticles/2+1:nparticles)=interval_a + (phasespace(4,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!                particle_qm(nparticles/2+1:nparticles)=sll_p_mass/sll_p_proton_mass
!
!
!                !Load fast electrons
!                idx_up=floor((nparticles*electron_ratio/2.0_f64))
!                particlespeed(1: idx_up  )=sll_normal_rnd(0.0_f64, 1.0_f64,phasespace(1,1:idx_up))
!                particleposition(1:idx_up)=interval_a + (phasespace(2,1:idx_up)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!                particle_qm(1:idx_up)=-1.0_f64
!
!                !Load slow H-
!                idx_low=idx_up+1
!                idx_up=nparticles/2
!                particlespeed(idx_low:idx_up  )= sll_normal_rnd(0.0_f64, 1.0_f64,phasespace(5,idx_low:idx_up))
!                !particleposition(idx_low:idx_up)=interval_a + (phasespace(6,idx_low:idx_up)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!                particle_qm(idx_low:idx_up)=-sll_p_mass/sll_p_proton_mass
!                particleposition(idx_low:idx_up)=interval_b
!                particleweight(idx_low:idx_up)=0
!
!                call sll_s_pic1d_ensure_boundary_conditions(particleposition, particlespeed)
!
!
!
!
!
!!                !Load fast electrons
!!                particlespeed(1:nparticles/2)=sll_normal_rnd(0.0_f64, 1.0_f64,phasespace(1,:))
!!
!!                particleposition(1:nparticles/2)=interval_a + (phasespace(2,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!!                particle_qm(1:nparticles/2)=-1.0_f64
!!
!!                particleposition(1:nparticles/2)=interval_a + (phasespace(2,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!!
!!                !Load Slow H ions
!!                particlespeed(nparticles/2+1:nparticles)=sll_normal_rnd(0.0_f64, &
    !!                                        sqrt(sll_p_mass/sll_p_proton_mass ),phasespace(3,:))
!!
!!                particleposition(nparticles/2+1:nparticles)=interval_a + (phasespace(4,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!!                particle_qm(nparticles/2+1:nparticles)=sll_p_mass/sll_p_proton_mass
!!
!!
!!                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
!!
!!                !particleposition=steadyparticleposition
!!
!!                particleweight=(2.0_f64/(nparticles*coll_size*1.0_f64))
!            endif
!            !######################################################################
!        case(sll_p_pic1d_testcase_ionbeam_electrons)
!            !Generate random numbers
!            SLL_CLEAR_ALLOCATE(phasespace(1:3,1:nparticles),ierr)
!            call sll_s_i8_sobol_generate ( 3_f64, nparticles, coll_rank*nparticles, phasespace )
!                particle_qm=-1.0_f64
!
!            if (sll_v_enable_deltaf .eqv. .TRUE.) then
!                print *, "DELTA-F NOT IMPLEMENTED YET"
!                stop
!            else
!                call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, numberof_streams)
!                steadyparticleposition=interval_a + (phasespace(2,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!                particleposition=interval_a + (phasespace(3,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
!
!                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
!                steadyparticleposition=sll_f_pic1d_ensure_periodicity(steadyparticleposition,  interval_a, interval_b)
!                !particleposition=steadyparticleposition
!                particleweight=(1.0_f64/(nparticles*coll_size*1.0_f64))
!            endif
!
!            !######################################################################
!        case(sll_p_pic1d_testcase_quiet)
!            SLL_CLEAR_ALLOCATE(phasespace(1:2,1:nparticles),ierr)
!            call sll_s_i8_sobol_generate ( 2_f64, nparticles, coll_rank*nparticles, phasespace )
!                particle_qm=-1.0_f64
!
!            if (sll_v_enable_deltaf .eqv. .TRUE.) then
!                print *, "DELTA-F NOT IMPLEMENTED YET"
!                stop
!            else
!                call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, numberof_streams)
!                particleposition=interval_a + phasespace(2,:)*(interval_b-interval_a)
!                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
!                steadyparticleposition =particleposition
!                particleweight=(1.0_f64/(nparticles*coll_size*1.0_f64))
!
!            endif
!            SLL_DEALLOCATE_ARRAY(phasespace,ierr)
!    end select
!
!
!
!
!
!
!
!
!    !particlespeed=(phasespace(1,:)-0.5_f64)*20.0_f64
!    !particleweight=particleweight*sll_normal_prb_kernel(0.0_f64, 1.0_f64, particlespeed)/(1.0_f64)
!    !particleweight=particleweight/sum(particleweight)
!
!
!    !particleweight=1.0_f64
!    !    particleweight=sll_normal_landaudamp_prb_kernel(mu, sigma, 0.01_f64, sll_p_pi*2*0.1, particlespeed)/&
    !        !                                sll_normal_prb_ kernel(mu, sigma, particlespeed)
!
!    !call random_number(steadyparticleposition);
!    !call random_number(particleposition);
!    !steadyparticleposition=interval_a + steadyparticleposition*(interval_b-interval_a)
!
!    !particleposition=interval_a + particleposition*(interval_b-interval_a)
!
!
!
!    !        !Initzialize absolute velocity as maxwellian
!    !  do i=1,nparticles
!    !      print *, phasespace(1,i)
!    !         call maxwell_cdf_inv( phasespace(1,i) , maxwellian_a, particlespeed(i) )
!    !  enddo
!    !For a quiet start determine velocity direction
!    !  where (phasespace(3,:)>0.5_f64)
!    !      particlespeed=-particlespeed
!    !  end where
!
!
!    !Manipulate generated numbers for maximum negative correlation
!    !phasespace(1 , nparticles/2 +1:nparticles)=1-phasespace(1, 1:nparticles/2)
!    !Initzialize absolute velocity as Gaussian
!    !particlespeed=gaussian_from_rnd( maxwellian_a, maxwellian_a**2, phasespace(1,:))
!
!
!    !phasespace(1,1:nparticles-1:2)= phasespace(1,1:nparticles/2)
!    !phasespace(1,2:nparticles:2)= phasespace(2,1:nparticles/2)
!    !particlespeed=gaussian_from_rnd( 0.0_f64, 20.0_f64**2, phasespace(1,:) )
!
!    !
!    !    do i=2,nparticles
!    !        call sll_s_normal_cdf_inv( phasespace(1,i) , 0.0_f64 , maxwellian_a , particlespeed(i) )
!    !    enddo
!    !call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, 2,  (/ -5.0_f64,5.0_f64 /))
!
!
!
!
!    ! particlespeed=phasespace(1,:) -0.5_f64
!    !particlespeed=birdsall_normal_1d(phasespace(1,:) , 12 )
!    !    particlespeed=particlespeed*1000000_f64
!    !stop
!    !Gaussian speed distribution
!    !do i=1,size(particlespeed)
!    !   particlespeed(i)=gaussianrnd(0, maxwellian_a  )
!    !end do
!
!    !do i=1,nparticles
!    !    call maxwell_pdf(abs(particlespeed(i)), maxwellian_a, particleweight(i))
!    !enddo
!
!    !sll_normal_prb_kernel(maxwellian_a, maxwellian_a**2, particlespeed)/(interval_b -interval_a)
!
!
!    if (size(particlespeed) <= 16) then
!        call sll_display(particlespeed,"(F8.4)")
!        call sll_display(particleposition,"(F8.4)")
!        call sll_display(steadyparticleposition,"(F8.4)")
!        print *, sum(particlespeed )
!    endif
!
!    SLL_ASSERT(minval(particleposition)>=interval_a)
!    SLL_ASSERT(maxval(particleposition)<=interval_b)
!
!end subroutine load_particles


  subroutine sll_s_initialize_intrinsic_mpi_random( collective )
    type(sll_t_collective_t), pointer :: collective

    sll_int32 :: coll_rank
    sll_int32 :: seed_size, idx
    sll_int32, allocatable :: seed(:)

    coll_rank = sll_f_get_collective_rank( collective )

    call random_seed(size=seed_size)  ! output: the size of the seed array
    allocate(seed(seed_size))
    do idx=1, seed_size
      seed(idx) = 10**9*coll_rank
    end do
    call random_seed(put=seed)

    call sll_s_collective_barrier(collective)

  end subroutine sll_s_initialize_intrinsic_mpi_random


  subroutine sll_s_load_particle_species( nparticles, &
      interval_a_user, interval_b_user, particle_species, pic1d_testcase )
    !    steadyparticleposition, &
    !        particleposition, particlespeed, particleweight, particleweight_constant, particle_qm)
    sll_int32 ,                  intent(in)    :: nparticles
    sll_real64,                  intent(in)    :: interval_a_user
    sll_real64,                  intent(in)    :: interval_b_user
    type(sll_t_particle_1d_group), intent(inout) :: particle_species(10)
    sll_int32,                   intent(in)    :: pic1d_testcase

    sll_int32 :: idx
    !sll_real64 :: x !funlandau,x

    !Ionbeam scenario
    sll_real64 :: electron_ratio, hplus_ratio, hminus_ratio
    sll_real64 :: electron_temp, hplus_temp, hminus_temp
    sll_int32 :: idx_up, idx_low

    integer :: ierr
    sll_real64 :: mu, sigma
    !integer :: i
    !real ( kind = 8 ) :: tmp
    !sll_real64 :: maxwellian_a=1.0_f64
    sll_real64, allocatable :: phasespace(:,:)

    !!sll_real64 :: landau_damping=0.01_f64
    !Plasma temperature in Kelvin 150000273.15
    !maxwellian_a=sqrt(sll_p_kb*150000273.15_f64/sll_p_mass)

    interval_a=interval_a_user
    interval_b=interval_b_user
    interval_length=interval_b-interval_a

    !print *, funlandau(particleposition(1))

    mu=0.0_f64
    sigma=1.0_f64

    !!!!!!!!!REMOVE THIS FROM HERE!!!! LATER WE HAVE ONE PIC MODULE
    coll_rank = sll_f_get_collective_rank( sll_v_world_collective )
    coll_size = sll_f_get_collective_size( sll_v_world_collective )
    call sll_s_collective_barrier(sll_v_world_collective)

    sll_v_num_species=0

    selectcase (pic1d_testcase)
        !######################################################################
        case(sll_p_pic1d_testcase_landau)
            sll_v_num_species=1
            SLL_ALLOCATE( particle_species(1)%particle(1:nparticles),ierr)
            SLL_ALLOCATE( particle_species(2)%particle(1:nparticles),ierr)

            !Generate random numbers
            SLL_CLEAR_ALLOCATE(phasespace(1:3,1:nparticles),ierr)
            !call sll_s_i8_sobol_generate ( 3_f64, nparticles, coll_rank*nparticles, phasespace )
            call random_number(phasespace(1:3,1:nparticles))
            particle_species(1)%qm=-1.0_f64 !Electrons
            particle_species(2)%qm=1.0_f64 !Ions

            particle_species(1)%particle%weight_const=0.0_f64
            if (sll_v_enable_deltaf .eqv. .TRUE.) then
                !Landau damping
                !Load deltaf and set weights

                call sll_pic1d_load_stream(phasespace(1,:) , particle_species(1)%particle%vx, 1)

                sampling_dist_v=>sll_pic1d_normalPDF
                initial_dist_x=>sll_pic_1d_landaudamp_PDF
                initial_dist_v=>sll_pic1d_normalPDF
                control_variate_v=>initial_dist_v

                !                control_variate_x=>sll_pic1d_constantPDFx
                !                sampling_dist_x=>sll_pic1d_abscosPDFlandau
                !                particleposition=sll_pic1d_abscosiCDF(interval_length,landau_mode,phasespace(3,:))


                !control_variate_x=>sll_pic1d_constantPDFx
                control_variate_x=>sll_pic_1d_landaudamp_PDF
                sampling_dist_x=>sll_pic_1d_landaudamp_PDF
                call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, phasespace(3,:) , &
                    particle_species(1)%particle%dx)


                particle_species(1)%particle%dx=sll_f_pic1d_ensure_periodicity(particle_species(1)%particle%dx, interval_a, interval_b)
                particle_species(1)%particle%weight=initial_dist_x(particle_species(1)%particle%dx)/sampling_dist_x(particle_species(1)%particle%dx)

                particle_species(1)%particle%weight_const=initial_dist_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx)/sampling_dist_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx)
                particle_species(1)%particle%weight=particle_species(1)%particle%weight_const                             &
                    -sll_f_control_variate_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx) /&
                    sampling_dist_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx)


                particle_species(1)%particle%weight_const= particle_species(1)%particle%weight_const/(coll_size*nparticles)
                particle_species(1)%particle%weight= particle_species(1)%particle%weight/(coll_size*nparticles)


                print *, "DELTA-F LOADING DONE"

                !project from 0-1 on actual interval
                particle_species(2)%particle%dx=interval_a + phasespace(2,:)*(interval_b-interval_a)
                particle_species(2)%particle%dx=sll_f_pic1d_ensure_periodicity(   particle_species(2)%particle%dx,  interval_a, interval_b)
            else
                particle_species(2)%particle%dx=interval_a + phasespace(2,:)*(interval_b-interval_a)
                particle_species(2)%particle%dx=sll_f_pic1d_ensure_periodicity(   particle_species(2)%particle%dx,  interval_a, interval_b)

                if  (landau_alpha/=0) then
                    call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, &
                        phasespace(3,:) ,particle_species(1)%particle%dx)
                else
                    particle_species(1)%particle%dx=phasespace(3,:)*(interval_b-interval_a) - interval_a
                endif
                particle_species(1)%particle%weight=(1.0_f64/(nparticles*coll_size*1.0_f64))
                particle_species(1)%particle%dx=sll_f_pic1d_ensure_periodicity(particle_species(1)%particle%dx,  interval_a, interval_b)

                call sll_pic1d_load_stream(phasespace(1,:) , particle_species(1)%particle%vx, numberof_streams)
            endif
            SLL_DEALLOCATE_ARRAY(phasespace,ierr)
            !######################################################################
        case(sll_p_pic1d_testcase_bumpontail)
            particle_species(1)%qm=-1.0_f64 !Electrons
            !Generate random numbers
            sll_v_num_species=1
            SLL_ALLOCATE( particle_species(1)%particle(1:nparticles),ierr)

            !Generate random numbers
            SLL_CLEAR_ALLOCATE(phasespace(1:2,1:nparticles),ierr)
            call sll_s_i8_sobol_generate ( 2_f64, nparticles, coll_rank*nparticles, phasespace )


            particle_species(1)%particle%weight_const=0.0_f64
            if (sll_v_enable_deltaf .eqv. .TRUE.) then

                call sll_pic1d_load_bumpontail_velocity(phasespace(1,:)  , particle_species(1)%particle%vx, &
                    bumpontail_a, bumpontail_v0  ,bumpontail_sigma)

                sampling_dist_v=>sll_pic_1d_bumpontail_PDF
                sampling_dist_x=>sll_pic_1d_landaudamp_PDF

                initial_dist_x=>sampling_dist_x
                initial_dist_v=>sampling_dist_v
                control_variate_v=>initial_dist_v
                control_variate_x=>initial_dist_x


                if  (landau_alpha/=0) then
                    call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, &
                        phasespace(2,:) ,particle_species(1)%particle%dx)
                else
                    particle_species(1)%particle%dx=phasespace(2,:)*(interval_b-interval_a) - interval_a
                endif
                particle_species(1)%particle%dx=sll_f_pic1d_ensure_periodicity(particle_species(1)%particle%dx,  interval_a, interval_b)

                particle_species(1)%particle%weight=initial_dist_x(particle_species(1)%particle%dx)/sampling_dist_x(particle_species(1)%particle%dx)

                particle_species(1)%particle%weight_const=initial_dist_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx)/sampling_dist_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx)
                particle_species(1)%particle%weight=particle_species(1)%particle%weight_const                             &
                    -sll_f_control_variate_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx) /&
                    sampling_dist_xv(particle_species(1)%particle%dx,   particle_species(1)%particle%vx)


                particle_species(1)%particle%weight_const= particle_species(1)%particle%weight_const/(coll_size*nparticles)
                particle_species(1)%particle%weight= particle_species(1)%particle%weight/(coll_size*nparticles)
                print *, "DELTA-F LOADING DONE"
            else

                call sll_pic1d_load_bumpontail_velocity(phasespace(1,:)  , particle_species(1)%particle%vx, &
                    bumpontail_a, bumpontail_v0  ,bumpontail_sigma)

                if  (landau_alpha/=0) then
                    call sll_pic1d_load_landau(landau_alpha,landau_mode, interval_a,interval_b, &
                        phasespace(2,:) ,particle_species(1)%particle%dx)
                else
                    particle_species(1)%particle%dx=phasespace(2,:)*(interval_b-interval_a) - interval_a
                endif
                particle_species(1)%particle%dx=sll_f_pic1d_ensure_periodicity(particle_species(1)%particle%dx,  interval_a, interval_b)


                particle_species(1)%particle%weight=(1.0_f64/(nparticles*coll_size*1.0_f64))
            endif
            SLL_DEALLOCATE_ARRAY(phasespace,ierr)

            !######################################################################
        case(sll_p_pic1d_testcase_ionbeam)
            !Generate random numbers
            !            SLL_CLEAR_ALLOCATE(phasespace(1:4,1:nparticles/2),ierr)
            !            call sll_s_i8_sobol_generate ( 6_f64, nparticles/2, coll_rank*(nparticles/2), phasespace )

            SLL_CLEAR_ALLOCATE(phasespace(1:6,1:nparticles/2),ierr)
            call sll_s_i8_sobol_generate ( 6_f64, nparticles/2, coll_rank*(nparticles/2), phasespace )
            sll_v_num_species=3

            if (sll_v_enable_deltaf .eqv. .TRUE.) then
                print *, "DELTA-F NOT IMPLEMENTED YET"
                stop
            else

                hplus_ratio=1.0_f64      !fixed

                electron_ratio=0.99_f64
                hminus_ratio=1.0_f64-electron_ratio

                hplus_temp= 0.0148_f64
                hminus_temp= 0.0148_f64
                electron_temp=1.0_f64


                !Load all particles at beginning and set the Hminus to zero.

                !Load Slow H+

                SLL_ALLOCATE( particle_species(1)%particle(1:nparticles/2),ierr)

                particle_species(1)%particle%vx=sll_normal_rnd(0.0_f64, &
                    hplus_temp,phasespace(3,:))

                particle_species(1)%particle%dx=interval_a + (phasespace(4,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
                particle_species(1)%qm=sll_p_mass/sll_p_proton_mass
                particle_species(1)%particle%weight=(2.0_f64/(nparticles*coll_size*1.0_f64))


                !Load fast electrons
                idx_up=floor(nparticles*electron_ratio/2.0_f64)
                SLL_ALLOCATE( particle_species(2)%particle(1:idx_up),ierr)

                particle_species(2)%particle%vx=sll_normal_rnd(0.0_f64, 1.0_f64,phasespace(1,1:idx_up))
                particle_species(2)%particle%dx=interval_a + (phasespace(2,1:idx_up)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
                particle_species(2)%qm=-1.0_f64
                particle_species(2)%particle%weight=(2.0_f64/(nparticles*coll_size*1.0_f64))

                !Load slow H-
                idx_low=idx_up+1
                idx_up=nparticles/2
                SLL_ALLOCATE( particle_species(3)%particle(1:idx_up-idx_low+1),ierr)

                particle_species(3)%particle%vx= sll_normal_rnd(0.0_f64, hminus_temp,phasespace(5,idx_low:idx_up))
                !particleposition(idx_low:idx_up)=interval_a + (phasespace(6,idx_low:idx_up)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
                particle_species(3)%qm=-sll_p_mass/sll_p_proton_mass
                particle_species(3)%particle%dx=interval_b
                particle_species(3)%particle%weight=0.0_f64

                do idx=1,sll_v_num_species
                    call sll_s_pic1d_ensure_boundary_conditions(  particle_species(idx)%particle%dx, particle_species(idx)%particle%vx,pic1d_testcase)
                enddo

            endif

            !######################################################################
            !        case(sll_p_pic1d_testcase_ionbeam_electrons)
            !            !Generate random numbers
            !            SLL_CLEAR_ALLOCATE(phasespace(1:3,1:nparticles),ierr)
            !            call sll_s_i8_sobol_generate ( 3_f64, nparticles, coll_rank*nparticles, phasespace )
            !                particle_qm=-1.0_f64
            !
            !            if (sll_v_enable_deltaf .eqv. .TRUE.) then
            !                print *, "DELTA-F NOT IMPLEMENTED YET"
            !                stop
            !            else
            !                call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, numberof_streams)
            !                steadyparticleposition=interval_a + (phasespace(2,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
            !                particleposition=interval_a + (phasespace(3,:)+1.0_f64)*(interval_b-interval_a)*plasma_size/2.0_f64
            !
            !                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
            !                steadyparticleposition=sll_f_pic1d_ensure_periodicity(steadyparticleposition,  interval_a, interval_b)
            !                !particleposition=steadyparticleposition
            !                particleweight=(1.0_f64/(nparticles*coll_size*1.0_f64))
            !            endif
            !
            !            !######################################################################
            !        case(sll_p_pic1d_testcase_quiet)
            !            SLL_CLEAR_ALLOCATE(phasespace(1:2,1:nparticles),ierr)
            !            call sll_s_i8_sobol_generate ( 2_f64, nparticles, coll_rank*nparticles, phasespace )
            !                particle_qm=-1.0_f64
            !
            !            if (sll_v_enable_deltaf .eqv. .TRUE.) then
            !                print *, "DELTA-F NOT IMPLEMENTED YET"
            !                stop
            !            else
            !                call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, numberof_streams)
            !                particleposition=interval_a + phasespace(2,:)*(interval_b-interval_a)
            !                particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
            !                steadyparticleposition =particleposition
            !                particleweight=(1.0_f64/(nparticles*coll_size*1.0_f64))
            !
            !            endif
            !            SLL_DEALLOCATE_ARRAY(phasespace,ierr)
    end select

    !particlespeed=(phasespace(1,:)-0.5_f64)*20.0_f64
    !particleweight=particleweight*sll_normal_prb_kernel(0.0_f64, 1.0_f64, particlespeed)/(1.0_f64)
    !particleweight=particleweight/sum(particleweight)

    !particleweight=1.0_f64
    !    particleweight=sll_normal_landaudamp_prb_kernel(mu, sigma, 0.01_f64, sll_p_pi*2*0.1, particlespeed)/&
        !                                sll_normal_prb_ kernel(mu, sigma, particlespeed)

    !call random_number(steadyparticleposition);
    !call random_number(particleposition);
    !steadyparticleposition=interval_a + steadyparticleposition*(interval_b-interval_a)

    !particleposition=interval_a + particleposition*(interval_b-interval_a)

    !        !Initzialize absolute velocity as maxwellian
    !  do i=1,nparticles
    !      print *, phasespace(1,i)
    !         call maxwell_cdf_inv( phasespace(1,i) , maxwellian_a, particlespeed(i) )
    !  enddo
    !For a quiet start determine velocity direction
    !  where (phasespace(3,:)>0.5_f64)
    !      particlespeed=-particlespeed
    !  end where

    !Manipulate generated numbers for maximum negative correlation
    !phasespace(1 , nparticles/2 +1:nparticles)=1-phasespace(1, 1:nparticles/2)
    !Initzialize absolute velocity as Gaussian
    !particlespeed=gaussian_from_rnd( maxwellian_a, maxwellian_a**2, phasespace(1,:))

    !phasespace(1,1:nparticles-1:2)= phasespace(1,1:nparticles/2)
    !phasespace(1,2:nparticles:2)= phasespace(2,1:nparticles/2)
    !particlespeed=gaussian_from_rnd( 0.0_f64, 20.0_f64**2, phasespace(1,:) )

    !
    !    do i=2,nparticles
    !        call sll_s_normal_cdf_inv( phasespace(1,i) , 0.0_f64 , maxwellian_a , particlespeed(i) )
    !    enddo
    !call sll_pic1d_load_stream(phasespace(1,:) , particlespeed, 2,  (/ -5.0_f64,5.0_f64 /))

    ! particlespeed=phasespace(1,:) -0.5_f64
    !particlespeed=birdsall_normal_1d(phasespace(1,:) , 12 )
    !    particlespeed=particlespeed*1000000_f64
    !stop
    !Gaussian speed distribution
    !do i=1,size(particlespeed)
    !   particlespeed(i)=gaussianrnd(0, maxwellian_a  )
    !end do

    !do i=1,nparticles
    !    call maxwell_pdf(abs(particlespeed(i)), maxwellian_a, particleweight(i))
    !enddo

    !sll_normal_prb_kernel(maxwellian_a, maxwellian_a**2, particlespeed)/(interval_b -interval_a)

    !
    !    if (size(particlespeed) <= 16) then
    !        call sll_display(particlespeed,"(F8.4)")
    !        call sll_display(particleposition,"(F8.4)")
    !        call sll_display(steadyparticleposition,"(F8.4)")
    !        print *, sum(particlespeed )
    !    endif
    !
    !    SLL_ASSERT(minval(particleposition)>=interval_a)
    !    SLL_ASSERT(maxval(particleposition)<=interval_b)

  end subroutine sll_s_load_particle_species


  subroutine sll_pic1d_load_landau( landau_alpha, landau_mode, &
      interval_a, interval_b, uniform_random, particleposition )
    sll_real64, intent(in)  :: landau_alpha
    sll_real64, intent(in)  :: landau_mode
    sll_real64, intent(in)  :: interval_a
    sll_real64, intent(in)  :: interval_b
    sll_real64, intent(in)  :: uniform_random(:)
    sll_real64, intent(out) :: particleposition(:)

    sll_int32  :: nparticles
    sll_real64 :: interval_length
    !sll_int32 :: i

    nparticles=size(uniform_random)
    SLL_ASSERT(size(particleposition)==nparticles)

    interval_length=interval_b-interval_a

    !        particleposition=uniform_random*interval_length+interval_a

    !Introduce Landau damping disturbance
    !
    !            do i=1, floor(nparticles*(1.0_f64-landau_alpha))
    !                    particleposition(i)=uniform_random(i)
    !
    !                    !particleposition(i)= phasespace(3,i)
    !            enddo

    !
    !            do i=floor(nparticles*(1.0_f64-landau_alpha)) +1, nparticles
    !                SLL_ASSERT(   )
    !                particleposition(i)= asin(interval_length*uniform_random(i)*landau_mode)/(landau_mode )&
        !                                        +interval_a
    !!                particleposition(i)= sin(uniform_random(i)*landau_mode)/(landau_mode*interval_length )&
        !!                                        +interval_a
    !            enddo

    !        particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)

    !Direct inverse sampling
    particleposition=sll_pic_1d_landaudamp_CDF(landau_alpha,landau_mode, &
        interval_length ,uniform_random)

    !particleposition=uniform_random+landau_alpha*sin(landau_mode*uniform_random);
    !particleposition=particleposition*(interval_length)

    particleposition=interval_a + particleposition
    particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)

  end subroutine sll_pic1d_load_landau


  subroutine sll_s_pic1d_ensure_boundary_conditions_species( particle_species, &
      pic1d_testcase )
    type(sll_t_particle_1d_group), intent(inout) :: particle_species(:)
    sll_int32                  , intent(in)    :: pic1d_testcase

    sll_int32 :: sll_v_num_species, jdx

    sll_v_num_species=size(particle_species)
    SLL_ASSERT(sll_v_num_species==size(particle_species))

    do jdx=1,sll_v_num_species
        call sll_s_pic1d_ensure_boundary_conditions(particle_species(jdx)%particle%dx,&
            particle_species(jdx)%particle%vx,pic1d_testcase )
    end do
  end subroutine sll_s_pic1d_ensure_boundary_conditions_species


  subroutine sll_s_pic1d_ensure_boundary_conditions( particle_position, &
      particlespeed, pic1d_testcase )
    sll_real64, intent(inout) :: particle_position(:)
    sll_real64, intent(inout) :: particlespeed(:)
    sll_int32 , intent(in)    :: pic1d_testcase

    selectcase (pic1d_testcase)
        case(sll_p_pic1d_testcase_ionbeam)
            !Deflect particles
            where (particle_position>interval_b)
                !particle_position=particle_position -2.0_f64*(particle_position- interval_b)
                particle_position=interval_b
                particlespeed=0.0_f64
            end where
            where (particle_position<=interval_a)
                particle_position=particle_position -2.0_f64*(particle_position- interval_a)
                particle_position=interval_a
                particlespeed=-particlespeed
            end where
        case default
            particle_position=sll_f_pic1d_ensure_periodicity( particle_position, &
                interval_a, interval_b)
    end select
  end subroutine sll_s_pic1d_ensure_boundary_conditions


  function sll_f_pic1d_ensure_periodicity( particle_position, &
        interval_a, interval_b) result( particle_position_out)
    sll_real64, intent(in) :: particle_position(:)
    sll_real64, intent(in) :: interval_a
    sll_real64, intent(in) :: interval_b

    sll_real64 :: particle_position_out(size(particle_position))
    sll_real64 :: interval_length
    !sll_int32 :: idx

    SLL_ASSERT(interval_a < interval_b)
    interval_length=interval_b - interval_a

    !        do idx=1, size(particle_position)
    !
    !            do while (particle_position_out(idx)<interval_a)
    !                    particle_position_out(idx)=particle_position_out(idx) +interval_length
    !            enddo
    !            do while (particle_position_out(idx)>=interval_b)
    !                    particle_position_out(idx)=particle_position_out(idx) -interval_length
    !            enddo
    !        enddo
    particle_position_out=particle_position- interval_a
    !        if (sll_pic1d_testcase==sll_p_pic1d_testcase_ionbeam)
    !            where (particle_position>interval_b) particle_position_out=interval_b*0.99999_f64

    do while (minval(particle_position_out)<0.0_f64)
        particle_position_out=particle_position_out +interval_length
    end do

    particle_position_out=mod(particle_position_out, interval_length)&
        +interval_a
    !    do i=1, nparticles
    !            particleposition(i)= (particleposition(i)-interval_a )*(1.0_f64+ 0.01_f64*cos((sll_p_kx/(interval_b-interval_a))*1.0_f64*(particleposition(i)-interval_a ))) + interval_a
    !    enddo
    !    particleposition=sll_f_pic1d_ensure_periodicity(particleposition,  interval_a, interval_b)
    !    if (minval(particle_position_out)==NaN) then
    !        print *, "NaN appeared"
    !    endif

    SLL_ASSERT(minval(particle_position_out)>=interval_a)
    SLL_ASSERT(maxval(particle_position_out)<interval_b)
    !            print *, maxval(particle_position_out)

  end function sll_f_pic1d_ensure_periodicity


  !The Original is 1+alpha*sin(2pi*x)
  function sll_pic_1d_landaudamp_CDF( alpha, mode, interval_length, uniform_random ) &
        result(inverseCDF)
    sll_real64, intent(in) :: alpha
    sll_real64, intent(in) :: mode
    sll_real64, intent(in) :: interval_length
    sll_real64, intent(in) :: uniform_random(:)

    sll_real64, dimension(size(uniform_random)):: inverseCDF, inverseCDF_new
    !sll_real  :: newtonstep
    sll_real64 :: numerror=1.0_f64
    sll_int32  :: idx

    !        if (uniform_random==0.0_f64 .OR.uniform_random==1.0_f64 ) then
    !             inverseCDF=uniform_random

    SLL_ASSERT(interval_length>0)

    !Newton method with fixpoint steps
    inverseCDF=uniform_random

    where  (uniform_random==0.0_f64 .OR.uniform_random==1.0_f64 )
        inverseCDF=interval_length*uniform_random
    end where

    idx=0
    do while (numerror>=1D-28 .AND. idx<=20 )
        idx=idx+1
        !        inverseCDF=inverseCDF - &
            !                        (uniform_random + (alpha/((sll_p_kx)*mode))*cos(sll_p_kx*mode*inverseCDF))/&
            !                        (alpha*sin(sll_p_kx*mode*inverseCDF))
        inverseCDF_new=(interval_length*uniform_random - (alpha/(mode))*sin(mode*inverseCDF))
        numerror=maxval((inverseCDF_new-inverseCDF)**2)/(maxval(inverseCDF_new)/maxval(inverseCDF))
        !print *, idx, "#Fixpoint Iter Error:" ,numerror
        inverseCDF=inverseCDF_new
    end do

  end function sll_pic_1d_landaudamp_CDF


  function sll_pic_1d_landaudamp_PDF1(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    y = (1.0_f64 + cos(landau_mode*x))/(interval_length)
  end function sll_pic_1d_landaudamp_PDF1


  function sll_pic_1d_landaudamp_PDF(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    y = (1.0_f64 + landau_alpha*cos(landau_mode*x))/(interval_length)
  end function sll_pic_1d_landaudamp_PDF


  function sll_pic_1d_bumpontail_PDF(v) result(y)
    sll_real64, intent(in) :: v(:)
    sll_real64 :: y(size(v))
    y = 1.0_f64/(1.0_f64+bumpontail_a)*( &
        exp(-0.5_f64*v**2)/sqrt(sll_p_kx) + &
        bumpontail_a*exp(-0.5_f64*((v-bumpontail_v0)/bumpontail_sigma)**2)&
        /sqrt(sll_p_kx)/bumpontail_sigma)
  end function sll_pic_1d_bumpontail_PDF


  function sll_pic_1d_landaudamp_PDFxv(x,v) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64, intent(in) :: v(:)
    sll_real64 :: y(size(x))
    y = (1.0_f64 + landau_alpha*cos(landau_mode*x))/(interval_length) * &
      exp(-0.5_f64*v**2)/sqrt(sll_p_kx)
  end function sll_pic_1d_landaudamp_PDFxv


  function sll_pic1d_normalPDF(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    y = 1.0_f64/sqrt(sll_p_kx)*exp(-0.5_f64*x**2)
  end function sll_pic1d_normalPDF


  function sll_pic_1d_cos_landau(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    y = (landau_alpha*cos(landau_mode*x))/(interval_length)
  end function sll_pic_1d_cos_landau


  function sll_pic1d_abscosPDF(L,k,x) result(y)
    sll_real64, intent(in) :: L
    sll_real64, intent(in) :: k
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    SLL_ASSERT(L>0)
    y = abs(cos(k*x))/(L*4.0_f64/sll_p_kx)
  end function sll_pic1d_abscosPDF


  function sll_pic1d_abscosPDFlandau(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    sll_real64 :: L,k
    L=interval_length
    k=landau_mode
    SLL_ASSERT(L>0)
    y = abs(cos(k*x))/(L*4.0_f64/sll_p_kx)
  end function sll_pic1d_abscosPDFlandau


  function sll_pic1d_constantPDFx(x) result(y)
    sll_real64, intent(in) :: x(:)
    sll_real64 :: y(size(x))
    y = 1.0_f64/(interval_length)
  end function sll_pic1d_constantPDFx


  function sll_pic1d_abscosiCDF(L,k,x) result(y)
    sll_real64, intent(in) :: L
    sll_real64, intent(in) :: k
    sll_real64, intent(in) :: x(:)

    sll_real64 :: y(size(x))

    SLL_ASSERT(L>0)

    where ( mod(x*k*(L/sll_p_kx),0.5_f64)<0.25_f64)  y=asin(  mod((L*4.0_f64/sll_p_kx)*x*k,1.0_f64) )/k
    where ( mod(x*k*(L/sll_p_kx),0.5_f64)>0.25_f64) &
        y=acos((0.5_f64-mod(2.0_f64*x*(L/sll_p_kx),1.0_f64/k)/2.0_f64*k )*4.0_f64)/k  + sll_p_pi/2.0_f64/k

    y = y+real(floor( x*k*(L/sll_p_kx)/0.5_f64),f64)*sll_p_pi/k

    !Matlab Code
    !f=@(x)
    !
    !g=@(x)  (sin(mod(k*x,pi/2))/k).*(mod(k*x,pi)<pi/2)...
    !+(1-cos(mod(k*x,pi/2)))/k.*(mod(k*x,pi)>=pi/2  ) ...
    !   + floor(k*x/(pi/2))/k
    !
    !ginv=@(x)asin(  mod((L*4/(2*pi))*x*k,1) )/k.*( mod(x*k*(L/2/pi),0.5)<0.25)...
    !            +  floor( x*k*(L/2/pi)/0.5)*pi/k...
    !+   (acos(   (0.5-mod(2*x*(L/2/pi),1/k)/2*k )*4)/k  + pi/2/k).*( mod(x*k*(L/2/pi),0.5)>0.25)...
    !
  end function sll_pic1d_abscosiCDF


  !< Normalizes the vector weights to norm, in a collective enviroment
  subroutine sll_pic1d_normalizeweights( weights, norm )
    sll_real64, intent(inout) :: weights(:)
    sll_real64, intent(in)    ::   norm

    sll_real64 :: sumweights

    if (norm/=0.0_f64) then
        sumweights=sum(weights)
        weights=(weights/(sumweights))*norm/coll_size
    else
        weights=(weights/(sumweights)-1.0_f64)/coll_size
    end if
  end subroutine sll_pic1d_normalizeweights


  !<Plasma Frequncy, density
  function sll_plasma_frequency(density) result(omega)
    sll_real64, intent(in) :: density
    sll_real64 :: omega

    omega = sqrt(density/(sll_p_mass*sll_p_epsilon_0) )*sll_p_charge
  end function sll_plasma_frequency


  !<Thermal velocity, temperature
  function sll_thermal_velocity(T) result(vth)
    sll_real64, intent(in) :: T
    sll_real64 :: vth
    vth = sqrt(sll_p_kb*T/sll_p_mass)
  end function sll_thermal_velocity


  function initial_dist_xv(x,v) result(p)
    sll_real64, intent(in) :: x(:)
    sll_real64, intent(in) :: v(:)
    sll_real64 :: p(size(x))
    p = initial_dist_x(x)*initial_dist_v(v)
  end function initial_dist_xv


  function sampling_dist_xv(x,v) result(p)
    sll_real64, intent(in) :: x(:)
    sll_real64, intent(in) :: v(:)
    sll_real64 :: p(size(x))
    p = sampling_dist_x(x)*sampling_dist_v(v)
  end function sampling_dist_xv


  function sll_f_control_variate_xv(x,v) result(p)
    sll_real64, intent(in) :: x(:)
    sll_real64, intent(in) :: v(:)
    sll_real64 :: p(size(x))
    p = control_variate_x(x)*control_variate_v(v)
  end function sll_f_control_variate_xv


  elemental sll_real64 function sll_local_maxwellian(v, vtherm, temperature)
    sll_real64, intent(in) :: v
    sll_real64, intent(in) :: vtherm
    sll_real64, intent(in) :: temperature

    sll_local_maxwellian = exp(-0.5_f64*(v-vtherm)**2/temperature)

  end function sll_local_maxwellian


  !>Note that this cannot be done in parallel
  function sll_rejection_sampling( probability_density, reject_const, nmark ) &
                    result(sample)
    procedure(sll_pdf_1d)  :: probability_density
    sll_real64, intent(in) :: reject_const
    sll_int32 , intent(in) :: nmark
    sll_real64 :: sample(nmark)

    sll_int32  :: idx!,jdx
    sll_real64 :: rndu
    sll_real64, dimension(1) :: rndg, f_rndg
    idx=1

    do while (idx<=nmark)
        call random_number(rndg)
        call random_number(rndu)
        f_rndg=probability_density(rndg)
        if (rndu <= f_rndg(1)/reject_const) then
            sample(idx)=rndg(1)
            idx=idx+1
        endif
    enddo

  end function sll_rejection_sampling


end module
