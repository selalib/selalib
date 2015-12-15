module sll_m_pic_1d_distribution

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_pic_visu, only: &
    sll_s_compute_df_cic

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_m_xdmf, only: &
    sll_s_xdmf_corect2d_nodes

  implicit none

  public :: &
    sll_t_pic1d_eulerian_distribution

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!==============================================================================

  type :: sll_t_pic1d_eulerian_distribution

    sll_real64, allocatable :: f(:,:)
    sll_real64, allocatable :: n(:), nu(:), nke(:)
!    sll_float64, allocatable :: x(:), v(:)
    sll_int32                :: nx, nv
    sll_real64               :: xmin, xmax
    sll_real64               :: vmin, vmax

  contains

    procedure :: initialize       => pic1d_ed__initialize
    procedure :: compute_f        => pic1d_ed__compute_f
    procedure :: compute_moments  => pic1d_ed__compute_moments
    procedure :: print_f          => pic1d_ed__print_f
    procedure :: print_moments    => pic1d_ed__print_moments
    
    procedure :: mean_velocity    => pic1d_ed__mean_velocity
    procedure :: mean_temperature => pic1d_ed__mean_temperature
    
  end type

!==============================================================================
contains
!==============================================================================

  subroutine pic1d_ed__initialize( self, xmin, xmax, nx, vmin, vmax, nv )
    class( sll_t_pic1d_eulerian_distribution ), intent( inout ) :: self
    sll_real64                          , intent( in    ) :: xmin
    sll_real64                          , intent( in    ) :: xmax
    sll_int32                           , intent( in    ) :: nx   
    sll_real64                          , intent( in    ) :: vmin
    sll_real64                          , intent( in    ) :: vmax
    sll_int32                           , intent( in    ) :: nv   

    sll_int32 :: ierr

    self%xmin = xmin
    self%xmax = xmax
    self%nx   = nx
    self%vmin = vmin
    self%vmax = vmax
    self%nv   = nv

    SLL_ALLOCATE( self%f  (nx,nv), ierr )
    SLL_ALLOCATE( self%n  (nx)   , ierr )
    SLL_ALLOCATE( self%nu (nx)   , ierr )
    SLL_ALLOCATE( self%nke(nx)   , ierr )

  end subroutine pic1d_ed__initialize

  !----------------------------------------------------------------------------
  subroutine pic1d_ed__compute_f( self, xp, vp, wp )
    class( sll_t_pic1d_eulerian_distribution ), intent( inout ) :: self
    sll_real64                          , intent( in    ) :: xp(:)
    sll_real64                          , intent( in    ) :: vp(:)
    sll_real64                          , intent( in    ) :: wp(:)
    
    ! Use cloud-in-cell (CIC) algorithm to compute 1D-1V distribution function
    call sll_s_compute_df_cic( xp, vp, wp, &
      self%xmin, self%xmax, self%nx, &
      self%vmin, self%vmax, self%nv, &
      self%f )

  end subroutine pic1d_ed__compute_f

  !----------------------------------------------------------------------------
  ! Computing the number of particle in each cell
  subroutine pic1d_ed__compute_moments( self )
    class( sll_t_pic1d_eulerian_distribution ), intent( inout ) :: self

    sll_int32  :: i
    sll_real64 :: dv
    sll_real64 :: v (self%nv)
    sll_real64 :: ke(self%nv)  ! v^2/2
    
    dv = (self%vmax-self%vmin)/(self%nv-1)
    do i=1,self%nv
      v(i) = self%vmin + dv*(i-1)
    end do
    ke = 0.5_f64*v**2
    
    do i=1,self%nx
      self%n  (i) = sum( self%f(i,:) ) * dv
      self%nu (i) = dot_product( v , self%f(i,:) ) * dv
      self%nke(i) = dot_product( ke, self%f(i,:) ) * dv
    end do

  end subroutine pic1d_ed__compute_moments
  
  !----------------------------------------------------------------------------
  ! Compute mean velocity in each cell
  function pic1d_ed__mean_velocity( self ) result( u )
    class( sll_t_pic1d_eulerian_distribution ), intent( in ) :: self
    sll_real64                                         :: u(self%nx)

    u = self%nu / self%n

  end function pic1d_ed__mean_velocity

  !----------------------------------------------------------------------------
  ! Compute mean temperature in each cell
  function pic1d_ed__mean_temperature( self ) result( T )
    class( sll_t_pic1d_eulerian_distribution ), intent( in ) :: self
    sll_real64                                         :: T(self%nx)

    
    T = (2.0_f64*self%nke - self%nu**2/self%n) / self%n
    !K*T/2=(kin_en-n*mean_velocity^2/2)/n
    !n*mean_velocity^2/2=nu^2/(2*n)
    !We fix K=1
    
  end function pic1d_ed__mean_temperature
  
  !----------------------------------------------------------------------------
  subroutine pic1d_ed__print_f( self, plot_name, iplot )
    class( sll_t_pic1d_eulerian_distribution ), intent( in ) :: self
    character( len=* )                  , intent( in ) :: plot_name
    sll_int32                           , intent( in ) :: iplot

    sll_real64       :: dx, dv
    character(len=4) :: fin

    dx = (self%xmax-self%xmin)/(self%nx-1)
    dv = (self%vmax-self%vmin)/(self%nv-1)
    call sll_s_int2string( iplot, fin )
    call sll_s_xdmf_corect2d_nodes( plot_name//'_'//fin, self%f, "f(x,v)", &
      self%xmin, dx, self%vmin, dv )

  end subroutine pic1d_ed__print_f
  
  !----------------------------------------------------------------------------
  subroutine pic1d_ed__print_moments( self, plot_name, iplot,root_path )
    class( sll_t_pic1d_eulerian_distribution ), intent( in ) :: self
    character( len=* )                  , intent( in ) :: plot_name
    sll_int32                           , intent( in ) :: iplot
    character(len=256)                  , intent( in ) :: root_path
    character(len=4)     :: fin
    integer , parameter  :: file_id=20 
    sll_int32            :: i   
    sll_real64           :: dx,x(self%nx)
    sll_real64                                         :: T(self%nx)
        
    dx = (self%xmax-self%xmin)/(self%nx-1)
    do i=1,self%nx
      x(i) = self%xmin + dx*(i-1)
    end do
    T=self%mean_temperature()                       
    call sll_s_int2string( iplot, fin )
    open(file_id, file = trim(root_path)//plot_name//'_'//fin//'.dat')
    write (file_id,*)"#distribution moments"
    write (file_id,*)  "#Time steps:", iplot
    write(file_id,*) "node_index                  " ,"x                    ","n                      ","u                   ","nke               ","Temperature            "
  
    do i=1,self%nx
      write(file_id,*)  i,x(i),self%n(i),self%nu(i),self%nke(i),T(i)
    end do
  close(file_id)
  
  
  end subroutine pic1d_ed__print_moments


end module sll_m_pic_1d_distribution

