module sll_m_particle_group_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_particle_group_base, &
    sll_t_species, &
    sll_f_index_1dto2d, &
    sll_f_index_1dto3d, &
    sll_f_index_2dto1d, &
    sll_f_index_3dto1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !============================================================================
  ! Particle species
  !============================================================================
  type :: sll_t_species

    character(len=64) :: name !< species name
    sll_real64        :: q    !< charge of a single particle
    sll_real64        :: m    !< mass   of a single particle

  contains
    procedure         :: q_over_m  !< charge over mass ratio
    procedure         :: init => initialize_species !< Constructor

  end type sll_t_species


  !============================================================================
  ! Particle group
  !============================================================================
  type, abstract :: sll_c_particle_group_base

    class( sll_t_species ), allocatable :: species
    sll_int32                           :: id
    sll_int32                           :: n_particles !< number of particles local to the processor
    sll_int32                           :: n_total_particles !< number of particles in total simulation    
    sll_int32                           :: n_weights !< number of weights per particle


  contains
    ! Getters ( for single particles )
    procedure( i_get_coords ), deferred :: get_x !< Get the physical coordinate of one particle   
    procedure( i_get_coords ), deferred :: get_v !< Get the velocity of one particle
    procedure( i_get_scalar ), deferred :: get_charge !< Get the charge of one particle
    procedure( i_get_scalar ), deferred :: get_mass !< Get the mass of one particle
    procedure( i_get_array  ), deferred :: get_weights !< Get the weights of one particle
    procedure                           :: get_box !< Assuming we have some elements/grid cells "box" for the representation of the fields, we get the (1d) index of this box
    procedure                           :: get_boxnd !< Assuming we have some elements/grid cells "box" for the representation of the fields, we get the nd index of this box (note this is for tensor product grids only)
    procedure                           :: get_xbox  !< Assuming we have some elements/grid cells "box" for the representation of the fields, we get the (normalized) offset of the position in the box.
    procedure                           :: get_patch
    procedure                           :: get_xpatch !< Get x in terms of the patch local coordinate system

    ! Getters ( for chunk of particles )
    ! All chunk-getters are calling the single versions by default
    ! We use pointers to avoid copying. Note that no setters are necessary since we use pointers to access the values.
    ! For chunks, only one chget_x is provided this gives the coordinates in the format that the positions are stored. 
    procedure                           :: chget_x
    procedure                           :: chget_v
    procedure                           :: chget_weights
    procedure                           :: chget_box
    procedure                           :: chget_boxnd
    procedure                           :: chget_patch
    
    ! Setters ( for single particles )
    procedure( i_set_coords ), deferred :: set_x !< Set the coordinate of one particle
    procedure( i_set_coords ), deferred :: set_v !< Set the velocity of one particle
    procedure( i_set_array  ), deferred :: set_weights !< Set the particle weight
    procedure( set_scalar ),   deferred :: set_common_weight !< Set the common particle weight

    procedure                           :: set_box
    procedure                           :: set_boxnd
    procedure                           :: set_patch

    procedure( empty ),        deferred :: free

!    ! Getters for the whole group
!    procedure( get_all_coords), deferred :: get_all_x
    !    procedure( get_all_coords), deferred :: get_all_v

    procedure                           :: print !< Prints out the particle array in some form (for debugging, by default nothing is printed out)

  end type sll_c_particle_group_base

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_int( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       ,   intent( in ) :: i
    sll_int32 :: r
  end function i_get_int
end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_intnd( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       ,   intent( in ) :: i
    sll_int32 :: r(3)
  end function i_get_intnd
end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_scalar( self, i , i_weight) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       ,   intent( in ) :: i
    sll_int32, optional             ,   intent( in ) :: i_weight
    sll_real64 :: r
   end function i_get_scalar
  end interface


  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_coords( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       ,   intent( in ) :: i
    sll_real64                                       :: r(3)
   end function i_get_coords
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_array( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       ,   intent( in ) :: i
    sll_real64                                       :: r(self%n_weights)
  end function i_get_array
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_coords( self, i, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_int32                         , intent( in    ) :: i
    sll_real64                        , intent( in    ) :: x(3)
   end subroutine i_set_coords
  end interface

!----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_scalar( self, i, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_int32                        , intent( in    ) :: i
    sll_real64                       , intent( in    ) :: x
  end subroutine i_set_scalar
  end interface

!----------------------------------------------------------------------------
  abstract interface
   subroutine set_scalar( self, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_real64                        , intent( in    ) :: x
  end subroutine set_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_array( self, i, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_int32                         , intent( in    ) :: i
    sll_real64                        , intent( in    ) :: x(self%n_weights)
  end subroutine i_set_array
  end interface

 !---------------------------------------------------------------------------!  
  abstract interface
     subroutine empty(self)
       import sll_c_particle_group_base
       class( sll_c_particle_group_base ), intent( inout ) :: self

     end subroutine empty
  end interface

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !----------------------------------------------------------------------------
  pure function q_over_m( self ) result( r )
    class( sll_t_species ), intent( in ) :: self
    sll_real64 :: r

    r = self%q / self%m
  end function q_over_m

  
  !----------------------------------------------------------------------------
  subroutine initialize_species(&
       self, &
       species_charge,     &
       species_mass        &
       )
    class(sll_t_species), intent ( out ) :: self
    sll_real64,           intent ( in )  :: species_charge
    sll_real64,           intent ( in )  :: species_mass


    self%q = species_charge
    self%m = species_mass

  end subroutine initialize_species


  
  !----------------------------------------------------------------------------!
  ! Single particle gets and setters for boxes and patches
  !----------------------------------------------------------------------------!
  
!< Assuming we have some elements/grid cells "box" for the representation of the fields, we get the nd index of this box, by default all are in box 1.
  pure function get_box( self, i ) result( r )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: i    !< particle number

    sll_int32 :: r

    r = 1

  end function get_box

  !< Assuming we have some elements/grid cells "box" for the representation of the fields, we get the nd index of this box (note this is for tensor product grids only), by default all are in box 1
  pure function get_boxnd( self, i ) result( r )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: i    !< particle number

    sll_int32 :: r(3)

    r = 1

  end function get_boxnd

  !< Assuming we have some elements/grid cells "box" for the representation of the fields, we get the (normalized) offset of the position in the box, by default we take the value of x 
  pure function get_xbox( self, i ) result( r )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,                       intent ( in ) :: i    !< particle number

    sll_real64 :: r(3)

    r = self%get_x( i )   

  end function get_xbox
  
  !> Get the number of the patch (per default we only have one patch and therefore always return a 1)
  pure function get_patch( self, i ) result( r )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: i    !< particle number

    sll_int32 :: r

    r = 1

  end function get_patch

  !> Get the coordinate in the patch-local system. Default: \a get_x
  pure function get_xpatch( self, i ) result( r )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,                       intent ( in ) :: i    !< particle number

    sll_real64 :: r(3)

    r = self%get_x( i )   

  end function get_xpatch


  !> Set the number of the box (by default nothing is done)
  subroutine set_box( self, i , box)
    class(sll_c_particle_group_base), intent ( in )   :: self !< particle group object
    sll_int32,                        intent ( in )   :: i    !< particle number
    sll_int32,                        intent (   out) :: box !< box number


  end subroutine set_box

  !> Set the number of the box from nd index (by default nothing is done)
  subroutine set_boxnd( self, i , box)
    class(sll_c_particle_group_base), intent ( in )   :: self !< particle group object
    sll_int32,                        intent ( in )   :: i    !< particle number
    sll_int32,                        intent (   out) :: box(:) !< box number


  end subroutine set_boxnd
  
  
  !> Set the number of the patch (by default nothing is done)
  subroutine set_patch( self, i , box)
    class(sll_c_particle_group_base), intent ( in )   :: self !< particle group object
    sll_int32,                        intent ( in )   :: i    !< particle number
    sll_int32,                        intent (   out) :: box !< box number


  end subroutine set_patch

  !----------------------------------------------------------------------------!
  ! Getters for chunks of particles
  !----------------------------------------------------------------------------!

  !< \a get_x or \a get_xbox for chunk of particle
  subroutine chget_x( self, chunk, x )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: chunk(1:2)    !< particle number range
    sll_int32, pointer,    intent( out ) :: x(:,:)     !< x of particle chunk

    sll_int32 :: i
    
    do i=chunk(1),chunk(2)
       x(:,i-chunk(1)+1) = self%get_x( i )
    end do

  end subroutine chget_x

  !< \a get_v for chunk of particle
  subroutine chget_v( self, chunk, v )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: chunk(1:2)    !< particle number range
    sll_int32, pointer,    intent( out ) :: v(:,:)     !< v of particle chunk

    sll_int32 :: i
    
    do i=chunk(1),chunk(2)
       v(:,i-chunk(1)+1) = self%get_v( i )
    end do

  end subroutine chget_v

  !< \a get_weights for chunk of particle
  subroutine chget_weights( self, chunk, weights )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: chunk(1:2)    !< particle number range
    sll_int32, pointer,    intent( out ) :: weights(:,:)     !< x of particle chunk

    sll_int32 :: i
    
    do i=chunk(1),chunk(2)
       weights(:,i-chunk(1)+1) = self%get_weights( i )
    end do

  end subroutine chget_weights
  
  !< \a get_box for chunk of particle
  subroutine chget_box( self, chunk, box )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: chunk(1:2)    !< particle number range
    sll_int32, pointer,    intent( out ) :: box(:) !< box number of particle chunk

    sll_int32 :: i
    
    do i=chunk(1),chunk(2)
       box(i-chunk(1)+1) = self%get_box( i )
    end do

  end subroutine chget_box

  !< \a get_boxnd for chunk of particle
  subroutine chget_boxnd( self, chunk, boxnd )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: chunk(1:2)    !< particle number range
    sll_int32, pointer,    intent( out ) :: boxnd(:,:)    !< box number of particle chunk

    sll_int32 :: i
    
    do i=chunk(1),chunk(2)
       boxnd(:,i-chunk(1)+1) = self%get_boxnd( i )
    end do

  end subroutine chget_boxnd
  
  !< \a get_patch for chunk of particle
  subroutine chget_patch( self, chunk, patch )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    sll_int32,            intent ( in ) :: chunk(1:2)    !< particle number range
    sll_int32, pointer,    intent( out ) :: patch(:) !< box number of particle chunk

    sll_int32 :: i
    
    do i=chunk(1),chunk(2)
       patch(i-chunk(1)+1) = self%get_patch( i )
    end do

  end subroutine chget_patch

  !> Dummy print function ( can be overwritten to print for debugging )
  subroutine print( self, filename )
    class(sll_c_particle_group_base), intent ( in ) :: self !< particle group object
    character(len=*), intent(in) :: filename

  end subroutine print

  
  !----------------------------------------------------------------------------!
  ! Index conversion helper functions
  !----------------------------------------------------------------------------!
  
  !> Helper function to compute the 2d index for a tensor product grid from the 1d index
  pure function sll_f_index_1dto2d( num_pts, ind1d) result( ind2d )
    sll_int32, intent( in    ) :: num_pts(2)
    sll_int32, intent( in    ) :: ind1d
    sll_int32                  :: ind2d(2)

    ind2d(2) = ( ind1d/num_pts(1) )
    ind2d(1) = ind1d - ind2d(2)*num_pts(1) 
    ind2d(2) = ind2d(2)+1

  end function sll_f_index_1dto2d

  !> Helper function to compute the 1d index for a tensor product grid from the 2d index
  pure function sll_f_index_2dto1d( num_pts, ind2d) result( ind1d )
    sll_int32, intent( in    ) :: num_pts(2)
    sll_int32, intent( in    ) :: ind2d(2)
    sll_int32                  :: ind1d


    ind1d = ind2d(1)+ (ind2d(2)-1) * num_pts(1)

  end function sll_f_index_2dto1d

  
  !> Helper function to compute the 2d index for a tensor product grid from the 1d index
  pure function sll_f_index_1dto3d( num_pts, ind1d) result( ind3d )
    sll_int32, intent( in    ) :: num_pts(3)
    sll_int32, intent( in    ) :: ind1d
    sll_int32                  :: ind3d(3)

    sll_int32 :: ind
    
    ind3d(3) = ( ind1d/(num_pts(1)*num_pts(2)) )
    ind = ind1d - ind3d(3)*(num_pts(1)*num_pts(2))
    ind3d(3) = ind3d(3) +1
    ind3d(2) = ( ind/num_pts(1) )
    ind3d(1) = ind - ind3d(2)*num_pts(1) 
    ind3d(2) = ind3d(2)+1

  end function sll_f_index_1dto3d

  !> Helper function to compute the 1d index for a tensor product grid from the 2d index
  pure function sll_f_index_3dto1d( num_pts, ind3d) result( ind1d )
    sll_int32, intent( in    ) :: num_pts(3)
    sll_int32, intent( in    ) :: ind3d(3)
    sll_int32                  :: ind1d


    ind1d = ind3d(1) + (ind3d(2)-1) * num_pts(1) + (ind3d(3)-1)* num_pts(1)*num_pts(2)

  end function sll_f_index_3dto1d


end module sll_m_particle_group_base
