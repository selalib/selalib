module geometry1d_module
  use used_precision
  implicit none
  private :: new_geometry1, new_geometry2
  public :: new
  type :: geometry1d
     real(wp) :: x0            ! coordinates of origin
     real(wp) :: dx            ! grid size
     integer     :: nx            ! number of grid points
     real(wp), dimension(:), pointer :: xgrid  ! coordinates of points
     logical :: bc   ! conditions au bord
  end type geometry1d
  interface new
     module procedure new_geometry1, new_geometry2
  end interface
  contains
    subroutine new_geometry1(geom,x0,nx,dx,bc,iflag)
      !initialize geometry from origin cordinates, space steps and number of points
      type(geometry1d) :: geom
      real(wp) :: x0,dx
      integer :: nx
      integer :: iflag ! indicateur d erreur en sortie
      ! variables locales
      integer :: i    ! indice de boucle
      integer :: err  ! indicateur d erreur
      logical  :: bc

      ! initialisation de l'indicateur d'erreur
      iflag = 0
      geom%bc = bc

      geom%x0=x0  ; geom%dx=dx ; geom%nx=nx

      ! allocation et initialisation des tableaux de coordonnees
      allocate(geom%xgrid(nx), stat=err)
      if (err.ne.0) then
         iflag = 10
         return
      end if
      do i=1,nx
         geom%xgrid(i)=geom%x0+(i-1)*geom%dx
      enddo
    end subroutine new_geometry1

    subroutine new_geometry2(geom,x0,x1,nx,bc,iflag)
      !initialize geometry from origin and end coordinates and number of points
      type(geometry1d) :: geom
      real(wp) :: x0,x1
      integer :: nx
      integer :: iflag ! indicateur d erreur en sortie
      ! variables locales
      integer :: i    ! indice de boucle
      integer :: err  ! indicateur d erreur
      logical :: bc

      ! initialisation de l'indicateur d'erreur
      iflag = 0
      geom%bc = bc ;  geom%x0=x0 ;  geom%nx=nx
      if (bc) then
         geom%dx=(x1-x0)/nx
      else
         geom%dx=(x1-x0)/(nx-1)
      end if

      ! allocation et initialisation du tableau des ordonnees
      allocate(geom%xgrid(nx), stat=err)
      if (err.ne.0) then
         iflag = 10
         return
      end if
      do i=1,nx
         geom%xgrid(i)=geom%x0+(i-1)*geom%dx
      enddo
    end subroutine new_geometry2      

    subroutine getlinw(this,x,inode,w,iflag)
      ! compute linear weights of point x in mesh geom
      type(geometry1d) :: this
      real(wp) :: x
      integer :: inode ! localisation of point (x) in mesh
      integer :: iflag
      real(wp), dimension(2) :: w  !weights
      ! variables locales
      real(wp) :: x1,x2,long
      ! INITIALISATIONS
      iflag = 0
      long = this%dx
      ! localize point (x) in mesh
      inode = int((x-this%x0)/this%dx)+1
      if ((inode.lt.1).or.(inode.gt.this%nx)) then
         iflag = 10
         return
      end if      

      ! compute weights
      w(1) = (this%xgrid(inode+1)-x)/long
      w(2) = (x-this%xgrid(inode))/long
    end subroutine getlinw

end module geometry1d_module
