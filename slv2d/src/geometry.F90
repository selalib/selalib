module geometry_module
  use used_precision
  implicit none
  public :: new
  !private :: new_geometry1, new_geometry2
  type :: geometry
     real(wp) :: x0,y0,x1,y1       ! coordinates of origin
     real(wp) :: dx, dy            ! grid size
     integer     :: nx, ny            ! number of grid points
     real(wp), dimension(:), pointer :: xgrid, ygrid  ! coordinates of points
     character(5) :: bc            ! boundary conditions
  end type geometry
  interface new
     module procedure new_geometry1, new_geometry2
  end interface
  contains
    subroutine new_geometry1(geom,x0,y0,nx,ny,dx,dy,iflag,bc)
      !initialize geometry from origin cordinates, space steps and number of points
      type(geometry) :: geom
      real(wp) :: x0,y0,dx,dy
      integer :: nx, ny
      integer :: iflag ! indicateur d erreur en sortie
      character(5) :: bc 
      ! variables locales
      integer :: i    ! indice de boucle
      integer :: err  ! indicateur d erreur

      ! initialisation de l'indicateur d'erreur
      iflag = 0

      geom%x0=x0
      geom%y0=y0
      geom%dx=dx
      geom%dy=dy
      geom%nx=nx
      geom%ny=ny
      geom%bc=bc

      ! allocation et initialisation des tableaux de coordonnees
      if (associated(geom%xgrid)) deallocate(geom%xgrid)
      allocate(geom%xgrid(nx), stat=err)
      if (err.ne.0) then
         iflag = 10
         return
      end if
      if (associated(geom%ygrid)) deallocate(geom%ygrid)
      allocate(geom%ygrid(ny), stat=err)
      if (err.ne.0) then
         iflag = 20
         return
      end if
      do i=1,nx
         geom%xgrid(i)=geom%x0+(i-1)*geom%dx
      enddo
      do i=1,ny
         geom%ygrid(i)=geom%y0+(i-1)*geom%dy
      enddo
    end subroutine new_geometry1

    subroutine new_geometry2(geom,x0,y0,x1,y1,nx,ny,iflag,bc)
      !initialize geometry from origin and end coordinates and number of points
      type(geometry), intent(inout) :: geom
      real(wp), intent(in) :: x0,y0,x1,y1
      integer, intent(in) :: nx, ny
      integer, intent(out) :: iflag ! indicateur d erreur en sortie
      character(len=*), intent(in) :: bc 
      ! variables locales
      integer :: i    ! indice de boucle
      integer :: err  ! indicateur d erreur

      ! initialisation de l'indicateur d'erreur
      iflag = 0

      geom%x0=x0;geom%x1=x1
      geom%y0=y0;geom%y1=y1
      if ((bc.eq."perxy").or.(bc.eq."perx")) then
         geom%dx=(x1-x0)/nx
      else
         geom%dx=(x1-x0)/(nx-1)
      endif
      if ((bc.eq."perxy").or.(bc.eq."pery")) then
         geom%dy=(y1-y0)/ny
      else
         geom%dy=(y1-y0)/(ny-1)
      endif
      geom%nx=nx
      geom%ny=ny
      geom%bc=bc

!      ! allocation et initialisation des tableaux de coordonnees
!      if (associated(geom%xgrid)) deallocate(geom%xgrid)
!      allocate(geom%xgrid(nx), stat=err)
!      if (err.ne.0) then
!         iflag = 10
!         return
!      end if
!      if (associated(geom%ygrid)) deallocate(geom%ygrid)
!      allocate(geom%ygrid(ny), stat=err)
!      if (err.ne.0) then
!         iflag = 20
!         return
!      end if
!      ! allocation et initialisation des tableaux de coordonnees
!      if (associated(geom%xgrid)) deallocate(geom%xgrid)
      allocate(geom%xgrid(nx), stat=err)
      if (err.ne.0) then
         iflag = 10
         return
      end if
!      if (associated(geom%ygrid)) deallocate(geom%ygrid)
      allocate(geom%ygrid(ny), stat=err)
      if (err.ne.0) then
         iflag = 20
         return
      end if
      do i=1,nx
         geom%xgrid(i)=geom%x0+(i-1)*geom%dx
      enddo
      do i=1,ny
         geom%ygrid(i)=geom%y0+(i-1)*geom%dy
      enddo
    end subroutine new_geometry2      

    subroutine getlinw(this,x,y,inode,jnode,w,iflag)
      ! compute linear weights of point x,y in mesh geom
      type(geometry) :: this
      real(wp) :: x,y
      integer :: inode, jnode  ! localisation of point (x,y) in mesh
      integer :: iflag
      real(wp), dimension(4) :: w  !weights
      ! variables locales
      real(wp) :: x1, y1, x2, y2, aire
      ! INITIALISATIONS
      iflag = 0
      aire = this%dx*this%dy
      ! localize point (x,y) in mesh
      inode = int((x-this%x0)/this%dx)+1
      if ((inode.lt.1).or.(inode.gt.this%nx)) then
         iflag = 10
         return
      end if      
      jnode = int((y-this%y0)/this%dy)+1
      if ((jnode.lt.1).or.(jnode.gt.this%ny)) then
         iflag = 20
         return
      end if
      ! compute weights
      w(1) = (this%xgrid(inode+1)-x)*(this%ygrid(jnode+1)-y)/aire
      w(2) = (x-this%xgrid(inode))*(this%ygrid(jnode+1)-y)/aire
      w(3) = (x-this%xgrid(inode))*(y-this%ygrid(jnode))/aire
      w(4) = (this%xgrid(inode+1)-x)*(y-this%ygrid(jnode))/aire
    end subroutine getlinw
end module geometry_module
