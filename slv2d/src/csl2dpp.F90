module csl2dpp_class
  use used_precision
  use geometry_module
  !use clock
  implicit none
  private
  public :: new, interpole
  type, public :: csl2dpp
     type (geometry) :: geom
  end type csl2dpp
  interface new
     module procedure new_csl2dpp
  end interface
  interface interpole
     module procedure interpole_csl2dpp,interpole_csl2dppdep
  end interface
contains
  subroutine new_csl2dpp(this,geom,iflag)
    type(csl2dpp), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag

    ! initialisation de geom
    this%geom = geom


 end subroutine new_csl2dpp

  subroutine interpole_csl2dpp(this,fin,fout,x,y) 
    type(csl2dpp), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille precedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destine a contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points dans les quels on veut 
    ! evaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout on trouve en sortie les valeurs de f(i,j) 
    ! dans les points x(i),y(i).
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr

    

  end subroutine interpole_csl2dpp

  subroutine interpole_csl2dppdep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(csl2dpp), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx, depy 
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    logical  :: aff

    !if (aff) then 
    !   call clck_temps(l_a)
    !end if

  end subroutine interpole_csl2dppdep




end module csl2dpp_class
