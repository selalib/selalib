
module map_function_module

   use sll_common_coordinate_transformations
   use sll_module_coordinate_transformations_2d

   class(sll_coordinate_transformation_2d_base),pointer :: tau

contains

   subroutine set_map_function(mytau)
   class(sll_coordinate_transformation_2d_base), target :: mytau

   tau => mytau

   end subroutine set_map_function


  ! fonction qui envoie le carré [0,1]x[0,1] sur le vrai domaine de calcul
  ! variables de référence: (u,v)
  ! variables physiques: (x,y)
  ! autres données calculées:
  ! jac, invjac, det: jacobienne, son inverse et déterminant de la jacobienne
  ! subroutine map(u,v,x,y,jac,invjac,det)
  ! pour l'instant on n'utilise pas la jacobienne
  subroutine map(u,v,x,y)
    implicit none
    real(8),intent(in) :: u,v
    real(8),intent(out) :: x,y
!    real(8) :: jac(2,2),invjac(2,2),det
    real(8),parameter :: pi=4*atan(1.d0)
    real(8) :: eta1, eta2

    !x=(1+u)*(1+v)*cos(pi*v)
    !y=(1+u)*sin(pi*v)
    
    eta1 = tau%mesh%eta1_min + u * tau%mesh%delta_eta1
    eta2 = tau%mesh%eta2_min + v * tau%mesh%delta_eta2
    x = tau%x1(eta1,eta2)
    y = tau%x2(eta1,eta2)

!    !x=u
!    !y=v
!    ! non utilisé
!    jac=0
!    jac(1,1)=1
!    jac(2,2)=1
!    invjac=jac
!    det=1
!
  end subroutine map

end module map_function_module

