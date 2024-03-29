module sll_m_map_function

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   use sll_m_cartesian_meshes, only: &
      sll_t_cartesian_mesh_2d

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base

   implicit none

   public :: &
      sll_s_map, &
      sll_s_set_map_function

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   class(sll_c_coordinate_transformation_2d_base), pointer :: tau

contains

   subroutine sll_s_set_map_function(mytau)
      class(sll_c_coordinate_transformation_2d_base), target :: mytau

      tau => mytau

   end subroutine sll_s_set_map_function

   ! fonction qui envoie le carré [0,1]x[0,1] sur le vrai domaine de calcul
   ! variables de référence: (u,v)
   ! variables physiques: (x,y)
   ! autres données calculées:
   ! jac, invjac, det: jacobienne, son inverse et déterminant de la jacobienne
   ! subroutine sll_s_map(u,v,x,y,jac,invjac,det)
   ! pour l'instant on n'utilise pas la jacobienne
   subroutine sll_s_map(u, v, x, y)
      implicit none
      real(8), intent(in) :: u, v
      real(8), intent(out) :: x, y
      real(8) :: eta1, eta2

      !x=(1+u)*(1+v)*cos(pi*v)
      !y=(1+u)*sin(pi*v)
      associate (mesh => tau%mesh)

         eta1 = mesh%eta1_min + u*(mesh%eta1_max - mesh%eta1_min)
         eta2 = mesh%eta2_min + v*(mesh%eta2_max - mesh%eta2_min)
         x = tau%x1(eta1, eta2)
         y = tau%x2(eta1, eta2)

      end associate

   end subroutine sll_s_map

end module sll_m_map_function

