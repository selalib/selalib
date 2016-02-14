!
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
module sll_m_maxwell_diga_hex_mesh

#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_m_hexagonal_meshes, only : &
       sll_hex_mesh_2d
use sll_m_triangle_dg_matrices, only : &
     ElementRef, &
     AssMatElem

implicit none
private

sll_real64, parameter :: xi = 0._f64
sll_real64, parameter :: c  = 1._f64

type, public :: maxwell_dg_hex_mesh

   sll_int32 :: degree
   sll_int32 :: n_ddl

   type(ElementRef) :: Elem

   sll_real64, dimension(:,:,:), pointer :: DxP
   sll_real64, dimension(:,:,:), pointer :: DyP

   sll_real64, dimension(:,:), pointer :: Ex
   sll_real64, dimension(:,:), pointer :: Ey
   sll_real64, dimension(:,:), pointer :: Bz
   sll_real64, dimension(:,:), pointer :: Jx
   sll_real64, dimension(:,:), pointer :: Jy
   sll_real64, dimension(:,:), pointer :: Ro

   sll_real64, dimension(:,:), pointer :: D_Ex
   sll_real64, dimension(:,:), pointer :: D_Ey
   sll_real64, dimension(:,:), pointer :: D_Bz

   sll_real64, dimension(:,:), pointer :: Po
   sll_real64, dimension(:,:), pointer :: D_Po

   sll_real64, dimension(:,:), pointer :: x_ddl
   sll_real64, dimension(:,:), pointer :: y_ddl

   sll_int32,  dimension(:,:), pointer :: ntri
   sll_int32,  dimension(:,:), pointer :: nvois
   sll_int32,  dimension(:,:), pointer :: nvoif

end type maxwell_dg_hex_mesh

sll_int32, private :: iplot

interface sll_create
   module procedure init_gd_solver_2d
end interface sll_create

interface sll_solve
   module procedure gd_solver_2d_te
end interface sll_solve

public :: initialize, solve

contains

subroutine init_gd_solver_2d(self, mesh, degree)

type(maxwell_dg_hex_mesh), intent(inout) :: self
type(sll_hex_mesh_2d), pointer, intent(in)   :: mesh
sll_int32, intent(in)                      :: degree
sll_int32 :: idl, error, iel
sll_real64 :: xs1, ys1, xs2, ys2, xs3, ys3
sll_real64 :: det
sll_real64, dimension(:), pointer :: xref, yref
sll_int32 :: k, n_ddl
sll_int32, dimension(:),   allocatable :: npoel1
sll_int32, dimension(:),   allocatable :: npoel2
sll_int32, dimension (:),  allocatable :: indc
sll_int32 :: i1, i2, i3, is, is1, is2, is3
sll_int32 :: jel1, jel2, jel3, nel1, nel2, nel3
sll_real64 :: x1, x2
sll_int32 :: i, j

self%degree  = degree
n_ddl = (degree+1)*(degree+2)/2
self%n_ddl = n_ddl
SLL_ALLOCATE(xref(n_ddl), error)
SLL_ALLOCATE(yref(n_ddl), error)

SLL_ALLOCATE(self%Elem%MassMat(n_ddl,n_ddl),error)
if (error.ne.0) print*, 'error in allocation of Elem%MassMat'
self%Elem%MassMat = 0._f64
SLL_ALLOCATE(self%Elem%MassMatInv(n_ddl,n_ddl),error)
if (error.ne.0) print*, 'error in allocation of Elem%MassMatInv'
self%Elem%MassMatInv = 0._f64
SLL_ALLOCATE(self%Elem%DxMat(n_ddl,n_ddl),error)
if (error.ne.0) print*, 'error in allocation of Elem%DxMat'
self%Elem%DxMat = 0._f64
SLL_ALLOCATE(self%Elem%DyMat(n_ddl,n_ddl),error)
if (error.ne.0) print*, 'error in allocation of Elem%DyMat'
self%Elem%DyMat = 0._f64
SLL_ALLOCATE(self%Elem%BndMassMat(degree+1,degree+1),error)
if (error.ne.0) print*, 'error in allocation of Elem%BndMassMat'
self%Elem%BndMassMat = 0._f64

call AssMatElem(self%Elem,degree,xref,yref)

SLL_CLEAR_ALLOCATE(self%Ex(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%Ey(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%Bz(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%Po(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%Jx(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%Jy(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%Ro(1:n_ddl,1:mesh%num_triangles),error)

SLL_CLEAR_ALLOCATE(self%D_Ex(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%D_Ey(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%D_Bz(1:n_ddl,1:mesh%num_triangles),error)
SLL_CLEAR_ALLOCATE(self%D_Po(1:n_ddl,1:mesh%num_triangles),error)

!Les matrices derivees
SLL_CLEAR_ALLOCATE(self%DxP(1:n_ddl,n_ddl,1:mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(self%DyP(1:n_ddl,n_ddl,1:mesh%num_triangles), error)

!Determinant de l'element et coordonnees des D.D.L.
SLL_ALLOCATE(self%x_ddl(n_ddl,mesh%num_triangles),error)
SLL_ALLOCATE(self%y_ddl(n_ddl,mesh%num_triangles),error)

SLL_ALLOCATE(self%ntri(3,mesh%num_triangles),error)

do iel = 1, mesh%num_triangles    !Boucle sur les elements

   x1 = mesh%center_cartesian_coord(1, iel)
   x2 = mesh%center_cartesian_coord(2, iel)

   call sll_s_get_cell_vertices_index( x1, x2, mesh, is1, is2, is3)

   !Coordonnees des sommets du triangle
   xs1 = mesh%global_to_x1(is1)
   ys1 = mesh%global_to_x2(is1)
   xs2 = mesh%global_to_x1(is2)
   ys2 = mesh%global_to_x2(is2)
   xs3 = mesh%global_to_x1(is3)
   ys3 = mesh%global_to_x2(is3)

   !Les triangles sont tous identiques det = delta*delta*sqrt(3)/2
   !Calcul du determinant = aire du triangle
   det = (xs2-xs1)*(ys3-ys1)-(xs3-xs1)*(ys2-ys1)

   if (det > 0) then
      self%ntri(:,iel) = [is1,is2,is3]
   else
      self%ntri(:,iel) = [is1,is3,is2]
   end if

   xs1 = mesh%global_to_x1(self%ntri(1,iel))
   ys1 = mesh%global_to_x2(self%ntri(1,iel))
   xs2 = mesh%global_to_x1(self%ntri(2,iel))
   ys2 = mesh%global_to_x2(self%ntri(2,iel))
   xs3 = mesh%global_to_x1(self%ntri(3,iel))
   ys3 = mesh%global_to_x2(self%ntri(3,iel))

   det = (xs2-xs1)*(ys3-ys1)-(xs3-xs1)*(ys2-ys1)

   !Positions absolues des D.D.L
   do idl = 1, n_ddl
      self%x_ddl(idl,iel) = xs1 + (xs2-xs1)*xref(idl) + (xs3-xs1)*yref(idl)
      self%y_ddl(idl,iel) = ys1 + (ys2-ys1)*xref(idl) + (ys3-ys1)*yref(idl)
   end do

   self%DxP(:,:,iel) = ((ys3-ys1)*self%Elem%DxMat+(ys1-ys2)*self%Elem%DyMat)/det
   self%DyP(:,:,iel) = ((xs1-xs3)*self%Elem%DxMat+(xs2-xs1)*self%Elem%DyMat)/det

end do

! ... recherche des elements ayant un sommet commun
!     creation du tableau npoel1(i+1)  contenant le nombre de 
!     triangles ayant le noeud i en commun

SLL_ALLOCATE(npoel1(mesh%num_pts_tot+1),error)

npoel1 = 0
do iel = 1,mesh%num_triangles

   npoel1(self%ntri(1,iel)+1) = npoel1(self%ntri(1,iel)+1)+1
   npoel1(self%ntri(2,iel)+1) = npoel1(self%ntri(2,iel)+1)+1
   npoel1(self%ntri(3,iel)+1) = npoel1(self%ntri(3,iel)+1)+1

end do

! ... le tableau npoel1 devient le tableau donnant l'adresse 
!     dans npoel2 du dernier element dans la suite des triangles
!     communs a un noeud

npoel1(1)=0
do i=3,mesh%num_pts_tot+1
   npoel1(i)=npoel1(i-1)+npoel1(i)
end do

! ... creation du tableau npoel2 contenant sequentiellement les 
!     numeros des triangles ayant un noeud en commun      
!     le premier triangle s'appuyant sur le noeud i est
!     adresse par "npoel1(i)+1" 
!     le nombre de triangles ayant le noeud i en commun est
!     "npoel1(i+1)-npoel1(i)"


SLL_ALLOCATE(npoel2(npoel1(mesh%num_pts_tot+1)),error)
SLL_ALLOCATE(indc(mesh%num_pts_tot),error)

indc   = 1  !Le tableau temporaire indc doit etre initialise a 1

do iel = 1,mesh%num_triangles

   do k = 1,3
      is = self%ntri(k,iel)
      npoel2(npoel1(is)+indc(is)) = iel
      indc(is) = indc(is)+1
   end do

end do
 
! --- Recherche des numeros des triangles voisins d'un triangle 

write(*,*)"*** Compute neighbors ***"

SLL_ALLOCATE(self%nvois(3,mesh%num_triangles),error); self%nvois = -1
SLL_ALLOCATE(self%nvoif(3,mesh%num_triangles),error); self%nvoif = -1

do iel=1,mesh%num_triangles

   ! ... numeros des 3 sommets du triangle
   is1 = self%ntri(1,iel)
   is2 = self%ntri(2,iel)
   is3 = self%ntri(3,iel)

   ! ... boucles imbriquees sur les elements pointant vers
   !     les 2 noeuds extremites de l'arete consideree
   !     Le voisin est le triangle commun (hormis iel)

   ! ... premiere arete (entre le sommet is1 et is2)

   nel1=npoel1(is1+1)-npoel1(is1) !nb de triangles communs a is1
   nel2=npoel1(is2+1)-npoel1(is2) !nb de triangles communs a is2

   loop1:do i1=1,nel1
     jel1=npoel2(npoel1(is1)+i1) !premier triangle is1
     if(jel1.ne.iel) then
       do i2=1,nel2
         jel2=npoel2(npoel1(is2)+i2)
         if(jel2 == jel1) then
           self%nvois(1,iel)  = jel1
           exit loop1
         end if
       end do
     end if
   end do loop1

   ! ... deuxieme arete (entre le sommet is2 et is3)

   nel2=npoel1(is2+1)-npoel1(is2)
   nel3=npoel1(is3+1)-npoel1(is3)

   loop2:do i2=1,nel2
     jel2=npoel2(npoel1(is2)+i2)
     if(jel2 /= iel) then
       do i3=1,nel3
         jel3=npoel2(npoel1(is3)+i3)
         if(jel3 == jel2) then
           self%nvois(2,iel)=jel2
           exit loop2
         end if
       end do
     end if
   end do loop2

   ! ... troisieme arete (entre le sommet is3 et is1)

   nel3=npoel1(is3+1)-npoel1(is3)
   nel1=npoel1(is1+1)-npoel1(is1)

   loop3:do i3=1,nel3
     jel3=npoel2(npoel1(is3)+i3)
     if(jel3 /= iel) then
       do i1=1,nel1
         jel1=npoel2(npoel1(is1)+i1)
         if(jel1 == jel3) then
           self%nvois(3,iel)=jel3
           exit loop3
         end if
       end do
     end if
   end do loop3

end do

!Calcul de nvoif : Numero local dans le voisin de la face j commune a l'elt i

self%nvoif = self%nvois

do i = 1, mesh%num_triangles
  do j = 1,3
    jel1 = self%nvois(j,i)
    if (jel1 > 0) then
      do k = 1,3
        jel2 = self%nvois(k,jel1)
        if (jel2 == i) then 
          self%nvoif(j,i) = k
          exit
        end if
      end do
    end if
  end do
end do

end subroutine init_gd_solver_2d

subroutine gd_solver_2d_TE(self, mesh)

type(maxwell_dg_hex_mesh),  intent(inout) :: self
type(sll_hex_mesh_2d), pointer, intent(in)    :: mesh

sll_int32 :: idl, jdl, i, j, l, jdv, i1, i2
sll_int32 :: iel, ifl, iev, ifv, ief
sll_int32 :: n_ddl

sll_real64 :: n1, n2
sll_real64 :: W_l(4), W_r(4)
sll_real64 :: A_m(4,4), A_p(4,4)
sll_real64 :: det
sll_real64 :: flux(self%degree+1,4)
sll_real64 :: xs(3), ys(3)
sll_real64 :: Esn(self%degree+1)

n_ddl = (self%degree+1)*(self%degree+2)/2

do iel = 1, mesh%num_triangles   !Boucle sur les elements

  self%D_Ex(:,iel) = c*c*( - matmul(self%DyP(:,:,iel),self%Bz(:,iel)) &
                   + xi*xi * matmul(self%DxP(:,:,iel),self%Po(:,iel)))
  self%D_Ey(:,iel) = c*c*( + matmul(self%DxP(:,:,iel),self%Bz(:,iel)) &
                   + xi*xi * matmul(self%DyP(:,:,iel),self%Po(:,iel)))
  self%D_Bz(:,iel) =    -    matmul(self%DyP(:,:,iel),self%Ex(:,iel)) &
                        +    matmul(self%DxP(:,:,iel),self%Ey(:,iel))
  self%D_Po(:,iel) = xi*xi * matmul(self%DxP(:,:,iel),self%Ex(:,iel)) &
                   + xi*xi * matmul(self%DyP(:,:,iel),self%Ey(:,iel))

  self%Ro(:,iel)   = matmul(self%DxP(:,:,iel),self%Ex(:,iel)) &
                   + matmul(self%DyP(:,:,iel),self%Ey(:,iel))

  xs(1) = mesh%global_to_x1(self%ntri(1,iel))
  ys(1) = mesh%global_to_x2(self%ntri(1,iel))
  xs(2) = mesh%global_to_x1(self%ntri(2,iel))
  ys(2) = mesh%global_to_x2(self%ntri(2,iel))
  xs(3) = mesh%global_to_x1(self%ntri(3,iel))
  ys(3) = mesh%global_to_x2(self%ntri(3,iel))

  do ifl = 1, 3   !Boucle sur les faces

    i1 = mod(ifl-1,3)+1   !Premier sommet
    i2 = mod(ifl  ,3)+1   !Deuxieme sommet
     
    n1  = ys(i2)-ys(i1)
    n2  = xs(i1)-xs(i2)

    A_p(1,:) = [ (n2*n2+xi*n1*n1)*c,     n2*n1*(xi-1)*c,  -n2*c*c, n1*c*c]
    A_p(2,:) = [     n2*n1*(xi-1)*c, (n1*n1+xi*n2*n2)*c,   n1*c*c, n2*c*c]
    A_p(3,:) = [                -n2,                 n1,        c,    0._f64]
    A_p(4,:) = [           n1*xi*xi,           n2*xi*xi,      0._f64,   xi*c]

    A_m(1,:) = [-(n2*n2+xi*n1*n1)*c,    -n2*n1*(xi-1)*c,  -n2*c*c, n1*c*c]
    A_m(2,:) = [    -n2*n1*(xi-1)*c,-(n1*n1+xi*n2*n2)*c,   n1*c*c, n2*c*c]
    A_m(3,:) = [                -n2,                 n1,       -c,    0._f64]
    A_m(4,:) = [           n1*xi*xi,           n2*xi*xi,      0._f64,  -xi*c]

    A_p = 0.5 * A_p
    A_m = 0.5 * A_m

    iev = self%nvois(ifl,iel)         !Numero du voisin
    ifv = self%nvoif(ifl,iel)         !Numero de face dans le voisin

    if ( iev > 0 ) then               !Cote interne

      do idl = 1, self%degree+1                !Boucle sur les D.D.L.
    
        jdl = iddl_local(ifl,idl,self%degree)   !indice local du DDL
     
        W_l(1) = self%Ex(jdl,iel)
        W_l(2) = self%Ey(jdl,iel)
        W_l(3) = self%Bz(jdl,iel)
        W_l(4) = self%Po(jdl,iel)

        jdv = iddl_voisin(ifv,idl,self%degree)  !indice dans le voisin 

        W_r(1) = self%Ex(jdv,iev)
        W_r(2) = self%Ey(jdv,iev)
        W_r(3) = self%Bz(jdv,iev)
        W_r(4) = self%Po(jdv,iev)

        flux(idl,:) = matmul(A_p,w_l)+matmul(A_m,w_r)
        Esn(idl) = 0.5*(n1*(W_l(1)+W_r(1))+n2*(W_l(2)+W_r(2)))

      end do

    else      !Cote frontiere

      ief = 1 !- iev

      do idl = 1, self%degree+1

            jdl = iddl_local(ifl,idl,self%degree)   !indice local du DDL

            W_l(1) = self%Ex(jdl,iel)
            W_l(2) = self%Ey(jdl,iel)
            W_l(3) = self%Bz(jdl,iel)
            W_l(4) = self%Po(jdl,iel)
     
            if (ief == 1) then          !Conducteur parfait E x n = 0, B . n = 0
               W_r(1) = - W_l(1)
               W_r(2) = - W_l(2)
               W_r(3) =   W_l(3)
               W_r(4) = - W_l(4) - 2*xi/c*(n1*W_l(1)+n2*W_l(2))
               flux(idl,:) = matmul(A_p,w_l)+matmul(A_m,w_r)
            else if (ief == 3) then
               flux(idl,:) = matmul(A_p,W_l)
            else
               flux = 0._f64
            end if

            Esn(idl) = n1*W_l(1)+n2*W_l(2)
     
         end do

      end if
      
      det = (xs(2)-xs(1))*(ys(3)-ys(1))-(xs(3)-xs(1))*(ys(2)-ys(1))

      flux(:,1) = matmul(self%Elem%BndMassMat,flux(:,1))/det
      flux(:,2) = matmul(self%Elem%BndMassMat,flux(:,2))/det
      flux(:,3) = matmul(self%Elem%BndMassMat,flux(:,3))/det
      flux(:,4) = matmul(self%Elem%BndMassMat,flux(:,4))/det

      Esn = matmul(self%Elem%BndMassMat,Esn)/det

      do i = 1, n_ddl
         do j = 1, self%degree+1
            l = iddl_local(ifl,j,self%degree)
            self%D_Ex(i,iel)=self%D_Ex(i,iel)-self%Elem%MassMatInv(i,l)*flux(j,1)
            self%D_Ey(i,iel)=self%D_Ey(i,iel)-self%Elem%MassMatInv(i,l)*flux(j,2)
            self%D_Bz(i,iel)=self%D_Bz(i,iel)-self%Elem%MassMatInv(i,l)*flux(j,3)
            self%D_Po(i,iel)=self%D_Po(i,iel)-self%Elem%MassMatInv(i,l)*flux(j,4)

            self%Ro(i,iel)=self%Ro(i,iel)-self%Elem%MassMatInv(i,l)*Esn(j)

         end do
      end do

   end do

end do

end subroutine gd_solver_2d_TE

!Function:iddl_local
sll_int32 function iddl_local(ifl,ik,degree)
sll_int32 :: ifl, ik, degree

if (ifl == 1) then
   iddl_local = ik 
else if (ifl == 2) then
   iddl_local = degree+ik
else if (ifl == 3) then
   iddl_local = 2*degree+ik
end if

if (iddl_local == 3*degree+1) then
   iddl_local = 1
end if

end function iddl_local

!Function:iddl_voisin
sll_int32 function iddl_voisin(ifv,ik, degree)
sll_int32 :: ifv, ik, idv, degree

if (ifv == 1) then
   idv = (degree+1)-(ik-1)
else if (ifv == 2) then
   idv = (2*degree+1)-(ik-1)
else if (ifv == 3) then
   idv = (3*degree+1)-(ik-1)
end if

if (idv == 3*degree+1) then
   iddl_voisin = 1
else
   iddl_voisin = idv
end if

end function iddl_voisin

end module sll_m_maxwell_diga_hex_mesh
