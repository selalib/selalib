module poisson
#include "sll_working_precision.h"
#include "sll_memory.h"
use zone

implicit none

sll_int32, private :: i, j, ii, jj

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Solveur de Poisson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine E_initial( tm )
 
type( tm_mesh_fields ) :: tm
sll_real64, dimension(:,:), allocatable :: div, lap, phi 
sll_real64, dimension(:), allocatable ::  rho
sll_real64 :: ddx, ddy
sll_int32 :: i, j, indice, info
sll_int32 :: error

! calcul de E verifiant div E = rho/e0

SLL_ALLOCATE( phi( 0:nx, 0:ny ) , error)
SLL_ALLOCATE( rho( (nx-1)*(ny-1) ), error )
SLL_ALLOCATE( div( 0:nx , 0:ny ), error )
SLL_ALLOCATE( lap( ny , (nx-1)*(ny-1) ), error )

! second membre :
! vecteur -(r22,r23,...,r2ny-1,r32,...,r3ny-1,....rnx-1ny-1)/e0

do i=1,nx-1
   do j=1,ny-1
     indice = (i-1)*(ny-1) + j
     rho(indice) = -tm%r1(i,j)/e0
   enddo
enddo

! entree de la matrice de -laplacien

ddx = hx(2)
ddy = hy(2)

lap=0.
do i=1,(nx-1)*(ny-1)
  lap(ny,i)=2/(ddx*ddx)+2/(ddy*ddy)
  lap(ny-1,i)=-1/(ddy*ddy)
  lap(1,i)=-1/(ddx*ddx)
end do
do i=0,nx-2
  lap(ny-1,i*(ny-1)+1)=0
end do
CALL DPBTRF('U',(nx-1)*(ny-1),ny-1,lap,ny,info) 
print*,'factorisation du laplacien',info

CALL DPBTRS('U',(nx-1)*(ny-1),ny-1,1,lap,ny,rho,(nx-1)*(ny-1),info) 
print*,'resolution du laplacien',info

phi=0.
do i=1,nx-1
   do j=1,ny-1
      indice = (i-1)*(ny-1) + j
      phi(i,j) = rho(indice)
   enddo
enddo

! reconstruction de E = grad phi
 
do i=0,nx-1
   do j=0,ny
      tm%ex(i,j)=(phi(i+1,j)-phi(i,j))/ddx
   enddo
enddo

do i=0,nx
   do j=0,ny-1
      tm%ey(i,j)=(phi(i,j+1)-phi(i,j))/ddy
   enddo
enddo

end subroutine E_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine poisson_clnulles( tm )

type(tm_mesh_fields) :: tm
sll_real64, dimension(4,4) :: Axelem, Ayelem, Melem
sll_real64, dimension((nx-1)*(ny-1),(nx-1)*(ny-1)) :: A, M
sll_real64, dimension(nx+1,(nx-1)*(ny-1)) :: mat
sll_real64, dimension((nx-1)*(ny-1)) :: b
sll_real64, dimension(0:nx,0:ny) :: phi
sll_real64 :: dum
sll_int32 :: Iel, info
sll_int32, dimension(4) :: Isom


!** Construction des matrices elementaires
dum = 1.d0/6.d0
Axelem(1,1)=2*dum ; Axelem(1,2)=-2*dum ; Axelem(1,3)=-dum ; Axelem(1,4)=dum ;
Axelem(2,1)=-2*dum ; Axelem(2,2)=2*dum ; Axelem(2,3)=dum ; Axelem(2,4)=-dum ;
Axelem(3,1)=-dum ; Axelem(3,2)=dum ; Axelem(3,3)=2*dum ; Axelem(3,4)=-2*dum ;
Axelem(4,1)=dum ; Axelem(4,2)=-dum ; Axelem(4,3)=-2*dum ; Axelem(4,4)=2*dum ;

Ayelem(1,1)=2*dum ; Ayelem(1,2)=dum ; Ayelem(1,3)=-dum ; Ayelem(1,4)=-2*dum ;
Ayelem(2,1)=dum ; Ayelem(2,2)=2*dum ; Ayelem(2,3)=-2*dum ; Ayelem(2,4)=-dum ;
Ayelem(3,1)=-dum ; Ayelem(3,2)=-2*dum ; Ayelem(3,3)=2*dum ; Ayelem(3,4)=dum ;
Ayelem(4,1)=-2*dum ; Ayelem(4,2)=-dum ; Ayelem(4,3)=dum ; Ayelem(4,4)=2*dum ;

dum = 1.d0/36.d0
Melem(1,1)=4*dum ; Melem(1,2)=2*dum ; Melem(1,3)=dum ; Melem(1,4)=2*dum ;
Melem(2,1)=2*dum ; Melem(2,2)=4*dum ; Melem(2,3)=2*dum ; Melem(2,4)=dum ;
Melem(3,1)=dum ; Melem(3,2)=2*dum ; Melem(3,3)=4*dum ; Melem(3,4)=2*dum ;
Melem(4,1)=2*dum ; Melem(4,2)=dum ; Melem(4,3)=2*dum ; Melem(4,4)=4*dum ;


A = 0.d0
B = 0.d0


!** Contribution des mailles interieures
do i=1,nx-2
   do j=1,ny-2
      Iel = i+(j-1)*(nx-1)
      Isom(1)=Iel; Isom(2)=Iel+1; Isom(3)=Iel+nx; Isom(4)=Iel+nx-1;
      do ii=1,4
         do jj=1,4
            A(Isom(ii),Isom(jj)) = A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * hy(j) / hx(i) &
                 & + Ayelem(ii,jj) * hx(i) / hy(j)
            M(Isom(ii),Isom(jj)) = M(Isom(ii),Isom(jj)) &
                 & + Melem(ii,jj) * hx(i) * hy(j)
         end do
      end do
   end do
end do

!** Contribution des mailles au sud et au nord
do i=1,nx-2
   Isom(3)=i+1; Isom(4)=i  !Sud
   do ii=3,4
      do jj=3,4
         A(Isom(ii),Isom(jj)) = A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * hy(0) / hx(i) &
                 & + Ayelem(ii,jj) * hx(i) / hy(0)
         M(Isom(ii),Isom(jj)) = M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * hx(i) * hy(0)
      end do
   end do
   Iel = (ny-2)*(nx-1)+i   !Nord
   Isom(1)=Iel; Isom(2)=Iel+1
   do ii=1,2
      do jj=1,2
         A(Isom(ii),Isom(jj)) = A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * hy(ny-1) / hx(i) &
                 & + Ayelem(ii,jj) * hx(i) / hy(ny-1)
         M(Isom(ii),Isom(jj)) = M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * hx(i) * hy(ny-1)
      end do
   end do
end do

!** Contribution des mailles a l'ouest et a l'est
do j=1,ny-2
   Isom(2)=1+(j-1)*(nx-1); Isom(3)=1+j*(nx-1) !Ouest
   do ii=2,3
      do jj=2,3
         A(Isom(ii),Isom(jj)) = A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * hy(j) / hx(0) &
                 & + Ayelem(ii,jj) * hx(0) / hy(j)
         M(Isom(ii),Isom(jj)) = M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * hx(0) * hy(j)
      end do
   end do
   Iel = j*(nx-1)                              !Est
   Isom(1)=Iel; Isom(4)=Iel+nx-1
   do ii=1,4,3
      do jj=1,4,3
         A(Isom(ii),Isom(jj)) = A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * hy(j) / hx(nx-1) &
                 & + Ayelem(ii,jj) * hx(nx-1) / hy(j)
         M(Isom(ii),Isom(jj)) = M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * hx(nx-1) * hy(j)
      end do
   end do
end do

!** Contribution des coins
Isom(3) = 1    !SW
A(1,1) = A(1,1) + Axelem(3,3) * hy(0) / hx(0) + Ayelem(3,3) * hx(0) / hy(0)
M(1,1) = M(1,1) + Melem(3,3) * hx(0) * hy(0)

Isom(4) = nx-1 !SE
A(nx-1,nx-1) = A(nx-1,nx-1) + Axelem(4,4) * hy(0) / hx(nx-1) &
     & + Ayelem(4,4) * hx(nx-1) / hy(0)
M(nx-1,nx-1) = M(nx-1,nx-1) + Melem(4,4) * hx(nx-1) * hy(0)

Isom(1) = (nx-1)*(ny-1)   !NE
A(Isom(1),Isom(1)) = A(Isom(1),Isom(1)) + Axelem(1,1) * hy(ny-1) / hx(nx-1) &
     & + Ayelem(1,1) * hx(nx-1) / hy(ny-1)
M(Isom(1),Isom(1)) = M(Isom(1),Isom(1)) + Melem(1,1) * hx(nx-1) * hy(ny-1)

Isom(2) = (nx-1)*(ny-2)+1 !NW
A(Isom(2),Isom(2)) = A(Isom(2),Isom(2)) + Axelem(2,2) * hy(ny-1) / hx(0) &
     & + Ayelem(2,2) * hx(0) / hy(ny-1)
M(Isom(2),Isom(2)) = M(Isom(2),Isom(2)) + Melem(2,2) * hx(0) * hy(ny-1)


!** Construction du second membre (rho a support compact --> projete)
do i=1,nx-1
   do j=1,ny-1
      b(i+(j-1)*(nx-1)) = tm%r1(i,j)/e0
   end do
end do

b = matmul(M,b)

!** Apres la construction des matrices de masse et de raideur,
!** il faut encore inverser le systeme matriciel  A*u=b

mat = 0.d0

mat(nx+1,1) = A(1,1)
do j=2,(nx-1)*(ny-1)
   mat(nx+1,j) = A(j,j)
   mat(nx,j)   = A(j-1,j)
end do
mat(3,nx) = A(2,nx)
mat(2,nx) = A(1,nx)
do j=nx+1,(nx-1)*(ny-1)
   mat(3,j) = A(j-nx+2,j)
   mat(2,j) = A(j-nx+1,j)
   mat(1,j) = A(j-nx,j)
end do
CALL DPBTRF('U',(nx-1)*(ny-1),nx,mat,nx+1,info)
print*,'factorisation pour Cholesky',info

CALL DPBTRS('U',(nx-1)*(ny-1),nx,1,mat,nx+1,b,(nx-1)*(ny-1),info) 
print*,'resolution par Cholesky',info

phi = 0.d0
do i=1,nx-1
   do j=1,ny-1
      phi(i,j) = b(i+(j-1)*(nx-1)) 
   end do
end do

!** Reconstruction du champ E

do i=0,nx-1
   do j=0,ny
      tm%ex(i,j) = - (phi(i+1,j)-phi(i,j)) / hx(i)
   end do
end do

do i=0,nx
   do j=0,ny-1
      tm%ey(i,j) = - (phi(i,j+1)-phi(i,j)) / hy(j)
   end do
end do

end subroutine poisson_clnulles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module poisson
